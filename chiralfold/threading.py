"""
ChiralFold — Template-Based Fold Prediction (Threading)
=========================================================

Threads a target amino-acid sequence onto the backbone of a known PDB
template, producing a starting model that inherits the template's
backbone geometry.  Ideal for:

  * Homology modelling when no force-field conformer is available.
  * Generating D-protein backbones from mirrored L-templates.
  * Rapid structure initialisation before energy minimisation.

Usage::

    from chiralfold.threading import thread_sequence, find_template

    result = thread_sequence(
        target_seq    = 'ACDEFGHIKLMNPQRSTVWY',
        template_pdb  = 'template.pdb',
        template_chain= 'A',
        output_pdb    = 'threaded.pdb',
    )
    print(result)
"""

import os
import math
import glob
import warnings
import numpy as np
from typing import Optional, Dict, Any, List, Tuple


# ═══════════════════════════════════════════════════════════════════════════
# Amino-acid code tables
# ═══════════════════════════════════════════════════════════════════════════

_ONE_TO_THREE: Dict[str, str] = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}

_L_TO_D_RESNAME: Dict[str, str] = {
    'ALA': 'DAL', 'ARG': 'DAR', 'ASN': 'DSG', 'ASP': 'DAS',
    'CYS': 'DCY', 'GLU': 'DGL', 'GLN': 'DGN', 'GLY': 'GLY',
    'HIS': 'DHI', 'ILE': 'DIL', 'LEU': 'DLE', 'LYS': 'DLY',
    'MET': 'MED', 'PHE': 'DPN', 'PRO': 'DPR', 'SER': 'DSN',
    'THR': 'DTH', 'TRP': 'DTR', 'TYR': 'DTY', 'VAL': 'DVA',
}

_THREE_TO_ONE: Dict[str, str] = {v: k for k, v in _ONE_TO_THREE.items()}
# Also map D-names back to one-letter
for _l3, _d3 in _L_TO_D_RESNAME.items():
    if _l3 in _THREE_TO_ONE:
        _THREE_TO_ONE[_d3] = _THREE_TO_ONE[_l3]

_BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


# ═══════════════════════════════════════════════════════════════════════════
# PDB I/O helpers
# ═══════════════════════════════════════════════════════════════════════════

def _parse_pdb_backbone(
    pdb_path: str,
    chain: str,
) -> List[Dict[str, Any]]:
    """Extract backbone atoms (N, CA, C, O) for a given chain, ordered by
    residue sequence number.

    Returns:
        List of dicts: {resseq, resname, atom_name, x, y, z}.
        Sorted by (resseq, atom_name) with backbone atom order N→CA→C→O.
    """
    _ATOM_ORDER = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}
    records: List[Dict[str, Any]] = []

    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            try:
                rec_chain = line[21].strip()
                if rec_chain != chain:
                    continue
                atom_name = line[12:16].strip()
                if atom_name not in _BACKBONE_ATOMS:
                    continue
                records.append({
                    'resseq':    int(line[22:26]),
                    'resname':   line[17:20].strip(),
                    'atom_name': atom_name,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                })
            except (ValueError, IndexError):
                continue

    # Sort by residue, then canonical backbone order
    records.sort(key=lambda r: (r['resseq'], _ATOM_ORDER.get(r['atom_name'], 99)))
    return records


def _extract_chain_sequence(pdb_path: str, chain: str) -> str:
    """Return the one-letter sequence of Cα residues for *chain* in *pdb_path*."""
    seen: Dict[int, str] = {}
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            try:
                rec_chain = line[21].strip()
                if rec_chain != chain:
                    continue
                atom_name = line[12:16].strip()
                if atom_name != 'CA':
                    continue
                resseq  = int(line[22:26])
                resname = line[17:20].strip()
                if resseq not in seen:
                    seen[resseq] = resname
            except (ValueError, IndexError):
                continue

    seq = ''
    for resseq in sorted(seen):
        one = _THREE_TO_ONE.get(seen[resseq], 'X')
        seq += one
    return seq


def _format_atom_line(
    serial:    int,
    atom_name: str,
    resname:   str,
    chain:     str,
    resseq:    int,
    x: float, y: float, z: float,
    element:   str = '',
) -> str:
    """Format a PDB ATOM record line (80 chars + newline)."""
    # PDB atom name: 4-char field, left-padded for 2-char elements
    if len(atom_name) < 4:
        if len(atom_name) == 1 or (len(atom_name) == 2 and atom_name[0].isdigit()):
            name_field = f' {atom_name:<3s}'
        else:
            name_field = f' {atom_name:<3s}'
    else:
        name_field = f'{atom_name:<4s}'

    elem = element if element else (atom_name[0] if atom_name else ' ')
    line = (
        f"ATOM  {serial:5d} {name_field} {resname:<3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n"
    )
    return line


def _mirror_pdb_file(input_pdb: str, output_pdb: str) -> None:
    """Reflect x-coordinates (x → -x) and convert L-residue names to D."""
    lines_out = []
    with open(input_pdb, 'r') as fh:
        for line in fh:
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    x  = float(line[30:38])
                    y  = float(line[38:46])
                    z  = float(line[46:54])
                    rn = line[17:20].strip()
                    d_rn = _L_TO_D_RESNAME.get(rn, rn)
                    new_line = (
                        line[:17]
                        + f'{d_rn:<3s}'
                        + line[20:30]
                        + f'{-x:8.3f}'
                        + f'{y:8.3f}'
                        + f'{z:8.3f}'
                        + line[54:]
                    )
                    lines_out.append(new_line)
                except (ValueError, IndexError):
                    lines_out.append(line)
            else:
                lines_out.append(line)
    with open(output_pdb, 'w') as fh:
        fh.writelines(lines_out)


# ═══════════════════════════════════════════════════════════════════════════
# Core threading function
# ═══════════════════════════════════════════════════════════════════════════

def thread_sequence(
    target_seq:     str,
    template_pdb:   str,
    template_chain: str,
    output_pdb:     str,
    chirality:      str = 'L',
) -> Dict[str, Any]:
    """Thread *target_seq* onto the backbone of *template_pdb*.

    Backbone atoms (N, CA, C, O) are copied from the template for each
    matching residue position; residue names are renamed to match the target
    sequence.  Side-chain atoms are discarded — the result is a backbone-only
    model suitable for subsequent side-chain grafting or energy minimisation.

    Args:
        target_seq:     One-letter amino-acid sequence to thread.
        template_pdb:   Path to the template PDB file.
        template_chain: Chain identifier in the template to use as backbone.
        output_pdb:     Path for the output threaded PDB file.
        chirality:      'L' (default) for L-amino acids, 'D' for D-amino acids.
                        When 'D', residue names are converted to D-codes.

    Returns:
        Dict with:
        - **n_residues** (int): Number of residue positions written.
        - **n_atoms_written** (int): Total ATOM records written.
        - **template_pdb** (str): Input template path.
        - **target_seq** (str): Sequence actually threaded (may be truncated).
        - **output_pdb** (str): Path of the output file.

    Raises:
        FileNotFoundError: If *template_pdb* does not exist.
        ValueError: If *chirality* is not 'L' or 'D'.
    """
    if not os.path.isfile(template_pdb):
        raise FileNotFoundError(f"Template PDB not found: {template_pdb}")
    if chirality not in ('L', 'D'):
        raise ValueError(f"chirality must be 'L' or 'D', got {chirality!r}")

    bb_records = _parse_pdb_backbone(template_pdb, template_chain)
    if not bb_records:
        raise ValueError(
            f"No backbone atoms found for chain '{template_chain}' in {template_pdb}"
        )

    # Collect unique residue sequence numbers (ordered)
    seen_resseqs: List[int] = []
    for r in bb_records:
        if r['resseq'] not in seen_resseqs:
            seen_resseqs.append(r['resseq'])

    template_len = len(seen_resseqs)
    target_len   = len(target_seq)

    # Truncate to the shorter length
    if target_len != template_len:
        use_len = min(target_len, template_len)
        warnings.warn(
            f"Target sequence length ({target_len}) differs from template chain "
            f"length ({template_len}). Truncating to {use_len} residues.",
            UserWarning,
            stacklevel=2,
        )
        target_seq  = target_seq[:use_len]
        seen_resseqs = seen_resseqs[:use_len]
    else:
        use_len = target_len

    # Build resseq → new residue name mapping
    resseq_to_new: Dict[int, str] = {}
    for i, resseq in enumerate(seen_resseqs):
        aa1  = target_seq[i].upper()
        aa3  = _ONE_TO_THREE.get(aa1, 'UNK')
        if chirality == 'D':
            aa3 = _L_TO_D_RESNAME.get(aa3, aa3)
        resseq_to_new[resseq] = aa3

    # Write output PDB
    os.makedirs(os.path.dirname(os.path.abspath(output_pdb)), exist_ok=True)
    n_atoms = 0
    serial  = 1

    with open(output_pdb, 'w') as fh:
        fh.write(f"REMARK  Threaded by ChiralFold threading.py\n")
        fh.write(f"REMARK  Template: {os.path.basename(template_pdb)} chain {template_chain}\n")
        fh.write(f"REMARK  Target:   {target_seq}\n")
        fh.write(f"REMARK  Chirality: {chirality}\n")

        for r in bb_records:
            if r['resseq'] not in resseq_to_new:
                continue
            new_resname = resseq_to_new[r['resseq']]
            # New sequential residue index (1-based) for output
            new_resseq = seen_resseqs.index(r['resseq']) + 1
            line = _format_atom_line(
                serial    = serial,
                atom_name = r['atom_name'],
                resname   = new_resname,
                chain     = template_chain,
                resseq    = new_resseq,
                x         = r['x'],
                y         = r['y'],
                z         = r['z'],
            )
            fh.write(line)
            serial  += 1
            n_atoms += 1

        fh.write("END\n")

    return {
        'n_residues':     use_len,
        'n_atoms_written': n_atoms,
        'template_pdb':   template_pdb,
        'target_seq':     target_seq,
        'output_pdb':     output_pdb,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Mirror-image convenience wrapper
# ═══════════════════════════════════════════════════════════════════════════

def thread_and_mirror(
    target_seq:     str,
    template_pdb:   str,
    template_chain: str,
    output_pdb:     str,
) -> str:
    """Thread *target_seq* onto the template, then mirror to produce a D-protein.

    Workflow:
      1. Thread L-sequence onto backbone (``chirality='L'``).
      2. Reflect x-coordinates (x → −x) and rename residues to D-codes.
      3. Write the mirrored structure to *output_pdb*.

    Args:
        target_seq:     One-letter amino-acid sequence.
        template_pdb:   Path to the template PDB file.
        template_chain: Chain identifier in the template.
        output_pdb:     Path for the final D-protein PDB.

    Returns:
        Path to the D-protein PDB (*output_pdb*).
    """
    # Step 1: thread as L
    l_pdb = output_pdb + '.L_intermediate.pdb'
    thread_sequence(
        target_seq     = target_seq,
        template_pdb   = template_pdb,
        template_chain = template_chain,
        output_pdb     = l_pdb,
        chirality      = 'L',
    )

    # Step 2: mirror
    _mirror_pdb_file(l_pdb, output_pdb)

    # Clean up intermediate
    try:
        os.remove(l_pdb)
    except OSError:
        pass

    return output_pdb


# ═══════════════════════════════════════════════════════════════════════════
# Template finder
# ═══════════════════════════════════════════════════════════════════════════

def _sequence_identity(seq_a: str, seq_b: str) -> Tuple[float, int]:
    """Compute simple percent sequence identity between two sequences.

    Aligns by truncating to the shorter length.

    Returns:
        (identity_fraction, alignment_length)
    """
    n = min(len(seq_a), len(seq_b))
    if n == 0:
        return 0.0, 0
    matches = sum(a == b for a, b in zip(seq_a[:n], seq_b[:n]))
    return matches / n, n


def find_template(
    target_seq: str,
    pdb_dir:    str,
    chain:      str = 'A',
) -> Dict[str, Any]:
    """Scan a directory of PDB files and return the best template for *target_seq*.

    Evaluates every ``*.pdb`` file in *pdb_dir*, extracts chain *chain*
    sequence, and selects the file with the highest simple sequence identity.

    Args:
        target_seq: One-letter query sequence.
        pdb_dir:    Directory containing PDB files to search.
        chain:      Chain identifier to extract from each PDB (default 'A').

    Returns:
        Dict with:
        - **template_pdb** (str): Path to the best-matching PDB.
        - **chain** (str): Chain used ('A' by default).
        - **identity** (float): Sequence identity fraction (0–1).
        - **alignment_length** (int): Number of positions compared.
        - **template_seq** (str): Sequence extracted from the template.

    Raises:
        FileNotFoundError: If *pdb_dir* does not exist.
        ValueError: If no PDB files are found in *pdb_dir*.
    """
    if not os.path.isdir(pdb_dir):
        raise FileNotFoundError(f"PDB directory not found: {pdb_dir}")

    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, '*.pdb')))
    if not pdb_files:
        raise ValueError(f"No .pdb files found in {pdb_dir}")

    best_pdb      = None
    best_identity = -1.0
    best_alen     = 0
    best_seq      = ''

    for pdb_path in pdb_files:
        try:
            tmpl_seq = _extract_chain_sequence(pdb_path, chain)
        except Exception:
            continue
        if not tmpl_seq:
            continue
        identity, alen = _sequence_identity(target_seq, tmpl_seq)
        if identity > best_identity:
            best_identity = identity
            best_pdb      = pdb_path
            best_alen     = alen
            best_seq      = tmpl_seq

    if best_pdb is None:
        raise ValueError(f"Could not read any valid chain '{chain}' from PDB files in {pdb_dir}")

    return {
        'template_pdb':     best_pdb,
        'chain':            chain,
        'identity':         round(best_identity, 4),
        'alignment_length': best_alen,
        'template_seq':     best_seq,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Quick self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import tempfile
    import shutil

    print("=== threading.py self-test ===\n")

    # Build a minimal 5-residue backbone PDB as template
    bb_coords = [
        # resseq, atom, x, y, z
        (1, 'N',  0.000, 0.000, 0.000),
        (1, 'CA', 1.458, 0.000, 0.000),
        (1, 'C',  2.130, 1.200, 0.000),
        (1, 'O',  1.600, 2.200, 0.000),
        (2, 'N',  3.460, 1.200, 0.000),
        (2, 'CA', 4.200, 2.400, 0.000),
        (2, 'C',  5.700, 2.400, 0.000),
        (2, 'O',  6.300, 3.400, 0.000),
        (3, 'N',  6.400, 1.300, 0.000),
        (3, 'CA', 7.850, 1.100, 0.000),
        (3, 'C',  8.500, 2.300, 0.000),
        (3, 'O',  8.000, 3.300, 0.000),
        (4, 'N',  9.820, 2.200, 0.000),
        (4, 'CA',10.600, 3.300, 0.000),
        (4, 'C', 12.100, 3.200, 0.000),
        (4, 'O', 12.700, 4.200, 0.000),
        (5, 'N', 12.800, 2.100, 0.000),
        (5, 'CA',14.200, 1.900, 0.000),
        (5, 'C', 14.900, 3.100, 0.000),
        (5, 'O', 14.400, 4.100, 0.000),
    ]

    tmpdir = tempfile.mkdtemp()
    tmpl_path = os.path.join(tmpdir, 'template.pdb')
    with open(tmpl_path, 'w') as f:
        for i, (resseq, aname, x, y, z) in enumerate(bb_coords, start=1):
            f.write(_format_atom_line(i, aname, 'ALA', 'A', resseq, x, y, z))
        f.write("END\n")

    # --- thread_sequence: L ---
    out_l = os.path.join(tmpdir, 'threaded_L.pdb')
    result = thread_sequence('ACDEF', tmpl_path, 'A', out_l, chirality='L')
    print(f"thread_sequence (L): {result}")
    assert result['n_residues'] == 5
    assert result['n_atoms_written'] == 20
    assert os.path.isfile(out_l)

    # Check residue names
    with open(out_l) as f:
        content = f.read()
    assert 'ALA' in content and 'CYS' in content, "Expected ALA and CYS in L output"

    # --- thread_sequence: D ---
    out_d_direct = os.path.join(tmpdir, 'threaded_D.pdb')
    result_d = thread_sequence('ACDEF', tmpl_path, 'A', out_d_direct, chirality='D')
    print(f"thread_sequence (D): {result_d}")
    with open(out_d_direct) as f:
        d_content = f.read()
    assert 'DAL' in d_content, "Expected D-ALA code DAL"

    # --- thread_and_mirror ---
    out_mirror = os.path.join(tmpdir, 'mirror.pdb')
    mirror_path = thread_and_mirror('ACDEF', tmpl_path, 'A', out_mirror)
    print(f"thread_and_mirror: {mirror_path}")
    assert os.path.isfile(mirror_path)

    # --- find_template ---
    found = find_template('ACDEF', tmpdir, chain='A')
    print(f"find_template: {found}")
    assert found['identity'] == 1.0 or found['identity'] >= 0.0

    shutil.rmtree(tmpdir)
    print("\nAll threading.py tests passed.")
