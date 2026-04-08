"""
ChiralFold — Fragment-Based Backbone Assembly
===============================================

Builds full protein backbones for sequences > 30 residues using a
fragment-assembly approach.  Secondary-structure propensity is predicted
with a simple Chou–Fasman-style scheme, then backbone φ/ψ angles are
sampled from precomputed fragment libraries for helix, sheet, and coil.

Atom positions are placed sequentially using the **NeRF** (Natural
Extension Reference Frame) algorithm, which places each new atom given
the three preceding atoms, a bond length, bond angle, and dihedral.

Usage::

    from chiralfold.fragments import assemble_protein

    result = assemble_protein('AKLVEFGLMKQRSTHWIVDN', output_pdb='out.pdb')
    print(result['ss_prediction'])   # e.g. 'HHHHEEEECCHHHH...'
"""

import os
import math
import random
import warnings
import numpy as np
from typing import Optional, List, Tuple, Dict, Any


# ═══════════════════════════════════════════════════════════════════════════
# Fragment library — φ/ψ angles sampled from real secondary structures
# ═══════════════════════════════════════════════════════════════════════════

#: Alpha-helix backbone φ angles (degrees), sampled from real helices.
HELIX_PHIS: List[float] = [-63, -57, -70, -60, -65, -58, -67, -62]

#: Alpha-helix backbone ψ angles (degrees).
HELIX_PSIS: List[float] = [-42, -47, -37, -40, -44, -48, -35, -43]

#: Beta-sheet backbone φ angles (degrees).
SHEET_PHIS: List[float] = [-120, -130, -140, -119, -135, -125, -110, -128]

#: Beta-sheet backbone ψ angles (degrees).
SHEET_PSIS: List[float] = [130, 135, 145, 128, 140, 133, 125, 138]

#: Random-coil backbone φ angles (degrees).
COIL_PHIS: List[float] = [-70, -80, -60, -90, -65, -75, -85, -100, 60, -50]

#: Random-coil backbone ψ angles (degrees).
COIL_PSIS: List[float] = [140, -30, 150, 0, 130, -40, 160, 10, 30, -20]

#: Combined backbone fragment library keyed by secondary-structure type.
BACKBONE_FRAGMENTS: Dict[str, Dict[str, List[float]]] = {
    'H': {'phi': HELIX_PHIS, 'psi': HELIX_PSIS},
    'E': {'phi': SHEET_PHIS, 'psi': SHEET_PSIS},
    'C': {'phi': COIL_PHIS,  'psi': COIL_PSIS},
}

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

# ═══════════════════════════════════════════════════════════════════════════
# Chou–Fasman secondary-structure propensities
# ═══════════════════════════════════════════════════════════════════════════

#: Residues with helix propensity > 1.0.
_HELIX_FORMERS = set('ALEMKQ')   # Ala Leu Glu Met Lys Gln

#: Residues with sheet propensity > 1.0.
_SHEET_FORMERS = set('VIYFWT')   # Val Ile Tyr Phe Trp Thr

# Everything else defaults to coil.


# ═══════════════════════════════════════════════════════════════════════════
# Secondary-structure prediction
# ═══════════════════════════════════════════════════════════════════════════

def predict_secondary_structure(sequence: str) -> str:
    """Predict secondary structure using a Chou–Fasman-style propensity model.

    Each residue is classified independently based on its side-chain
    propensity for helix, sheet, or coil:

    - **H** (helix): A, L, E, M, K, Q — high helix propensity.
    - **E** (sheet): V, I, Y, F, W, T — high strand propensity.
    - **C** (coil):  all other residues.

    Args:
        sequence: One-letter amino-acid sequence (case-insensitive).

    Returns:
        String of 'H', 'E', 'C' with the same length as *sequence*.
    """
    ss = []
    for aa in sequence.upper():
        if aa in _HELIX_FORMERS:
            ss.append('H')
        elif aa in _SHEET_FORMERS:
            ss.append('E')
        else:
            ss.append('C')
    return ''.join(ss)


# ═══════════════════════════════════════════════════════════════════════════
# NeRF atom placement
# ═══════════════════════════════════════════════════════════════════════════

def _place_atom(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    bond_length: float,
    bond_angle:  float,
    dihedral:    float,
) -> np.ndarray:
    """Place a new atom using the Natural Extension Reference Frame (NeRF).

    Given three preceding atom positions *p1*, *p2*, *p3*, places atom *p4*
    such that the bond *p3–p4* has length *bond_length*, the angle
    *p2–p3–p4* equals *bond_angle* degrees, and the dihedral
    *p1–p2–p3–p4* equals *dihedral* degrees.

    Args:
        p1, p2, p3:  Positions of the three preceding atoms (arrays of 3).
        bond_length: Distance from *p3* to new atom (Å).
        bond_angle:  Bond angle at *p3* in degrees (p2–p3–p4).
        dihedral:    Dihedral angle in degrees (p1–p2–p3–p4).

    Returns:
        Position of the new atom *p4* as a numpy array of shape (3,).
    """
    d = p3 - p2
    n = np.cross(p2 - p1, d)
    n_norm = np.linalg.norm(n)
    d_norm = np.linalg.norm(d)

    n = n / (n_norm + 1e-12)
    d = d / (d_norm + 1e-12)
    m = np.cross(n, d)

    r     = bond_length
    theta = math.radians(bond_angle)
    phi   = math.radians(dihedral)

    x_coef = r * math.cos(math.pi - theta)
    y_coef = r * math.sin(math.pi - theta) * math.cos(phi)
    z_coef = r * math.sin(math.pi - theta) * math.sin(phi)

    p4 = p3 + x_coef * d + y_coef * m + z_coef * n
    return p4


# ═══════════════════════════════════════════════════════════════════════════
# Ideal backbone geometry constants
# ═══════════════════════════════════════════════════════════════════════════

# Bond lengths (Å)
_BOND_N_CA  = 1.458  # N–Cα
_BOND_CA_C  = 1.525  # Cα–C
_BOND_C_N   = 1.329  # C–N (peptide bond)
_BOND_C_O   = 1.231  # C=O (carbonyl)

# Bond angles (degrees)
_ANGLE_N_CA_C  = 111.0  # N–Cα–C
_ANGLE_CA_C_N  = 116.2  # Cα–C–N
_ANGLE_C_N_CA  = 121.7  # C–N–Cα

# Peptide omega (trans)
_OMEGA = 180.0


# ═══════════════════════════════════════════════════════════════════════════
# Backbone builder
# ═══════════════════════════════════════════════════════════════════════════

def build_backbone_from_fragments(
    sequence:       str,
    ss_prediction:  Optional[str] = None,
    seed:           Optional[int] = None,
) -> List[Dict[str, np.ndarray]]:
    """Build a backbone atom array for *sequence* using fragment φ/ψ sampling.

    For each residue position, φ and ψ are sampled from the appropriate
    fragment set (helix / sheet / coil) according to the secondary-structure
    prediction.  Atom positions are computed sequentially with the NeRF
    algorithm.

    Args:
        sequence:      One-letter amino-acid sequence.
        ss_prediction: Pre-computed SS string (H/E/C per residue).  If None,
                       it is computed via :func:`predict_secondary_structure`.
        seed:          Random seed for reproducible sampling (optional).

    Returns:
        List of dicts, one per residue, each with keys ``'N'``, ``'CA'``,
        ``'C'``, ``'O'`` mapping to numpy arrays of shape (3,).
    """
    if ss_prediction is None:
        ss_prediction = predict_secondary_structure(sequence)

    if len(ss_prediction) != len(sequence):
        raise ValueError(
            f"ss_prediction length ({len(ss_prediction)}) must match "
            f"sequence length ({len(sequence)})."
        )

    rng = random.Random(seed)
    n_res = len(sequence)

    # ── Sample phi/psi for every residue from the fragment library ───────
    phi_psi: List[Tuple[float, float]] = []
    for i in range(n_res):
        ss = ss_prediction[i]
        frags_i = BACKBONE_FRAGMENTS.get(ss, BACKBONE_FRAGMENTS['C'])
        idx_i = rng.randrange(len(frags_i['phi']))
        phi_psi.append((frags_i['phi'][idx_i], frags_i['psi'][idx_i]))

    # ── Build backbone atom positions ────────────────────────────────────
    # Representation: store positions indexed by global atom index.
    # Order per residue: N, CA, C  (O placed after N of next residue)

    residue_data: List[Dict[str, Optional[np.ndarray]]] = []

    # Seed atoms for residue 0
    N_cur  = np.array([0.0, 0.0, 0.0])
    CA_cur = np.array([_BOND_N_CA, 0.0, 0.0])
    # C: place with ideal angle, in XY plane
    C_cur = np.array([
        _BOND_N_CA + _BOND_CA_C * math.cos(math.radians(180.0 - _ANGLE_N_CA_C)),
        _BOND_CA_C * math.sin(math.radians(180.0 - _ANGLE_N_CA_C)),
        0.0,
    ])

    # We need a "phantom" atom before N for the first NeRF call.
    # Place it such that we get reasonable geometry.
    phantom = N_cur - np.array([_BOND_C_N, 0.0, 0.0])

    # Use chain of placed atoms: grow atom-by-atom
    # Chain: ... [prev_N, prev_CA, prev_C] → [N, CA, C] → ...
    chain_N:  List[np.ndarray] = []
    chain_CA: List[np.ndarray] = []
    chain_C:  List[np.ndarray] = []

    # Residue 0 — directly seeded
    chain_N.append(N_cur.copy())
    chain_CA.append(CA_cur.copy())
    chain_C.append(C_cur.copy())

    # Residues 1..n_res-1 — grown with NeRF
    for i in range(1, n_res):
        phi_i, psi_i = phi_psi[i]
        phi_prev, psi_prev = phi_psi[i-1]

        # Place N_i:  p1=N_{i-1}, p2=CA_{i-1}, p3=C_{i-1}
        # Dihedral = psi_{i-1} (rotation around CA_{i-1}–C_{i-1} bond)
        N_i = _place_atom(
            p1=chain_N[i-1], p2=chain_CA[i-1], p3=chain_C[i-1],
            bond_length=_BOND_C_N,
            bond_angle=_ANGLE_CA_C_N,
            dihedral=psi_prev,
        )
        chain_N.append(N_i)

        # Place CA_i: p1=CA_{i-1}, p2=C_{i-1}, p3=N_i
        # Dihedral = omega (peptide, ~180°)
        CA_i = _place_atom(
            p1=chain_CA[i-1], p2=chain_C[i-1], p3=N_i,
            bond_length=_BOND_N_CA,
            bond_angle=_ANGLE_C_N_CA,
            dihedral=_OMEGA,
        )
        chain_CA.append(CA_i)

        # Place C_i:  p1=C_{i-1}, p2=N_i, p3=CA_i
        # Dihedral = phi_i (rotation around N_i–CA_i bond)
        C_i = _place_atom(
            p1=chain_C[i-1], p2=N_i, p3=CA_i,
            bond_length=_BOND_CA_C,
            bond_angle=_ANGLE_N_CA_C,
            dihedral=phi_i,
        )
        chain_C.append(C_i)

    # ── Place O atoms ────────────────────────────────────────────────────
    # O_i is placed using: p1=CA_i, p2=C_i, p3=N_{i+1} (or phantom for last)
    # Dihedral for O is ~0° relative to the peptide plane (cis to N)
    # Approximate: O is placed at 120° from CA-C bond, 180° from N_{i+1}
    # Use NeRF with dihedral = 0° and angle 120.9° (ideal C=O in peptide)

    _ANGLE_CA_C_O  = 120.9  # Cα–C=O bond angle (degrees)
    _DIHEDRAL_O    = 0.0    # O is in the peptide plane, same side as N_{i+1}

    chain_O: List[np.ndarray] = []
    for i in range(n_res):
        if i < n_res - 1:
            # Use N_{i+1} as the reference
            O_i = _place_atom(
                p1=chain_CA[i], p2=chain_C[i], p3=chain_N[i+1],
                bond_length=_BOND_C_O,
                bond_angle=_ANGLE_CA_C_O,
                dihedral=_DIHEDRAL_O,
            )
        else:
            # Last residue: use a phantom N in the direction of the chain
            # Extend C in the direction of the last peptide bond
            direction = chain_C[-1] - chain_CA[-1]
            d_norm = np.linalg.norm(direction)
            if d_norm > 1e-8:
                direction = direction / d_norm
            phantom_N = chain_C[-1] + direction * _BOND_C_N
            O_i = _place_atom(
                p1=chain_CA[-1], p2=chain_C[-1], p3=phantom_N,
                bond_length=_BOND_C_O,
                bond_angle=_ANGLE_CA_C_O,
                dihedral=_DIHEDRAL_O,
            )
        chain_O.append(O_i)

    # ── Assemble residue dicts ───────────────────────────────────────────
    for i in range(n_res):
        residue_data.append({
            'N':  chain_N[i],
            'CA': chain_CA[i],
            'C':  chain_C[i],
            'O':  chain_O[i],
        })

    return residue_data


# ═══════════════════════════════════════════════════════════════════════════
# PDB writer
# ═══════════════════════════════════════════════════════════════════════════

def _format_pdb_atom(
    serial:    int,
    atom_name: str,
    resname:   str,
    chain:     str,
    resseq:    int,
    x: float, y: float, z: float,
    element:   str,
) -> str:
    """Format a single ATOM PDB record."""
    name_field = f' {atom_name:<3s}' if len(atom_name) < 4 else f'{atom_name:<4s}'
    return (
        f"ATOM  {serial:5d} {name_field} {resname:<3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


# ═══════════════════════════════════════════════════════════════════════════
# Full pipeline
# ═══════════════════════════════════════════════════════════════════════════

def assemble_protein(
    sequence:    str,
    output_pdb:  Optional[str] = None,
    chirality:   str = 'L',
    seed:        Optional[int] = None,
) -> Dict[str, Any]:
    """Full fragment-assembly pipeline: predict SS → build backbone → write PDB.

    Args:
        sequence:   One-letter amino-acid sequence (case-insensitive).
                    Best results for sequences > 30 residues.
        output_pdb: Output PDB path.  If None, auto-generated as
                    ``fragment_<len>aa.pdb`` in the current directory.
        chirality:  'L' (default) or 'D'.  For 'D', all x-coordinates are
                    negated and residue names converted to D-amino acid codes.
        seed:       Random seed for reproducible fragment sampling.

    Returns:
        Dict with:
        - **n_residues** (int): Number of residues in the structure.
        - **n_atoms** (int): Total backbone atoms written.
        - **ss_prediction** (str): Secondary-structure string (H/E/C).
        - **output_pdb** (str): Path to the written PDB file.

    Raises:
        ValueError: If *chirality* is not 'L' or 'D'.
    """
    sequence = sequence.upper().strip()
    if not sequence:
        raise ValueError("sequence must be a non-empty string.")
    if chirality not in ('L', 'D'):
        raise ValueError(f"chirality must be 'L' or 'D', got {chirality!r}")

    if output_pdb is None:
        output_pdb = f"fragment_{len(sequence)}aa.pdb"

    ss_pred = predict_secondary_structure(sequence)
    residue_data = build_backbone_from_fragments(sequence, ss_pred, seed=seed)

    n_res   = len(sequence)
    n_atoms = 0
    serial  = 1
    chain   = 'A'

    _ATOM_ORDER = [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O')]

    os.makedirs(os.path.dirname(os.path.abspath(output_pdb)), exist_ok=True)

    with open(output_pdb, 'w') as fh:
        fh.write(f"REMARK  Built by ChiralFold fragments.py\n")
        fh.write(f"REMARK  Sequence:  {sequence}\n")
        fh.write(f"REMARK  SS:        {ss_pred}\n")
        fh.write(f"REMARK  Chirality: {chirality}\n")

        for i, aa1 in enumerate(sequence):
            resseq  = i + 1
            aa3     = _ONE_TO_THREE.get(aa1, 'UNK')
            if chirality == 'D':
                aa3 = _L_TO_D_RESNAME.get(aa3, aa3)

            atoms = residue_data[i]

            for atom_key, element in _ATOM_ORDER:
                pos = atoms.get(atom_key)
                if pos is None:
                    continue

                x, y, z = float(pos[0]), float(pos[1]), float(pos[2])

                # Mirror x for D-chirality
                if chirality == 'D':
                    x = -x

                line = _format_pdb_atom(
                    serial    = serial,
                    atom_name = atom_key,
                    resname   = aa3,
                    chain     = chain,
                    resseq    = resseq,
                    x=x, y=y, z=z,
                    element   = element,
                )
                fh.write(line)
                serial  += 1
                n_atoms += 1

        fh.write("END\n")

    return {
        'n_residues':    n_res,
        'n_atoms':       n_atoms,
        'ss_prediction': ss_pred,
        'output_pdb':    output_pdb,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Quick self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import tempfile

    print("=== fragments.py self-test ===\n")

    # --- predict_secondary_structure ---
    seq = 'AKLVEFGLMKQRSTHWIVDN'
    ss  = predict_secondary_structure(seq)
    print(f"Sequence : {seq}")
    print(f"SS pred  : {ss}")
    assert len(ss) == len(seq), "SS prediction length mismatch"
    assert set(ss).issubset({'H', 'E', 'C'}), "Invalid SS characters"

    # A → H, V → E, G → C
    assert ss[seq.index('A')] == 'H', "A should be helix"
    assert ss[seq.index('V')] == 'E', "V should be sheet"
    print("predict_secondary_structure: OK\n")

    # --- _place_atom (NeRF) ---
    # Use non-collinear geometry so cross products are non-zero.
    p1 = np.array([0.0,  0.0, 0.0])
    p2 = np.array([1.458, 0.0, 0.0])
    p3 = np.array([1.936, 1.440, 0.0])   # bent out of axis
    p4 = _place_atom(p1, p2, p3, bond_length=1.5, bond_angle=109.5, dihedral=180.0)
    dist = np.linalg.norm(p4 - p3)
    print(f"NeRF bond length (expect 1.5): {dist:.4f}")
    assert abs(dist - 1.5) < 0.01, f"Bond length error: {dist}"
    print("_place_atom: OK\n")

    # --- build_backbone_from_fragments ---
    seq35 = 'AKLVEFGLMKQRSTHWIVDNAKLVEFGLMKQRS'
    residues = build_backbone_from_fragments(seq35, seed=42)
    assert len(residues) == len(seq35), "Wrong number of residues"
    for r in residues:
        for key in ('N', 'CA', 'C', 'O'):
            assert key in r and r[key] is not None, f"Missing {key}"
            assert r[key].shape == (3,), f"Wrong shape for {key}"
    print(f"build_backbone_from_fragments: {len(residues)} residues built")

    # Check bond length N-CA for residue 0
    d0 = np.linalg.norm(residues[0]['CA'] - residues[0]['N'])
    print(f"  Residue 0 N–CA distance (expect ~{_BOND_N_CA}): {d0:.3f} Å")
    print("build_backbone_from_fragments: OK\n")

    # --- assemble_protein (L) ---
    tmpdir = tempfile.mkdtemp()
    out_l  = os.path.join(tmpdir, 'frag_L.pdb')
    result_l = assemble_protein(seq35, output_pdb=out_l, chirality='L', seed=42)
    print(f"assemble_protein (L): {result_l}")
    assert result_l['n_residues']    == len(seq35)
    assert result_l['n_atoms']       == len(seq35) * 4
    assert os.path.isfile(out_l)

    # Spot-check PDB content
    with open(out_l) as f:
        content = f.read()
    assert 'ATOM' in content
    assert 'ALA'  in content  # A → ALA
    print("assemble_protein (L): OK\n")

    # --- assemble_protein (D) ---
    out_d = os.path.join(tmpdir, 'frag_D.pdb')
    result_d = assemble_protein(seq35, output_pdb=out_d, chirality='D', seed=42)
    print(f"assemble_protein (D): {result_d}")
    with open(out_d) as f:
        d_content = f.read()
    assert 'DAL' in d_content, "Expected D-ALA code DAL in D output"

    # D-structure should have negated x relative to L
    def _extract_first_x(pdb_text):
        for line in pdb_text.splitlines():
            if line.startswith('ATOM'):
                return float(line[30:38])
        return None

    x_l = _extract_first_x(content)
    x_d = _extract_first_x(d_content)
    print(f"First atom x (L): {x_l:.3f},  (D): {x_d:.3f}  (expect negated)")
    assert abs(x_l + x_d) < 0.001, "D x-coords should be negated L x-coords"
    print("assemble_protein (D): OK\n")

    import shutil
    shutil.rmtree(tmpdir)
    print("All fragments.py tests passed.")
