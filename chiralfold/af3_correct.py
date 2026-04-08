"""
ChiralFold — AlphaFold 3 Chirality Correction Pipeline
========================================================

Detects and corrects chirality violations in AlphaFold 3 outputs.
AF3's diffusion architecture produces ~51% chirality violations on
D-peptides (Childs et al. 2025). This module fixes them.

Correction method:
  For each Cα with wrong chirality, reflect its position across the
  plane defined by its three backbone neighbors (N, C, Cβ). This
  inverts the stereocenter while minimally perturbing the backbone.

Reference:
  Childs, Zhou & Donald (2025). "Has AlphaFold 3 Solved the Protein
  Folding Problem for D-Peptides?" bioRxiv 2025.03.14.643307
"""

from __future__ import annotations

import math
import os
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# D-amino acid residue names (HETATM records in PDB)
D_RESNAMES: Set[str] = {
    "DAL", "DAR", "DSG", "DAS", "DCY", "DGL", "DGN",
    "DHI", "DIL", "DLE", "DLY", "MED", "DPN", "DPR",
    "DSN", "DTH", "DTR", "DTY", "DVA",
}

# Standard L-amino acid residue names (ATOM records)
L_RESNAMES: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "MSE", "SEP", "TPO", "PTR", "HSD", "HSE", "HSP",
    "HIE", "HID", "HIP", "CYX", "CYM",
}

GLYCINE_RESNAMES: Set[str] = {"GLY"}

# Signed-volume convention (matches auditor.py):
#   _signed_volume(ca, n, c, cb) = (n-ca) . ((c-ca) x (cb-ca))
#   positive -> L-chirality (S-configuration, standard amino acids)
#   negative -> D-chirality (R-configuration, mirror-image)
_SV_THRESHOLD = 0.0


# ─────────────────────────────────────────────────────────────────────────────
# PDB Parsing (lightweight, mirrors auditor.py style)
# ─────────────────────────────────────────────────────────────────────────────

class _Atom:
    """Lightweight representation of a parsed PDB atom record."""

    __slots__ = [
        "record", "serial", "name", "altloc", "resname",
        "chain", "resseq", "icode",
        "x", "y", "z", "element",
        "line_idx",  # original line index for rewriting
    ]

    def __init__(
        self,
        record: str, serial: int, name: str, altloc: str,
        resname: str, chain: str, resseq: int, icode: str,
        x: float, y: float, z: float, element: str,
        line_idx: int = -1,
    ):
        self.record   = record
        self.serial   = serial
        self.name     = name
        self.altloc   = altloc
        self.resname  = resname
        self.chain    = chain
        self.resseq   = resseq
        self.icode    = icode
        self.x        = x
        self.y        = y
        self.z        = z
        self.element  = element
        self.line_idx = line_idx

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    @xyz.setter
    def xyz(self, coords: np.ndarray) -> None:
        self.x, self.y, self.z = float(coords[0]), float(coords[1]), float(coords[2])


def _parse_pdb_full(pdb_path: str):
    """
    Parse ATOM and HETATM records from *pdb_path*.

    Returns:
        (lines, atoms) where lines is the raw list of all PDB lines,
        and atoms is a list of _Atom objects with line_idx set.
    """
    atoms: List[_Atom] = []
    lines: List[str] = []
    seen_altloc: Dict[Tuple, str] = {}

    with open(pdb_path) as fh:
        all_lines = fh.readlines()

    lines = list(all_lines)

    for i, line in enumerate(all_lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            continue

        resname = line[17:20].strip()
        if resname in ("HOH", "WAT", "DOD"):
            continue

        altloc = line[16]
        if altloc not in (" ", "A", "1", ""):
            continue

        try:
            record  = line[0:6].strip()
            serial  = int(line[6:11])
            name    = line[12:16].strip()
            chain   = line[21] if len(line) > 21 else " "
            resseq  = int(line[22:26])
            icode   = line[26] if len(line) > 26 else " "
            x       = float(line[30:38])
            y       = float(line[38:46])
            z       = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else ""
        except (ValueError, IndexError):
            continue

        key = (chain, resseq, icode, name)
        if key in seen_altloc:
            continue
        seen_altloc[key] = altloc

        atoms.append(_Atom(
            record=record, serial=serial, name=name, altloc=altloc,
            resname=resname, chain=chain, resseq=resseq, icode=icode,
            x=x, y=y, z=z, element=element, line_idx=i,
        ))

    return lines, atoms


def _group_by_residue(atoms: List[_Atom]) -> Dict[Tuple, List[_Atom]]:
    """Group atoms into residues keyed by (chain, resseq, icode, resname)."""
    groups: Dict[Tuple, List[_Atom]] = defaultdict(list)
    for a in atoms:
        key = (a.chain, a.resseq, a.icode, a.resname)
        groups[key].append(a)
    return dict(sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])))


# ─────────────────────────────────────────────────────────────────────────────
# Chirality geometry
# ─────────────────────────────────────────────────────────────────────────────

def _signed_volume(ca: np.ndarray, n: np.ndarray,
                   c: np.ndarray, cb: np.ndarray) -> float:
    """
    Compute the signed volume of the tetrahedron (Cα; N, C, Cβ).

    Sign convention (consistent with auditor.py):
      positive -> L-chirality (S-configuration)
      negative -> D-chirality (R-configuration)

    Uses the scalar triple product: (N-Cα) · ((C-Cα) × (Cβ-Cα))
    """
    v1 = n  - ca
    v2 = c  - ca
    v3 = cb - ca
    return float(np.dot(v1, np.cross(v2, v3)))


def _estimate_cb(ca: np.ndarray, n: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Estimate an ideal Cβ position from backbone atoms (for Gly or missing Cβ).

    Places a pseudo-Cβ 1.52 Å from Cα, perpendicular to the N-CA-C plane,
    on the same side as it would be for an L-amino acid.
    """
    v1 = n - ca
    v2 = c - ca
    v1 = v1 / (np.linalg.norm(v1) + 1e-10)
    v2 = v2 / (np.linalg.norm(v2) + 1e-10)
    # Bisector direction (points roughly toward where Cβ would be)
    bisector = v1 + v2
    norm_bisector = np.linalg.norm(bisector)
    if norm_bisector < 1e-10:
        bisector = np.cross(v1, v2)
    else:
        bisector = bisector / norm_bisector

    # Normal to backbone plane
    normal = np.cross(v1, v2)
    norm_n = np.linalg.norm(normal)
    if norm_n < 1e-10:
        normal = np.array([0.0, 0.0, 1.0])
    else:
        normal = normal / norm_n

    # Cβ direction: combine bisector and normal
    cb_dir = -bisector + normal
    cb_dir = cb_dir / (np.linalg.norm(cb_dir) + 1e-10)
    return ca + 1.52 * cb_dir


def _reflect_across_plane(point: np.ndarray,
                           p1: np.ndarray,
                           p2: np.ndarray,
                           p3: np.ndarray) -> np.ndarray:
    """
    Reflect *point* across the plane defined by three points p1, p2, p3.

    Algorithm:
      1. Compute plane normal from cross product of (p2-p1) and (p3-p1).
      2. Project the point onto the plane.
      3. Reflect: new_pos = point - 2 * d * normal
         where d is the signed distance from point to the plane.
    """
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    norm_len = np.linalg.norm(normal)
    if norm_len < 1e-10:
        return point.copy()
    normal = normal / norm_len
    d = np.dot(point - p1, normal)
    return point - 2.0 * d * normal


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def detect_chirality_violations(pdb_path: str) -> Dict:
    """
    Parse a PDB file and detect chirality violations at Cα stereocenters.

    For ATOM records (standard L-amino acids): expects L-chirality
    (negative signed volume).
    For HETATM records with D-residue names (DAL, DTR, etc.): expects
    D-chirality (positive signed volume).
    Glycine and proline are skipped.

    Args:
        pdb_path: Path to the PDB file to analyse.

    Returns:
        dict with keys:
          - n_residues (int): Total residues parsed.
          - n_checked (int): Residues where chirality was evaluated.
          - n_violations (int): Residues with incorrect chirality.
          - pct_violations (float): Percentage of checked residues violated.
          - violations (list[dict]): Per-violation details with resnum,
            resname, chain, expected_chirality, found_chirality, signed_volume.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    _, atoms = _parse_pdb_full(pdb_path)
    residues = _group_by_residue(atoms)

    n_residues  = len(residues)
    n_checked   = 0
    violations  = []

    for (chain, resseq, icode, resname), res_atoms in residues.items():
        # Skip glycine (achiral) and proline (ring complicates geometry)
        if resname in GLYCINE_RESNAMES or resname in {"PRO", "DPR"}:
            continue

        # Determine expected chirality from residue name and record type
        if resname in D_RESNAMES:
            expected = "D"
        elif resname in L_RESNAMES:
            expected = "L"
        else:
            # Unknown residue — try to infer from HETATM vs ATOM record
            records = {a.record for a in res_atoms}
            if "HETATM" in records:
                expected = "D"
            else:
                expected = "L"

        # Find backbone atoms
        atom_by_name = {a.name: a for a in res_atoms}
        ca_atom = atom_by_name.get("CA")
        n_atom  = atom_by_name.get("N")
        c_atom  = atom_by_name.get("C")
        cb_atom = atom_by_name.get("CB")

        if ca_atom is None or n_atom is None or c_atom is None:
            continue  # Cannot evaluate without backbone atoms

        ca_pos = ca_atom.xyz
        n_pos  = n_atom.xyz
        c_pos  = c_atom.xyz

        if cb_atom is not None:
            cb_pos = cb_atom.xyz
        else:
            cb_pos = _estimate_cb(ca_pos, n_pos, c_pos)

        sv = _signed_volume(ca_pos, n_pos, c_pos, cb_pos)
        n_checked += 1

        # Determine found chirality
        # positive SV = L (matches auditor.py convention)
        if sv > _SV_THRESHOLD:
            found = "L"
        else:
            found = "D"

        if found != expected:
            violations.append({
                "resnum":             resseq,
                "resname":            resname,
                "chain":              chain,
                "expected_chirality": expected,
                "found_chirality":    found,
                "signed_volume":      sv,
            })

    n_violations  = len(violations)
    pct_violations = (n_violations / n_checked * 100.0) if n_checked > 0 else 0.0

    return {
        "n_residues":    n_residues,
        "n_checked":     n_checked,
        "n_violations":  n_violations,
        "pct_violations": pct_violations,
        "violations":    violations,
    }


def _format_coord(v: float) -> str:
    """Format a coordinate for a PDB ATOM line (8.3f, 8 chars wide)."""
    return f"{v:8.3f}"


def _rewrite_atom_coords(line: str, x: float, y: float, z: float) -> str:
    """
    Rewrite the X, Y, Z fields of a PDB ATOM/HETATM line.

    PDB format: columns 31-38 = X, 39-46 = Y, 47-54 = Z (1-indexed).
    In 0-indexed Python: [30:38], [38:46], [46:54].
    """
    if len(line) < 54:
        return line
    new_line = (
        line[:30]
        + f"{x:8.3f}"
        + f"{y:8.3f}"
        + f"{z:8.3f}"
        + line[54:]
    )
    return new_line


def correct_chirality(pdb_path: str, output_path: str) -> Dict:
    """
    Correct chirality violations in a PDB file by reflecting offending Cα atoms.

    For each Cα with incorrect chirality:
      1. Identify its N, C, Cβ neighbors.
      2. Reflect Cα across the plane defined by N, C, Cβ — this inverts
         the stereocenter while minimally displacing the backbone.
      3. If an H atom is present on Cα (HA), also reflect it symmetrically.
      4. Write all corrected coordinates to output_path.

    Args:
        pdb_path:    Path to the input PDB file (AF3 output).
        output_path: Path to write the corrected PDB file.

    Returns:
        dict with keys:
          - n_corrected (int): Number of stereocenters corrected.
          - violations_before (int): Violations detected before correction.
          - violations_after (int): Violations remaining after correction.
          - output_path (str): Path to the output file.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    # First detect violations
    before_report = detect_chirality_violations(pdb_path)
    violations = before_report["violations"]

    if not violations:
        # No corrections needed — copy file as-is
        lines, _ = _parse_pdb_full(pdb_path)
        with open(output_path, "w") as fh:
            fh.writelines(lines)
        return {
            "n_corrected":       0,
            "violations_before": 0,
            "violations_after":  0,
            "output_path":       output_path,
        }

    # Parse full atom list with line indices for rewriting
    lines, atoms = _parse_pdb_full(pdb_path)
    residues = _group_by_residue(atoms)

    # Build a set of violation keys for quick lookup
    viol_keys = {
        (v["chain"], v["resnum"]): v for v in violations
    }

    # Work on a mutable copy of lines
    new_lines = list(lines)
    n_corrected = 0

    for (chain, resseq, icode, resname), res_atoms in residues.items():
        if (chain, resseq) not in viol_keys:
            continue

        atom_by_name = {a.name: a for a in res_atoms}
        ca_atom = atom_by_name.get("CA")
        n_atom  = atom_by_name.get("N")
        c_atom  = atom_by_name.get("C")
        cb_atom = atom_by_name.get("CB")
        ha_atom = atom_by_name.get("HA")

        if ca_atom is None or n_atom is None or c_atom is None:
            continue

        ca_pos = ca_atom.xyz
        n_pos  = n_atom.xyz
        c_pos  = c_atom.xyz

        if cb_atom is not None:
            cb_pos = cb_atom.xyz
        else:
            cb_pos = _estimate_cb(ca_pos, n_pos, c_pos)

        # Reflect Cα across the plane N, C, Cβ
        ca_new = _reflect_across_plane(ca_pos, n_pos, c_pos, cb_pos)

        # Update Cα in the line buffer
        if ca_atom.line_idx >= 0:
            new_lines[ca_atom.line_idx] = _rewrite_atom_coords(
                new_lines[ca_atom.line_idx], ca_new[0], ca_new[1], ca_new[2]
            )

        # Also adjust HA (the Cα-bound hydrogen) by reflecting it the same way
        if ha_atom is not None and ha_atom.line_idx >= 0:
            ha_new = _reflect_across_plane(ha_atom.xyz, n_pos, c_pos, cb_pos)
            new_lines[ha_atom.line_idx] = _rewrite_atom_coords(
                new_lines[ha_atom.line_idx], ha_new[0], ha_new[1], ha_new[2]
            )

        n_corrected += 1

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as fh:
        fh.writelines(new_lines)

    # Verify corrections
    after_report = detect_chirality_violations(output_path)

    return {
        "n_corrected":       n_corrected,
        "violations_before": before_report["n_violations"],
        "violations_after":  after_report["n_violations"],
        "output_path":       output_path,
    }


def correct_af3_output(pdb_path: str, output_path: Optional[str] = None) -> Dict:
    """
    Convenience wrapper: detect, correct, and verify AF3 chirality violations.

    Detects chirality violations in the AF3 PDB output, corrects them by
    reflecting offending Cα atoms, then re-detects to confirm the fix.

    Args:
        pdb_path:    Path to the AF3-generated PDB file.
        output_path: Where to write the corrected PDB. If None, derives a
                     path by inserting '_corrected' before the extension.

    Returns:
        dict with full report:
          - input_path (str)
          - output_path (str)
          - before (dict): detect_chirality_violations result before correction
          - correction (dict): correct_chirality result (n_corrected, etc.)
          - after (dict): detect_chirality_violations result after correction
          - success (bool): True if all violations were resolved
          - summary (str): Human-readable one-line summary
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if output_path is None:
        base, ext = os.path.splitext(pdb_path)
        output_path = f"{base}_corrected{ext}"

    # Step 1: Detect
    before = detect_chirality_violations(pdb_path)

    if before["n_violations"] == 0:
        # Nothing to do — just copy
        import shutil
        shutil.copy2(pdb_path, output_path)
        after = before
        correction = {
            "n_corrected":       0,
            "violations_before": 0,
            "violations_after":  0,
            "output_path":       output_path,
        }
    else:
        # Step 2: Correct
        correction = correct_chirality(pdb_path, output_path)
        # Step 3: Re-detect to verify
        after = detect_chirality_violations(output_path)

    success = after["n_violations"] == 0

    summary = (
        f"Corrected {correction['n_corrected']} of "
        f"{before['n_violations']} violations "
        f"({before['pct_violations']:.1f}% → "
        f"{after['pct_violations']:.1f}%)"
    )

    return {
        "input_path":  pdb_path,
        "output_path": output_path,
        "before":      before,
        "correction":  correction,
        "after":       after,
        "success":     success,
        "summary":     summary,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI / test block
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    import tempfile
    import textwrap

    # ── Minimal synthetic PDB for self-test ──────────────────────────────────
    # A single alanine residue with an intentionally wrong (D) Cα geometry
    # encoded as an L-residue (ATOM ALA). The signed volume will be positive,
    # indicating a violation.
    SYNTHETIC_PDB = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1      -1.200   2.095   0.000  1.00  0.00           O
        ATOM      5  CB  ALA A   1       0.000  -0.500  -1.200  1.00  0.00           C
        END
    """)
    # CB is at -z: signed volume is negative -> D geometry, but residue is ALA (L)
    # -> should detect 1 violation

    # A correct L-Ala: CB at +z gives positive signed volume (L-chirality)
    CORRECT_PDB = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1      -1.200   2.095   0.000  1.00  0.00           O
        ATOM      5  CB  ALA A   1       0.000  -0.500   1.200  1.00  0.00           C
        END
    """)

    print("=" * 60)
    print("ChiralFold — AF3 Correction Self-Test")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Test detect on correct L-ala
        correct_file = os.path.join(tmpdir, "correct.pdb")
        with open(correct_file, "w") as f:
            f.write(CORRECT_PDB)

        report = detect_chirality_violations(correct_file)
        print(f"\n[Correct L-Ala] n_violations = {report['n_violations']} "
              f"(expected 0)")

        # Test detect on intentionally wrong D-like Cα with L label
        wrong_file = os.path.join(tmpdir, "wrong.pdb")
        with open(wrong_file, "w") as f:
            f.write(SYNTHETIC_PDB)

        report2 = detect_chirality_violations(wrong_file)
        print(f"[Wrong L-Ala]   n_violations = {report2['n_violations']} "
              f"(expected 1 if geometry is inverted)")
        if report2["violations"]:
            v = report2["violations"][0]
            print(f"  Residue {v['resnum']} {v['resname']} chain {v['chain']}: "
                  f"expected {v['expected_chirality']}, "
                  f"found {v['found_chirality']}, "
                  f"SV = {v['signed_volume']:.4f}")

        # Test correction pipeline
        out_file = os.path.join(tmpdir, "corrected.pdb")
        full_report = correct_af3_output(wrong_file, out_file)
        print(f"\n[correct_af3_output] {full_report['summary']}")
        print(f"  Success: {full_report['success']}")

    # If a PDB path is passed as argument, run on it
    if len(sys.argv) > 1:
        pdb = sys.argv[1]
        out = sys.argv[2] if len(sys.argv) > 2 else None
        print(f"\nProcessing: {pdb}")
        result = correct_af3_output(pdb, out)
        print(result["summary"])
        for k, v in result["before"].items():
            if k != "violations":
                print(f"  Before {k}: {v}")
        for k, v in result["after"].items():
            if k != "violations":
                print(f"  After  {k}: {v}")
