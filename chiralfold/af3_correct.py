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

import os
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ._pdbio import read_lines_and_records

# Module-level: do NOT call warnings.filterwarnings globally — that suppresses
# warnings for all downstream user code.


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
        "chain", "resseq", "icode", "model",
        "x", "y", "z", "element",
        "line_idx",  # original line index for rewriting
    ]

    def __init__(
        self,
        record: str, serial: int, name: str, altloc: str,
        resname: str, chain: str, resseq: int, icode: str,
        x: float, y: float, z: float, element: str,
        line_idx: int = -1, model: int = 1,
    ):
        self.record   = record
        self.serial   = serial
        self.name     = name
        self.altloc   = altloc
        self.resname  = resname
        self.chain    = chain
        self.resseq   = resseq
        self.icode    = icode
        self.model    = model
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

    Reads **every** model: a multi-model file is corrected model by model, and
    each residue is keyed by its model number so stereocentres in models 2..n
    are neither skipped nor conflated with model 1. Alternate locations are
    resolved to one self-consistent conformation per residue by
    :mod:`chiralfold._pdbio`.

    Returns:
        (lines, atoms) where lines is the raw list of all PDB lines,
        and atoms is a list of _Atom objects with line_idx set.
    """
    lines, records = read_lines_and_records(pdb_path, first_model_only=False)
    atoms = [
        _Atom(
            record=r.record, serial=r.serial, name=r.name, altloc=r.altloc,
            resname=r.resname, chain=r.chain, resseq=r.resseq, icode=r.icode,
            x=r.x, y=r.y, z=r.z, element=r.element,
            line_idx=r.line_idx, model=r.model,
        )
        for r in records
    ]
    return lines, atoms


def _group_by_residue(atoms: List[_Atom]) -> Dict[Tuple, List[_Atom]]:
    """Group atoms into residues keyed by (model, chain, resseq, icode, resname).

    The insertion code and model number are part of the key: keying on
    ``(chain, resseq)`` alone merged residues ``47`` and ``47A``, and merged the
    same residue across NMR models.
    """
    groups: Dict[Tuple, List[_Atom]] = defaultdict(list)
    for a in atoms:
        key = (a.model, a.chain, a.resseq, a.icode, a.resname)
        groups[key].append(a)
    return dict(sorted(groups.items(), key=lambda kv: kv[0][:4]))


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
    Estimate an **L-configured** Cβ position from backbone atoms alone.

    Places a pseudo-Cβ 1.52 Å from Cα, off the N–Cα–C plane on the side an
    L-amino acid's Cβ would occupy.

    .. warning::
       Because the construction fixes the side of the plane, the signed volume
       computed with this pseudo-Cβ is **positive by construction** and carries
       no information about the deposited stereochemistry. It must therefore
       never be used to *classify* a residue — doing so labels every Cβ-less
       residue "L" and reports a false violation for every D-residue with an
       unmodelled side chain. :func:`detect_chirality_violations` instead uses
       Hα when Cβ is absent, and otherwise reports the residue as
       unassignable. This helper is kept for callers that want an idealised
       L-Cβ position (e.g. rebuilding).
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

#: Hα names recognised as the fourth Cα substituent, in preference order.
_HA_NAMES: Tuple[str, ...] = ("HA", "HA2", "1HA", "2HA", "HA1")

#: |V| below which the Cα tetrahedron is treated as planar and its sign as
#: noise. Matches ``chiralfold.auditor.CHIRALITY_PLANAR_EPS``.
CHIRALITY_PLANAR_EPS = 0.05


def _residue_substituents(
    res_atoms: List["_Atom"],
) -> Tuple[Optional[Tuple["_Atom", "_Atom", "_Atom", "_Atom"]], float, str]:
    """Pick the four atoms defining a residue's Cα tetrahedron.

    Returns ``((ca, n, c, fourth), sign, basis)``. *sign* is ``+1`` when the
    fourth substituent is Cβ and ``-1`` when it is Hα, whose tetrahedron carries
    the opposite convention (Hα sits on the far side of the N/C/Cβ plane).
    *basis* is ``"cb"``, ``"ha"``, ``"missing_backbone"`` or
    ``"no_fourth_substituent"``.
    """
    atom_by_name = {a.name: a for a in res_atoms}
    ca_atom = atom_by_name.get("CA")
    n_atom = atom_by_name.get("N")
    c_atom = atom_by_name.get("C")
    if ca_atom is None or n_atom is None or c_atom is None:
        return None, 0.0, "missing_backbone"

    cb_atom = atom_by_name.get("CB")
    if cb_atom is not None:
        return (ca_atom, n_atom, c_atom, cb_atom), 1.0, "cb"

    ha_atom = next(
        (atom_by_name[nm] for nm in _HA_NAMES if nm in atom_by_name), None
    )
    if ha_atom is None:
        return None, 0.0, "no_fourth_substituent"
    return (ca_atom, n_atom, c_atom, ha_atom), -1.0, "ha"


def _residue_chirality(res_atoms: List["_Atom"]) -> Tuple[Optional[str], float, str]:
    """Classify one residue's Cα stereochemistry from its coordinates.

    Returns ``(found, signed_volume, basis)`` where *found* is ``"L"``, ``"D"``
    or ``None`` when the residue cannot be assigned, and *basis* is one of
    ``"cb"``, ``"ha"``, ``"missing_backbone"``, ``"no_fourth_substituent"`` or
    ``"planar"``.
    """
    quad, sign, basis = _residue_substituents(res_atoms)
    if quad is None:
        return None, 0.0, basis

    ca, n, c, fourth = quad
    sv = sign * _signed_volume(ca.xyz, n.xyz, c.xyz, fourth.xyz)
    if abs(sv) < CHIRALITY_PLANAR_EPS:
        return None, sv, "planar"
    return ("L" if sv > _SV_THRESHOLD else "D"), sv, basis


def _signed_volumes_batch(
    ca: np.ndarray, n: np.ndarray, c: np.ndarray, fourth: np.ndarray
) -> np.ndarray:
    """Row-wise ``(N−Cα)·[(C−Cα)×(X−Cα)]`` over ``(m, 3)`` arrays."""
    v1 = n - ca
    v2 = c - ca
    v3 = fourth - ca
    return np.einsum("ij,ij->i", v1, np.cross(v2, v3))


def _expected_chirality(resname: str, res_atoms: List["_Atom"]) -> str:
    """Expected handedness from the residue name, falling back to record type."""
    if resname in D_RESNAMES:
        return "D"
    if resname in L_RESNAMES:
        return "L"
    # Unknown residue — try to infer from HETATM vs ATOM record
    return "D" if any(a.record == "HETATM" for a in res_atoms) else "L"


def detect_chirality_violations(pdb_path: str) -> Dict:
    """
    Parse a PDB file and detect chirality violations at Cα stereocenters.

    For ATOM records (standard L-amino acids): expects L-chirality
    (positive signed volume). For HETATM records with D-residue names
    (DAL, DTR, etc.): expects D-chirality (negative signed volume).
    Glycine and proline are skipped.

    A residue is only classified when a fourth Cα substituent is present:
    Cβ if modelled, otherwise Hα (with the opposite sign convention). Residues
    with neither, and residues whose tetrahedron is degenerate
    (``|V| < CHIRALITY_PLANAR_EPS``), are reported as *unassignable* rather than
    forced into a class. Earlier releases substituted an idealised L-Cβ built
    from the backbone, which is positive by construction and therefore reported
    a false violation for every D-residue with an unmodelled side chain.

    Args:
        pdb_path: Path to the PDB file to analyse.

    Returns:
        dict with keys:
          - n_residues (int): Total residues parsed (all models).
          - n_models (int): Number of models present.
          - n_checked (int): Residues where chirality was evaluated.
          - n_unassignable (int): Residues skipped for lack of a fourth
            substituent or a degenerate tetrahedron.
          - n_violations (int): Residues with incorrect chirality.
          - pct_violations (float): Percentage of checked residues violated.
          - violations (list[dict]): Per-violation details with resnum,
            resname, chain, icode, model, expected_chirality,
            found_chirality, signed_volume, basis.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    _, atoms = _parse_pdb_full(pdb_path)
    return _detect_from_atoms(atoms)


def _detect_from_atoms(atoms: List["_Atom"]) -> Dict:
    """Detection core shared by the public entry points (no file I/O).

    Signed volumes for every candidate residue are evaluated in one batched
    call: a per-residue ``np.cross``/``np.dot`` pair costs more in numpy
    dispatch than in arithmetic, and this is the inner loop of a PDB-wide scan.
    """
    residues = _group_by_residue(atoms)

    n_residues = len(residues)
    n_models = len({a.model for a in atoms}) if atoms else 0
    n_unassignable = 0

    # ── Gather the tetrahedron of every evaluable residue ─────────────────
    keys: List[Tuple] = []
    bases: List[str] = []
    signs: List[float] = []
    expectations: List[str] = []
    coords: List[Tuple] = []

    for key, res_atoms in residues.items():
        resname = key[4]
        # Skip glycine (achiral) and proline (ring complicates geometry)
        if resname in GLYCINE_RESNAMES or resname in {"PRO", "DPR"}:
            continue

        quad, sign, basis = _residue_substituents(res_atoms)
        if basis == "missing_backbone":
            continue  # Cannot evaluate without backbone atoms
        if quad is None:
            n_unassignable += 1
            continue

        ca, n_at, c_at, fourth = quad
        keys.append(key)
        bases.append(basis)
        signs.append(sign)
        expectations.append(_expected_chirality(resname, res_atoms))
        coords.append((
            (ca.x, ca.y, ca.z), (n_at.x, n_at.y, n_at.z),
            (c_at.x, c_at.y, c_at.z), (fourth.x, fourth.y, fourth.z),
        ))

    violations: List[Dict] = []
    n_checked = 0

    if coords:
        arr = np.asarray(coords, dtype=float)   # (m, 4, 3)
        sv = np.asarray(signs) * _signed_volumes_batch(
            arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
        )
        assignable = np.abs(sv) >= CHIRALITY_PLANAR_EPS
        n_unassignable += int((~assignable).sum())
        n_checked = int(assignable.sum())

        found_l = sv > _SV_THRESHOLD
        for i in np.flatnonzero(assignable):
            expected = expectations[i]
            found = "L" if found_l[i] else "D"
            if found == expected:
                continue
            model, chain, resseq, icode, resname = keys[i]
            violations.append({
                "resnum":             resseq,
                "resname":            resname,
                "chain":              chain,
                "icode":              icode,
                "model":              model,
                "expected_chirality": expected,
                "found_chirality":    found,
                "signed_volume":      float(sv[i]),
                "basis":              bases[i],
            })

    n_violations  = len(violations)
    pct_violations = (n_violations / n_checked * 100.0) if n_checked > 0 else 0.0

    return {
        "n_residues":    n_residues,
        "n_models":      n_models,
        "n_checked":     n_checked,
        "n_unassignable": n_unassignable,
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

    # Parse once, then detect from the in-memory atoms. Earlier releases parsed
    # the input up to four times per correction (detect, re-parse for lines,
    # verify, and again in correct_af3_output).
    lines, atoms = _parse_pdb_full(pdb_path)
    before_report = _detect_from_atoms(atoms)
    violations = before_report["violations"]

    if not violations:
        # No corrections needed — copy file as-is
        with open(output_path, "w") as fh:
            fh.writelines(lines)
        return {
            "n_corrected":       0,
            "violations_before": 0,
            "violations_after":  0,
            "output_path":       output_path,
        }

    residues = _group_by_residue(atoms)

    # Violation lookup keyed on the full residue identity. Keying on
    # (chain, resnum) alone made a violation at residue 47 also "correct"
    # residue 47A, and a violation in model 1 also flip models 2..n.
    viol_keys = {
        (v["model"], v["chain"], v["resnum"], v["icode"]) for v in violations
    }

    # Work on a mutable copy of lines
    new_lines = list(lines)
    n_corrected = 0
    n_skipped = 0

    for (model, chain, resseq, icode, _resname), res_atoms in residues.items():
        if (model, chain, resseq, icode) not in viol_keys:
            continue

        atom_by_name = {a.name: a for a in res_atoms}
        ca_atom = atom_by_name.get("CA")
        n_atom  = atom_by_name.get("N")
        c_atom  = atom_by_name.get("C")
        cb_atom = atom_by_name.get("CB")

        if ca_atom is None or n_atom is None or c_atom is None:
            continue

        ca_pos = ca_atom.xyz
        n_pos  = n_atom.xyz
        c_pos  = c_atom.xyz

        if cb_atom is not None:
            mirror_plane_third = cb_atom.xyz
        else:
            # Without Cβ the reflection plane would have to be invented, and
            # reflecting across an invented plane does not invert the real
            # stereocentre. Detection only flags such residues via Hα, so
            # reflect across N, C, Hα instead when that is what was measured.
            ha_atom = next(
                (atom_by_name[nm] for nm in _HA_NAMES if nm in atom_by_name), None
            )
            if ha_atom is None:
                n_skipped += 1
                continue
            mirror_plane_third = ha_atom.xyz

        # Reflect Cα across the plane through N, C and the third substituent
        ca_new = _reflect_across_plane(ca_pos, n_pos, c_pos, mirror_plane_third)

        # Update Cα in the line buffer
        if ca_atom.line_idx >= 0:
            new_lines[ca_atom.line_idx] = _rewrite_atom_coords(
                new_lines[ca_atom.line_idx], ca_new[0], ca_new[1], ca_new[2]
            )

        # Also adjust the Cα-bound hydrogens by reflecting them the same way, so
        # the local geometry stays consistent. All Hα variants are handled, not
        # just a literal "HA" (Gly-style HA2/HA3 were previously left behind).
        if cb_atom is not None:
            for nm in _HA_NAMES:
                ha = atom_by_name.get(nm)
                if ha is None or ha.line_idx < 0:
                    continue
                ha_new = _reflect_across_plane(
                    ha.xyz, n_pos, c_pos, mirror_plane_third
                )
                new_lines[ha.line_idx] = _rewrite_atom_coords(
                    new_lines[ha.line_idx], ha_new[0], ha_new[1], ha_new[2]
                )

        n_corrected += 1

    # Write output
    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.writelines(new_lines)

    # Verify corrections by re-reading what was actually written
    after_report = detect_chirality_violations(output_path)

    return {
        "n_corrected":       n_corrected,
        "n_skipped":         n_skipped,
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

    # Step 1: Detect (single parse, reused below)
    _, atoms = _parse_pdb_full(pdb_path)
    before = _detect_from_atoms(atoms)

    if before["n_violations"] == 0:
        # Nothing to do — just copy
        import shutil
        shutil.copy2(pdb_path, output_path)
        after = before
        correction = {
            "n_corrected":       0,
            "n_skipped":         0,
            "violations_before": 0,
            "violations_after":  0,
            "output_path":       output_path,
        }
    else:
        # Step 2: Correct (re-detects internally, then verifies from the file it
        # wrote — so the reported "after" numbers describe the output on disk)
        correction = correct_chirality(pdb_path, output_path)
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
