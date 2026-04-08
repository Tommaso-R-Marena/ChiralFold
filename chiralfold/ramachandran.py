"""
ChiralFold — Ramachandran Scoring
===================================

Provides per-residue φ/ψ classification and conformer-ensemble filtering
based on Ramachandran plot regions.  Works for both L- and D-amino acids,
glycine, and proline.

Ramachandran region definitions
--------------------------------
All regions defined as (φ_min, φ_max, ψ_min, ψ_max) tuples and stored in
``_REGIONS``.  For D-amino acids the signs of φ and ψ are both mirrored,
so the same logic applies after negation.

References
----------
Lovell et al. (2003) Proteins 50:437–450.  Approximate boundaries used;
for production use replace the polygon lookup with MolProbity tables.
"""

from __future__ import annotations

import math
import warnings
from typing import List, Optional, Sequence, Tuple

import numpy as np

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# Region definitions for L-amino acids (φ, ψ in degrees)
# Each entry: (label, φ_min, φ_max, ψ_min, ψ_max)
# ─────────────────────────────────────────────────────────────────────────────

# General (non-Gly, non-Pro) L-amino acids
_GENERAL_FAVORED: List[Tuple[float, float, float, float]] = [
    # α-helix core
    (-80.0, -40.0, -60.0, -20.0),
    # β-sheet core
    (-160.0, -80.0,  80.0, 180.0),
    (-160.0, -80.0, -180.0, -150.0),
    # left-handed α (rare but allowed)
    (40.0,  80.0,  20.0,  60.0),
    # PPII-like
    (-100.0, -50.0, 110.0, 170.0),
]

_GENERAL_ALLOWED: List[Tuple[float, float, float, float]] = [
    # Extended α-helix
    (-100.0, -30.0, -80.0,  10.0),
    # Extended β-sheet
    (-180.0, -55.0,  60.0, 180.0),
    (-180.0, -55.0, -180.0, -120.0),
    # ε region
    (-130.0, -80.0,  40.0,  80.0),
    # Left-handed extended
    (20.0, 100.0, -30.0,  80.0),
    # γ-turn
    (50.0, 100.0, -50.0,  20.0),
]

# Glycine — symmetric, larger allowed area
_GLY_FAVORED: List[Tuple[float, float, float, float]] = [
    (-80.0, -40.0, -60.0, -20.0),
    (-160.0, -80.0,  80.0, 180.0),
    (-160.0, -80.0, -180.0, -150.0),
    (40.0,   80.0,  20.0,   60.0),
    (-100.0, -50.0, 110.0, 170.0),
    # Mirror (glycine is achiral)
    (40.0,   80.0,  20.0,   60.0),
    (80.0,  160.0, -180.0, -80.0),
    (80.0,  160.0, 150.0, 180.0),
    (20.0,   80.0,  -60.0,  20.0),
    (50.0,  100.0, -110.0, -170.0),
]

_GLY_ALLOWED: List[Tuple[float, float, float, float]] = [
    # All four quadrants broadly allowed
    (-180.0, 180.0, -180.0, 180.0),
]

# Proline — φ constrained near -60°
_PRO_FAVORED: List[Tuple[float, float, float, float]] = [
    (-80.0, -40.0, -60.0, -20.0),
    (-80.0, -40.0, 110.0, 180.0),
    (-80.0, -40.0, -180.0, -160.0),
]

_PRO_ALLOWED: List[Tuple[float, float, float, float]] = [
    (-100.0, -30.0, -80.0, 10.0),
    (-100.0, -30.0,  80.0, 180.0),
    (-100.0, -30.0, -180.0, -130.0),
]


def _in_any_box(
    phi: float,
    psi: float,
    boxes: List[Tuple[float, float, float, float]],
) -> bool:
    """Return True if (phi, psi) lies inside any of the rectangular boxes."""
    for phi_min, phi_max, psi_min, psi_max in boxes:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def score_ramachandran(
    phi: float,
    psi: float,
    residue_type: str = "general",
) -> str:
    """
    Classify a (φ, ψ) point as 'favored', 'allowed', or 'outlier'.

    Args:
        phi: φ dihedral in degrees (−180 to +180).
        psi: ψ dihedral in degrees (−180 to +180).
        residue_type: One of 'general', 'glycine', 'proline',
                      'D-general', 'D-proline'.  D-amino acids have
                      their φ/ψ negated before lookup against the L-tables.

    Returns:
        'favored', 'allowed', or 'outlier'.

    Examples
    --------
    >>> score_ramachandran(-60, -45)            # α-helix
    'favored'
    >>> score_ramachandran(-120, 130)           # β-sheet
    'favored'
    >>> score_ramachandran(60, 45, 'D-general') # D α-helix
    'favored'
    """
    rtype = residue_type.lower()

    # Mirror coordinates for D-amino acids
    if rtype.startswith("d-"):
        phi = -phi
        psi = -psi
        rtype = rtype[2:]  # 'general' or 'proline'

    if rtype in ("glycine", "gly", "g"):
        favored_boxes = _GLY_FAVORED
        allowed_boxes = _GLY_ALLOWED
    elif rtype in ("proline", "pro", "p"):
        favored_boxes = _PRO_FAVORED
        allowed_boxes = _PRO_ALLOWED
    else:
        favored_boxes = _GENERAL_FAVORED
        allowed_boxes = _GENERAL_ALLOWED

    if _in_any_box(phi, psi, favored_boxes):
        return "favored"
    if _in_any_box(phi, psi, allowed_boxes):
        return "allowed"
    return "outlier"


# ─────────────────────────────────────────────────────────────────────────────
# Dihedral helpers (shared with auditor.py via direct import)
# ─────────────────────────────────────────────────────────────────────────────

def _dihedral_deg(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
) -> float:
    """
    Compute the dihedral angle (p1-p2-p3-p4) in degrees using the
    BioPython/IUPAC convention (−180 to +180).

    Positive values correspond to clockwise rotation when viewed
    along the p2→p3 bond axis.
    """
    b0 = -(p2 - p1)   # negated for BioPython/IUPAC sign convention
    b1 = p3 - p2
    b2 = p4 - p3

    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-8:
        return float("nan")
    b1_u = b1 / b1_norm

    # Components of b0 and b2 orthogonal to b1
    v = b0 - np.dot(b0, b1_u) * b1_u
    w = b2 - np.dot(b2, b1_u) * b1_u

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_u, v), w)

    return math.degrees(math.atan2(y, x))


def _get_backbone_coords_from_mol(mol, conf_id: int = 0):
    """
    Extract per-residue backbone atom coordinates from an RDKit molecule.

    Returns a list of dicts, each with keys:
        'res_idx'  (int), 'N', 'CA', 'C', 'O'  (np.ndarray or None)

    This is a heuristic that identifies backbone atoms by connectivity
    pattern: N bonded to CA, CA bonded to C, C bonded to N (next residue)
    and to O (carbonyl).

    Relies on the same peptide-bond detection logic as geometry.py.
    """
    from .geometry import find_peptide_bonds

    conf = mol.GetConformer(conf_id)

    def _pos(idx):
        p = conf.GetAtomPosition(idx)
        return np.array([p.x, p.y, p.z])

    peptide_bonds = find_peptide_bonds(mol)
    if not peptide_bonds:
        return []

    # Collect unique backbone C, N, CA atoms
    # Each peptide bond tuple: (ca1_idx, c_idx, n_idx, ca2_idx)
    residue_data: dict = {}

    for ca1, c, n, ca2 in peptide_bonds:
        # Residue i: ca1, c
        for ca_idx in (ca1, ca2):
            if ca_idx not in residue_data:
                residue_data[ca_idx] = {"CA": _pos(ca_idx), "C": None,
                                         "N": None, "O": None}

        # CA1's carbonyl C
        residue_data[ca1]["C"] = _pos(c)

        # CA2's N
        residue_data[ca2]["N"] = _pos(n)

        # Find O bonded to C
        c_atom = mol.GetAtomWithIdx(c)
        for nb in c_atom.GetNeighbors():
            if nb.GetSymbol() == "O":
                bnd = mol.GetBondBetweenAtoms(c, nb.GetIdx())
                if bnd and bnd.GetBondTypeAsDouble() >= 1.5:
                    residue_data[ca1]["O"] = _pos(nb.GetIdx())
                    break

    # Build ordered list by CA atom index (proxy for sequence order)
    result = []
    for i, (ca_idx, d) in enumerate(sorted(residue_data.items())):
        result.append({
            "res_idx": i,
            "ca_atom_idx": ca_idx,
            "N": d["N"],
            "CA": d["CA"],
            "C": d["C"],
            "O": d["O"],
        })

    return result


def _compute_phi_psi_from_residues(
    residues: List[dict],
) -> List[Tuple[Optional[float], Optional[float]]]:
    """
    Given an ordered list of residue backbone dicts (from
    _get_backbone_coords_from_mol), compute (φ, ψ) pairs.

    Returns a list of (phi, psi) tuples; None where undefined
    (first/last residue may lack one angle).
    """
    angles = []
    n = len(residues)
    for i, res in enumerate(residues):
        phi = psi = None

        # φ = C(i-1) - N(i) - CA(i) - C(i)
        if i > 0:
            c_prev = residues[i - 1].get("C")
            n_curr = res.get("N")
            ca_curr = res.get("CA")
            c_curr = res.get("C")
            if all(v is not None for v in [c_prev, n_curr, ca_curr, c_curr]):
                phi = _dihedral_deg(c_prev, n_curr, ca_curr, c_curr)

        # ψ = N(i) - CA(i) - C(i) - N(i+1)
        if i < n - 1:
            n_curr = res.get("N")
            ca_curr = res.get("CA")
            c_curr = res.get("C")
            n_next = residues[i + 1].get("N")
            if all(v is not None for v in [n_curr, ca_curr, c_curr, n_next]):
                psi = _dihedral_deg(n_curr, ca_curr, c_curr, n_next)

        angles.append((phi, psi))

    return angles


# ─────────────────────────────────────────────────────────────────────────────
# Conformer filtering
# ─────────────────────────────────────────────────────────────────────────────

def filter_conformers_by_ramachandran(
    mol,
    residues: List[str],
    conf_ids: Sequence[int],
    min_pct_favored: float = 50.0,
) -> List[int]:
    """
    Filter a conformer ensemble to keep only conformers that pass a
    Ramachandran quality threshold.

    Args:
        mol: RDKit Mol with explicit hydrogens and multiple 3D conformers.
        residues: Per-residue type strings, one per residue.  Each entry
                  should be one of: 'general', 'glycine', 'proline',
                  'D-general', 'D-proline'.  Use one-letter codes or
                  three-letter codes and they will be normalised.
        conf_ids: Sequence of conformer IDs to evaluate.
        min_pct_favored: Minimum percentage of favored residues required
                         to keep a conformer (default 50%).

    Returns:
        Sorted list of conformer IDs that pass the threshold.

    Notes
    -----
    The function uses ``geometry.find_peptide_bonds`` to identify backbone
    atoms and computes φ/ψ angles directly from 3D coordinates — no SMILES
    topology needed.  Residues without both angles computable (i.e. terminal
    residues) are excluded from the denominator.
    """
    # Normalise residue types
    _ONE_GLY = {"G", "GLY"}
    _ONE_PRO = {"P", "PRO"}

    def _normalise(r: str, is_d: bool) -> str:
        ru = r.strip().upper()
        if ru in _ONE_GLY:
            return "glycine"
        if ru in _ONE_PRO:
            return "D-proline" if is_d else "proline"
        return "D-general" if is_d else "general"

    # Detect if residue strings already carry 'D-' prefix
    normalised_residues = []
    for r in residues:
        is_d = r.upper().startswith("D-") or r.upper().startswith("D_")
        base = r.lstrip("Dd-_")
        normalised_residues.append(_normalise(base, is_d))

    passing = []
    for cid in conf_ids:
        if cid >= mol.GetNumConformers():
            continue

        try:
            backbone = _get_backbone_coords_from_mol(mol, conf_id=cid)
        except Exception:
            continue

        if not backbone:
            continue

        phi_psi_list = _compute_phi_psi_from_residues(backbone)

        n_favored = 0
        n_total = 0

        for idx, (phi, psi) in enumerate(phi_psi_list):
            if phi is None or psi is None:
                continue
            if math.isnan(phi) or math.isnan(psi):
                continue

            rtype = (normalised_residues[idx]
                     if idx < len(normalised_residues)
                     else "general")

            result = score_ramachandran(phi, psi, rtype)
            n_total += 1
            if result == "favored":
                n_favored += 1

        if n_total == 0:
            # No computable angles — keep the conformer conservatively
            passing.append(cid)
            continue

        pct_favored = 100.0 * n_favored / n_total
        if pct_favored >= min_pct_favored:
            passing.append(cid)

    return sorted(passing)
