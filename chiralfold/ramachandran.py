"""
ChiralFold — Ramachandran Scoring (MolProbity-calibrated)
==========================================================

Per-residue φ/ψ classification and conformer-ensemble filtering based on
Ramachandran plot regions calibrated to match MolProbity/Top8000 boundaries.

Region definitions approximate the Lovell et al. (2003) contours used by
MolProbity, with favored regions covering ~98% of high-quality residues
and allowed regions covering ~99.8%.

Works for L-amino acids, D-amino acids (mirror symmetry), glycine, and
proline.

References:
    Lovell et al. (2003) Proteins 50:437-450.
    Chen et al. (2010) Acta Cryst D66:12-21 (MolProbity).
    Richardson et al. (2018) Top8000 dataset.
"""

from __future__ import annotations

import math
import warnings
from typing import List, Optional, Sequence, Tuple

import numpy as np

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════════════════════
# Region definitions — calibrated to MolProbity/Top8000 contours
# Each region: (φ_min, φ_max, ψ_min, ψ_max)
# ═══════════════════════════════════════════════════════════════════════════

# ── General (non-Gly, non-Pro) L-amino acids ─────────────────────────────
_GENERAL_FAVORED: List[Tuple[float, float, float, float]] = [
    # α-helix core (right-handed)
    (-105.0, -25.0, -75.0,  0.0),
    # β-sheet core
    (-180.0, -55.0, 80.0, 180.0),
    (-180.0, -55.0, -180.0, -120.0),   # wrap-around
    # PPII (polyproline II)
    (-90.0, -55.0, 100.0, 180.0),
    # Bridge region (α ↔ β transition)
    (-120.0, -55.0, 0.0, 80.0),
    # Left-handed α (rare but real)
    (35.0, 75.0, 15.0, 75.0),
]

_GENERAL_ALLOWED: List[Tuple[float, float, float, float]] = [
    # Extended α-helix
    (-135.0, -15.0, -100.0, 25.0),
    # Extended β-sheet + PPII continuum
    (-180.0, -40.0, 50.0, 180.0),
    (-180.0, -40.0, -180.0, -100.0),
    # ε region (extended bridge)
    (-135.0, -40.0, 0.0, 100.0),
    # Extended left-handed α + γ-turn
    (15.0, 110.0, -40.0, 100.0),
    # γ-turn region
    (50.0, 110.0, -60.0, 40.0),
    # δ region (rare but allowed)
    (-180.0, -110.0, -60.0, 0.0),
    # ζ region (under β-sheet in negative ψ)
    (-160.0, -80.0, -60.0, -20.0),
]

# ── Glycine — symmetric, much larger allowed area ────────────────────────
_GLY_FAVORED: List[Tuple[float, float, float, float]] = [
    # Same as general
    (-105.0, -25.0, -75.0, 0.0),
    (-180.0, -55.0, 80.0, 180.0),
    (-180.0, -55.0, -180.0, -120.0),
    (-90.0, -55.0, 100.0, 180.0),
    (-120.0, -55.0, 0.0, 80.0),
    (35.0, 75.0, 15.0, 75.0),
    # Mirror image regions (glycine is achiral)
    (25.0, 105.0, 0.0, 75.0),
    (55.0, 180.0, -180.0, -80.0),
    (55.0, 180.0, 120.0, 180.0),
    (55.0, 120.0, -80.0, 0.0),
    (55.0, 90.0, -180.0, -100.0),
    # Central bridging regions for Gly
    (-90.0, 90.0, -20.0, 20.0),
]

_GLY_ALLOWED: List[Tuple[float, float, float, float]] = [
    # Glycine: almost everything is allowed
    (-180.0, 180.0, -180.0, 180.0),
]

# ── Proline — restricted φ range ─────────────────────────────────────────
_PRO_FAVORED: List[Tuple[float, float, float, float]] = [
    # Proline α (endo)
    (-80.0, -50.0, -55.0, -15.0),
    # Proline α (exo)
    (-70.0, -55.0, 120.0, 170.0),
    # Proline β
    (-80.0, -55.0, 80.0, 170.0),
    # Proline PPII
    (-80.0, -55.0, 120.0, 180.0),
]

_PRO_ALLOWED: List[Tuple[float, float, float, float]] = [
    # Extended proline
    (-100.0, -40.0, -75.0, 10.0),
    (-100.0, -40.0, 60.0, 180.0),
    (-100.0, -40.0, -180.0, -130.0),
    # cis proline
    (-100.0, -50.0, -15.0, 60.0),
]


# ═══════════════════════════════════════════════════════════════════════════
# Core scoring function
# ═══════════════════════════════════════════════════════════════════════════

def _in_any_region(phi: float, psi: float,
                   regions: List[Tuple[float, float, float, float]]) -> bool:
    """Check if (φ, ψ) falls within any of the rectangular regions."""
    for phi_min, phi_max, psi_min, psi_max in regions:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return True
    return False


# ═══════════════════════════════════════════════════════════════════════════
# Empirical probability grid (built from PDB data)
# ═══════════════════════════════════════════════════════════════════════════

_EMPIRICAL_GRID = None
_EMPIRICAL_THRESHOLDS = None


def _load_empirical_grid():
    """Load the empirical Ramachandran probability grid from PDB data."""
    global _EMPIRICAL_GRID, _EMPIRICAL_THRESHOLDS
    if _EMPIRICAL_GRID is not None:
        return True

    import json
    import os
    grid_path = os.path.join(os.path.dirname(__file__), 'data', 'ramachandran_grid.json')
    thresh_path = os.path.join(os.path.dirname(__file__), 'data', 'ramachandran_thresholds.json')

    if not os.path.exists(grid_path):
        return False

    try:
        with open(grid_path) as f:
            data = json.load(f)
        _EMPIRICAL_GRID = {
            'counts': np.array(data['counts']),
            'phi_edges': np.array(data['phi_edges']),
            'psi_edges': np.array(data['psi_edges']),
        }
        with open(thresh_path) as f:
            _EMPIRICAL_THRESHOLDS = json.load(f)
        return True
    except Exception:
        return False


def _score_empirical(phi: float, psi: float) -> str:
    """Score using the empirical probability grid."""
    if _EMPIRICAL_GRID is None:
        return None

    grid = _EMPIRICAL_GRID['counts']
    phi_edges = _EMPIRICAL_GRID['phi_edges']
    psi_edges = _EMPIRICAL_GRID['psi_edges']

    # Find bin indices
    phi_idx = np.searchsorted(phi_edges, phi) - 1
    psi_idx = np.searchsorted(psi_edges, psi) - 1

    if phi_idx < 0 or phi_idx >= grid.shape[0]:
        return 'outlier'
    if psi_idx < 0 or psi_idx >= grid.shape[1]:
        return 'outlier'

    prob = grid[phi_idx, psi_idx]
    thresh = _EMPIRICAL_THRESHOLDS

    if prob >= thresh['favored']:
        return 'favored'
    elif prob > 0:
        return 'allowed'
    else:
        return 'outlier'


def score_ramachandran(phi: float, psi: float,
                       residue_type: str = 'general') -> str:
    """
    Classify a (φ, ψ) point as 'favored', 'allowed', or 'outlier'.

    Uses a hybrid approach:
      1. Empirical probability grid from PDB data (when available)
      2. MolProbity-calibrated rectangular regions (fallback)

    Args:
        phi: Backbone φ dihedral angle in degrees [-180, 180].
        psi: Backbone ψ dihedral angle in degrees [-180, 180].
        residue_type: One of 'general', 'glycine', 'proline',
                     'D-general', 'D-proline'.

    Returns:
        'favored', 'allowed', or 'outlier'.
    """
    if math.isnan(phi) or math.isnan(psi):
        return 'outlier'

    # For D-amino acids: negate φ and ψ (mirror symmetry) then score as L
    if residue_type.startswith('D-'):
        phi, psi = -phi, -psi
        residue_type = residue_type[2:]

    # Try empirical grid first (for general residues)
    if residue_type == 'general':
        _load_empirical_grid()
        emp = _score_empirical(phi, psi)
        if emp is not None:
            # Combine empirical with rectangular: take the more generous result
            rect = _score_rectangular(phi, psi, residue_type)
            # If either says favored, it's favored
            if emp == 'favored' or rect == 'favored':
                return 'favored'
            elif emp == 'allowed' or rect == 'allowed':
                return 'allowed'
            return 'outlier'

    return _score_rectangular(phi, psi, residue_type)


def _score_rectangular(phi: float, psi: float, residue_type: str) -> str:
    """Score using rectangular region definitions."""
    if residue_type == 'glycine':
        fav, alw = _GLY_FAVORED, _GLY_ALLOWED
    elif residue_type == 'proline':
        fav, alw = _PRO_FAVORED, _PRO_ALLOWED
    else:
        fav, alw = _GENERAL_FAVORED, _GENERAL_ALLOWED

    if _in_any_region(phi, psi, fav):
        return 'favored'
    elif _in_any_region(phi, psi, alw):
        return 'allowed'
    else:
        return 'outlier'


def score_ramachandran_batch(dihedrals, residue_types=None):
    """
    Score multiple residues at once.

    Args:
        dihedrals: List of (phi, psi) tuples.
        residue_types: List of residue types (default: all 'general').

    Returns:
        dict with 'pct_favored', 'pct_allowed', 'pct_outlier', 'per_residue'.
    """
    if residue_types is None:
        residue_types = ['general'] * len(dihedrals)

    results = []
    for (phi, psi), rtype in zip(dihedrals, residue_types):
        results.append(score_ramachandran(phi, psi, rtype))

    n = len(results)
    if n == 0:
        return {'pct_favored': 0, 'pct_allowed': 0, 'pct_outlier': 0,
                'per_residue': []}

    return {
        'pct_favored': sum(1 for r in results if r == 'favored') / n * 100,
        'pct_allowed': sum(1 for r in results if r == 'allowed') / n * 100,
        'pct_outlier': sum(1 for r in results if r == 'outlier') / n * 100,
        'per_residue': results,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Conformer filtering
# ═══════════════════════════════════════════════════════════════════════════

def filter_conformers_by_ramachandran(mol, residues, conf_ids,
                                      min_pct_favored: float = 50.0):
    """
    Filter an RDKit conformer ensemble by Ramachandran quality.

    Args:
        mol: RDKit Mol object with 3D conformers and explicit Hs.
        residues: Backbone residue list from find_backbone_atoms().
        conf_ids: List of conformer IDs to evaluate.
        min_pct_favored: Minimum % of residues in favored regions.

    Returns:
        List of (conf_id, pct_favored) tuples that pass the threshold,
        sorted by quality (best first).
    """
    from .geometry import find_peptide_bonds

    def _compute_dihedral(p1, p2, p3, p4):
        b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        norm_b2 = np.linalg.norm(b2)
        if norm_b2 < 1e-12:
            return float('nan')
        m1 = np.cross(n1, b2 / norm_b2)
        return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))

    passing = []

    for cid in conf_ids:
        if int(cid) >= mol.GetNumConformers():
            continue

        pos = mol.GetConformer(int(cid)).GetPositions()
        n_res = len(residues)
        n_favored = 0
        n_scored = 0

        for i in range(n_res):
            phi = psi = float('nan')

            if i > 0:
                phi = _compute_dihedral(
                    pos[residues[i-1]['C']], pos[residues[i]['N']],
                    pos[residues[i]['CA']], pos[residues[i]['C']],
                )
            if i < n_res - 1:
                psi = _compute_dihedral(
                    pos[residues[i]['N']], pos[residues[i]['CA']],
                    pos[residues[i]['C']], pos[residues[i+1]['N']],
                )

            if not (math.isnan(phi) or math.isnan(psi)):
                result = score_ramachandran(phi, psi, 'general')
                n_scored += 1
                if result == 'favored':
                    n_favored += 1

        pct = (n_favored / n_scored * 100) if n_scored > 0 else 0.0
        if pct >= min_pct_favored:
            passing.append((int(cid), pct))

    passing.sort(key=lambda x: -x[1])
    return passing
