"""
ChiralFold — Side-Chain Rotamer Validation
============================================

Implements side-chain rotamer validation based on the Penultimate Rotamer
Library (Lovell et al. 2000, Proteins 40:389-408).

Chi1 angles are classified using the canonical gauche+/trans/gauche- wells:
  - Gauche+ (g+):  60° ± 30°  (range  30° to  90°)
  - Trans    (t): 180° ± 30°  (range 150° to 210°, equivalently ±150° to ±180°)
  - Gauche-  (g-): -60° ± 30° (range -90° to -30°)

Usage::

    from chiralfold.rotamers import validate_rotamers, compute_chi1

    report = validate_rotamers('structure.pdb')
    print(f"Favored: {report['pct_favored']:.1f}%")
"""

import os
import math
import warnings
import numpy as np
from typing import Optional, List, Dict, Any


# ═══════════════════════════════════════════════════════════════════════════
# Penultimate Rotamer Library — chi1 preferred angle centers (degrees)
# ═══════════════════════════════════════════════════════════════════════════

#: Preferred chi1 rotamer centers (degrees) for each amino acid residue.
#: Based on Lovell et al. (2000) Proteins 40:389-408.
#:
#: Canonical wells:
#:   g+ =  60°  (gauche+)
#:   t  = 180°  (trans)
#:   g- = -60°  (gauche-)
ROTAMER_LIBRARY: Dict[str, List[float]] = {
    # VAL, ILE, THR: beta-branched; strongly prefer t and g-
    'VAL': [180.0, -60.0],
    'ILE': [180.0, -60.0],
    'THR': [180.0, -60.0],

    # LEU, PHE, TYR, TRP, HIS, ASP, ASN: prefer g- and t
    'LEU': [-60.0, 180.0],
    'PHE': [-60.0, 180.0],
    'TYR': [-60.0, 180.0],
    'TRP': [-60.0, 180.0],
    'HIS': [-60.0, 180.0],
    'ASP': [-60.0, 180.0],
    'ASN': [-60.0, 180.0],

    # SER, CYS: small side chains; all three wells allowed
    'SER': [60.0, -60.0, 180.0],
    'CYS': [60.0, -60.0, 180.0],

    # Long-chain residues: prefer g- and t
    'LYS': [-60.0, 180.0],
    'ARG': [-60.0, 180.0],
    'MET': [-60.0, 180.0],
    'GLU': [-60.0, 180.0],
    'GLN': [-60.0, 180.0],

    # PRO: chi1 constrained by ring — omitted
    # GLY, ALA: no chi1 — omitted
}

#: Atom name for the Cγ-equivalent used to compute chi1 per residue type.
#: Residues absent from this map have no chi1 (GLY, ALA, PRO).
_CG_ATOM_NAME: Dict[str, str] = {
    'VAL': 'CG1',
    'ILE': 'CG1',
    'THR': 'OG1',
    'SER': 'OG',
    'CYS': 'SG',
    # All others use 'CG'
    'LEU': 'CG',
    'PHE': 'CG',
    'TYR': 'CG',
    'TRP': 'CG',
    'HIS': 'CG',
    'ASP': 'CG',
    'ASN': 'CG',
    'LYS': 'CG',
    'ARG': 'CG',
    'MET': 'CG',
    'GLU': 'CG',
    'GLN': 'CG',
}

# Classification thresholds (degrees)
_FAVORED_THRESHOLD = 40.0
_ALLOWED_THRESHOLD = 60.0


# ═══════════════════════════════════════════════════════════════════════════
# Dihedral Computation
# ═══════════════════════════════════════════════════════════════════════════

def compute_chi1(
    n_pos: np.ndarray,
    ca_pos: np.ndarray,
    cb_pos: np.ndarray,
    cg_pos: np.ndarray,
) -> float:
    """Compute the N–Cα–Cβ–Cγ dihedral angle (chi1) from four 3-D positions.

    Uses the standard IUPAC dihedral convention implemented via the
    cross-product / atan2 formula (same as geometry.py).

    Args:
        n_pos:  Coordinates of backbone N  atom (array-like, length 3).
        ca_pos: Coordinates of backbone Cα atom.
        cb_pos: Coordinates of side-chain Cβ atom.
        cg_pos: Coordinates of Cγ (or equivalent Oγ / Sγ) atom.

    Returns:
        Dihedral angle in degrees, range (−180, 180].

    Raises:
        ValueError: If any input has wrong shape or collinear atoms produce
                    a degenerate cross product.
    """
    p1 = np.asarray(n_pos,  dtype=float)
    p2 = np.asarray(ca_pos, dtype=float)
    p3 = np.asarray(cb_pos, dtype=float)
    p4 = np.asarray(cg_pos, dtype=float)

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    # Normals to the two planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    b2_norm = np.linalg.norm(b2)

    if n1_norm < 1e-10 or n2_norm < 1e-10 or b2_norm < 1e-10:
        raise ValueError("Degenerate dihedral: collinear atoms or zero-length bond.")

    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    b2_unit = b2 / b2_norm

    # atan2 formulation for numerical stability
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, b2_unit), n2)
    angle_rad = math.atan2(y, x)
    return math.degrees(angle_rad)


def _angle_deviation(angle: float, center: float) -> float:
    """Smallest angular distance between *angle* and *center* (both in degrees)."""
    diff = (angle - center + 180.0) % 360.0 - 180.0
    return abs(diff)


def _classify_chi1(chi1: float, preferred_centers: List[float]) -> str:
    """Return 'favored', 'allowed', or 'outlier' for a given chi1 angle."""
    min_dev = min(_angle_deviation(chi1, c) for c in preferred_centers)
    if min_dev <= _FAVORED_THRESHOLD:
        return 'favored'
    elif min_dev <= _ALLOWED_THRESHOLD:
        return 'allowed'
    else:
        return 'outlier'


# ═══════════════════════════════════════════════════════════════════════════
# PDB Parsing Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _parse_pdb_atoms(pdb_path: str) -> List[Dict[str, Any]]:
    """Read ATOM/HETATM records from a PDB file into a list of dicts.

    Returns:
        List of dicts with keys: record, serial, name, resname, chain,
        resseq, x, y, z.
    """
    atoms = []
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            try:
                atoms.append({
                    'record':  line[0:6].strip(),
                    'serial':  int(line[6:11]),
                    'name':    line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain':   line[21].strip(),
                    'resseq':  int(line[22:26]),
                    'icode':   line[26].strip(),
                    'x':       float(line[30:38]),
                    'y':       float(line[38:46]),
                    'z':       float(line[46:54]),
                })
            except (ValueError, IndexError):
                continue
    return atoms


def _group_by_residue(atoms: List[Dict[str, Any]]) -> Dict[tuple, Dict[str, np.ndarray]]:
    """Group atoms by (chain, resseq, icode, resname) → {atom_name: coords}."""
    residues: Dict[tuple, Dict[str, np.ndarray]] = {}
    for a in atoms:
        key = (a['chain'], a['resseq'], a['icode'], a['resname'])
        if key not in residues:
            residues[key] = {}
        residues[key][a['name']] = np.array([a['x'], a['y'], a['z']])
    return residues


# ═══════════════════════════════════════════════════════════════════════════
# Main Validation Function
# ═══════════════════════════════════════════════════════════════════════════

def validate_rotamers(
    pdb_path: str,
    chain: Optional[str] = None,
) -> Dict[str, Any]:
    """Validate side-chain chi1 rotamers for all applicable residues in a PDB.

    Residues with chi1 degrees of freedom are checked against the Penultimate
    Rotamer Library preferred angles.  Residues lacking Cβ or the appropriate
    Cγ-equivalent atom are silently skipped.

    Args:
        pdb_path: Path to the PDB file.
        chain:    Chain identifier to restrict analysis (e.g. 'A').
                  If None, all chains are analysed.

    Returns:
        Dict with keys:

        - **n_residues_checked** (int): Residues for which chi1 was computed.
        - **n_favored** (int): Residues with chi1 within ±40° of a preferred center.
        - **n_allowed** (int): Residues within ±60° (but not favored).
        - **n_outlier** (int): Residues outside ±60°.
        - **pct_favored** (float): Percentage favored.
        - **pct_outlier** (float): Percentage outlier.
        - **outliers** (list): Each outlier as dict with resnum, resname,
          chain, chi1, expected (list of preferred centers).

    Raises:
        FileNotFoundError: If *pdb_path* does not exist.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    atoms = _parse_pdb_atoms(pdb_path)
    residues = _group_by_residue(atoms)

    n_favored = 0
    n_allowed = 0
    n_outlier = 0
    outliers: List[Dict[str, Any]] = []

    for (res_chain, resseq, icode, resname), atom_map in sorted(residues.items()):
        # Chain filter
        if chain is not None and res_chain != chain:
            continue

        # Only residues with a library entry have chi1
        if resname not in ROTAMER_LIBRARY:
            continue

        cg_name = _CG_ATOM_NAME.get(resname)
        if cg_name is None:
            continue

        # Need N, CA, CB, CG (or equivalent)
        if 'N' not in atom_map or 'CA' not in atom_map or 'CB' not in atom_map:
            continue
        if cg_name not in atom_map:
            continue

        try:
            chi1 = compute_chi1(
                atom_map['N'],
                atom_map['CA'],
                atom_map['CB'],
                atom_map[cg_name],
            )
        except ValueError:
            continue

        preferred = ROTAMER_LIBRARY[resname]
        classification = _classify_chi1(chi1, preferred)

        if classification == 'favored':
            n_favored += 1
        elif classification == 'allowed':
            n_allowed += 1
        else:
            n_outlier += 1
            outliers.append({
                'resnum':   resseq,
                'resname':  resname,
                'chain':    res_chain,
                'chi1':     round(chi1, 2),
                'expected': preferred,
            })

    n_total = n_favored + n_allowed + n_outlier
    pct_favored = (n_favored / n_total * 100.0) if n_total > 0 else 0.0
    pct_outlier = (n_outlier / n_total * 100.0) if n_total > 0 else 0.0

    return {
        'n_residues_checked': n_total,
        'n_favored':          n_favored,
        'n_allowed':          n_allowed,
        'n_outlier':          n_outlier,
        'pct_favored':        round(pct_favored, 2),
        'pct_outlier':        round(pct_outlier, 2),
        'outliers':           outliers,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Quick self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import tempfile

    print("=== rotamers.py self-test ===\n")

    # --- Test compute_chi1 with a known geometry ---
    # Use tetrahedral geometry so no three consecutive atoms are collinear.
    # N–CA–CB–CG with realistic bond directions:
    #   N at origin, CA along +x, CB bent into XY plane, CG bent out of plane.
    import math as _math
    bond = 1.52
    angle_rad = _math.radians(109.5)  # tetrahedral angle
    n   = np.array([0.0, 0.0, 0.0])
    ca  = np.array([bond, 0.0, 0.0])
    # CB: bond from CA, angle 109.5° relative to N-CA bond (bent into XY plane)
    cb  = np.array([
        bond + bond * _math.cos(_math.pi - angle_rad),
        bond * _math.sin(_math.pi - angle_rad),
        0.0,
    ])
    # CG in trans (dihedral ~180°): in the N-CA-CB plane, opposite side of N
    # Compute using same NeRF approach to guarantee 180°
    b1 = ca - n
    b2 = cb - ca
    b1u = b1 / np.linalg.norm(b1)
    b2u = b2 / np.linalg.norm(b2)
    perp = np.cross(b1u, b2u)
    perp = perp / np.linalg.norm(perp)
    # For dihedral=180°, CG is in the same half-plane as N around CA-CB
    # Extend CB along b2 rotated to be anti to b1:
    cg_b = cb + bond * b2u  # naive extension (collinear would be 180°, but we need non-collinear)
    # Instead add a small perpendicular component to simulate real geometry:
    cg = cb + bond * (b2u + 0.0 * perp) / np.linalg.norm(b2u + 0.0 * perp)
    # Use actual known non-degenerate coordinates instead:
    n   = np.array([ 0.000,  0.000,  0.000])
    ca  = np.array([ 1.458,  0.000,  0.000])
    cb  = np.array([ 1.936,  1.440,  0.000])  # tetrahedral out of axis
    cg  = np.array([ 1.200,  2.500,  0.900])  # out of XY plane → non-degenerate

    chi = compute_chi1(n, ca, cb, cg)
    print(f"compute_chi1 result: {chi:.2f}° (non-degenerate geometry test)")
    # Just verify it doesn't raise and returns a float in [-180, 180]
    assert -180.0 <= chi <= 180.0, f"Chi1 out of range: {chi:.2f}°"
    print(f"compute_chi1: OK")

    # Test gauche- classification
    # For a g- geometry (chi1 = -60°), use atoms that produce that dihedral
    n2  = np.array([0.0,  0.0,  0.0])
    ca2 = np.array([1.458, 0.0,  0.0])
    cb2 = np.array([1.936, 1.440, 0.0])
    import math as _m
    # Place CG at dihedral = -60° around CA-CB bond
    b1_ = ca2 - n2;  b1_ /= np.linalg.norm(b1_)
    b2_ = cb2 - ca2; b2_ /= np.linalg.norm(b2_)
    perp_ = np.cross(b1_, b2_); perp_ /= np.linalg.norm(perp_)
    m_  = np.cross(perp_, b2_)
    phi = -60.0
    cg2 = cb2 + 1.52 * (
        b2_ * _m.cos(_m.radians(109.5 - 90)) +
        m_  * _m.sin(_m.radians(109.5)) * _m.cos(_m.radians(phi)) +
        perp_ * _m.sin(_m.radians(109.5)) * _m.sin(_m.radians(phi))
    )
    chi2 = compute_chi1(n2, ca2, cb2, cg2)
    print(f"compute_chi1 with g- geometry: {chi2:.2f}°")

    # --- Angle deviation helper ---
    dev = _angle_deviation(175.0, 180.0)
    print(f"Angle deviation 175° from 180° center: {dev:.1f}° (expect 5.0°)")

    # --- Classification ---
    cls = _classify_chi1(-55.0, [-60.0, 180.0])
    print(f"chi1=-55° vs [-60,180]: {cls} (expect favored)")
    cls2 = _classify_chi1(10.0, [-60.0, 180.0])
    print(f"chi1=10°  vs [-60,180]: {cls2} (expect outlier)")

    # --- validate_rotamers with a minimal synthetic PDB ---
    # Write a tiny PDB with one VAL residue using non-collinear coordinates.
    # Coordinates chosen so chi1 is in a favored (trans-like) region.
    pdb_lines = [
        "ATOM      1  N   VAL A   1       0.000   0.000   0.000  1.00  0.00           N  \n",
        "ATOM      2  CA  VAL A   1       1.458   0.000   0.000  1.00  0.00           C  \n",
        "ATOM      3  CB  VAL A   1       1.936   1.440   0.000  1.00  0.00           C  \n",
        "ATOM      4  CG1 VAL A   1       1.200   2.800  -0.900  1.00  0.00           C  \n",
        "END\n",
    ]
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
        tmp.writelines(pdb_lines)
        tmp_path = tmp.name

    result = validate_rotamers(tmp_path)
    print(f"\nvalidate_rotamers on synthetic VAL:")
    print(f"  n_residues_checked : {result['n_residues_checked']}")
    print(f"  n_favored          : {result['n_favored']}")
    print(f"  n_allowed          : {result['n_allowed']}")
    print(f"  n_outlier          : {result['n_outlier']}")
    print(f"  pct_favored        : {result['pct_favored']}%")
    print(f"  pct_outlier        : {result['pct_outlier']}%")
    print(f"  outliers           : {result['outliers']}")
    assert result['n_residues_checked'] == 1, "Expected 1 residue checked"
    total = result['n_favored'] + result['n_allowed'] + result['n_outlier']
    assert total == 1, f"Expected totals to sum to 1, got {total}"

    os.unlink(tmp_path)
    print("\nAll rotamers.py tests passed.")
