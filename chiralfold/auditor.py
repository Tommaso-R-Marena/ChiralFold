"""
ChiralFold — PDB Structure Quality Auditor (v3)
=================================================

Validates ANY protein structure — pure L, pure D, or mixed L/D — against
six complementary quality criteria and returns a comprehensive report dict.

Quick start
-----------
>>> from chiralfold.auditor import audit_pdb
>>> report = audit_pdb("results/1UBQ.pdb")
>>> print(report["overall_score"])

Checks performed
----------------
1. **Cα Chirality**      — Signed-volume test at each Cα tetrahedron.
2. **Bond Geometry**     — Backbone bond lengths and angles vs ideal values.
3. **Ramachandran**      — φ/ψ region classification (favored/allowed/outlier).
4. **Peptide Planarity** — ω dihedral deviation from ±180°.
5. **Clash Detection**   — Non-bonded atom pairs closer than vdW sum − 0.4 Å.
6. **Summary Score**     — Composite 0–100 quality score.

All parsing is done directly from PDB ATOM/HETATM records with no external
dependencies beyond NumPy.  RDKit is *not* required for audit_pdb().
"""

from __future__ import annotations

import math
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Ideal backbone geometry (Engh & Huber 1991)
IDEAL_BOND_LENGTHS: Dict[str, float] = {
    "N-CA":  1.458,
    "CA-C":  1.525,
    "C-N":   1.329,   # peptide bond (partial double-bond)
    "C-O":   1.231,   # carbonyl
}

IDEAL_BOND_ANGLES: Dict[str, float] = {
    "N-CA-C":  111.0,
    "CA-C-N":  116.2,
    "C-N-CA":  121.7,
}

# 3-sigma tolerances (approximate)
BL_SIGMA = 0.02   # Å  — typical ESD for bond lengths
BA_SIGMA = 2.5    # °  — typical ESD for bond angles

# Van der Waals radii (Bondi 1964, reduced set)
VDW_RADII: Dict[str, float] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "H": 1.20,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
}
VDW_DEFAULT = 1.70

# Bonded atom-name pairs to skip during clash detection
_BACKBONE_BONDS: Set[Tuple[str, str]] = {
    frozenset(["N", "CA"]),
    frozenset(["CA", "C"]),
    frozenset(["C", "O"]),
    frozenset(["C", "OXT"]),
}

# D-amino acid residue names (from pdb_pipeline.py mapping)
D_RESNAMES: Set[str] = {
    "DAL", "DAR", "DSG", "DAS", "DCY", "DGL", "DGN",
    "DHI", "DIL", "DLE", "DLY", "MED", "DPN", "DPR",
    "DSN", "DTH", "DTR", "DTY", "DVA",
}

# Standard L-amino acid residue names
L_RESNAMES: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    # Common variants
    "MSE", "SEP", "TPO", "PTR", "HSD", "HSE", "HSP",
    "HIE", "HID", "HIP", "CYX", "CYM",
}

GLYCINE_RESNAMES: Set[str] = {"GLY"}
PROLINE_RESNAMES: Set[str] = {"PRO", "DPR"}


# ─────────────────────────────────────────────────────────────────────────────
# PDB parsing
# ─────────────────────────────────────────────────────────────────────────────

class _Atom:
    """Lightweight representation of a parsed PDB atom record."""

    __slots__ = [
        "record", "serial", "name", "altloc", "resname",
        "chain", "resseq", "icode",
        "x", "y", "z", "element",
    ]

    def __init__(
        self,
        record: str, serial: int, name: str, altloc: str,
        resname: str, chain: str, resseq: int, icode: str,
        x: float, y: float, z: float, element: str,
    ):
        self.record  = record
        self.serial  = serial
        self.name    = name
        self.altloc  = altloc
        self.resname = resname
        self.chain   = chain
        self.resseq  = resseq
        self.icode   = icode
        self.x       = x
        self.y       = y
        self.z       = z
        self.element = element

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    @property
    def element_upper(self) -> str:
        if self.element:
            return self.element.upper().strip()
        # Infer from atom name (column 13-14, first non-digit char)
        for ch in self.name.strip():
            if ch.isalpha():
                return ch.upper()
        return "C"


def _parse_pdb(pdb_path: str) -> List[_Atom]:
    """
    Parse ATOM and HETATM records from *pdb_path*.

    Rules applied:
    - HOH / WAT / DOD water residues are skipped.
    - Only the first alternative location (altloc ≤ 'A' or blank) is kept.
    - Lines shorter than 54 characters are skipped.
    """
    atoms: List[_Atom] = []
    seen_altloc: Dict[Tuple, str] = {}  # (chain, resseq, icode, name) → first altloc

    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if len(line) < 54:
                continue

            resname = line[17:20].strip()
            if resname in ("HOH", "WAT", "DOD"):
                continue

            altloc = line[16]
            if altloc not in (" ", "A", "1", ""):
                continue  # keep blank and first alternate only

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

            # Deduplicate by (chain, resseq, icode, name) keeping first seen
            key = (chain, resseq, icode, name)
            if key in seen_altloc:
                continue
            seen_altloc[key] = altloc

            atoms.append(_Atom(
                record=record, serial=serial, name=name, altloc=altloc,
                resname=resname, chain=chain, resseq=resseq, icode=icode,
                x=x, y=y, z=z, element=element,
            ))

    return atoms


def _group_by_residue(
    atoms: List[_Atom],
) -> Dict[Tuple, List[_Atom]]:
    """
    Group atoms into residues keyed by (chain, resseq, icode, resname).

    Returns an ordered dict sorted by (chain, resseq, icode).
    """
    groups: Dict[Tuple, List[_Atom]] = defaultdict(list)
    for a in atoms:
        key = (a.chain, a.resseq, a.icode, a.resname)
        groups[key].append(a)

    ordered = dict(sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])))
    return ordered


# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

def _vec_norm(v: np.ndarray) -> float:
    n = np.linalg.norm(v)
    return float(n)


def _bond_length(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(b - a))


def _bond_angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle at vertex b (a-b-c) in degrees."""
    v1 = a - b
    v2 = c - b
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return float("nan")
    cos_a = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return math.degrees(math.acos(cos_a))


def _dihedral_deg(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
) -> float:
    """
    Dihedral angle (p1–p2–p3–p4) in degrees, range −180 to +180.

    Uses the same sign convention as BioPython calc_dihedral:
    positive values correspond to clockwise rotation when viewed
    along the p2→p3 axis.
    """
    b0 = -(p2 - p1)    # negated for BioPython/IUPAC convention
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


def _signed_volume(
    center: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
) -> float:
    """
    Signed volume of the tetrahedron (center, a, b, c).

    Positive → the substituents a, b, c form a right-handed triple
    when viewed from *center*.
    """
    va = a - center
    vb = b - center
    vc = c - center
    return float(np.dot(va, np.cross(vb, vc)))


# ─────────────────────────────────────────────────────────────────────────────
# Check 1: Cα Chirality
# ─────────────────────────────────────────────────────────────────────────────

def _check_chirality(
    residue_groups: Dict[Tuple, List[_Atom]],
) -> dict:
    """
    For each residue, compute the signed volume of the Cα tetrahedron
    (N, C, Cβ, Hα — or infer sign from N/C/Cβ only when H is absent).

    Expected configuration:
    - L-amino acids (ATOM records, L_RESNAMES) → positive signed volume
      (S-configuration, left-handed triple N→C→Cβ at Cα)
    - D-amino acids (HETATM / D_RESNAMES) → negative signed volume
      (R-configuration, mirror image)
    - Glycine → achiral, skip
    """
    n_correct  = 0
    n_wrong    = 0
    n_glycine  = 0
    n_error    = 0
    violations = []

    for key, atoms in residue_groups.items():
        chain, resseq, icode, resname = key

        # Classify residue
        if resname in GLYCINE_RESNAMES:
            n_glycine += 1
            continue

        is_d_residue = resname in D_RESNAMES
        is_l_residue = resname in L_RESNAMES

        # Build name→xyz lookup
        atom_map: Dict[str, np.ndarray] = {a.name: a.xyz for a in atoms}

        # Require Cα
        if "CA" not in atom_map:
            n_error += 1
            continue

        ca = atom_map["CA"]

        # Collect substituents: N, C (backbone carbonyl), CB
        n_pos  = atom_map.get("N")
        c_pos  = atom_map.get("C")
        cb_pos = atom_map.get("CB")

        if n_pos is None or c_pos is None or cb_pos is None:
            # Cannot compute — only CA found or terminal residue
            n_error += 1
            continue

        # Look for Hα (HA, HA2, HA3 …)
        h_pos: Optional[np.ndarray] = None
        for hname in ("HA", "HA2", "1HA", "2HA"):
            if hname in atom_map:
                h_pos = atom_map[hname]
                break

        if h_pos is not None:
            vol = _signed_volume(ca, n_pos, c_pos, cb_pos)
            # With all four substituents use N, C, Cβ, Hα
            vol4 = _signed_volume(ca, n_pos, cb_pos, h_pos)
            # Fall back to three-substituent version (ignore H sign ambiguity)
            signed_vol = vol4
        else:
            # Three-substituent version: N, C, Cβ
            signed_vol = _signed_volume(ca, n_pos, c_pos, cb_pos)

        if abs(signed_vol) < 0.05:
            # Essentially planar — cannot assign chirality
            n_error += 1
            continue

        # L-configuration convention: N→C→CB forms a left-handed triple
        # At Cα the signed volume (N, C, CB) > 0 for L; < 0 for D
        # (This matches the standard CIP S for L-amino acids in PDB frame)
        is_l_geometry = signed_vol > 0
        is_d_geometry = signed_vol < 0

        if is_d_residue:
            expected_d = True
        elif is_l_residue:
            expected_d = False
        else:
            # Unknown residue type — compare to ATOM/HETATM record
            expected_d = (atoms[0].record == "HETATM")

        if expected_d:
            correct = is_d_geometry
        else:
            correct = is_l_geometry

        if correct:
            n_correct += 1
        else:
            n_wrong += 1
            violations.append({
                "chain":    chain,
                "resseq":   resseq,
                "resname":  resname,
                "expected": "D" if expected_d else "L",
                "observed": "D" if is_d_geometry else "L",
                "signed_volume": round(signed_vol, 4),
            })

    total = n_correct + n_wrong
    pct_correct = 100.0 * n_correct / total if total > 0 else 100.0

    return {
        "n_correct":  n_correct,
        "n_wrong":    n_wrong,
        "n_glycine":  n_glycine,
        "n_error":    n_error,
        "pct_correct": round(pct_correct, 2),
        "violations": violations,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Check 2: Bond Geometry
# ─────────────────────────────────────────────────────────────────────────────

def _check_bond_geometry(
    residue_groups: Dict[Tuple, List[_Atom]],
    ordered_keys: List[Tuple],
) -> dict:
    """
    Measure backbone bond lengths and angles; flag deviations > 3σ.

    Returns bl_rmsd, ba_rmsd, and a list of outlier records.
    """
    bl_deviations: List[float] = []
    ba_deviations: List[float] = []
    outliers: List[dict] = []

    # Convert to list of per-residue atom maps
    res_maps: List[Dict[str, np.ndarray]] = []
    res_labels: List[str] = []
    for key in ordered_keys:
        amap = {a.name: a.xyz for a in residue_groups[key]}
        res_maps.append(amap)
        chain, resseq, icode, resname = key
        res_labels.append(f"{chain}:{resname}{resseq}{icode.strip()}")

    n = len(res_maps)
    for i, amap in enumerate(res_maps):
        label = res_labels[i]

        ca = amap.get("CA")
        n_  = amap.get("N")
        c   = amap.get("C")
        o   = amap.get("O")

        # Intra-residue bond lengths
        if n_ is not None and ca is not None:
            bl = _bond_length(n_, ca)
            dev = bl - IDEAL_BOND_LENGTHS["N-CA"]
            bl_deviations.append(dev)
            if abs(dev) > 3 * BL_SIGMA:
                outliers.append({
                    "residue": label, "type": "bond_length",
                    "bond": "N-CA", "value": round(bl, 4),
                    "ideal": IDEAL_BOND_LENGTHS["N-CA"],
                    "deviation": round(dev, 4),
                    "sigma": round(abs(dev) / BL_SIGMA, 1),
                })

        if ca is not None and c is not None:
            bl = _bond_length(ca, c)
            dev = bl - IDEAL_BOND_LENGTHS["CA-C"]
            bl_deviations.append(dev)
            if abs(dev) > 3 * BL_SIGMA:
                outliers.append({
                    "residue": label, "type": "bond_length",
                    "bond": "CA-C", "value": round(bl, 4),
                    "ideal": IDEAL_BOND_LENGTHS["CA-C"],
                    "deviation": round(dev, 4),
                    "sigma": round(abs(dev) / BL_SIGMA, 1),
                })

        if c is not None and o is not None:
            bl = _bond_length(c, o)
            dev = bl - IDEAL_BOND_LENGTHS["C-O"]
            bl_deviations.append(dev)
            if abs(dev) > 3 * BL_SIGMA:
                outliers.append({
                    "residue": label, "type": "bond_length",
                    "bond": "C=O", "value": round(bl, 4),
                    "ideal": IDEAL_BOND_LENGTHS["C-O"],
                    "deviation": round(dev, 4),
                    "sigma": round(abs(dev) / BL_SIGMA, 1),
                })

        # Intra-residue bond angle N-CA-C
        if n_ is not None and ca is not None and c is not None:
            ba = _bond_angle_deg(n_, ca, c)
            if not math.isnan(ba):
                dev = ba - IDEAL_BOND_ANGLES["N-CA-C"]
                ba_deviations.append(dev)
                if abs(dev) > 3 * BA_SIGMA:
                    outliers.append({
                        "residue": label, "type": "bond_angle",
                        "angle": "N-CA-C", "value": round(ba, 2),
                        "ideal": IDEAL_BOND_ANGLES["N-CA-C"],
                        "deviation": round(dev, 2),
                        "sigma": round(abs(dev) / BA_SIGMA, 1),
                    })

        # Inter-residue (peptide bond) geometry
        if i + 1 < n:
            next_map = res_maps[i + 1]
            n_next = next_map.get("N")
            ca_next = next_map.get("CA")

            if c is not None and n_next is not None:
                bl = _bond_length(c, n_next)
                dev = bl - IDEAL_BOND_LENGTHS["C-N"]
                bl_deviations.append(dev)
                if abs(dev) > 3 * BL_SIGMA:
                    outliers.append({
                        "residue": f"{label}→{res_labels[i+1]}",
                        "type": "bond_length",
                        "bond": "C-N(peptide)", "value": round(bl, 4),
                        "ideal": IDEAL_BOND_LENGTHS["C-N"],
                        "deviation": round(dev, 4),
                        "sigma": round(abs(dev) / BL_SIGMA, 1),
                    })

            # CA-C-N angle
            if ca is not None and c is not None and n_next is not None:
                ba = _bond_angle_deg(ca, c, n_next)
                if not math.isnan(ba):
                    dev = ba - IDEAL_BOND_ANGLES["CA-C-N"]
                    ba_deviations.append(dev)
                    if abs(dev) > 3 * BA_SIGMA:
                        outliers.append({
                            "residue": f"{label}→{res_labels[i+1]}",
                            "type": "bond_angle",
                            "angle": "CA-C-N", "value": round(ba, 2),
                            "ideal": IDEAL_BOND_ANGLES["CA-C-N"],
                            "deviation": round(dev, 2),
                            "sigma": round(abs(dev) / BA_SIGMA, 1),
                        })

            # C-N-CA angle
            if c is not None and n_next is not None and ca_next is not None:
                ba = _bond_angle_deg(c, n_next, ca_next)
                if not math.isnan(ba):
                    dev = ba - IDEAL_BOND_ANGLES["C-N-CA"]
                    ba_deviations.append(dev)
                    if abs(dev) > 3 * BA_SIGMA:
                        outliers.append({
                            "residue": f"{label}→{res_labels[i+1]}",
                            "type": "bond_angle",
                            "angle": "C-N-CA", "value": round(ba, 2),
                            "ideal": IDEAL_BOND_ANGLES["C-N-CA"],
                            "deviation": round(dev, 2),
                            "sigma": round(abs(dev) / BA_SIGMA, 1),
                        })

    bl_rmsd = float(np.sqrt(np.mean(np.array(bl_deviations) ** 2))) if bl_deviations else 0.0
    ba_rmsd = float(np.sqrt(np.mean(np.array(ba_deviations) ** 2))) if ba_deviations else 0.0

    return {
        "bl_rmsd":  round(bl_rmsd, 5),
        "ba_rmsd":  round(ba_rmsd, 3),
        "n_bonds_checked": len(bl_deviations),
        "n_angles_checked": len(ba_deviations),
        "outliers": outliers,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Check 3: Ramachandran
# ─────────────────────────────────────────────────────────────────────────────

def _get_rtype(resname: str) -> str:
    """Map a PDB residue name to a Ramachandran region type string."""
    if resname in GLYCINE_RESNAMES:
        return "glycine"
    if resname in PROLINE_RESNAMES:
        if resname == "DPR":
            return "D-proline"
        return "proline"
    if resname in D_RESNAMES:
        return "D-general"
    return "general"


def _check_ramachandran(
    residue_groups: Dict[Tuple, List[_Atom]],
    ordered_keys: List[Tuple],
) -> dict:
    """
    Compute φ/ψ for each residue and classify using score_ramachandran().

    Returns pct_favored, pct_allowed, pct_outlier, and outlier records.
    """
    from .ramachandran import score_ramachandran

    res_maps: List[Dict[str, np.ndarray]] = []
    res_labels: List[str] = []
    res_types: List[str] = []

    for key in ordered_keys:
        amap = {a.name: a.xyz for a in residue_groups[key]}
        res_maps.append(amap)
        chain, resseq, icode, resname = key
        res_labels.append(f"{chain}:{resname}{resseq}{icode.strip()}")
        res_types.append(_get_rtype(resname))

    results = []
    n = len(res_maps)

    for i in range(n):
        amap = res_maps[i]
        ca = amap.get("CA")
        n_  = amap.get("N")
        c   = amap.get("C")

        if ca is None or n_ is None or c is None:
            continue

        phi = psi = None

        # φ = C(i-1) - N(i) - CA(i) - C(i)
        if i > 0:
            c_prev = res_maps[i - 1].get("C")
            if c_prev is not None:
                phi = _dihedral_deg(c_prev, n_, ca, c)

        # ψ = N(i) - CA(i) - C(i) - N(i+1)
        if i < n - 1:
            n_next = res_maps[i + 1].get("N")
            if n_next is not None:
                psi = _dihedral_deg(n_, ca, c, n_next)

        if phi is None or psi is None:
            continue
        if math.isnan(phi) or math.isnan(psi):
            continue

        region = score_ramachandran(phi, psi, res_types[i])
        results.append({
            "label":  res_labels[i],
            "rtype":  res_types[i],
            "phi":    round(phi, 2),
            "psi":    round(psi, 2),
            "region": region,
        })

    n_total   = len(results)
    n_favored = sum(1 for r in results if r["region"] == "favored")
    n_allowed = sum(1 for r in results if r["region"] == "allowed")
    n_outlier = sum(1 for r in results if r["region"] == "outlier")

    pct_favored = 100.0 * n_favored / n_total if n_total else 0.0
    pct_allowed = 100.0 * n_allowed / n_total if n_total else 0.0
    pct_outlier = 100.0 * n_outlier / n_total if n_total else 0.0

    outlier_records = [r for r in results if r["region"] == "outlier"]

    return {
        "n_evaluated":  n_total,
        "n_favored":    n_favored,
        "n_allowed":    n_allowed,
        "n_outlier":    n_outlier,
        "pct_favored":  round(pct_favored, 2),
        "pct_allowed":  round(pct_allowed, 2),
        "pct_outlier":  round(pct_outlier, 2),
        "outliers":     outlier_records,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Check 4: Peptide Planarity (ω dihedral)
# ─────────────────────────────────────────────────────────────────────────────

def _check_planarity(
    residue_groups: Dict[Tuple, List[_Atom]],
    ordered_keys: List[Tuple],
) -> dict:
    """
    Measure ω dihedral (CA_i–C_i–N_{i+1}–CA_{i+1}) for each peptide bond.

    Flags bonds with |ω − 180°| > 6° (trans) or |ω| > 6° (cis).
    Reports % within 6°, mean deviation, and outlier records.
    """
    res_maps: List[Dict[str, np.ndarray]] = []
    res_labels: List[str] = []

    for key in ordered_keys:
        amap = {a.name: a.xyz for a in residue_groups[key]}
        res_maps.append(amap)
        chain, resseq, icode, resname = key
        res_labels.append(f"{chain}:{resname}{resseq}{icode.strip()}")

    deviations: List[float] = []
    omega_values: List[float] = []
    outliers: List[dict] = []
    n = len(res_maps)

    for i in range(n - 1):
        ca_i  = res_maps[i].get("CA")
        c_i   = res_maps[i].get("C")
        n_j   = res_maps[i + 1].get("N")
        ca_j  = res_maps[i + 1].get("CA")

        if any(v is None for v in [ca_i, c_i, n_j, ca_j]):
            continue

        omega = _dihedral_deg(ca_i, c_i, n_j, ca_j)
        if math.isnan(omega):
            continue

        omega_values.append(omega)

        # Deviation from nearest standard value (trans = ±180, cis = 0)
        dev_trans = min(abs(omega - 180.0), abs(omega + 180.0))
        dev_cis   = abs(omega)
        dev = min(dev_trans, dev_cis)
        deviations.append(dev)

        if dev > 6.0:
            bond_type = "cis" if dev_cis < dev_trans else "trans"
            outliers.append({
                "peptide_bond": f"{res_labels[i]}→{res_labels[i+1]}",
                "omega": round(omega, 2),
                "deviation": round(dev, 2),
                "type": bond_type,
            })

    n_total = len(deviations)
    n_within = sum(1 for d in deviations if d <= 6.0)
    pct_within = 100.0 * n_within / n_total if n_total else 100.0
    mean_dev   = float(np.mean(deviations)) if deviations else 0.0

    return {
        "n_bonds_checked":  n_total,
        "n_within_6deg":    n_within,
        "pct_within_6deg":  round(pct_within, 2),
        "mean_deviation":   round(mean_dev, 3),
        "outliers":         outliers,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Check 5: Clash Detection
# ─────────────────────────────────────────────────────────────────────────────

def _vdw_radius(atom: _Atom) -> float:
    """Return the van der Waals radius for an atom."""
    elem = atom.element_upper
    return VDW_RADII.get(elem, VDW_DEFAULT)


def _are_bonded_or_angled(a: _Atom, b: _Atom) -> bool:
    """
    Return True if atoms a and b are 1-2 (bonded) or 1-3 (angle partners)
    and should be excluded from clash detection.

    Rules:
    - Same residue: exclude if distance < 2.5 Å (covers bonds + angle partners).
    - Adjacent residues (|resseq difference| ≤ 1): exclude the C(i)–N(i+1)
      peptide bond, plus 1-3 partners C(i)–CA(i+1) and N(i+1)–CA(i) which
      are connected through the C–N–CA or CA–C–N angle.
    """
    same_res = (a.chain == b.chain
                and a.resseq == b.resseq
                and a.icode == b.icode)

    if same_res:
        dist = _bond_length(a.xyz, b.xyz)
        return dist < 2.6   # 1-2 and 1-3 in same residue

    adjacent = (a.chain == b.chain and abs(a.resseq - b.resseq) <= 1)
    if adjacent:
        # Peptide bond atoms: C(i)–N(i+1) is 1-2
        # C(i)–CA(i+1) and O(i)–N(i+1) and C(i)–H(i+1)N are 1-3
        # Exclude by distance threshold for adjacent residues
        dist = _bond_length(a.xyz, b.xyz)
        return dist < 2.7   # peptide bond ~1.33 Å; 1-3 ~2.4–2.5 Å

    return False


# Keep the old name for backward compatibility within the module
_are_bonded = _are_bonded_or_angled


def _check_clashes(atoms: List[_Atom]) -> dict:
    """
    Detect steric clashes between non-bonded heavy atoms.

    Two atoms clash when their distance < (rvdw_A + rvdw_B - 0.4) Å.
    Hydrogens are included if present.

    Clash score = clashes per 1000 atoms (MolProbity convention).
    """
    # Sort atoms by chain/residue so nearby atoms are adjacent in list
    heavy = [a for a in atoms if a.element_upper not in ("H", "D", "")]
    # Include H atoms if present (they can clash too in high-resolution structures)
    all_atoms_for_check = atoms  # include H

    n_atoms = len(atoms)
    if n_atoms < 2:
        return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}

    clashes: List[dict] = []
    # We need an efficient spatial search; for simplicity use O(n²) with
    # an early cutoff at 5 Å for the bounding box.

    # Build a simple grid for efficiency
    coords = np.array([[a.x, a.y, a.z] for a in all_atoms_for_check])
    radii  = np.array([_vdw_radius(a) for a in all_atoms_for_check])
    max_clash_dist = float(np.max(radii) * 2)

    n = len(all_atoms_for_check)
    seen: Set[Tuple[int, int]] = set()

    for i in range(n):
        ai = all_atoms_for_check[i]
        xi = coords[i]

        for j in range(i + 1, n):
            aj = all_atoms_for_check[j]

            # Skip bonded and 1-3 angle partners
            if _are_bonded(ai, aj):
                continue

            xj = coords[j]
            diff = xi - xj
            # Quick rejection: sum of max radii
            sq_dist = float(np.dot(diff, diff))
            clash_limit = radii[i] + radii[j] - 0.4
            if sq_dist >= clash_limit * clash_limit:
                continue

            dist = math.sqrt(sq_dist)
            overlap = (radii[i] + radii[j]) - dist
            if overlap >= 0.4:
                pair_key = (min(ai.serial, aj.serial), max(ai.serial, aj.serial))
                if pair_key in seen:
                    continue
                seen.add(pair_key)
                clashes.append({
                    "atom1": f"{ai.chain}:{ai.resname}{ai.resseq}.{ai.name}",
                    "atom2": f"{aj.chain}:{aj.resname}{aj.resseq}.{aj.name}",
                    "distance":  round(dist, 3),
                    "overlap":   round(overlap, 3),
                    "vdw_sum":   round(radii[i] + radii[j], 3),
                })

    clashes.sort(key=lambda c: -c["overlap"])
    n_clashes   = len(clashes)
    clash_score = 1000.0 * n_clashes / n_atoms if n_atoms > 0 else 0.0

    return {
        "n_clashes":    n_clashes,
        "clash_score":  round(clash_score, 2),
        "worst_clashes": clashes[:20],   # top 20 by overlap magnitude
    }


# ─────────────────────────────────────────────────────────────────────────────
# Overall quality score
# ─────────────────────────────────────────────────────────────────────────────

def _compute_overall_score(
    chirality: dict,
    bond_geo: dict,
    rama: dict,
    planarity: dict,
    clashes: dict,
) -> float:
    """
    Compute a composite 0–100 quality score.

    Component weights (sum to 100):
      Ramachandran favored       30 pts
      Chirality correctness      25 pts
      Peptide planarity          20 pts
      Clash score                15 pts
      Bond geometry              10 pts
    """
    score = 0.0

    # Ramachandran (30 pts) — 100% favored → 30 pts; 0% → 0 pts
    pct_fav = rama.get("pct_favored", 0.0)
    # Interpolate: <50% → 0 pts; 98% → ~30 pts
    rama_score = max(0.0, min(30.0, (pct_fav - 50.0) / 50.0 * 30.0))
    score += rama_score

    # Chirality (25 pts)
    pct_chir = chirality.get("pct_correct", 100.0)
    score += 25.0 * pct_chir / 100.0

    # Planarity (20 pts) — 100% within 6° → 20 pts
    pct_plan = planarity.get("pct_within_6deg", 100.0)
    score += 20.0 * pct_plan / 100.0

    # Clashes (15 pts) — clash_score 0 → 15 pts; ≥20 → 0 pts
    cs = clashes.get("clash_score", 0.0)
    clash_pts = max(0.0, 15.0 * (1.0 - cs / 20.0))
    score += clash_pts

    # Bond geometry (10 pts) — bl_rmsd < 0.01 Å → 10 pts; ≥0.05 → 0 pts
    bl_rmsd = bond_geo.get("bl_rmsd", 0.0)
    geo_pts = max(0.0, 10.0 * (1.0 - bl_rmsd / 0.05))
    score += geo_pts

    return round(min(100.0, max(0.0, score)), 1)


# ─────────────────────────────────────────────────────────────────────────────
# Main public function
# ─────────────────────────────────────────────────────────────────────────────

def audit_pdb(pdb_path: str) -> dict:
    """
    Run a comprehensive quality audit on a PDB structure file.

    Validates Cα chirality, backbone bond geometry, Ramachandran
    angles, peptide planarity, and steric clashes.  Works for pure
    L-proteins, pure D-peptides, and mixed L/D structures.

    Args:
        pdb_path: Path to a PDB file (ATOM and/or HETATM records).

    Returns:
        A dict with the following top-level keys:

        ``n_residues`` (int)
            Total number of residue groups found (excluding water).

        ``n_atoms`` (int)
            Total number of atoms parsed.

        ``chains`` (list[str])
            Unique chain IDs present in the file.

        ``chirality`` (dict)
            ``n_correct``, ``n_wrong``, ``n_glycine``, ``pct_correct``,
            ``violations`` list.

        ``bond_geometry`` (dict)
            ``bl_rmsd`` (Å), ``ba_rmsd`` (°), ``outliers`` list.

        ``ramachandran`` (dict)
            ``pct_favored``, ``pct_allowed``, ``pct_outlier``,
            ``outliers`` list.

        ``planarity`` (dict)
            ``pct_within_6deg``, ``mean_deviation``, ``outliers`` list.

        ``clashes`` (dict)
            ``n_clashes``, ``clash_score`` (per 1000 atoms),
            ``worst_clashes`` list.

        ``overall_score`` (float)
            Composite 0–100 quality score.

    Raises:
        FileNotFoundError: If *pdb_path* does not exist.
        ValueError: If the file contains no parseable atom records.

    Examples
    --------
    >>> report = audit_pdb("results/1UBQ.pdb")
    >>> print(f"Score: {report['overall_score']}")
    >>> print(f"Ramachandran favored: {report['ramachandran']['pct_favored']}%")
    >>> for v in report['chirality']['violations']:
    ...     print(v)
    """
    import os
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    atoms = _parse_pdb(pdb_path)
    if not atoms:
        raise ValueError(f"No parseable ATOM/HETATM records found in: {pdb_path}")

    residue_groups = _group_by_residue(atoms)
    ordered_keys   = list(residue_groups.keys())

    chains = sorted({key[0] for key in ordered_keys})
    n_residues = len(ordered_keys)
    n_atoms    = len(atoms)

    # Run all checks
    chirality   = _check_chirality(residue_groups)
    bond_geo    = _check_bond_geometry(residue_groups, ordered_keys)
    rama        = _check_ramachandran(residue_groups, ordered_keys)
    planarity   = _check_planarity(residue_groups, ordered_keys)
    clashes     = _check_clashes(atoms)
    overall     = _compute_overall_score(chirality, bond_geo, rama, planarity, clashes)

    return {
        "n_residues":     n_residues,
        "n_atoms":        n_atoms,
        "chains":         chains,
        "chirality":      chirality,
        "bond_geometry":  bond_geo,
        "ramachandran":   rama,
        "planarity":      planarity,
        "clashes":        clashes,
        "overall_score":  overall,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: pretty-print a report
# ─────────────────────────────────────────────────────────────────────────────

def format_report(report: dict) -> str:
    """
    Return a human-readable summary of an audit_pdb() report.

    Args:
        report: dict returned by audit_pdb().

    Returns:
        Multi-line string suitable for printing.
    """
    lines = [
        "═" * 60,
        "ChiralFold PDB Auditor — Quality Report",
        "═" * 60,
        f"  Residues : {report['n_residues']}",
        f"  Atoms    : {report['n_atoms']}",
        f"  Chains   : {', '.join(report['chains'])}",
        f"  Score    : {report['overall_score']} / 100",
        "",
        "── Cα Chirality ──────────────────────────────────",
        (f"  Correct  : {report['chirality']['n_correct']}"
         f"  Wrong : {report['chirality']['n_wrong']}"
         f"  Gly : {report['chirality']['n_glycine']}"
         f"  ({report['chirality']['pct_correct']}%)"),
    ]
    for v in report["chirality"]["violations"][:5]:
        lines.append(
            f"  VIOLATION: {v['chain']}:{v['resname']}{v['resseq']}"
            f" expected {v['expected']}, observed {v['observed']}"
        )

    lines += [
        "",
        "── Bond Geometry ─────────────────────────────────",
        f"  BL RMSD  : {report['bond_geometry']['bl_rmsd']:.4f} Å",
        f"  BA RMSD  : {report['bond_geometry']['ba_rmsd']:.2f}°",
        f"  Outliers : {len(report['bond_geometry']['outliers'])}",
        "",
        "── Ramachandran ──────────────────────────────────",
        (f"  Favored  : {report['ramachandran']['pct_favored']:.1f}%"
         f"  Allowed : {report['ramachandran']['pct_allowed']:.1f}%"
         f"  Outlier : {report['ramachandran']['pct_outlier']:.1f}%"),
        "",
        "── Peptide Planarity ─────────────────────────────",
        (f"  Within 6°: {report['planarity']['pct_within_6deg']:.1f}%"
         f"  Mean dev : {report['planarity']['mean_deviation']:.2f}°"
         f"  Outliers : {len(report['planarity']['outliers'])}"),
        "",
        "── Clash Detection ───────────────────────────────",
        (f"  Clashes  : {report['clashes']['n_clashes']}"
         f"  Score : {report['clashes']['clash_score']:.1f} / 1000 atoms"),
        "═" * 60,
    ]
    return "\n".join(lines)
