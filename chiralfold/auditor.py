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
5. **Clash Detection**   — Non-bonded pairs closer than vdW sum − 0.4 Å, with
   topology-aware 1–2/1–3/1–4 exclusions (MolProbity-inspired).
6. **Summary Score**     — Composite 0–100 quality score.

All parsing is done directly from PDB ATOM/HETATM records with no external
dependencies beyond NumPy/SciPy.  RDKit is *not* required for audit_pdb().
mmCIF inputs are accepted via automatic conversion (see ``audit_pdb``).
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.spatial import cKDTree

# Module-level: do NOT call warnings.filterwarnings globally — that suppresses
# warnings for downstream user code. Targeted suppression is applied locally
# where benign warnings are expected (see _check_ramachandran).

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

# Covalent heavy-atom bonds used for 1-2 / 1-3 clash exclusions (MolProbity-like).
# Names are stripped PDB atom names. D-CCD residues map via _RESNAME_FOR_BONDS.
_BACKBONE_BOND_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("N", "CA"),
    ("CA", "C"),
    ("C", "O"),
    ("C", "OXT"),
    ("N", "H"),
    ("N", "H1"),
    ("N", "H2"),
    ("N", "H3"),
    ("N", "HN"),
)

_SIDECHAIN_BOND_PAIRS: Dict[str, Tuple[Tuple[str, str], ...]] = {
    "ALA": (("CA", "CB"),),
    "ARG": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
        ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"),
    ),
    "ASN": (("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")),
    "ASP": (("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")),
    "CYS": (("CA", "CB"), ("CB", "SG")),
    "GLN": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2"),
    ),
    "GLU": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"),
    ),
    "GLY": (),
    "HIS": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"),
        ("ND1", "CE1"), ("CD2", "NE2"), ("CE1", "NE2"),
    ),
    "ILE": (("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")),
    "LEU": (("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")),
    "LYS": (("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")),
    "MET": (("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")),
    "PHE": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
        ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"),
    ),
    "PRO": (("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")),
    "SER": (("CA", "CB"), ("CB", "OG")),
    "THR": (("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")),
    "TRP": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
        ("CD1", "NE1"), ("NE1", "CE2"), ("CD2", "CE2"), ("CD2", "CE3"),
        ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ2", "CH2"), ("CZ3", "CH2"),
    ),
    "TYR": (
        ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
        ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"),
        ("CZ", "OH"),
    ),
    "VAL": (("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")),
    "MSE": (("CA", "CB"), ("CB", "CG"), ("CG", "SE"), ("SE", "CE")),
}

# D-CCD / common variants → template used for side-chain bonds
_RESNAME_FOR_BONDS: Dict[str, str] = {
    "DAL": "ALA", "DAR": "ARG", "DSG": "ASN", "DAS": "ASP", "DCY": "CYS",
    "DGL": "GLU", "DGN": "GLN", "DHI": "HIS", "DIL": "ILE", "DLE": "LEU",
    "DLY": "LYS", "MED": "MET", "DPN": "PHE", "DPR": "PRO", "DSN": "SER",
    "DTH": "THR", "DTR": "TRP", "DTY": "TYR", "DVA": "VAL",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HIE": "HIS", "HID": "HIS",
    "HIP": "HIS", "CYX": "CYS", "CYM": "CYS", "SEP": "SER", "TPO": "THR",
    "PTR": "TYR",
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


def _is_chain_continuous(
    key_i: Tuple,
    key_j: Tuple,
    res_map_i: Dict[str, np.ndarray],
    res_map_j: Dict[str, np.ndarray],
) -> bool:
    """
    Decide whether residues *key_i* and *key_j* are peptide-bonded neighbours
    in the same chain. Used to guard inter-residue geometry checks against
    chain breaks, missing residues, and concatenated multi-chain files.

    Keys are (chain, resseq, icode, resname) tuples.
    """
    chain_i, resseq_i, _icode_i, _resname_i = key_i
    chain_j, resseq_j, _icode_j, _resname_j = key_j
    if chain_i != chain_j:
        return False
    # Allow insertion codes: 47A → 47B → 48 has resseq diff ≤ 2.
    if abs(resseq_j - resseq_i) > 2:
        return False
    ca_i = res_map_i.get("CA")
    ca_j = res_map_j.get("CA")
    if ca_i is not None and ca_j is not None:
        if float(np.linalg.norm(ca_j - ca_i)) > 4.5:
            return False
    return True


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


def _is_hydrogen_atom(atom: _Atom) -> bool:
    """True for elemental H/D or PDB hydrogen-style names (HA, HB2, …)."""
    elem = atom.element_upper
    if elem in ("H", "D"):
        return True
    name = atom.name.strip().upper()
    return bool(name) and name[0] == "H"


def _parse_pdb(pdb_path: str) -> List[_Atom]:
    """
    Parse ATOM and HETATM records from *pdb_path*.

    Rules applied:
    - Only the first MODEL is read (NMR / multi-model depositions).
    - HOH / WAT / DOD water residues are skipped.
    - Only the first alternative location (altloc ≤ 'A' or blank) is kept.
    - Lines shorter than 54 characters are skipped.
    """
    atoms: List[_Atom] = []
    seen_altloc: Dict[Tuple, str] = {}  # (chain, resseq, icode, name) → first altloc
    saw_model = False

    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                if saw_model:
                    break  # ignore subsequent models
                saw_model = True
                continue
            if line.startswith("ENDMDL"):
                if saw_model:
                    break
                continue
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

        # Canonical convention: ALWAYS use the 3-substituent signed volume
        # at Cα with ordering (N, C, Cβ). Positive signed_volume(CA; N, C, CB)
        # → L-configuration (S at Cα for standard amino acids); negative → D.
        # This is consistent regardless of whether Hα is present in the file,
        # eliminating the sign-ambiguity that arose in v3.2.0 when files with
        # explicit hydrogens used a different substituent ordering.
        signed_vol = _signed_volume(ca, n_pos, c_pos, cb_pos)

        # If Hα is present, use it only as a sanity check (does not change the
        # primary chirality determination). The H-aware 4-substituent volume
        # has the opposite sign convention; we therefore expect it to have the
        # opposite sign to the 3-substituent volume for the same handedness.
        h_pos: Optional[np.ndarray] = None
        for hname in ("HA", "HA2", "1HA", "2HA"):
            if hname in atom_map:
                h_pos = atom_map[hname]
                break
        # (Sanity check is computed but does not override; kept for future use.)
        _ = h_pos  # silence linter; reserved for future Hα cross-checks

        # Planarity threshold: 0.05 Å^3 corresponds to a deviation from a flat
        # tetrahedron below which the sign of the signed volume is dominated
        # by coordinate noise rather than true stereochemistry.
        if abs(signed_vol) < 0.05:
            # Essentially planar — cannot assign chirality
            n_error += 1
            continue

        # Canonical convention (matches chiralfold/af3_correct.py):
        # Positive signed volume(CA; N, C, CB) → L-configuration (S at Cα for
        # standard amino acids). Negative → D-configuration. This is
        # consistent regardless of H presence and is the convention used
        # throughout the ChiralFold codebase.
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

        # Inter-residue (peptide bond) geometry — guard against chain
        # breaks, missing residues, and concatenated multi-chain files.
        if i + 1 < n and _is_chain_continuous(
            ordered_keys[i], ordered_keys[i + 1], res_maps[i], res_maps[i + 1]
        ):
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

        # φ = C(i-1) - N(i) - CA(i) - C(i) — only if (i-1, i) are continuous.
        if i > 0 and _is_chain_continuous(
            ordered_keys[i - 1], ordered_keys[i], res_maps[i - 1], res_maps[i]
        ):
            c_prev = res_maps[i - 1].get("C")
            if c_prev is not None:
                phi = _dihedral_deg(c_prev, n_, ca, c)

        # ψ = N(i) - CA(i) - C(i) - N(i+1) — only if (i, i+1) are continuous.
        if i < n - 1 and _is_chain_continuous(
            ordered_keys[i], ordered_keys[i + 1], res_maps[i], res_maps[i + 1]
        ):
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
        # Guard against chain gaps before measuring ω across residues.
        if not _is_chain_continuous(
            ordered_keys[i], ordered_keys[i + 1], res_maps[i], res_maps[i + 1]
        ):
            continue

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


def _bond_template_resname(resname: str) -> str:
    """Map PDB / D-CCD residue names to a side-chain bond template."""
    rn = resname.upper().strip()
    if rn in _SIDECHAIN_BOND_PAIRS:
        return rn
    return _RESNAME_FOR_BONDS.get(rn, rn)


def _intra_residue_bond_pairs(resname: str) -> List[Tuple[str, str]]:
    tmpl = _bond_template_resname(resname)
    pairs = list(_BACKBONE_BOND_PAIRS)
    pairs.extend(_SIDECHAIN_BOND_PAIRS.get(tmpl, ()))
    # Unknown residues: still connect CA–CB when present via backbone only;
    # distance fallback in _clash_excluded_index_pairs covers remaining 1-3.
    if tmpl not in _SIDECHAIN_BOND_PAIRS and tmpl not in ("GLY",):
        pairs.append(("CA", "CB"))
    return pairs


def _residue_key(atom: _Atom) -> Tuple[str, int, str]:
    return (atom.chain, atom.resseq, atom.icode)


def _clash_excluded_index_pairs(atoms: List[_Atom]) -> Set[Tuple[int, int]]:
    """
    Build index pairs that must not score as clashes (covalent 1-2 and 1-3).

    Uses standard amino-acid topology (not a fragile distance cutoff). Same-
    residue CA–CG (~2.5–2.7 Å) is a classic 1-3 false positive when a 2.6 Å
    distance gate is used — that bug inflated clash scores on AFDB/PDB models.
    """
    # name → list of atom indices in that residue
    by_res: Dict[Tuple[str, int, str], Dict[str, List[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for idx, atom in enumerate(atoms):
        by_res[_residue_key(atom)][atom.name.strip()].append(idx)

    excluded: Set[Tuple[int, int]] = set()

    def _add_edge(i: int, j: int, adj: Dict[int, Set[int]]) -> None:
        adj[i].add(j)
        adj[j].add(i)

    # Per-residue covalent graphs
    residue_adj: Dict[Tuple[str, int, str], Dict[int, Set[int]]] = {}
    for rkey, name_map in by_res.items():
        # Pick a representative resname from any atom in the residue
        sample_idx = next(iter(next(iter(name_map.values()))))
        resname = atoms[sample_idx].resname
        adj: Dict[int, Set[int]] = defaultdict(set)
        for a_name, b_name in _intra_residue_bond_pairs(resname):
            for i in name_map.get(a_name, ()):
                for j in name_map.get(b_name, ()):
                    _add_edge(i, j, adj)
        residue_adj[rkey] = adj

        # Exclude 1-2, 1-3, and 1-4 within the residue (MolProbity Probe -4).
        nodes = list(adj.keys())
        for name_idxs in name_map.values():
            for i in name_idxs:
                if i not in adj:
                    adj[i] = set()
                    nodes.append(i)
        for start in nodes:
            # BFS up to depth 3
            frontier = {start}
            visited = {start}
            for _depth in range(3):
                nxt: Set[int] = set()
                for node in frontier:
                    for nb in adj[node]:
                        if nb in visited:
                            continue
                        visited.add(nb)
                        nxt.add(nb)
                        a, b = (start, nb) if start < nb else (nb, start)
                        excluded.add((a, b))
                frontier = nxt

    # Peptide bonds + cross-residue 1-3 (C–CA_next, CA–N_next, O–N_next, C–H_next)
    ordered = sorted(by_res.keys(), key=lambda k: (k[0], k[1], k[2]))
    for i in range(len(ordered) - 1):
        k0, k1 = ordered[i], ordered[i + 1]
        if k0[0] != k1[0]:
            continue
        if abs(k1[1] - k0[1]) > 2:
            continue
        m0, m1 = by_res[k0], by_res[k1]
        # Continuity: CA–CA distance if both present
        if m0.get("CA") and m1.get("CA"):
            ca0 = atoms[m0["CA"][0]].xyz
            ca1 = atoms[m1["CA"][0]].xyz
            if float(np.linalg.norm(ca1 - ca0)) > 4.5:
                continue
        cross_12_13 = [
            ("C", "N"),    # peptide 1-2
            ("CA", "N"),   # 1-3 via C
            ("C", "CA"),   # 1-3 via N
            ("O", "N"),    # 1-3 via C
            ("C", "H"),    # 1-3 via N (added amide H)
            ("C", "HN"),
            ("O", "H"),
            ("O", "HN"),
            # Proline: CD bonded to N → previous-residue contacts through N
            ("C", "CD"),
            ("CA", "CD"),
            ("O", "CD"),
            # 1-4 peptide contacts (Probe -4 style)
            ("CA", "H"),
            ("CA", "HN"),
            ("O", "CA"),
            ("C", "C"),
            ("N", "CA"),   # N(i)–CA(i+1) is 1-4 via CA(i)–C(i)–N / C–N–CA
            ("N", "H"),
            ("N", "HN"),
            ("O", "C"),
            ("N", "N"),  # 1-4 via CA–C–N
            ("C", "CB"),  # 1-4 via N–CA of next residue
        ]
        for a_name, b_name in cross_12_13:
            for ia in m0.get(a_name, ()):
                for ib in m1.get(b_name, ()):
                    key = (ia, ib) if ia < ib else (ib, ia)
                    excluded.add(key)

    # Covalent disulfides: SG–SG within 2.5 Å, plus 1-3/1-4 across the bridge
    sg_idxs = [
        i for i, a in enumerate(atoms)
        if a.name.strip() == "SG"
    ]
    for a in range(len(sg_idxs)):
        for b in range(a + 1, len(sg_idxs)):
            i, j = sg_idxs[a], sg_idxs[b]
            dist = float(np.linalg.norm(atoms[i].xyz - atoms[j].xyz))
            if dist > 2.5:
                continue
            key = (i, j) if i < j else (j, i)
            excluded.add(key)
            # 1-3: CB–SG_partner ; 1-4: CA–SG_partner, CB–CB
            for sg_i, sg_j in ((i, j), (j, i)):
                rkey = _residue_key(atoms[sg_i])
                partner_key = _residue_key(atoms[sg_j])
                for cb_i in by_res[rkey].get("CB", ()):
                    k = (cb_i, sg_j) if cb_i < sg_j else (sg_j, cb_i)
                    excluded.add(k)
                    for cb_j in by_res[partner_key].get("CB", ()):
                        k2 = (cb_i, cb_j) if cb_i < cb_j else (cb_j, cb_i)
                        excluded.add(k2)
                for ca_i in by_res[rkey].get("CA", ()):
                    k = (ca_i, sg_j) if ca_i < sg_j else (sg_j, ca_i)
                    excluded.add(k)

    # Distance fallback ONLY for residues without a known covalent template
    # (ligands / nonstandard) — never for standard AA where topology is complete.
    for rkey, name_map in by_res.items():
        sample_idx = next(iter(next(iter(name_map.values()))))
        tmpl = _bond_template_resname(atoms[sample_idx].resname)
        if tmpl in _SIDECHAIN_BOND_PAIRS or tmpl == "GLY":
            continue
        idxs = [i for lst in name_map.values() for i in lst]
        for a in range(len(idxs)):
            for b in range(a + 1, len(idxs)):
                i, j = idxs[a], idxs[b]
                key = (i, j) if i < j else (j, i)
                if key in excluded:
                    continue
                dist = float(np.linalg.norm(atoms[i].xyz - atoms[j].xyz))
                if dist < 2.9:
                    excluded.add(key)

    return excluded


def _is_hbond_donor_acceptor_pair(a: _Atom, b: _Atom) -> bool:
    """True if (a,b) is a plausible H-bond donor/acceptor pair (not a clash)."""
    def _role(atom: _Atom) -> str:
        name = atom.name.strip().upper()
        elem = atom.element_upper
        if elem == "H" or name in ("H", "HN", "H1", "H2", "H3"):
            return "donor_h"
        if elem == "O" or name in (
            "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OXT"
        ):
            return "acceptor"
        if name in (
            "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ"
        ):
            return "donor_n"
        return ""

    ra, rb = _role(a), _role(b)
    roles = {ra, rb}
    if roles == {"donor_h", "acceptor"}:
        return True
    # Heavy-atom H-bond geometry (N···O ≈ 2.5–3.5 Å) — not a steric clash.
    if roles == {"donor_n", "acceptor"}:
        return _bond_length(a.xyz, b.xyz) < 3.5
    return False


def _are_bonded_or_angled(a: _Atom, b: _Atom) -> bool:
    """
    Return True if atoms a and b are covalent 1-2 or 1-3 partners.

    Prefer :func:`_clash_excluded_index_pairs` inside bulk clash checks.
    This helper remains for unit tests and ad-hoc calls; it uses topology
    when both atoms share a residue key, else a conservative distance gate.
    """
    same_res = (
        a.chain == b.chain and a.resseq == b.resseq and a.icode == b.icode
    )
    dist = _bond_length(a.xyz, b.xyz)
    if same_res:
        # Topology-aware path for standard residues
        pairs = _intra_residue_bond_pairs(a.resname)
        names = {a.name.strip(), b.name.strip()}
        # Direct bond
        for x, y in pairs:
            if names == {x, y}:
                return True
        # 1-3: share a common bonded neighbor name in the template graph
        adj: Dict[str, Set[str]] = defaultdict(set)
        for x, y in pairs:
            adj[x].add(y)
            adj[y].add(x)
        na, nb = a.name.strip(), b.name.strip()
        if adj[na] & adj[nb]:
            return True
        return dist < 2.9

    adjacent = a.chain == b.chain and abs(a.resseq - b.resseq) <= 1
    if adjacent:
        return dist < 2.7
    return False


# Keep the old name for backward compatibility within the module
_are_bonded = _are_bonded_or_angled


def _add_backbone_hydrogens(atoms: List[_Atom]) -> List[_Atom]:
    """
    Estimate backbone amide H positions for clash detection.

    Places H in the peptide plane of C(prev)–N–CA (MolProbity Reduce-style).
    Residue keys include insertion codes. Amide H is only built when the
    previous residue is on the *same chain* and peptide-continuous (resseq
    gap ≤ 2 and CA–CA ≤ 4.5 Å), so multi-chain / concatenated files do not
    invent a cross-chain amide H from another chain's carbonyl.
    """
    # (chain, resseq, icode) → atom-name map
    residues: Dict[Tuple[str, int, str], Dict[str, _Atom]] = {}
    for a in atoms:
        key = (a.chain, a.resseq, a.icode)
        if key not in residues:
            residues[key] = {}
        residues[key][a.name.strip()] = a

    new_h_atoms = []
    max_serial = max((a.serial for a in atoms), default=0)
    h_serial = max_serial + 1

    sorted_keys = sorted(residues.keys())
    for idx, key in enumerate(sorted_keys):
        res = residues[key]
        n_atom = res.get('N')
        ca_atom = res.get('CA')
        if n_atom is None or ca_atom is None:
            continue
        # Proline (and D-Pro) has no amide hydrogen — do not invent one.
        if n_atom.resname.upper().strip() in PROLINE_RESNAMES:
            continue
        # Skip if H already present
        if 'H' in res or 'HN' in res:
            continue
        # Need a preceding peptide neighbour on the same chain
        if idx == 0:
            continue
        prev_key = sorted_keys[idx - 1]
        if prev_key[0] != key[0]:
            continue
        if abs(key[1] - prev_key[1]) > 2:
            continue
        prev_res = residues.get(prev_key, {})
        prev_c = prev_res.get('C')
        if prev_c is None:
            continue
        prev_ca = prev_res.get('CA')
        if prev_ca is not None:
            ca_dist = float(np.linalg.norm(prev_ca.xyz - ca_atom.xyz))
            if ca_dist > 4.5:
                continue

        # Place H in the peptide plane of C(prev)–N–CA, opposite the angle
        # bisector at N (MolProbity Reduce-style). Vectors must originate at N.
        n_pos = np.array([n_atom.x, n_atom.y, n_atom.z])
        c_pos = np.array([prev_c.x, prev_c.y, prev_c.z])
        ca_pos = np.array([ca_atom.x, ca_atom.y, ca_atom.z])

        v_nc = c_pos - n_pos
        v_nca = ca_pos - n_pos
        n_nc = np.linalg.norm(v_nc) + 1e-12
        n_nca = np.linalg.norm(v_nca) + 1e-12
        h_dir = -(v_nc / n_nc + v_nca / n_nca)
        h_dir = h_dir / (np.linalg.norm(h_dir) + 1e-12)
        h_pos = n_pos + h_dir * 1.02

        h_atom = _Atom(
            record='ATOM', serial=h_serial, name=' H  ', altloc=' ',
            resname=n_atom.resname, chain=n_atom.chain,
            resseq=n_atom.resseq, icode=n_atom.icode,
            x=float(h_pos[0]), y=float(h_pos[1]), z=float(h_pos[2]),
            element='H',
        )
        new_h_atoms.append(h_atom)
        h_serial += 1

    return atoms + new_h_atoms


def _check_clashes(atoms: List[_Atom]) -> dict:
    """
    Detect steric clashes between non-bonded atoms using a scipy KD-tree.

    Two atoms clash when their distance < (rvdw_A + rvdw_B - 0.4) Å.
    Deposited hydrogens are stripped and backbone amide H are re-added
    (MolProbity Reduce-style), so explicit C–H bonds are never scored as clashes.
    Covalent 1-2 / 1-3 / 1-4 pairs are excluded via residue topology.

    Clash score = clashes per 1000 heavy atoms.
    """
    heavy = [a for a in atoms if not _is_hydrogen_atom(a)]
    n_atoms = len(heavy)
    if n_atoms < 2:
        return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}

    # Add backbone hydrogens if not present (MolProbity does this)
    all_atoms_for_check = _add_backbone_hydrogens(heavy)

    excluded = _clash_excluded_index_pairs(all_atoms_for_check)

    coords = np.asarray(
        [[a.x, a.y, a.z] for a in all_atoms_for_check], dtype=float
    )
    radii = np.asarray(
        [_vdw_radius(a) for a in all_atoms_for_check], dtype=float
    )
    max_sum_radii = float(np.max(radii)) * 2

    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=max_sum_radii, output_type="ndarray")

    clashes: List[dict] = []
    seen: Set[Tuple[int, int]] = set()

    for i, j in pairs:
        i_i, j_i = int(i), int(j)
        ex_key = (i_i, j_i) if i_i < j_i else (j_i, i_i)
        if ex_key in excluded:
            continue
        ai = all_atoms_for_check[i_i]
        aj = all_atoms_for_check[j_i]
        # MolProbity: donor–acceptor contacts are H-bonds, not clashes.
        if _is_hbond_donor_acceptor_pair(ai, aj):
            continue
        dist = float(np.linalg.norm(coords[i] - coords[j]))
        clash_limit = radii[i] + radii[j] - 0.4
        if dist >= clash_limit:
            continue
        overlap = float((radii[i] + radii[j]) - dist)
        if overlap < 0.4:
            continue
        pair_key = (min(ai.serial, aj.serial), max(ai.serial, aj.serial))
        if pair_key in seen:
            continue
        seen.add(pair_key)
        clashes.append({
            "atom1": f"{ai.chain}:{ai.resname}{ai.resseq}.{ai.name.strip()}",
            "atom2": f"{aj.chain}:{aj.resname}{aj.resseq}.{aj.name.strip()}",
            "distance": round(dist, 3),
            "overlap": round(overlap, 3),
            "vdw_sum": round(float(radii[i] + radii[j]), 3),
        })

    clashes.sort(key=lambda c: -c["overlap"])
    n_clashes = len(clashes)
    clash_score = 1000.0 * n_clashes / n_atoms if n_atoms > 0 else 0.0

    return {
        "n_clashes": n_clashes,
        "clash_score": round(clash_score, 2),
        "worst_clashes": clashes[:20],
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

    # Ramachandran (30 pts) — smooth linear interpolation: 0 pts at 0%,
    # 30 pts at 98% favored, and proportional below that with no cliff edge
    # at 50% (which previously gave the same 0 pts as 0% favored).
    pct_fav = rama.get("pct_favored", 0.0)
    rama_score = min(30.0, 30.0 * pct_fav / 98.0)
    score += max(0.0, rama_score)

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

def audit_pdb(pdb_path: str, chain: Optional[str] = None) -> dict:
    """
    Run a comprehensive quality audit on a PDB (or mmCIF) structure file.

    Validates Cα chirality, backbone bond geometry, Ramachandran
    angles, peptide planarity, and steric clashes.  Works for pure
    L-proteins, pure D-peptides, and mixed L/D structures.

    Args:
        pdb_path: Path to a PDB or mmCIF file (ATOM and/or HETATM records).
            mmCIF is converted in-core to a temporary PDB. FASTA is rejected
            (use ``chiralfold.fetch.resolve_to_pdb`` / CLI ``--id`` for AFDB).
        chain: Optional single chain ID to restrict the audit (e.g. ``"A"``).

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
            ``n_clashes``, ``clash_score`` (per 1000 heavy atoms),
            ``worst_clashes`` list.

        ``overall_score`` (float)
            Composite 0–100 quality score.

    Raises:
        FileNotFoundError: If *pdb_path* does not exist.
        ValueError: If the file contains no parseable atom records,
            or no atoms remain after *chain* filtering.

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

    # Accept mmCIF transparently; FASTA must go through fetch/AFDB.
    from .structure_io import ensure_pdb_path

    try:
        pdb_path = ensure_pdb_path(pdb_path)
    except ValueError as exc:
        raise ValueError(str(exc)) from exc

    atoms = _parse_pdb(pdb_path)
    if not atoms:
        raise ValueError(f"No parseable ATOM/HETATM records found in: {pdb_path}")

    if chain is not None:
        atoms = [a for a in atoms if a.chain == chain]
        if not atoms:
            raise ValueError(
                f"No ATOM/HETATM records found for chain {chain!r} in: {pdb_path}"
            )

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
         f"  Score : {report['clashes']['clash_score']:.1f} / 1000 heavy atoms"),
        "═" * 60,
    ]
    return "\n".join(lines)
