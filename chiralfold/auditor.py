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
   topology-aware 1–2/1–3/1–4 exclusions, Probe-style HN/HA hydrogens, and
   geometry-gated H-bond suppression (MolProbity-inspired).
6. **Summary Score**     — Composite 0–100 quality score.

All parsing is done directly from PDB ATOM/HETATM records with no external
dependencies beyond NumPy/SciPy.  RDKit is *not* required for audit_pdb().
mmCIF inputs are accepted via automatic conversion (see ``audit_pdb``).
"""

from __future__ import annotations

import math
from collections import defaultdict
from functools import lru_cache
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.spatial import cKDTree

from ._pdbio import read_atom_records

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
    # Probe uses ~1.17 Å for H (slightly tighter than Bondi 1.20)
    "H": 1.17,
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
    # Cα hydrogens (Reduce / Probe style)
    ("CA", "HA"),
    ("CA", "HA2"),
    ("CA", "HA3"),
)

# Ideal bond lengths for placed hydrogens (Å)
_BOND_LEN_NH = 1.02
_BOND_LEN_CH = 1.09

# Polar (H-bond-capable) hydrogen names. HA/HA2/HA3 are carbon-bound and excluded.
_POLAR_H_NAMES: Set[str] = {
    "H", "HN", "H1", "H2", "H3",
    "HG", "HG1",
    "HD1", "HD2",
    "HE", "HE1", "HE2",
    "HH", "HH11", "HH12", "HH21", "HH22",
    "HZ", "HZ1", "HZ2", "HZ3",
}

# Probe-inspired H-bond geometry gates (suppress clash only if satisfied)
_HBOND_H_ACC_MAX = 2.50   # Å  H···acceptor
_HBOND_ANGLE_MIN = 120.0  # °   donor–H–acceptor
_HBOND_NO_MIN = 2.50      # Å  N···O heavy-atom window
_HBOND_NO_MAX = 3.50

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
    """Lightweight representation of a parsed PDB atom record.

    ``sname`` (the stripped, upper-cased atom name) is computed in ``__init__``
    and the remaining derived quantities (``xyz``, ``element_upper``) are
    memoised: the auditor touches them O(10) times per atom, and rebuilding a
    3-vector or re-stripping a name on every access dominated the profile.
    The array returned by :attr:`xyz` is marked read-only so that an accidental
    in-place write fails loudly instead of silently corrupting the structure.
    """

    __slots__ = [
        "record", "serial", "name", "altloc", "resname",
        "chain", "resseq", "icode",
        "x", "y", "z", "element",
        "sname", "_xyz", "_elem", "_is_h",
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
        self.sname   = name.strip().upper()
        self._xyz    = None
        self._elem   = None
        self._is_h   = None

    @property
    def xyz(self) -> np.ndarray:
        v = self._xyz
        if v is None:
            v = np.array([self.x, self.y, self.z], dtype=float)
            v.setflags(write=False)
            self._xyz = v
        return v

    @property
    def element_upper(self) -> str:
        e = self._elem
        if e is not None:
            return e
        if self.element:
            e = self.element.upper().strip()
        else:
            # Infer from atom name (column 13-14, first non-digit char)
            e = "C"
            for ch in self.name.strip():
                if ch.isalpha():
                    e = ch.upper()
                    break
        self._elem = e
        return e


def _is_hydrogen_atom(atom: _Atom) -> bool:
    """True for elemental H/D or PDB hydrogen-style names (HA, HB2, …)."""
    cached = atom._is_h
    if cached is not None:
        return cached
    elem = atom.element_upper
    if elem in ("H", "D"):
        result = True
    else:
        name = atom.sname
        result = bool(name) and name[0] == "H"
    atom._is_h = result
    return result


def _parse_pdb(pdb_path: str) -> List[_Atom]:
    """
    Parse ATOM and HETATM records from *pdb_path*.

    Rules applied:
    - Only the first MODEL is read (NMR / multi-model depositions).
    - HOH / WAT / DOD water residues are skipped.
    - Alternate locations are resolved **per residue**: blank-altloc atoms are
      always kept, and among the lettered alternates of a residue the label
      with the highest summed occupancy wins (ties broken by label). Prior
      releases hard-coded the accept-list ``{' ', 'A', '1'}`` per *atom*, which
      (a) silently dropped every atom of a residue modelled only as altloc
      ``B``/``C`` and (b) could keep the minority conformer, or even mix atoms
      from different alternates within one residue, producing spurious
      intra-residue geometry and clash signals.
    - Lines shorter than 54 characters are skipped.
    """
    return [
        _Atom(
            record=r.record, serial=r.serial, name=r.name, altloc=r.altloc,
            resname=r.resname, chain=r.chain, resseq=r.resseq, icode=r.icode,
            x=r.x, y=r.y, z=r.z, element=r.element,
        )
        for r in read_atom_records(pdb_path)
    ]


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
# Batched geometry kernels
#
# The scalar helpers above are retained as the documented single-tetrahedron /
# single-dihedral API (and are used by tests and downstream callers). The audit
# itself evaluates every residue at once through the ``*_v`` kernels below:
# ``np.cross``/``np.linalg.norm`` on 3-vectors carry ~10 µs of dispatch overhead
# each, which dominated the auditor profile at ~30 numpy calls per residue.
# Each kernel reproduces its scalar counterpart's arithmetic elementwise, so
# results agree to floating-point round-off.
# ─────────────────────────────────────────────────────────────────────────────

def _dots_v(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Row-wise dot product of two (m, 3) arrays."""
    return np.einsum("ij,ij->i", u, v)


def _norms_v(u: np.ndarray) -> np.ndarray:
    """Row-wise Euclidean norm of an (m, 3) array."""
    return np.sqrt(np.einsum("ij,ij->i", u, u))


def _bond_lengths_v(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise |b − a| for (m, 3) arrays."""
    return _norms_v(b - a)


def _bond_angles_v(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Row-wise angle a–b–c in degrees (vertex at *b*); NaN when degenerate."""
    out = np.full(a.shape[0], np.nan)
    v1 = a - b
    v2 = c - b
    n1 = _norms_v(v1)
    n2 = _norms_v(v2)
    ok = (n1 > 1e-8) & (n2 > 1e-8)
    if not ok.any():
        return out
    cos_a = np.clip(_dots_v(v1[ok], v2[ok]) / (n1[ok] * n2[ok]), -1.0, 1.0)
    out[ok] = np.degrees(np.arccos(cos_a))
    return out


def _dihedrals_v(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray
) -> np.ndarray:
    """Row-wise dihedral p1–p2–p3–p4 in degrees; NaN when p2≡p3."""
    out = np.full(p1.shape[0], np.nan)
    b0 = -(p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    n1 = _norms_v(b1)
    ok = n1 > 1e-8
    if not ok.any():
        return out
    b1u = b1[ok] / n1[ok][:, None]
    b0o = b0[ok]
    b2o = b2[ok]
    v = b0o - _dots_v(b0o, b1u)[:, None] * b1u
    w = b2o - _dots_v(b2o, b1u)[:, None] * b1u
    out[ok] = np.degrees(np.arctan2(_dots_v(np.cross(b1u, v), w), _dots_v(v, w)))
    return out


def _signed_volumes_v(
    center: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray
) -> np.ndarray:
    """Row-wise signed tetrahedron volume, same convention as ``_signed_volume``."""
    return _dots_v(a - center, np.cross(b - center, c - center))


# ─────────────────────────────────────────────────────────────────────────────
# Shared per-residue coordinate table
# ─────────────────────────────────────────────────────────────────────────────

#: Backbone atoms every geometry check needs, in a fixed slot order.
_TRACKED_ATOMS: Tuple[str, ...] = ("N", "CA", "C", "O", "CB")
_SLOT_N, _SLOT_CA, _SLOT_C, _SLOT_O, _SLOT_CB = range(len(_TRACKED_ATOMS))


class _ResidueTable:
    """Column-oriented view of the tracked backbone atoms of every residue.

    Built **once** per audit and shared by the chirality, bond-geometry,
    Ramachandran and planarity checks. Previously each of those four checks
    rebuilt its own ``{atom_name: xyz}`` dict per residue — four full passes
    over the structure and four ``np.array`` allocations per atom.

    Attributes:
        pos: ``(5, n_res, 3)`` coordinates, NaN where the atom is absent.
        has: ``(5, n_res)`` presence mask, indexed by the ``_SLOT_*`` constants.
        continuous: ``(n_res - 1,)`` bool — residue *i* and *i+1* are
            peptide-bonded neighbours in the same chain (same predicate as
            :func:`_is_chain_continuous`).
    """

    __slots__ = ("keys", "labels", "rtypes", "resnames", "records",
                 "pos", "has", "n", "continuous")

    def __init__(
        self,
        residue_groups: Dict[Tuple, List[_Atom]],
        ordered_keys: List[Tuple],
    ):
        n = len(ordered_keys)
        self.keys = ordered_keys
        self.n = n
        self.pos = np.full((len(_TRACKED_ATOMS), n, 3), np.nan, dtype=float)
        self.has = np.zeros((len(_TRACKED_ATOMS), n), dtype=bool)
        slot_of = {name: k for k, name in enumerate(_TRACKED_ATOMS)}

        labels: List[str] = []
        rtypes: List[str] = []
        resnames: List[str] = []
        records: List[str] = []
        pos = self.pos
        has = self.has

        for i, key in enumerate(ordered_keys):
            chain, resseq, icode, resname = key
            labels.append(f"{chain}:{resname}{resseq}{icode.strip()}")
            rtypes.append(_get_rtype(resname))
            resnames.append(resname)
            group = residue_groups[key]
            records.append(group[0].record if group else "ATOM")
            for a in group:
                k = slot_of.get(a.sname)
                if k is not None and not has[k, i]:
                    pos[k, i, 0] = a.x
                    pos[k, i, 1] = a.y
                    pos[k, i, 2] = a.z
                    has[k, i] = True

        self.labels = labels
        self.rtypes = rtypes
        self.resnames = resnames
        self.records = records
        self.continuous = self._compute_continuity()

    def _compute_continuity(self) -> np.ndarray:
        n = self.n
        if n < 2:
            return np.zeros(max(n - 1, 0), dtype=bool)
        chains = np.array([k[0] for k in self.keys], dtype=object)
        resseq = np.array([int(k[1]) for k in self.keys], dtype=np.int64)
        same_chain = np.array(
            [chains[i] == chains[i + 1] for i in range(n - 1)], dtype=bool
        )
        near_seq = np.abs(resseq[1:] - resseq[:-1]) <= 2
        ca = self.pos[_SLOT_CA]
        has_ca = self.has[_SLOT_CA]
        both_ca = has_ca[:-1] & has_ca[1:]
        span_ok = np.ones(n - 1, dtype=bool)
        if both_ca.any():
            d = _bond_lengths_v(ca[:-1][both_ca], ca[1:][both_ca])
            span_ok[both_ca] = d <= 4.5
        return same_chain & near_seq & span_ok

    def next_pair_mask(self) -> np.ndarray:
        """Bool mask over residues *i* with a continuous successor *i+1*."""
        mask = np.zeros(self.n, dtype=bool)
        if self.n > 1:
            mask[:-1] = self.continuous
        return mask


def _as_table(source, ordered_keys: Optional[List[Tuple]] = None) -> _ResidueTable:
    """Coerce *source* to a :class:`_ResidueTable`.

    Accepts an existing table (returned unchanged) or the legacy
    ``(residue_groups[, ordered_keys])`` form used by the individual
    ``_check_*`` helpers before they were consolidated onto a shared table.
    """
    if isinstance(source, _ResidueTable):
        return source
    keys = list(source.keys()) if ordered_keys is None else list(ordered_keys)
    return _ResidueTable(source, keys)


# ─────────────────────────────────────────────────────────────────────────────
# Check 1: Cα Chirality
# ─────────────────────────────────────────────────────────────────────────────

#: Signed-volume magnitude (Å³) below which the Cα tetrahedron is treated as
#: planar and its sign as noise rather than stereochemistry.
CHIRALITY_PLANAR_EPS = 0.05


def _check_chirality(source, ordered_keys: Optional[List[Tuple]] = None) -> dict:
    """
    For each residue, compute the signed volume of the Cα tetrahedron
    (N, C, Cβ at Cα) and compare its sign to the residue label.

    Expected configuration:
    - L-amino acids (ATOM records, L_RESNAMES) → positive signed volume
      (S-configuration, left-handed triple N→C→Cβ at Cα)
    - D-amino acids (HETATM / D_RESNAMES) → negative signed volume
      (R-configuration, mirror image)
    - Glycine → achiral, skip

    Canonical convention (matches :mod:`chiralfold.af3_correct`): positive
    ``signed_volume(CA; N, C, CB)`` → L-configuration, negative → D. This holds
    regardless of whether Hα is present in the file, which removes the
    sign-ambiguity that arose in v3.2.0 for files with explicit hydrogens.

    Accepts a prebuilt :class:`_ResidueTable` (the audit path) or the legacy
    ``residue_groups`` mapping plus optional ordered key list.
    """
    table = _as_table(source, ordered_keys)
    n = table.n
    if n == 0:
        return {
            "n_correct": 0, "n_wrong": 0, "n_glycine": 0, "n_error": 0,
            "pct_correct": 100.0, "violations": [],
        }

    resnames = table.resnames
    is_gly = np.fromiter(
        (rn in GLYCINE_RESNAMES for rn in resnames), dtype=bool, count=n
    )
    complete = (
        table.has[_SLOT_CA] & table.has[_SLOT_N]
        & table.has[_SLOT_C] & table.has[_SLOT_CB]
    )
    evaluable = (~is_gly) & complete

    signed_vol = np.zeros(n, dtype=float)
    if evaluable.any():
        signed_vol[evaluable] = _signed_volumes_v(
            table.pos[_SLOT_CA][evaluable],
            table.pos[_SLOT_N][evaluable],
            table.pos[_SLOT_C][evaluable],
            table.pos[_SLOT_CB][evaluable],
        )

    # Planar tetrahedra carry no reliable sign — count with the incomplete ones.
    planar = evaluable & (np.abs(signed_vol) < CHIRALITY_PLANAR_EPS)
    scored = evaluable & ~planar
    n_glycine = int(is_gly.sum())
    n_error = int(((~is_gly) & ~complete).sum() + planar.sum())

    expected_d = np.fromiter(
        (
            True if rn in D_RESNAMES
            else False if rn in L_RESNAMES
            # Unknown residue type — fall back to the ATOM/HETATM record.
            else rec == "HETATM"
            for rn, rec in zip(resnames, table.records)
        ),
        dtype=bool, count=n,
    )
    observed_d = signed_vol < 0
    correct = scored & (observed_d == expected_d)
    wrong = scored & (observed_d != expected_d)

    violations: List[dict] = []
    for i in np.flatnonzero(wrong):
        chain, resseq, _icode, resname = table.keys[i]
        violations.append({
            "chain":    chain,
            "resseq":   resseq,
            "resname":  resname,
            "expected": "D" if expected_d[i] else "L",
            "observed": "D" if observed_d[i] else "L",
            "signed_volume": round(float(signed_vol[i]), 4),
        })

    n_correct = int(correct.sum())
    n_wrong = len(violations)
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

#: Bond-length columns measured per residue, in report order.
_BL_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("N-CA", "N-CA"),
    ("CA-C", "CA-C"),
    ("C=O", "C-O"),
    ("C-N(peptide)", "C-N"),
)

#: Bond-angle columns measured per residue, in report order.
_BA_COLUMNS: Tuple[str, ...] = ("N-CA-C", "CA-C-N", "C-N-CA")


def _check_bond_geometry(source, ordered_keys: Optional[List[Tuple]] = None) -> dict:
    """
    Measure backbone bond lengths and angles; flag deviations > 3σ.

    Four bond lengths (N–Cα, Cα–C, C=O, peptide C–N) and three angles
    (N–Cα–C, Cα–C–N, C–N–Cα) are evaluated for every residue at once. The
    inter-residue terms are gated on :attr:`_ResidueTable.continuous`, which
    guards against chain breaks, missing residues, and concatenated
    multi-chain files.

    Returns bl_rmsd, ba_rmsd, and a list of outlier records.
    """
    table = _as_table(source, ordered_keys)
    n = table.n
    if n == 0:
        return {
            "bl_rmsd": 0.0, "ba_rmsd": 0.0,
            "n_bonds_checked": 0, "n_angles_checked": 0, "outliers": [],
        }

    pos, has = table.pos, table.has
    nxt = table.next_pair_mask()

    # Shifted views of residue i+1 (row i holds residue i+1's atom).
    def _shift(arr: np.ndarray) -> np.ndarray:
        out = np.full_like(arr, np.nan)
        out[:-1] = arr[1:]
        return out

    n_next_pos = _shift(pos[_SLOT_N])
    ca_next_pos = _shift(pos[_SLOT_CA])
    has_n_next = np.zeros(n, dtype=bool)
    has_ca_next = np.zeros(n, dtype=bool)
    has_n_next[:-1] = has[_SLOT_N][1:]
    has_ca_next[:-1] = has[_SLOT_CA][1:]

    bl = np.full((n, 4), np.nan)
    bl_mask = np.zeros((n, 4), dtype=bool)
    ba = np.full((n, 3), np.nan)
    ba_mask = np.zeros((n, 3), dtype=bool)

    def _fill_bl(col: int, mask: np.ndarray, a: np.ndarray, b: np.ndarray) -> None:
        if mask.any():
            bl[mask, col] = _bond_lengths_v(a[mask], b[mask])
            bl_mask[:, col] = mask

    def _fill_ba(
        col: int, mask: np.ndarray,
        a: np.ndarray, b: np.ndarray, c: np.ndarray,
    ) -> None:
        if not mask.any():
            return
        vals = _bond_angles_v(a[mask], b[mask], c[mask])
        ba[mask, col] = vals
        # Degenerate angles come back NaN and are dropped, matching the
        # scalar implementation's ``if not math.isnan(ba)`` guard.
        keep = mask.copy()
        keep[mask] = ~np.isnan(vals)
        ba_mask[:, col] = keep

    _fill_bl(0, has[_SLOT_N] & has[_SLOT_CA], pos[_SLOT_N], pos[_SLOT_CA])
    _fill_bl(1, has[_SLOT_CA] & has[_SLOT_C], pos[_SLOT_CA], pos[_SLOT_C])
    _fill_bl(2, has[_SLOT_C] & has[_SLOT_O], pos[_SLOT_C], pos[_SLOT_O])
    _fill_bl(3, nxt & has[_SLOT_C] & has_n_next, pos[_SLOT_C], n_next_pos)

    _fill_ba(
        0, has[_SLOT_N] & has[_SLOT_CA] & has[_SLOT_C],
        pos[_SLOT_N], pos[_SLOT_CA], pos[_SLOT_C],
    )
    _fill_ba(
        1, nxt & has[_SLOT_CA] & has[_SLOT_C] & has_n_next,
        pos[_SLOT_CA], pos[_SLOT_C], n_next_pos,
    )
    _fill_ba(
        2, nxt & has[_SLOT_C] & has_n_next & has_ca_next,
        pos[_SLOT_C], n_next_pos, ca_next_pos,
    )

    bl_ideal = np.array([IDEAL_BOND_LENGTHS[k] for _, k in _BL_COLUMNS])
    ba_ideal = np.array([IDEAL_BOND_ANGLES[k] for k in _BA_COLUMNS])
    # Zero out the unmeasured cells so downstream comparisons never touch NaN.
    bl_dev = np.where(bl_mask, bl - bl_ideal[None, :], 0.0)
    ba_dev = np.where(ba_mask, ba - ba_ideal[None, :], 0.0)

    # Row-major ravel reproduces the per-residue interleaving of the previous
    # implementation, so the RMSD summation order (and thus its last bits) is
    # unchanged.
    bl_flat = bl_dev[bl_mask.nonzero()] if bl_mask.any() else np.empty(0)
    ba_flat = ba_dev[ba_mask.nonzero()] if ba_mask.any() else np.empty(0)
    bl_rmsd = float(np.sqrt(np.mean(bl_flat ** 2))) if bl_flat.size else 0.0
    ba_rmsd = float(np.sqrt(np.mean(ba_flat ** 2))) if ba_flat.size else 0.0

    bl_out = bl_mask & (np.abs(bl_dev) > 3 * BL_SIGMA)
    ba_out = ba_mask & (np.abs(ba_dev) > 3 * BA_SIGMA)

    outliers: List[dict] = []
    rows = np.flatnonzero(bl_out.any(axis=1) | ba_out.any(axis=1))
    labels = table.labels
    for i in rows:
        label = labels[i]
        pair_label = f"{label}→{labels[i + 1]}" if i + 1 < n else label
        # Emission order per residue: intra bonds, intra angle, then the
        # inter-residue bond and its two angles.
        for col, is_bond in ((0, True), (1, True), (2, True), (0, False),
                             (3, True), (1, False), (2, False)):
            if is_bond:
                if not bl_out[i, col]:
                    continue
                name, ideal_key = _BL_COLUMNS[col]
                dev = float(bl_dev[i, col])
                outliers.append({
                    "residue": pair_label if col == 3 else label,
                    "type": "bond_length",
                    "bond": name, "value": round(float(bl[i, col]), 4),
                    "ideal": IDEAL_BOND_LENGTHS[ideal_key],
                    "deviation": round(dev, 4),
                    "sigma": round(abs(dev) / BL_SIGMA, 1),
                })
            else:
                if not ba_out[i, col]:
                    continue
                name = _BA_COLUMNS[col]
                dev = float(ba_dev[i, col])
                outliers.append({
                    "residue": label if col == 0 else pair_label,
                    "type": "bond_angle",
                    "angle": name, "value": round(float(ba[i, col]), 2),
                    "ideal": IDEAL_BOND_ANGLES[name],
                    "deviation": round(dev, 2),
                    "sigma": round(abs(dev) / BA_SIGMA, 1),
                })

    return {
        "bl_rmsd":  round(bl_rmsd, 5),
        "ba_rmsd":  round(ba_rmsd, 3),
        "n_bonds_checked": int(bl_mask.sum()),
        "n_angles_checked": int(ba_mask.sum()),
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


def backbone_dihedrals(table: _ResidueTable) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (phi, psi, valid) arrays over the residues of *table*.

    φ = C(i−1)–N(i)–Cα(i)–C(i) and ψ = N(i)–Cα(i)–C(i)–N(i+1), each gated on
    peptide continuity with the relevant neighbour. ``valid`` marks residues
    where both angles are defined and finite.
    """
    n = table.n
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    if n == 0:
        return phi, psi, np.zeros(0, dtype=bool)

    pos, has = table.pos, table.has
    core = has[_SLOT_CA] & has[_SLOT_N] & has[_SLOT_C]

    prev_ok = np.zeros(n, dtype=bool)
    if n > 1:
        prev_ok[1:] = table.continuous
    prev_has_c = np.zeros(n, dtype=bool)
    prev_has_c[1:] = has[_SLOT_C][:-1]
    c_prev = np.full((n, 3), np.nan)
    c_prev[1:] = pos[_SLOT_C][:-1]

    next_ok = table.next_pair_mask()
    next_has_n = np.zeros(n, dtype=bool)
    next_has_n[:-1] = has[_SLOT_N][1:]
    n_next = np.full((n, 3), np.nan)
    n_next[:-1] = pos[_SLOT_N][1:]

    phi_mask = core & prev_ok & prev_has_c
    if phi_mask.any():
        phi[phi_mask] = _dihedrals_v(
            c_prev[phi_mask], pos[_SLOT_N][phi_mask],
            pos[_SLOT_CA][phi_mask], pos[_SLOT_C][phi_mask],
        )
    psi_mask = core & next_ok & next_has_n
    if psi_mask.any():
        psi[psi_mask] = _dihedrals_v(
            pos[_SLOT_N][psi_mask], pos[_SLOT_CA][psi_mask],
            pos[_SLOT_C][psi_mask], n_next[psi_mask],
        )

    valid = phi_mask & psi_mask & ~np.isnan(phi) & ~np.isnan(psi)
    return phi, psi, valid


def _check_ramachandran(source, ordered_keys: Optional[List[Tuple]] = None) -> dict:
    """
    Compute φ/ψ for each residue and classify with the Ramachandran scorer.

    Returns pct_favored, pct_allowed, pct_outlier, and outlier records.
    """
    from .ramachandran import classify_regions

    table = _as_table(source, ordered_keys)
    phi, psi, valid = backbone_dihedrals(table)
    idx = np.flatnonzero(valid)
    n_total = int(idx.size)

    if n_total == 0:
        return {
            "n_evaluated": 0, "n_favored": 0, "n_allowed": 0, "n_outlier": 0,
            "pct_favored": 0.0, "pct_allowed": 0.0, "pct_outlier": 0.0,
            "outliers": [],
        }

    rtypes = [table.rtypes[i] for i in idx]
    regions = classify_regions(phi[idx], psi[idx], rtypes)

    outlier_records: List[dict] = []
    n_favored = n_allowed = n_outlier = 0
    for k, i in enumerate(idx):
        region = regions[k]
        if region == "favored":
            n_favored += 1
            continue
        if region == "allowed":
            n_allowed += 1
            continue
        n_outlier += 1
        outlier_records.append({
            "label":  table.labels[i],
            "rtype":  rtypes[k],
            "phi":    round(float(phi[i]), 2),
            "psi":    round(float(psi[i]), 2),
            "region": region,
        })

    return {
        "n_evaluated":  n_total,
        "n_favored":    n_favored,
        "n_allowed":    n_allowed,
        "n_outlier":    n_outlier,
        "pct_favored":  round(100.0 * n_favored / n_total, 2),
        "pct_allowed":  round(100.0 * n_allowed / n_total, 2),
        "pct_outlier":  round(100.0 * n_outlier / n_total, 2),
        "outliers":     outlier_records,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Check 4: Peptide Planarity (ω dihedral)
# ─────────────────────────────────────────────────────────────────────────────

def _check_planarity(source, ordered_keys: Optional[List[Tuple]] = None) -> dict:
    """
    Measure ω dihedral (CA_i–C_i–N_{i+1}–CA_{i+1}) for each peptide bond.

    Flags bonds with |ω − 180°| > 6° (trans) or |ω| > 6° (cis).
    Reports % within 6°, mean deviation, and outlier records.
    """
    table = _as_table(source, ordered_keys)
    n = table.n
    if n < 2:
        return {
            "n_bonds_checked": 0, "n_within_6deg": 0,
            "pct_within_6deg": 100.0, "mean_deviation": 0.0, "outliers": [],
        }

    pos, has = table.pos, table.has
    mask = (
        table.continuous
        & has[_SLOT_CA][:-1] & has[_SLOT_C][:-1]
        & has[_SLOT_N][1:] & has[_SLOT_CA][1:]
    )
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return {
            "n_bonds_checked": 0, "n_within_6deg": 0,
            "pct_within_6deg": 100.0, "mean_deviation": 0.0, "outliers": [],
        }

    omega = _dihedrals_v(
        pos[_SLOT_CA][idx], pos[_SLOT_C][idx],
        pos[_SLOT_N][idx + 1], pos[_SLOT_CA][idx + 1],
    )
    finite = ~np.isnan(omega)
    idx = idx[finite]
    omega = omega[finite]

    dev_trans = np.minimum(np.abs(omega - 180.0), np.abs(omega + 180.0))
    dev_cis = np.abs(omega)
    dev = np.minimum(dev_trans, dev_cis)

    n_total = int(dev.size)
    n_within = int((dev <= 6.0).sum())
    pct_within = 100.0 * n_within / n_total if n_total else 100.0
    mean_dev = float(np.mean(dev)) if n_total else 0.0

    outliers: List[dict] = []
    labels = table.labels
    for k in np.flatnonzero(dev > 6.0):
        i = int(idx[k])
        outliers.append({
            "peptide_bond": f"{labels[i]}→{labels[i + 1]}",
            "omega": round(float(omega[k]), 2),
            "deviation": round(float(dev[k]), 2),
            "type": "cis" if dev_cis[k] < dev_trans[k] else "trans",
        })

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


@lru_cache(maxsize=512)
def _intra_residue_bond_pairs_cached(resname: str) -> Tuple[Tuple[str, str], ...]:
    tmpl = _bond_template_resname(resname)
    pairs = list(_BACKBONE_BOND_PAIRS)
    pairs.extend(_SIDECHAIN_BOND_PAIRS.get(tmpl, ()))
    # Unknown residues: still connect CA–CB when present via backbone only;
    # distance fallback in _clash_excluded_index_pairs covers remaining 1-3.
    if tmpl not in _SIDECHAIN_BOND_PAIRS and tmpl not in ("GLY",):
        pairs.append(("CA", "CB"))
    return tuple(pairs)


def _intra_residue_bond_pairs(resname: str) -> List[Tuple[str, str]]:
    """Covalent name pairs inside a residue (memoised per residue name)."""
    return list(_intra_residue_bond_pairs_cached(resname))


def _residue_key(atom: _Atom) -> Tuple[str, int, str]:
    return (atom.chain, atom.resseq, atom.icode)


#: Graph distance up to which same-residue pairs are always excluded
#: (1-2 / 1-3 / 1-4 in MolProbity Probe's ``-4`` convention).
_INTRA_EXCLUDE_DEPTH = 3
#: Extra depth applied when at least one partner is a hydrogen (1-5).
_INTRA_EXCLUDE_DEPTH_H = 4


@lru_cache(maxsize=4096)
def _intra_residue_excluded_name_pairs(
    resname: str,
    signature: frozenset,
) -> Tuple[Tuple[str, str], ...]:
    """Same-residue atom-name pairs that must not score as clashes.

    Args:
        resname: PDB residue name (mapped to a side-chain bond template).
        signature: ``frozenset`` of ``(atom_name, is_hydrogen)`` for the atom
            names present in this residue. Making the hydrogen flag part of the key
            keeps the result exact for oddities such as an ``H``-prefixed
            heavy-atom name or a deuterium whose name does not start with ``H``.

    Returns:
        Tuple of ``(name_a, name_b)`` pairs, each appearing once.

    A pair is excluded when its covalent graph distance is ≤ 3, or exactly 4
    with at least one hydrogen partner — identical to the breadth-first rule the
    previous per-residue implementation applied, but computed once per distinct
    residue signature.
    """
    is_h = dict(signature)
    names = sorted(is_h)
    present = set(names)

    adj: Dict[str, Set[str]] = {nm: set() for nm in names}
    for a_name, b_name in _intra_residue_bond_pairs_cached(resname):
        if a_name in present and b_name in present:
            adj[a_name].add(b_name)
            adj[b_name].add(a_name)

    out: Set[Tuple[str, str]] = set()
    for start in names:
        start_is_h = is_h[start]
        max_depth = _INTRA_EXCLUDE_DEPTH_H if start_is_h else _INTRA_EXCLUDE_DEPTH
        frontier = {start}
        visited = {start}
        for depth in range(max_depth):
            nxt: Set[str] = set()
            for node in frontier:
                for nb in adj[node]:
                    if nb in visited:
                        continue
                    visited.add(nb)
                    nxt.add(nb)
                    # Depth 3 (a 1-5 contact) needs a hydrogen partner.
                    if depth == _INTRA_EXCLUDE_DEPTH and not (
                        start_is_h or is_h[nb]
                    ):
                        continue
                    out.add((start, nb) if start < nb else (nb, start))
            frontier = nxt
    return tuple(sorted(out))


#: Cross-residue name pairs excluded across a peptide bond (i → i+1).
#: Hoisted to module scope: this was rebuilt for every adjacent residue pair.
_CROSS_RESIDUE_EXCLUDED_PAIRS: Tuple[Tuple[str, str], ...] = (
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
    # Cα hydrogens on residue i+1 (C–N–CA–HA* = 1-4)
    ("C", "HA"),
    ("C", "HA2"),
    ("C", "HA3"),
    # Peptide 1-5 involving model-built HA (O/CA/N ··· HA*)
    ("O", "HA"),
    ("O", "HA2"),
    ("O", "HA3"),
    ("CA", "HA"),
    ("CA", "HA2"),
    ("CA", "HA3"),
    ("N", "HA"),
    ("N", "HA2"),
    ("N", "HA3"),
    # X–Pro: HA(i) ··· CD(i+1) is 1-5 via C–N
    ("HA", "CD"),
    ("HA2", "CD"),
    ("HA3", "CD"),
)


def _clash_excluded_index_pairs(atoms: List[_Atom]) -> Set[Tuple[int, int]]:
    """
    Build index pairs that must not score as clashes (covalent 1-2/1-3/1-4).

    Uses standard amino-acid topology (not a fragile distance cutoff). Same-
    residue CA–CG (~2.5–2.7 Å) is a classic 1-3 false positive when a 2.6 Å
    distance gate is used — that bug inflated clash scores on AFDB/PDB models.

    Returns a set of ``(i, j)`` index pairs with ``i < j``. The audit itself
    uses :func:`_excluded_pair_keys`, which returns the same information as a
    sorted int64 array so that candidate filtering stays vectorised.
    """
    n = len(atoms)
    keys = _excluded_pair_keys(atoms)
    return {(int(k // n), int(k % n)) for k in keys}


def _excluded_pair_keys(atoms: List[_Atom]) -> np.ndarray:
    """Excluded atom-index pairs as a sorted, unique array of ``i * n + j`` keys.

    Packing pairs into a single int64 lets :func:`_check_clashes` reject
    covalently-connected candidates with one ``searchsorted`` over the whole
    candidate array instead of a Python-level set lookup per candidate pair.
    """
    n_atoms = len(atoms)
    # name → list of atom indices in that residue
    by_res: Dict[Tuple[str, int, str], Dict[str, List[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for idx, atom in enumerate(atoms):
        by_res[_residue_key(atom)][atom.sname].append(idx)

    # Excluded pairs are accumulated as packed ``min * n_atoms + max`` integers
    # in a plain set: an int key avoids the tuple allocation the previous
    # implementation paid for each of the O(20) exclusions per residue, and
    # ``add`` below binds ``set.add`` directly so the inner loops make no
    # Python-level function calls.
    excluded: Set[int] = set()
    add = excluded.add

    # Exclude 1-2, 1-3, and 1-4 within each residue (MolProbity Probe -4). For
    # pairs involving a hydrogen, also exclude 1-5: model-built H positions make
    # same-residue HA···CD1-style contacts noisy otherwise.
    #
    # The topology is a property of the residue *name* and the set of atom names
    # actually present, not of the individual residue, so the graph traversal is
    # memoised per (resname, atom-name signature). A structure with hundreds of
    # residues typically has only a few dozen distinct signatures, which turns
    # the per-residue BFS into a dictionary hit.
    for rkey, name_map in by_res.items():
        sample_idx = next(iter(next(iter(name_map.values()))))
        resname = atoms[sample_idx].resname
        signature = frozenset(
            (nm, _is_hydrogen_atom(atoms[idxs[0]]))
            for nm, idxs in name_map.items()
        )
        for a_name, b_name in _intra_residue_excluded_name_pairs(resname, signature):
            for i in name_map[a_name]:
                for j in name_map[b_name]:
                    add(i * n_atoms + j if i < j else j * n_atoms + i)

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
        for a_name, b_name in _CROSS_RESIDUE_EXCLUDED_PAIRS:
            for ia in m0.get(a_name, ()):
                for ib in m1.get(b_name, ()):
                    add(ia * n_atoms + ib if ia < ib else ib * n_atoms + ia)

    # Covalent disulfides: SG–SG within 2.5 Å, plus 1-3/1-4 across the bridge
    sg_idxs = [i for i, a in enumerate(atoms) if a.sname == "SG"]
    for a in range(len(sg_idxs)):
        for b in range(a + 1, len(sg_idxs)):
            i, j = sg_idxs[a], sg_idxs[b]
            dist = float(np.linalg.norm(atoms[i].xyz - atoms[j].xyz))
            if dist > 2.5:
                continue
            add(i * n_atoms + j if i < j else j * n_atoms + i)
            # 1-3: CB–SG_partner ; 1-4: CA–SG_partner, CB–CB
            for sg_i, sg_j in ((i, j), (j, i)):
                rkey = _residue_key(atoms[sg_i])
                partner_key = _residue_key(atoms[sg_j])
                for cb_i in by_res[rkey].get("CB", ()):
                    add(cb_i * n_atoms + sg_j if cb_i < sg_j
                        else sg_j * n_atoms + cb_i)
                    for cb_j in by_res[partner_key].get("CB", ()):
                        add(cb_i * n_atoms + cb_j if cb_i < cb_j
                            else cb_j * n_atoms + cb_i)
                for ca_i in by_res[rkey].get("CA", ()):
                    add(ca_i * n_atoms + sg_j if ca_i < sg_j
                        else sg_j * n_atoms + ca_i)

    # Distance fallback ONLY for residues without a known covalent template
    # (ligands / nonstandard) — never for standard AA where topology is complete.
    for rkey, name_map in by_res.items():
        sample_idx = next(iter(next(iter(name_map.values()))))
        tmpl = _bond_template_resname(atoms[sample_idx].resname)
        if tmpl in _SIDECHAIN_BOND_PAIRS or tmpl == "GLY":
            continue
        idxs = sorted(i for lst in name_map.values() for i in lst)
        if len(idxs) < 2:
            continue
        sub = np.array([[atoms[i].x, atoms[i].y, atoms[i].z] for i in idxs])
        d2 = np.einsum(
            "ijk,ijk->ij",
            sub[:, None, :] - sub[None, :, :],
            sub[:, None, :] - sub[None, :, :],
        )
        close_a, close_b = np.nonzero(np.triu(d2 < 2.9 ** 2, k=1))
        for a, b in zip(close_a, close_b):
            i, j = idxs[a], idxs[b]
            add(i * n_atoms + j if i < j else j * n_atoms + i)

    if not excluded:
        return np.empty(0, dtype=np.int64)
    return np.fromiter(sorted(excluded), dtype=np.int64, count=len(excluded))


def _angle_degrees(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle ABC in degrees (vertex at *b*)."""
    ba = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    bc = np.asarray(c, dtype=float) - np.asarray(b, dtype=float)
    na = float(np.linalg.norm(ba)) + 1e-12
    nc = float(np.linalg.norm(bc)) + 1e-12
    cosang = float(np.dot(ba, bc) / (na * nc))
    cosang = max(-1.0, min(1.0, cosang))
    return float(math.degrees(math.acos(cosang)))


def _is_acceptor_atom(atom: _Atom) -> bool:
    name = atom.name.strip().upper()
    elem = atom.element_upper
    if elem == "O" or name in (
        "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OXT"
    ):
        return True
    # Occasional N acceptors (HIS NE2/ND1) — treat named N as acceptor only when
    # not acting as the donor in the pair (handled by caller roles).
    return False


def _is_donor_n_atom(atom: _Atom) -> bool:
    name = atom.name.strip().upper()
    return name in (
        "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ"
    )


def _is_polar_h_atom(atom: _Atom) -> bool:
    """True for H-bond-capable hydrogens (not Cα HA/HA2/HA3)."""
    if not _is_hydrogen_atom(atom):
        return False
    return atom.name.strip().upper() in _POLAR_H_NAMES


def _is_hbond_donor_acceptor_pair(
    a: _Atom,
    b: _Atom,
    parent_n_by_h_index: Optional[Dict[int, _Atom]] = None,
    atom_index: Optional[Dict[int, int]] = None,
) -> bool:
    """
    Probe-inspired H-bond filter: suppress a clash only when geometry is good.

    - Carbon-bound HA/HA2/HA3 never count as H-bond donors.
    - Polar H···acceptor: H-bond if H···Acc ≤ 2.5 Å and donor–H–Acc ≥ 120°.
    - Heavy N···O: H-bond only in the 2.5–3.5 Å window (side-chain NH without H).
    """
    parent_n_by_h_index = parent_n_by_h_index or {}

    def _idx(atom: _Atom) -> Optional[int]:
        if atom_index is None:
            return None
        return atom_index.get(id(atom))

    # Polar H ··· acceptor
    h_atom: Optional[_Atom] = None
    acc: Optional[_Atom] = None
    if _is_polar_h_atom(a) and _is_acceptor_atom(b):
        h_atom, acc = a, b
    elif _is_polar_h_atom(b) and _is_acceptor_atom(a):
        h_atom, acc = b, a

    if h_atom is not None and acc is not None:
        dist = float(np.linalg.norm(h_atom.xyz - acc.xyz))
        if dist > _HBOND_H_ACC_MAX:
            return False
        h_i = _idx(h_atom)
        parent = parent_n_by_h_index.get(h_i) if h_i is not None else None
        if parent is None:
            # No parent → require only distance (still stricter than blanket skip)
            return True
        ang = _angle_degrees(parent.xyz, h_atom.xyz, acc.xyz)
        return ang >= _HBOND_ANGLE_MIN

    # Heavy-atom N···O (no polar H involved in this pair)
    n_atom: Optional[_Atom] = None
    o_atom: Optional[_Atom] = None
    if _is_donor_n_atom(a) and _is_acceptor_atom(b):
        n_atom, o_atom = a, b
    elif _is_donor_n_atom(b) and _is_acceptor_atom(a):
        n_atom, o_atom = b, a
    if n_atom is not None and o_atom is not None:
        dist = float(np.linalg.norm(n_atom.xyz - o_atom.xyz))
        return _HBOND_NO_MIN <= dist <= _HBOND_NO_MAX

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


def _place_tetrahedral_h(
    center: np.ndarray,
    neighbor_positions: List[np.ndarray],
    bond_length: float,
) -> np.ndarray:
    """Place H opposite the sum of unit vectors from *center* to neighbours."""
    c = np.asarray(center, dtype=float)
    acc = np.zeros(3, dtype=float)
    for npos in neighbor_positions:
        v = np.asarray(npos, dtype=float) - c
        acc = acc + v / (float(np.linalg.norm(v)) + 1e-12)
    direction = -acc / (float(np.linalg.norm(acc)) + 1e-12)
    return c + direction * bond_length


def _place_gly_ha2_ha3(
    n_pos: np.ndarray,
    ca_pos: np.ndarray,
    c_pos: np.ndarray,
    bond_length: float = _BOND_LEN_CH,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Place Gly HA2/HA3 in a tetrahedral arrangement about CA.

    Uses the N–CA–C bisector and the peptide-plane normal (Reduce-style).
    """
    ca = np.asarray(ca_pos, dtype=float)
    v_n = np.asarray(n_pos, dtype=float) - ca
    v_c = np.asarray(c_pos, dtype=float) - ca
    v_n = v_n / (float(np.linalg.norm(v_n)) + 1e-12)
    v_c = v_c / (float(np.linalg.norm(v_c)) + 1e-12)
    bisector = v_n + v_c
    bisector = bisector / (float(np.linalg.norm(bisector)) + 1e-12)
    normal = np.cross(v_n, v_c)
    nrm = float(np.linalg.norm(normal))
    if nrm < 1e-8:
        # Degenerate N–CA–C: pick an axis not parallel to the bisector
        normal = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(bisector, normal))) > 0.9:
            normal = np.array([0.0, 1.0, 0.0])
        normal = normal - bisector * float(np.dot(normal, bisector))
        nrm = float(np.linalg.norm(normal))
    normal = normal / (nrm + 1e-12)
    # ± half of tetrahedral angle from the direction opposite N/C
    half = math.radians(54.75)
    ref = -bisector
    ha2_dir = ref * math.cos(half) + normal * math.sin(half)
    ha3_dir = ref * math.cos(half) - normal * math.sin(half)
    ha2_dir = ha2_dir / (float(np.linalg.norm(ha2_dir)) + 1e-12)
    ha3_dir = ha3_dir / (float(np.linalg.norm(ha3_dir)) + 1e-12)
    return ca + ha2_dir * bond_length, ca + ha3_dir * bond_length


def _place_tetrahedral_h_batch(
    center: np.ndarray,
    neighbor_positions: Tuple[np.ndarray, ...],
    bond_length: float,
) -> np.ndarray:
    """Batched :func:`_place_tetrahedral_h` over (m, 3) coordinate arrays."""
    acc = np.zeros_like(center)
    for npos in neighbor_positions:
        v = npos - center
        acc = acc + v / (_norms_v(v)[:, None] + 1e-12)
    direction = -acc / (_norms_v(acc)[:, None] + 1e-12)
    return center + direction * bond_length


def _place_amide_h_batch(
    n_pos: np.ndarray,
    prev_c_pos: np.ndarray,
    ca_pos: np.ndarray,
    bond_length: float = _BOND_LEN_NH,
) -> np.ndarray:
    """Batched backbone amide HN placement (bisector of N→C(i−1) and N→Cα)."""
    v_nc = prev_c_pos - n_pos
    v_nca = ca_pos - n_pos
    h_dir = -(
        v_nc / (_norms_v(v_nc)[:, None] + 1e-12)
        + v_nca / (_norms_v(v_nca)[:, None] + 1e-12)
    )
    h_dir = h_dir / (_norms_v(h_dir)[:, None] + 1e-12)
    return n_pos + h_dir * bond_length


def _place_gly_ha2_ha3_batch(
    n_pos: np.ndarray,
    ca_pos: np.ndarray,
    c_pos: np.ndarray,
    bond_length: float = _BOND_LEN_CH,
) -> Tuple[np.ndarray, np.ndarray]:
    """Batched :func:`_place_gly_ha2_ha3`.

    Rows with a degenerate (collinear) N–Cα–C triad fall back to the scalar
    implementation, which resolves the ambiguity by choosing an axis that is not
    parallel to the bisector.
    """
    v_n = n_pos - ca_pos
    v_c = c_pos - ca_pos
    v_n = v_n / (_norms_v(v_n)[:, None] + 1e-12)
    v_c = v_c / (_norms_v(v_c)[:, None] + 1e-12)
    bisector = v_n + v_c
    bisector = bisector / (_norms_v(bisector)[:, None] + 1e-12)
    normal = np.cross(v_n, v_c)
    nrm = _norms_v(normal)
    normal = normal / (nrm[:, None] + 1e-12)

    half = math.radians(54.75)
    ref = -bisector
    cos_h, sin_h = math.cos(half), math.sin(half)
    ha2_dir = ref * cos_h + normal * sin_h
    ha3_dir = ref * cos_h - normal * sin_h
    ha2_dir = ha2_dir / (_norms_v(ha2_dir)[:, None] + 1e-12)
    ha3_dir = ha3_dir / (_norms_v(ha3_dir)[:, None] + 1e-12)
    ha2 = ca_pos + ha2_dir * bond_length
    ha3 = ca_pos + ha3_dir * bond_length

    for row in np.flatnonzero(nrm < 1e-8):
        ha2[row], ha3[row] = _place_gly_ha2_ha3(
            n_pos[row], ca_pos[row], c_pos[row], bond_length
        )
    return ha2, ha3


def _make_h_atom(
    serial: int,
    name: str,
    template: _Atom,
    pos: np.ndarray,
) -> _Atom:
    # PDB atom name: right-align in 4 columns when length ≤ 3 (standard H style)
    nm = name.strip().upper()
    if len(nm) == 1:
        pdb_name = f" {nm}  "
    elif len(nm) == 2:
        pdb_name = f" {nm} "
    elif len(nm) == 3:
        pdb_name = f" {nm}"
    else:
        pdb_name = nm[:4]
    return _Atom(
        record="ATOM",
        serial=serial,
        name=pdb_name,
        altloc=" ",
        resname=template.resname,
        chain=template.chain,
        resseq=template.resseq,
        icode=template.icode,
        x=float(pos[0]),
        y=float(pos[1]),
        z=float(pos[2]),
        element="H",
    )


def _add_backbone_hydrogens(atoms: List[_Atom]) -> List[_Atom]:
    """
    Place hydrogens for MolProbity/Probe-style clash detection.

    1. Backbone amide HN (not proline), only for same-chain peptide-continuous
       neighbours (resseq gap ≤ 2, CA–CA ≤ 4.5 Å; insertion codes included).
    2. Cα HA for residues with CB (tetrahedral from N/C/CB).
    3. Glycine HA2/HA3 (tetrahedral from N/C and the peptide-plane normal).

    Placement is planned in a first Python pass and then evaluated in three
    batched geometry calls, so the number of numpy dispatches is O(1) in the
    number of residues rather than O(6) per residue. Atom serials are still
    assigned in the original interleaved order (HN then Cα-H, residue by
    residue) so downstream serial-keyed bookkeeping is unchanged.
    """
    residues: Dict[Tuple[str, int, str], Dict[str, _Atom]] = {}
    for a in atoms:
        key = (a.chain, a.resseq, a.icode)
        res = residues.get(key)
        if res is None:
            res = residues[key] = {}
        res[a.sname] = a

    sorted_keys = sorted(residues.keys())

    # ── Pass 1: decide what to place, collecting geometry inputs ──────────
    # plan entries: ("H", hn_row) | ("GLY", gly_row) | ("HA", ha_row)
    plan: List[Tuple[str, int, _Atom]] = []
    hn_n: List[np.ndarray] = []
    hn_prev_c: List[np.ndarray] = []
    hn_ca: List[np.ndarray] = []
    gly_n: List[np.ndarray] = []
    gly_ca: List[np.ndarray] = []
    gly_c: List[np.ndarray] = []
    ha_ca: List[np.ndarray] = []
    ha_n: List[np.ndarray] = []
    ha_c: List[np.ndarray] = []
    ha_cb: List[np.ndarray] = []

    for idx, key in enumerate(sorted_keys):
        res = residues[key]
        n_atom = res.get("N")
        ca_atom = res.get("CA")
        c_atom = res.get("C")
        if ca_atom is None:
            continue

        resname_u = (n_atom or ca_atom).resname.upper().strip()
        tmpl = _bond_template_resname(resname_u)

        # ── Amide HN ──────────────────────────────────────────────────
        if (
            n_atom is not None
            and resname_u not in PROLINE_RESNAMES
            and "H" not in res
            and "HN" not in res
            and idx > 0
        ):
            prev_key = sorted_keys[idx - 1]
            if prev_key[0] == key[0] and abs(key[1] - prev_key[1]) <= 2:
                prev_res = residues.get(prev_key, {})
                prev_c = prev_res.get("C")
                prev_ca = prev_res.get("CA")
                ok = prev_c is not None
                if ok and prev_ca is not None:
                    if float(np.linalg.norm(prev_ca.xyz - ca_atom.xyz)) > 4.5:
                        ok = False
                if ok:
                    plan.append(("H", len(hn_n), n_atom))
                    hn_n.append(n_atom.xyz)
                    hn_prev_c.append(prev_c.xyz)
                    hn_ca.append(ca_atom.xyz)

        # ── Cα hydrogens ──────────────────────────────────────────────
        if n_atom is None or c_atom is None:
            continue

        if tmpl == "GLY" or resname_u in GLYCINE_RESNAMES:
            if "HA2" in res or "HA3" in res or "HA" in res:
                continue
            plan.append(("GLY", len(gly_n), ca_atom))
            gly_n.append(n_atom.xyz)
            gly_ca.append(ca_atom.xyz)
            gly_c.append(c_atom.xyz)
        else:
            cb_atom = res.get("CB")
            if cb_atom is None or "HA" in res:
                continue
            plan.append(("HA", len(ha_ca), ca_atom))
            ha_ca.append(ca_atom.xyz)
            ha_n.append(n_atom.xyz)
            ha_c.append(c_atom.xyz)
            ha_cb.append(cb_atom.xyz)

    if not plan:
        return list(atoms)

    # ── Pass 2: batched geometry ──────────────────────────────────────────
    hn_pos = (
        _place_amide_h_batch(
            np.asarray(hn_n), np.asarray(hn_prev_c), np.asarray(hn_ca)
        )
        if hn_n else np.empty((0, 3))
    )
    if gly_n:
        gly_ha2, gly_ha3 = _place_gly_ha2_ha3_batch(
            np.asarray(gly_n), np.asarray(gly_ca), np.asarray(gly_c)
        )
    else:
        gly_ha2 = gly_ha3 = np.empty((0, 3))
    ha_pos = (
        _place_tetrahedral_h_batch(
            np.asarray(ha_ca),
            (np.asarray(ha_n), np.asarray(ha_c), np.asarray(ha_cb)),
            _BOND_LEN_CH,
        )
        if ha_ca else np.empty((0, 3))
    )

    # ── Pass 3: materialise atoms in the original order ───────────────────
    new_h_atoms: List[_Atom] = []
    h_serial = max((a.serial for a in atoms), default=0) + 1
    for kind, row, template in plan:
        if kind == "H":
            new_h_atoms.append(_make_h_atom(h_serial, "H", template, hn_pos[row]))
            h_serial += 1
        elif kind == "GLY":
            new_h_atoms.append(
                _make_h_atom(h_serial, "HA2", template, gly_ha2[row])
            )
            h_serial += 1
            new_h_atoms.append(
                _make_h_atom(h_serial, "HA3", template, gly_ha3[row])
            )
            h_serial += 1
        else:
            new_h_atoms.append(_make_h_atom(h_serial, "HA", template, ha_pos[row]))
            h_serial += 1

    return atoms + new_h_atoms


#: Minimum vdW overlap (Å) that counts as a clash (MolProbity Probe default).
CLASH_OVERLAP_CUTOFF = 0.4


def _check_clashes(atoms: List[_Atom]) -> dict:
    """
    Detect steric clashes between non-bonded atoms using a scipy KD-tree.

    Two atoms clash when their distance < (rvdw_A + rvdw_B − 0.4) Å.
    Deposited hydrogens are stripped; amide HN and Cα HA/HA2/HA3 are
    re-added (MolProbity Reduce-style). Covalent 1-2 / 1-3 / 1-4 pairs are
    excluded via residue topology. H-bonds suppress clashes only when
    Probe-like distance/angle geometry is satisfied.

    Clash score = clashes per 1000 heavy atoms.

    The candidate set is pruned in increasing order of cost: the KD-tree
    neighbour radius is the largest distance that can possibly overlap
    (``2·max(rvdw) − 0.4``, not ``2·max(rvdw)``, which cut candidate pairs by
    roughly a quarter), then a vectorised overlap test, then the packed
    covalent-exclusion lookup, and only then the per-pair H-bond geometry test.
    Reordering is safe because the H-bond rule can only *suppress* a clash.
    """
    heavy = [a for a in atoms if not _is_hydrogen_atom(a)]
    n_atoms = len(heavy)
    if n_atoms < 2:
        return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}

    all_atoms_for_check = _add_backbone_hydrogens(heavy)
    n_all = len(all_atoms_for_check)

    coords = np.empty((n_all, 3), dtype=float)
    radii = np.empty(n_all, dtype=float)
    for i, a in enumerate(all_atoms_for_check):
        coords[i, 0] = a.x
        coords[i, 1] = a.y
        coords[i, 2] = a.z
        radii[i] = VDW_RADII.get(a.element_upper, VDW_DEFAULT)

    # A pair can only clash if dist < r_i + r_j − 0.4 ≤ 2·max(r) − 0.4.
    query_radius = float(np.max(radii)) * 2.0 - CLASH_OVERLAP_CUTOFF
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=query_radius, output_type="ndarray")
    if pairs.size == 0:
        return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}

    lo = np.minimum(pairs[:, 0], pairs[:, 1]).astype(np.int64)
    hi = np.maximum(pairs[:, 0], pairs[:, 1]).astype(np.int64)

    # 1. Vectorised overlap test.
    delta = coords[hi] - coords[lo]
    dist = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    vdw_sum = radii[lo] + radii[hi]
    overlap = vdw_sum - dist
    keep = overlap > CLASH_OVERLAP_CUTOFF
    if not keep.any():
        return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}
    lo, hi, dist, vdw_sum, overlap = (
        lo[keep], hi[keep], dist[keep], vdw_sum[keep], overlap[keep]
    )

    # 2. Drop covalent 1-2 / 1-3 / 1-4 (and H 1-5) pairs in one pass.
    excl = _excluded_pair_keys(all_atoms_for_check)
    if excl.size:
        packed = lo * n_all + hi
        slot = np.searchsorted(excl, packed)
        np.clip(slot, 0, excl.size - 1, out=slot)
        bonded = excl[slot] == packed
        keep = ~bonded
        if not keep.any():
            return {"n_clashes": 0, "clash_score": 0.0, "worst_clashes": []}
        lo, hi, dist, vdw_sum, overlap = (
            lo[keep], hi[keep], dist[keep], vdw_sum[keep], overlap[keep]
        )

    # 3. Probe-style H-bond suppression on the surviving candidates only.
    parent_n_by_h_index, atom_index = _hbond_context(all_atoms_for_check)

    clashes: List[dict] = []
    seen: Set[Tuple[int, int]] = set()
    for k in range(lo.size):
        i_i = int(lo[k])
        j_i = int(hi[k])
        ai = all_atoms_for_check[i_i]
        aj = all_atoms_for_check[j_i]
        if _is_hbond_donor_acceptor_pair(
            ai, aj,
            parent_n_by_h_index=parent_n_by_h_index,
            atom_index=atom_index,
        ):
            continue
        pair_key = (min(ai.serial, aj.serial), max(ai.serial, aj.serial))
        if pair_key in seen:
            continue
        seen.add(pair_key)
        clashes.append({
            "atom1": f"{ai.chain}:{ai.resname}{ai.resseq}.{ai.sname}",
            "atom2": f"{aj.chain}:{aj.resname}{aj.resseq}.{aj.sname}",
            "distance": round(float(dist[k]), 3),
            "overlap": round(float(overlap[k]), 3),
            "vdw_sum": round(float(vdw_sum[k]), 3),
        })

    clashes.sort(key=lambda c: -c["overlap"])
    n_clashes = len(clashes)
    clash_score = 1000.0 * n_clashes / n_atoms if n_atoms > 0 else 0.0

    return {
        "n_clashes": n_clashes,
        "clash_score": round(clash_score, 2),
        "worst_clashes": clashes[:20],
    }


def _hbond_context(
    all_atoms: List[_Atom],
) -> Tuple[Dict[int, _Atom], Dict[int, int]]:
    """Build the (parent-N per polar H, id→index) maps the H-bond test needs."""
    by_res_n: Dict[Tuple[str, int, str], _Atom] = {}
    for a in all_atoms:
        if a.sname == "N":
            by_res_n[_residue_key(a)] = a
    parent_n_by_h_index: Dict[int, _Atom] = {}
    atom_index: Dict[int, int] = {}
    for i, a in enumerate(all_atoms):
        atom_index[id(a)] = i
        if _is_polar_h_atom(a):
            parent = by_res_n.get(_residue_key(a))
            if parent is not None:
                parent_n_by_h_index[i] = parent
    return parent_n_by_h_index, atom_index


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

    # One shared coordinate table feeds all four backbone checks.
    table = _ResidueTable(residue_groups, ordered_keys)

    # Run all checks
    chirality   = _check_chirality(table)
    bond_geo    = _check_bond_geometry(table)
    rama        = _check_ramachandran(table)
    planarity   = _check_planarity(table)
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
