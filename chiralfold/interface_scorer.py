"""
ChiralFold — Binding Interface Scorer
=======================================

Computes geometric and energetic metrics at the binder:target interface.
Designed for evaluating D-peptide therapeutic candidates against L-protein
receptors, but works for any receptor:ligand pair.

Metrics computed
----------------
- Buried Surface Area (BSA): Approximate from interface atom count × ~10 Å²
- Shape Complementarity (SC): Fraction of interface pairs in the 3.5–5.0 Å
  complementary range vs. clashing (<3.5 Å) or no-contact (>5.0 Å).
- Hydrogen bonds: N/O…N/O donor-acceptor pairs within 2.5–3.5 Å.
- Salt bridges: Oppositely charged residue pairs (K/R/H vs D/E) within 4.0 Å
  Cα–Cα distance.
- Hydrophobic contacts: F/W/Y/L/I/V/A/M Cα pairs within 6.0 Å.
- Interface residues: Residues with ≥1 heavy atom within 5 Å of the partner.
- Interface score: Weighted composite 0–100.
"""

from __future__ import annotations

import os
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Contact distance thresholds (Å)
INTERFACE_CUTOFF    = 5.0    # any atom closer than this is "at interface"
HBOND_MIN           = 2.5    # hydrogen bond minimum N/O...N/O distance
HBOND_MAX           = 3.5    # hydrogen bond maximum
HBOND_DONORS        = {"N", "O"}   # simplified: any N or O can donate/accept
CLASH_CUTOFF        = 3.5    # pairs closer than this are clashing
COMPLEMENTARY_MAX   = 5.0    # pairs in 3.5–5.0 Å are complementary
SALT_BRIDGE_CUTOFF  = 4.0    # Cα–Cα cutoff for salt bridge detection (Å)
HYDROPHOBIC_CUTOFF  = 6.0    # Cα–Cα cutoff for hydrophobic contacts (Å)

# Residue sets
POSITIVE_RESIDUES   = {"LYS", "ARG", "HIS", "HIE", "HID", "HIP", "HSE", "HSD", "HSP",
                        "DLY", "DAR", "DHI"}
NEGATIVE_RESIDUES   = {"ASP", "GLU", "DAS", "DGL"}
HYDROPHOBIC_RESIDUES = {"PHE", "TRP", "TYR", "LEU", "ILE", "VAL", "ALA", "MET",
                         "DPN", "DTR", "DTY", "DLE", "DIL", "DVA", "DAL", "MED"}

# Approximate solvent-accessible area per interface atom (Å²)
SA_PER_ATOM = 10.0

# Score weights for composite interface_score
_WEIGHTS = {
    "hbonds":            3.0,
    "salt_bridges":      5.0,
    "hydrophobic":       0.5,
    "shape_complement":  30.0,
    "bsa_bonus":         0.002,  # per Å²
}
_SCORE_CAP = 100.0


# ─────────────────────────────────────────────────────────────────────────────
# PDB parsing (lightweight, no external dependencies)
# ─────────────────────────────────────────────────────────────────────────────

class _Atom:
    """Lightweight PDB atom record."""

    __slots__ = [
        "record", "serial", "name", "resname",
        "chain", "resseq", "x", "y", "z", "element",
    ]

    def __init__(self, record, serial, name, resname,
                 chain, resseq, x, y, z, element):
        self.record  = record
        self.serial  = serial
        self.name    = name
        self.resname = resname
        self.chain   = chain
        self.resseq  = resseq
        self.x       = x
        self.y       = y
        self.z       = z
        self.element = element

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    @property
    def element_sym(self) -> str:
        if self.element:
            return self.element.upper().strip()
        for ch in self.name.strip():
            if ch.isalpha():
                return ch.upper()
        return "C"


def _parse_atoms(pdb_path: str,
                 chains: Optional[List[str]] = None) -> List[_Atom]:
    """
    Parse ATOM/HETATM records from *pdb_path*, optionally filtered by chain(s).

    Skips water molecules, keeps only the first altloc position.
    """
    atoms: List[_Atom] = []
    seen: Dict[Tuple, bool] = {}

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
                continue

            try:
                record  = line[0:6].strip()
                serial  = int(line[6:11])
                name    = line[12:16].strip()
                chain   = line[21] if len(line) > 21 else " "
                resseq  = int(line[22:26])
                x       = float(line[30:38])
                y       = float(line[38:46])
                z       = float(line[46:54])
                element = line[76:78].strip() if len(line) >= 78 else ""
            except (ValueError, IndexError):
                continue

            if chains is not None and chain not in chains:
                continue

            # Deduplicate
            key = (chain, resseq, name)
            if key in seen:
                continue
            seen[key] = True

            atoms.append(_Atom(
                record=record, serial=serial, name=name, resname=resname,
                chain=chain, resseq=resseq, x=x, y=y, z=z, element=element,
            ))

    return atoms


def _infer_chains(pdb_path: str) -> List[str]:
    """Return all unique chain IDs present in the PDB file."""
    chains = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                chains.add(line[21])
    return sorted(chains)


# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

def _dist(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def _build_coords_array(atoms: List[_Atom]) -> np.ndarray:
    return np.array([a.xyz for a in atoms])


def _find_interface_pairs(
    rec_atoms: List[_Atom],
    lig_atoms: List[_Atom],
    cutoff: float = INTERFACE_CUTOFF,
) -> List[Tuple[_Atom, _Atom, float]]:
    """
    Find all receptor-ligand atom pairs within *cutoff* Å.

    Returns list of (rec_atom, lig_atom, distance) tuples.
    """
    if not rec_atoms or not lig_atoms:
        return []

    rec_coords = _build_coords_array(rec_atoms)  # (N_rec, 3)
    lig_coords = _build_coords_array(lig_atoms)  # (N_lig, 3)

    # Pairwise distances via broadcasting
    diff = rec_coords[:, np.newaxis, :] - lig_coords[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=-1))  # (N_rec, N_lig)

    pairs = []
    rec_i, lig_j = np.where(dists <= cutoff)
    for ri, lj in zip(rec_i, lig_j):
        pairs.append((rec_atoms[ri], lig_atoms[lj], float(dists[ri, lj])))

    return pairs


# ─────────────────────────────────────────────────────────────────────────────
# Metric computations
# ─────────────────────────────────────────────────────────────────────────────

def _compute_bsa(pairs: List[Tuple[_Atom, _Atom, float]]) -> float:
    """
    Approximate Buried Surface Area (Å²).

    Counts unique receptor and ligand atoms at the interface and
    multiplies by SA_PER_ATOM per atom.
    """
    rec_atoms_at_iface = {id(r) for r, l, d in pairs}
    lig_atoms_at_iface = {id(l) for r, l, d in pairs}
    n_interface_atoms  = len(rec_atoms_at_iface) + len(lig_atoms_at_iface)
    return n_interface_atoms * SA_PER_ATOM


def _compute_shape_complementarity(
    pairs: List[Tuple[_Atom, _Atom, float]]
) -> Dict:
    """
    Shape complementarity based on pairwise distance distribution.

    Categories:
      - clashing:      d < CLASH_CUTOFF (3.5 Å)
      - complementary: CLASH_CUTOFF ≤ d ≤ COMPLEMENTARY_MAX (3.5–5.0 Å)
      - no_contact:    d > COMPLEMENTARY_MAX (>5.0 Å)  [should be 0 since we filtered]

    Returns:
      n_clashing, n_complementary, n_no_contact, sc_fraction
    """
    if not pairs:
        return {"n_clashing": 0, "n_complementary": 0, "n_no_contact": 0,
                "sc_fraction": 0.0}

    n_clash  = sum(1 for r, l, d in pairs if d < CLASH_CUTOFF)
    n_comp   = sum(1 for r, l, d in pairs if CLASH_CUTOFF <= d <= COMPLEMENTARY_MAX)
    n_no     = sum(1 for r, l, d in pairs if d > COMPLEMENTARY_MAX)
    total    = len(pairs)
    sc_frac  = n_comp / total if total > 0 else 0.0

    return {
        "n_clashing":      n_clash,
        "n_complementary": n_comp,
        "n_no_contact":    n_no,
        "sc_fraction":     sc_frac,
    }


def _compute_hbonds(pairs: List[Tuple[_Atom, _Atom, float]]) -> int:
    """
    Count potential hydrogen bonds: N or O on one side, N or O on the other,
    with distance in [HBOND_MIN, HBOND_MAX] = [2.5, 3.5] Å.
    """
    count = 0
    for rec_atom, lig_atom, d in pairs:
        if HBOND_MIN <= d <= HBOND_MAX:
            rec_elem = rec_atom.element_sym
            lig_elem = lig_atom.element_sym
            if rec_elem in HBOND_DONORS and lig_elem in HBOND_DONORS:
                count += 1
    return count


def _ca_pairs(
    rec_atoms: List[_Atom],
    lig_atoms: List[_Atom],
    cutoff: float,
) -> List[Tuple[_Atom, _Atom, float]]:
    """
    Return all Cα–Cα pairs between receptor and ligand within *cutoff* Å.
    """
    rec_ca = [a for a in rec_atoms if a.name == "CA"]
    lig_ca = [a for a in lig_atoms if a.name == "CA"]

    if not rec_ca or not lig_ca:
        return []

    pairs = []
    rec_coords = _build_coords_array(rec_ca)
    lig_coords = _build_coords_array(lig_ca)
    diff  = rec_coords[:, np.newaxis, :] - lig_coords[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=-1))
    ri_arr, lj_arr = np.where(dists <= cutoff)
    for ri, lj in zip(ri_arr, lj_arr):
        pairs.append((rec_ca[ri], lig_ca[lj], float(dists[ri, lj])))

    return pairs


def _compute_salt_bridges(
    rec_atoms: List[_Atom],
    lig_atoms: List[_Atom],
) -> int:
    """
    Count salt bridges: oppositely charged residue Cα pairs within
    SALT_BRIDGE_CUTOFF (4.0 Å).
    """
    ca_pairs = _ca_pairs(rec_atoms, lig_atoms, SALT_BRIDGE_CUTOFF)
    count = 0
    for rec_ca, lig_ca, d in ca_pairs:
        rec_pos = rec_ca.resname.upper() in POSITIVE_RESIDUES
        rec_neg = rec_ca.resname.upper() in NEGATIVE_RESIDUES
        lig_pos = lig_ca.resname.upper() in POSITIVE_RESIDUES
        lig_neg = lig_ca.resname.upper() in NEGATIVE_RESIDUES
        if (rec_pos and lig_neg) or (rec_neg and lig_pos):
            count += 1
    return count


def _compute_hydrophobic(
    rec_atoms: List[_Atom],
    lig_atoms: List[_Atom],
) -> int:
    """
    Count hydrophobic contacts: Cα pairs of hydrophobic residues within
    HYDROPHOBIC_CUTOFF (6.0 Å).
    """
    ca_pairs = _ca_pairs(rec_atoms, lig_atoms, HYDROPHOBIC_CUTOFF)
    count = 0
    for rec_ca, lig_ca, d in ca_pairs:
        if (rec_ca.resname.upper() in HYDROPHOBIC_RESIDUES and
                lig_ca.resname.upper() in HYDROPHOBIC_RESIDUES):
            count += 1
    return count


def _interface_residues(
    rec_atoms: List[_Atom],
    lig_atoms: List[_Atom],
    pairs: List[Tuple[_Atom, _Atom, float]],
) -> Dict:
    """
    Identify residues at the interface (any atom within INTERFACE_CUTOFF).

    Returns dicts of receptor and ligand interface residue identifiers.
    """
    rec_iface_atoms = {id(r) for r, l, d in pairs}
    lig_iface_atoms = {id(l) for r, l, d in pairs}

    rec_residues = sorted({
        (a.chain, a.resseq, a.resname)
        for a in rec_atoms
        if id(a) in rec_iface_atoms
    })
    lig_residues = sorted({
        (a.chain, a.resseq, a.resname)
        for a in lig_atoms
        if id(a) in lig_iface_atoms
    })

    return {
        "receptor": [{"chain": c, "resnum": r, "resname": n}
                     for c, r, n in rec_residues],
        "ligand":   [{"chain": c, "resnum": r, "resname": n}
                     for c, r, n in lig_residues],
    }


def _composite_score(
    bsa: float,
    sc: Dict,
    hbonds: int,
    salt_bridges: int,
    hydrophobic: int,
) -> float:
    """
    Compute a weighted composite interface score (0–100).

    Higher is better. Score components:
      - hbonds:           3 pts each
      - salt_bridges:     5 pts each
      - hydrophobic:      0.5 pts each
      - shape_complement: sc_fraction × 30
      - bsa_bonus:        bsa × 0.002
    """
    raw = (
        hbonds       * _WEIGHTS["hbonds"]
        + salt_bridges  * _WEIGHTS["salt_bridges"]
        + hydrophobic   * _WEIGHTS["hydrophobic"]
        + sc["sc_fraction"] * _WEIGHTS["shape_complement"]
        + bsa           * _WEIGHTS["bsa_bonus"]
    )
    return round(min(raw, _SCORE_CAP), 2)


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def score_interface(
    receptor_pdb: str,
    ligand_pdb: str,
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
) -> Dict:
    """
    Score the binding interface between a receptor and ligand PDB.

    Computes geometric and pseudo-energetic metrics at the interface.
    Both receptor and ligand can be the same PDB file (multi-chain complex)
    or separate PDB files.

    Args:
        receptor_pdb:   Path to receptor PDB (typically the protein target).
        ligand_pdb:     Path to ligand PDB (typically the D-peptide binder).
        receptor_chain: Single chain ID to use from receptor (None = all chains).
        ligand_chain:   Single chain ID to use from ligand (None = all chains).

    Returns:
        dict with keys:
          - bsa (float): Buried surface area estimate (Å²)
          - shape_complementarity (dict): sc_fraction, n_clashing, n_complementary
          - hbonds (int): Number of putative hydrogen bonds
          - salt_bridges (int): Number of putative salt bridges
          - hydrophobic_contacts (int): Number of hydrophobic Cα contacts
          - interface_residues (dict): {receptor: [...], ligand: [...]}
          - n_interface_pairs (int): Total atom pairs within 5 Å
          - interface_score (float): Weighted composite score 0–100
          - receptor_pdb (str): Input receptor path
          - ligand_pdb (str): Input ligand path
    """
    for path in (receptor_pdb, ligand_pdb):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"PDB not found: {path}")

    rec_chains = [receptor_chain] if receptor_chain else None
    lig_chains = [ligand_chain]   if ligand_chain   else None

    rec_atoms = _parse_atoms(receptor_pdb, chains=rec_chains)
    lig_atoms = _parse_atoms(ligand_pdb,   chains=lig_chains)

    if not rec_atoms:
        raise ValueError(
            f"No atoms found in receptor PDB "
            f"(chain={receptor_chain}): {receptor_pdb}"
        )
    if not lig_atoms:
        raise ValueError(
            f"No atoms found in ligand PDB "
            f"(chain={ligand_chain}): {ligand_pdb}"
        )

    # Core: find all interface atom pairs within 5 Å
    pairs = _find_interface_pairs(rec_atoms, lig_atoms, cutoff=INTERFACE_CUTOFF)

    bsa            = _compute_bsa(pairs)
    sc             = _compute_shape_complementarity(pairs)
    hbonds         = _compute_hbonds(pairs)
    salt_bridges   = _compute_salt_bridges(rec_atoms, lig_atoms)
    hydrophobic    = _compute_hydrophobic(rec_atoms, lig_atoms)
    iface_residues = _interface_residues(rec_atoms, lig_atoms, pairs)
    score          = _composite_score(bsa, sc, hbonds, salt_bridges, hydrophobic)

    return {
        "bsa":                    round(bsa, 1),
        "shape_complementarity":  sc,
        "hbonds":                 hbonds,
        "salt_bridges":           salt_bridges,
        "hydrophobic_contacts":   hydrophobic,
        "interface_residues":     iface_residues,
        "n_interface_pairs":      len(pairs),
        "interface_score":        score,
        "receptor_pdb":           receptor_pdb,
        "ligand_pdb":             ligand_pdb,
    }


def compare_interfaces(
    complexes: List[Tuple[str, str, str]],
) -> List[Dict]:
    """
    Score and compare multiple receptor:ligand complexes.

    Args:
        complexes: List of (receptor_pdb, ligand_pdb, label) tuples.

    Returns:
        List of dicts, each containing the label and all interface metrics,
        sorted by interface_score descending.
    """
    results = []
    for receptor_pdb, ligand_pdb, label in complexes:
        try:
            metrics = score_interface(receptor_pdb, ligand_pdb)
            metrics["label"] = label
        except Exception as exc:
            metrics = {
                "label":           label,
                "receptor_pdb":    receptor_pdb,
                "ligand_pdb":      ligand_pdb,
                "interface_score": 0.0,
                "error":           str(exc),
            }
        results.append(metrics)

    results.sort(key=lambda d: d.get("interface_score", 0.0), reverse=True)
    return results


def format_comparison_table(results: List[Dict]) -> str:
    """
    Format a compare_interfaces result as a human-readable table.

    Args:
        results: Output from compare_interfaces().

    Returns:
        Multi-line string with aligned columns.
    """
    if not results:
        return "(no results)"

    header = (
        f"{'Rank':<5} {'Label':<30} {'Score':>7} {'BSA(Å²)':>9} "
        f"{'HBonds':>7} {'SaltBr':>7} {'Hydro':>6} {'SC':>6}"
    )
    sep    = "-" * len(header)
    lines  = [header, sep]

    for i, r in enumerate(results, 1):
        if "error" in r:
            lines.append(f"{i:<5} {r['label']:<30}  ERROR: {r['error']}")
            continue
        sc_frac = r.get("shape_complementarity", {}).get("sc_fraction", 0.0)
        lines.append(
            f"{i:<5} {r.get('label', ''):<30} "
            f"{r.get('interface_score', 0.0):>7.1f} "
            f"{r.get('bsa', 0.0):>9.0f} "
            f"{r.get('hbonds', 0):>7d} "
            f"{r.get('salt_bridges', 0):>7d} "
            f"{r.get('hydrophobic_contacts', 0):>6d} "
            f"{sc_frac:>6.3f}"
        )

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# CLI / test block
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    import tempfile
    import textwrap

    # Minimal two-chain PDB: chain A = receptor (L-ALA), chain B = ligand (D-ALA)
    # Placed ~4.5 Å apart — should produce an interface with complementary contacts.
    TWO_CHAIN_PDB = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1       2.100   1.200   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1       1.800   2.300   0.000  1.00  0.00           O
        ATOM      5  CB  ALA A   1       1.900  -0.800   1.200  1.00  0.00           C
        HETATM    6  N   DAL B   1       5.000   0.000   0.000  1.00  0.00           N
        HETATM    7  CA  DAL B   1       6.458   0.000   0.000  1.00  0.00           C
        HETATM    8  C   DAL B   1       7.100   1.200   0.000  1.00  0.00           C
        HETATM    9  O   DAL B   1       6.800   2.300   0.000  1.00  0.00           O
        HETATM   10  CB  DAL B   1       6.900  -0.800  -1.200  1.00  0.00           C
        END
    """)

    print("=" * 60)
    print("ChiralFold — Interface Scorer Self-Test")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb = os.path.join(tmpdir, "complex.pdb")
        with open(pdb, "w") as f:
            f.write(TWO_CHAIN_PDB)

        print("\n[score_interface] receptor chain A, ligand chain B")
        result = score_interface(pdb, pdb, receptor_chain="A", ligand_chain="B")
        for k, v in result.items():
            if k not in ("interface_residues",):
                print(f"  {k}: {v}")
        print(f"  interface_residues.receptor: "
              f"{result['interface_residues']['receptor']}")
        print(f"  interface_residues.ligand:   "
              f"{result['interface_residues']['ligand']}")

        print("\n[compare_interfaces] two identical entries with different labels")
        comp = compare_interfaces([
            (pdb, pdb, "complex_A"),
            (pdb, pdb, "complex_B"),
        ])
        # Override chains for comparison since both are full-file
        print(format_comparison_table(comp))

    if len(sys.argv) >= 3:
        rec = sys.argv[1]
        lig = sys.argv[2]
        rec_ch = sys.argv[3] if len(sys.argv) > 3 else None
        lig_ch = sys.argv[4] if len(sys.argv) > 4 else None
        print(f"\nScoring: {rec} (receptor) vs {lig} (ligand)")
        res = score_interface(rec, lig, receptor_chain=rec_ch, ligand_chain=lig_ch)
        for k, v in res.items():
            print(f"  {k}: {v}")
