"""Tests for the mmCIF -> PDB chain-ID remapping in the Ramachandran benchmark.

These guard the fix for a bug in ``_cif_to_pdb`` where multi-character mmCIF
``auth_asym_id`` values (e.g. "AA", "AB", "BA", "BB") were truncated to a single
character (``[:1]``), collapsing distinct chains onto the same PDB chain ID and
corrupting per-chain residue identity for large / modern structures.
"""

import os
import sys

# Allow importing the standalone benchmark script from ``benchmarks/``.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_BENCH_DIR = os.path.join(_REPO_ROOT, "benchmarks")
if _BENCH_DIR not in sys.path:
    sys.path.insert(0, _BENCH_DIR)

from expand_ramachandran_benchmark import (  # noqa: E402
    _cif_to_pdb,
    PDB_CHAIN_CHARS,
)


def _make_cif(chain_ids):
    """Build a minimal mmCIF ``_atom_site`` loop, one GLY CA per chain.

    Each chain gets a single C-alpha atom at a distinct coordinate so that the
    converter has unambiguous, well-formed input.
    """
    header = [
        "data_TEST",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_atom_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.pdbx_PDB_model_num",
    ]
    rows = []
    for n, ch in enumerate(chain_ids, start=1):
        x = float(n)
        # group id serial element atomname comp labelasym labelseq authseq
        # authasym authcomp authatom x y z occ b icode model
        rows.append(
            f"ATOM {n} C CA GLY {ch} 1 1 {ch} GLY CA "
            f"{x:.3f} {x:.3f} {x:.3f} 1.00 10.00 ? 1"
        )
    return "\n".join(header + rows) + "\n#\n"


def _parse_atom_records(pdb_text):
    """Return list of (chain, resseq, icode) tuples from PDB ATOM/HETATM lines."""
    recs = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            resseq = line[22:26].strip()
            icode = line[26]
            recs.append((chain, resseq, icode))
    return recs


def test_multichar_chains_remapped_uniquely():
    """Four multi-character chains map to four distinct single-char IDs."""
    cif = _make_cif(["AA", "AB", "BA", "BB"])
    pdb_text = _cif_to_pdb(cif)
    assert pdb_text is not None

    recs = _parse_atom_records(pdb_text)
    assert len(recs) == 4, f"expected 4 ATOM records, got {len(recs)}"

    chains = [r[0] for r in recs]
    # No collisions: four distinct, single-character chain IDs.
    assert len(set(chains)) == 4, f"chain collision: {chains}"
    assert all(len(c) == 1 for c in chains)

    # No two ATOM records share the same (chain, resseq) combination.
    chain_res = [(r[0], r[1]) for r in recs]
    assert len(set(chain_res)) == len(chain_res), (
        f"duplicate (chain, resseq) pairs: {chain_res}"
    )

    # Deterministic mapping: first four targets are A, B, C, D.
    assert chains == PDB_CHAIN_CHARS[:4] == ["A", "B", "C", "D"]


def test_exactly_36_chains_ok():
    """A structure with exactly 36 distinct chains is representable."""
    cif = _make_cif([f"X{i}" for i in range(36)])
    pdb_text = _cif_to_pdb(cif)
    assert pdb_text is not None
    chains = [r[0] for r in _parse_atom_records(pdb_text)]
    assert len(set(chains)) == 36


def test_more_than_36_chains_returns_none():
    """A structure with 37 distinct chains cannot be represented -> None."""
    cif = _make_cif([f"Y{i}" for i in range(37)])
    pdb_text = _cif_to_pdb(cif)
    assert pdb_text is None
