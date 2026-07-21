"""Regression tests for MolProbity-like clash exclusions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chiralfold import auditor as auditor_mod
from chiralfold.auditor import audit_pdb


TOY = (
    Path(__file__).resolve().parents[1]
    / "chiralfold"
    / "data"
    / "examples"
    / "toy_ubiquitin_fragment.pdb"
)


def test_leu_ca_cg_is_excluded_as_1_3():
    """CA–CG (~2.6 Å) is covalent 1-3 via CB — must never score as a clash."""
    atoms = auditor_mod._parse_pdb(str(TOY))
    leus = {}
    for a in atoms:
        if a.resname == "LEU" and a.name.strip() in ("CA", "CG", "CB"):
            leus.setdefault(a.resseq, {})[a.name.strip()] = a
    assert leus, "toy fragment should contain LEU"
    for resseq, m in leus.items():
        if "CA" in m and "CG" in m:
            assert auditor_mod._are_bonded_or_angled(m["CA"], m["CG"]) is True


def test_amide_hydrogen_placement_not_on_prev_ca():
    atoms = auditor_mod._parse_pdb(str(TOY))
    with_h = auditor_mod._add_backbone_hydrogens(atoms)
    # No proline should receive a fake amide H
    for a in with_h:
        if a.element_upper == "H" and a.resname in auditor_mod.PROLINE_RESNAMES:
            pytest.fail(f"Proline {a.resseq} must not get an amide H")
    # Previous CA to amide H should be a normal 1-4 distance (~2.0–2.8 Å)
    by = {}
    for a in with_h:
        by.setdefault((a.chain, a.resseq), {})[a.name.strip()] = a
    keys = sorted(by)
    checked = 0
    for (c0, r0), (c1, r1) in zip(keys, keys[1:]):
        if c0 != c1 or r1 - r0 != 1:
            continue
        if "CA" not in by[(c0, r0)] or "H" not in by[(c1, r1)]:
            continue
        dist = float(
            np.linalg.norm(by[(c0, r0)]["CA"].xyz - by[(c1, r1)]["H"].xyz)
        )
        assert dist > 2.0, f"amide H too close to prev CA: {dist}"
        checked += 1
    assert checked >= 1


def test_toy_ubiquitin_clashscore_is_reasonable():
    """After topology exclusions, a clean fragment must not look clash-ridden."""
    report = audit_pdb(str(TOY))
    assert report["clashes"]["clash_score"] < 50.0
    assert report["chirality"]["n_wrong"] == 0


def test_deposited_hydrogens_do_not_inflate_clashscore(tmp_path: Path):
    """Explicit C–H bonds (~1.09 Å) must not be scored as clashes."""
    # Minimal ALA with deposited HA (would be a ~1.8 Å "clash" if kept)
    pdb = tmp_path / "ala_with_h.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.200   1.300   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.700   2.400   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.000  -1.200   0.000  1.00  0.00           C\n"
        "ATOM      6  HA  ALA A   1       1.800   0.300   0.900  1.00  0.00           H\n"
        "END\n"
    )
    report = audit_pdb(str(pdb))
    assert report["clashes"]["n_clashes"] == 0
    assert report["clashes"]["clash_score"] == 0.0


def test_only_first_model_is_parsed(tmp_path: Path):
    pdb = tmp_path / "two_models.pdb"
    pdb.write_text(
        "MODEL        1\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.200   1.300   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.700   2.400   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.000  -1.200   0.000  1.00  0.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      6  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      7  CA  ALA A   1      11.458  10.000  10.000  1.00  0.00           C\n"
        "ATOM      8  C   ALA A   1      12.200  11.300  10.000  1.00  0.00           C\n"
        "ATOM      9  O   ALA A   1      11.700  12.400  10.000  1.00  0.00           O\n"
        "ATOM     10  CB  ALA A   1      12.000   8.800  10.000  1.00  0.00           C\n"
        "ENDMDL\n"
    )
    atoms = auditor_mod._parse_pdb(str(pdb))
    assert len(atoms) == 5
    assert max(a.x for a in atoms) < 5.0


def test_multichain_does_not_build_cross_chain_amide_h(tmp_path: Path):
    """Chain B residue 1 must not get an H from chain A's carbonyl."""
    pdb = tmp_path / "two_chains.pdb"
    # Two well-separated ALA residues on different chains, same resseq.
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.009  -0.771   1.200  1.00  0.00           C\n"
        "ATOM      6  N   ALA B   1      20.000  20.000  20.000  1.00  0.00           N\n"
        "ATOM      7  CA  ALA B   1      21.458  20.000  20.000  1.00  0.00           C\n"
        "ATOM      8  C   ALA B   1      22.009  21.420  20.000  1.00  0.00           C\n"
        "ATOM      9  O   ALA B   1      21.251  22.390  20.000  1.00  0.00           O\n"
        "ATOM     10  CB  ALA B   1      22.009  19.229  21.200  1.00  0.00           C\n"
        "END\n"
    )
    atoms = auditor_mod._parse_pdb(str(pdb))
    with_h = auditor_mod._add_backbone_hydrogens(atoms)
    h_atoms = [a for a in with_h if a.element_upper == "H"]
    assert h_atoms == [], f"unexpected amide H(s): {[(a.chain, a.resseq) for a in h_atoms]}"


def test_disulfide_sg_sg_is_excluded(tmp_path: Path):
    """Covalent SG–SG (~2.0 Å) must not count as a clash."""
    pdb = tmp_path / "cys_ss.pdb"
    # Two CYS with SG placed 2.05 Å apart (typical disulfide)
    pdb.write_text(
        "ATOM      1  N   CYS A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  CYS A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   CYS A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   CYS A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  CYS A   1       2.000  -0.800   1.100  1.00  0.00           C\n"
        "ATOM      6  SG  CYS A   1       1.500  -0.500   2.700  1.00  0.00           S\n"
        "ATOM      7  N   CYS A   2       3.400   1.700   0.000  1.00  0.00           N\n"
        "ATOM      8  CA  CYS A   2       4.100   2.900   0.200  1.00  0.00           C\n"
        "ATOM      9  C   CYS A   2       5.600   2.700   0.200  1.00  0.00           C\n"
        "ATOM     10  O   CYS A   2       6.200   1.700   0.000  1.00  0.00           O\n"
        "ATOM     11  CB  CYS A   2       3.600   3.500   1.400  1.00  0.00           C\n"
        "ATOM     12  SG  CYS A   2       2.200  -0.300   4.400  1.00  0.00           S\n"
        "END\n"
    )
    atoms = auditor_mod._parse_pdb(str(pdb))
    sgs = [a for a in atoms if a.name.strip() == "SG"]
    assert len(sgs) == 2
    dist = float(np.linalg.norm(sgs[0].xyz - sgs[1].xyz))
    assert dist < 2.5, f"fixture SG–SG should be disulfide-like, got {dist}"
    excluded = auditor_mod._clash_excluded_index_pairs(atoms)
    sg_idxs = [i for i, a in enumerate(atoms) if a.name.strip() == "SG"]
    i, j = sg_idxs[0], sg_idxs[1]
    key = (i, j) if i < j else (j, i)
    assert key in excluded
    report = audit_pdb(str(pdb))
    for c in report["clashes"]["worst_clashes"]:
        names = {c["atom1"].split(".")[-1].strip(), c["atom2"].split(".")[-1].strip()}
        assert names != {"SG"}, f"disulfide scored as clash: {c}"


def test_audit_chain_filter_and_mmcifs(tmp_path: Path):
    pdb = tmp_path / "ab.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.009  -0.771   1.200  1.00  0.00           C\n"
        "ATOM      6  N   ALA B   1      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      7  CA  ALA B   1      11.458  10.000  10.000  1.00  0.00           C\n"
        "ATOM      8  C   ALA B   1      12.009  11.420  10.000  1.00  0.00           C\n"
        "ATOM      9  O   ALA B   1      11.251  12.390  10.000  1.00  0.00           O\n"
        "ATOM     10  CB  ALA B   1      12.009   9.229  11.200  1.00  0.00           C\n"
        "END\n"
    )
    only_a = audit_pdb(str(pdb), chain="A")
    assert only_a["chains"] == ["A"]
    assert only_a["n_residues"] == 1
    with pytest.raises(ValueError, match="chain 'Z'"):
        audit_pdb(str(pdb), chain="Z")

    # Minimal mmCIF round-trip through audit_pdb (same schema as fetch tests)
    cif = tmp_path / "ala.cif"
    cif.write_text(
        "data_TEST\n"
        "#\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 N N ALA A 1 1 A ALA N 0.000 0.000 0.000 1.00 0.00 ? 1\n"
        "ATOM 2 C CA ALA A 1 1 A ALA CA 1.458 0.000 0.000 1.00 0.00 ? 1\n"
        "ATOM 3 C C ALA A 1 1 A ALA C 2.009 1.420 0.000 1.00 0.00 ? 1\n"
        "ATOM 4 O O ALA A 1 1 A ALA O 1.251 2.390 0.000 1.00 0.00 ? 1\n"
        "ATOM 5 C CB ALA A 1 1 A ALA CB 2.009 -0.771 -1.200 1.00 0.00 ? 1\n"
        "#\n"
    )
    cif_report = audit_pdb(str(cif))
    assert cif_report["n_residues"] == 1
    assert cif_report["chirality"]["n_correct"] == 1


def test_d_leu_ca_cg_topology_exclusion():
    """D-CCD LEU (DLE) must use the LEU side-chain template for 1–3 exclusion."""
    assert auditor_mod._bond_template_resname("DLE") == "LEU"
    pairs = auditor_mod._intra_residue_bond_pairs("DLE")
    assert ("CA", "CB") in pairs
    assert ("CB", "CG") in pairs
