from pathlib import Path
import math

import numpy as np
import pytest

from chiralfold import auditor as auditor_mod
from chiralfold.auditor import audit_pdb, format_report


L_ALA_CORRECT = """\
ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -1.200   2.095   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.000  -0.500   1.200  1.00  0.00           C
END
"""


def test_audit_pdb_minimal_correct_l_ala(tmp_path):
    pdb_path = tmp_path / "l_ala.pdb"
    pdb_path.write_text(L_ALA_CORRECT)

    report = audit_pdb(str(pdb_path))

    assert set(report) >= {
        "n_residues",
        "n_atoms",
        "chains",
        "chirality",
        "bond_geometry",
        "ramachandran",
        "planarity",
        "clashes",
        "overall_score",
    }
    assert set(report["chirality"]) >= {
        "n_correct",
        "n_wrong",
        "n_glycine",
        "n_error",
        "pct_correct",
        "violations",
    }
    assert report["n_residues"] == 1
    assert report["chirality"]["n_correct"] == 1
    assert report["chirality"]["n_wrong"] == 0


def test_audit_pdb_detects_fixture_inverted_l_ala_if_available():
    fixture = Path(__file__).parent / "fixtures" / "af3_correction" / "synthetic_l_ala_inverted.pdb"
    if not fixture.exists():
        pytest.skip("AF3 correction inverted fixture not available")

    report = audit_pdb(str(fixture))

    assert report["chirality"]["n_wrong"] == 1
    assert report["chirality"]["violations"][0]["expected"] == "L"
    assert report["chirality"]["violations"][0]["observed"] == "D"


def test_audit_pdb_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        audit_pdb(str(tmp_path / "missing.pdb"))


def test_audit_pdb_empty_file_raises_value_error(tmp_path):
    empty = tmp_path / "empty.pdb"
    empty.write_text("REMARK no atoms\n")

    with pytest.raises(ValueError, match="No parseable"):
        audit_pdb(str(empty))


def test_format_report_includes_violation_details():
    fixture = Path(__file__).parent / "fixtures" / "af3_correction" / "synthetic_l_ala_inverted.pdb"
    report = audit_pdb(str(fixture))

    text = format_report(report)

    assert "ChiralFold PDB Auditor" in text
    assert "VIOLATION" in text
    assert "expected L, observed D" in text


def test_auditor_parse_filters_water_altloc_duplicates_and_bad_records(tmp_path):
    pdb = tmp_path / "parse_edges.pdb"
    pdb.write_text(
        """\
ATOM      1  N   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  CA AALA A   2       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA BALA A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      4  CA AALA A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM  broken
END
"""
    )

    atoms = auditor_mod._parse_pdb(str(pdb))
    groups = auditor_mod._group_by_residue(atoms)

    assert len(atoms) == 1
    assert atoms[0].altloc == "A"
    assert list(groups) == [("A", 2, " ", "ALA")]


def test_auditor_geometry_and_chain_continuity_helpers():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, 0.0])
    c = np.array([1.0, 0.0, 0.0])

    assert auditor_mod._vec_norm(c) == 1.0
    assert auditor_mod._bond_length(a, c) == 1.0
    assert math.isnan(auditor_mod._bond_angle_deg(a, b, c))
    assert math.isnan(auditor_mod._dihedral_deg(a, b, b, np.array([2.0, 0.0, 0.0])))

    key_a = ("A", 1, " ", "ALA")
    key_b = ("A", 2, " ", "ALA")
    key_gap = ("A", 5, " ", "ALA")
    key_other_chain = ("B", 2, " ", "ALA")
    res_a = {"CA": np.array([0.0, 0.0, 0.0])}
    res_near = {"CA": np.array([3.8, 0.0, 0.0])}
    res_far = {"CA": np.array([6.0, 0.0, 0.0])}

    assert auditor_mod._is_chain_continuous(key_a, key_b, res_a, res_near) is True
    assert auditor_mod._is_chain_continuous(key_a, key_other_chain, res_a, res_near) is False
    assert auditor_mod._is_chain_continuous(key_a, key_gap, res_a, res_near) is False
    assert auditor_mod._is_chain_continuous(key_a, key_b, res_a, res_far) is False


def _audit_atom(serial, name, resname, chain, resseq, x, y, z, record="ATOM", element="C"):
    return auditor_mod._Atom(record, serial, name, " ", resname, chain, resseq, " ", x, y, z, element)


def test_check_chirality_covers_skip_planar_and_unknown_paths():
    groups = {
        ("A", 1, " ", "GLY"): [
            _audit_atom(1, "N", "GLY", "A", 1, 1.0, 0.0, 0.0, element="N"),
            _audit_atom(2, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0),
            _audit_atom(3, "C", "GLY", "A", 1, -1.0, 0.0, 0.0),
        ],
        ("A", 2, " ", "ALA"): [_audit_atom(4, "N", "ALA", "A", 2, 0.0, 0.0, 0.0)],
        ("A", 3, " ", "ALA"): [
            _audit_atom(5, "N", "ALA", "A", 3, 1.0, 0.0, 0.0, element="N"),
            _audit_atom(6, "CA", "ALA", "A", 3, 0.0, 0.0, 0.0),
            _audit_atom(7, "C", "ALA", "A", 3, -1.0, 0.0, 0.0),
        ],
        ("A", 4, " ", "ALA"): [
            _audit_atom(8, "N", "ALA", "A", 4, 1.0, 0.0, 0.0, element="N"),
            _audit_atom(9, "CA", "ALA", "A", 4, 0.0, 0.0, 0.0),
            _audit_atom(10, "C", "ALA", "A", 4, -1.0, 0.0, 0.0),
            _audit_atom(11, "CB", "ALA", "A", 4, 0.0, 1.0, 0.0),
        ],
        ("B", 1, " ", "UNK"): [
            _audit_atom(12, "N", "UNK", "B", 1, 1.201, 0.847, 0.0, "HETATM", "N"),
            _audit_atom(13, "CA", "UNK", "B", 1, 0.0, 0.0, 0.0, "HETATM"),
            _audit_atom(14, "C", "UNK", "B", 1, -1.25, 0.881, 0.0, "HETATM"),
            _audit_atom(15, "CB", "UNK", "B", 1, 0.0, -0.5, -1.2, "HETATM"),
        ],
    }

    report = auditor_mod._check_chirality(groups)

    assert report["n_glycine"] == 1
    assert report["n_error"] == 3
    assert report["n_correct"] == 1


def test_check_bond_geometry_outliers_and_planarity_outlier():
    groups = {
        ("A", 1, " ", "ALA"): [
            _audit_atom(1, "N", "ALA", "A", 1, -2.0, 0.0, 0.0, element="N"),
            _audit_atom(2, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _audit_atom(3, "C", "ALA", "A", 1, 1.0, 0.0, 0.0),
            _audit_atom(4, "O", "ALA", "A", 1, 3.0, 0.0, 0.0, element="O"),
        ],
        ("A", 2, " ", "ALA"): [
            _audit_atom(5, "N", "ALA", "A", 2, 1.0, 1.0, 0.0, element="N"),
            _audit_atom(6, "CA", "ALA", "A", 2, 1.0, 1.0, 1.0),
            _audit_atom(7, "C", "ALA", "A", 2, 2.0, 1.0, 1.0),
        ],
    }
    keys = list(groups)

    bond_geo = auditor_mod._check_bond_geometry(groups, keys)
    planarity = auditor_mod._check_planarity(groups, keys)

    assert bond_geo["outliers"]
    assert planarity["outliers"]


def test_clash_helpers_cover_single_atom_hydrogen_and_default_radius():
    atom = _audit_atom(1, "XX", "UNK", "A", 1, 0.0, 0.0, 0.0, element="")
    assert auditor_mod._vdw_radius(atom) == auditor_mod.VDW_DEFAULT
    assert auditor_mod._check_clashes([atom]) == {
        "n_clashes": 0,
        "clash_score": 0.0,
        "worst_clashes": [],
    }

    atoms = [
        _audit_atom(1, "C", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _audit_atom(2, "N", "ALA", "A", 2, 1.3, 0.0, 0.0, element="N"),
        _audit_atom(3, "CA", "ALA", "A", 2, 2.5, 0.0, 0.0),
    ]
    with_h = auditor_mod._add_backbone_hydrogens(atoms)
    assert len(with_h) == 4
    assert with_h[-1].element == "H"

    same_res_a = _audit_atom(4, "CA", "ALA", "A", 3, 0.0, 0.0, 0.0)
    same_res_b = _audit_atom(5, "CB", "ALA", "A", 3, 1.0, 0.0, 0.0)
    far_atom = _audit_atom(6, "CB", "ALA", "A", 5, 10.0, 0.0, 0.0)
    assert auditor_mod._are_bonded_or_angled(same_res_a, same_res_b) is True
    assert auditor_mod._are_bonded_or_angled(same_res_a, far_atom) is False
