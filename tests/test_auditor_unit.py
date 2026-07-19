from pathlib import Path

import pytest

from chiralfold.auditor import audit_pdb


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
