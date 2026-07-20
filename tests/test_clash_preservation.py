"""Prove mirror isometry preserves clashscore; AF3 correction clears violations."""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from chiralfold.af3_correct import correct_af3_output, detect_chirality_violations
from chiralfold.auditor import audit_pdb
from chiralfold.pdb_pipeline import mirror_pdb

TOY = (
    Path(__file__).resolve().parents[1]
    / "chiralfold"
    / "data"
    / "examples"
    / "toy_ubiquitin_fragment.pdb"
)
INVERTED = Path(__file__).parent / "fixtures" / "af3_correction" / "synthetic_l_ala_inverted.pdb"


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_mirror_preserves_clashscore(tmp_path: Path, axis: str):
    """Global reflection is an isometry — clashscore must be identical."""
    src = tmp_path / "toy.pdb"
    src.write_text(TOY.read_text())
    before = audit_pdb(str(src))["clashes"]

    out = tmp_path / f"mirror_{axis}.pdb"
    mirror_pdb(str(src), str(out), axis=axis, rename_residues=True)
    after = audit_pdb(str(out))["clashes"]

    assert after["n_clashes"] == before["n_clashes"]
    assert after["clash_score"] == before["clash_score"]


def test_mirror_preserves_clashscore_without_rename(tmp_path: Path):
    src = tmp_path / "toy.pdb"
    src.write_text(TOY.read_text())
    before = audit_pdb(str(src))["clashes"]
    out = tmp_path / "mirror.pdb"
    mirror_pdb(str(src), str(out), axis="x", rename_residues=False)
    after = audit_pdb(str(out))["clashes"]
    assert after["n_clashes"] == before["n_clashes"]
    assert after["clash_score"] == before["clash_score"]


def test_af3_correction_zero_residual_violations(tmp_path: Path):
    src = tmp_path / "inv.pdb"
    src.write_text(INVERTED.read_text())
    before = detect_chirality_violations(str(src))
    assert before["n_violations"] >= 1

    out = tmp_path / "fixed.pdb"
    result = correct_af3_output(str(src), str(out))
    assert result["after"]["n_violations"] == 0

    # Tiny synthetic systems typically have zero clashes; ensure correction
    # does not create a worse clashscore than a large spike.
    b_clash = audit_pdb(str(src))["clashes"]["clash_score"]
    a_clash = audit_pdb(str(out))["clashes"]["clash_score"]
    assert a_clash <= b_clash + 50.0  # generous bound for 1-residue toy


def test_double_mirror_restores_geometry(tmp_path: Path):
    """Mirror twice returns to the original chirality class."""
    src = tmp_path / "toy.pdb"
    src.write_text(TOY.read_text())
    mid = tmp_path / "m1.pdb"
    back = tmp_path / "m2.pdb"
    mirror_pdb(str(src), str(mid), axis="x", rename_residues=False)
    mirror_pdb(str(mid), str(back), axis="x", rename_residues=False)
    a0 = audit_pdb(str(src))
    a2 = audit_pdb(str(back))
    assert a0["clashes"]["n_clashes"] == a2["clashes"]["n_clashes"]
    assert a0["chirality"]["n_wrong"] == a2["chirality"]["n_wrong"]
