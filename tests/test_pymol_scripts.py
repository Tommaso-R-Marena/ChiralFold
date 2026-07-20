from pathlib import Path

from chiralfold.af3_correct import correct_af3_output, detect_chirality_violations
from chiralfold.viz import (
    write_af3_correction_session,
    write_chirality_audit_session,
    write_mirror_comparison_session,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "af3_correction"


def test_write_chirality_audit_session(tmp_path):
    pdb = FIXTURE_DIR / "synthetic_l_ala_inverted.pdb"
    report = detect_chirality_violations(str(pdb))
    out = tmp_path / "audit.pml"

    write_chirality_audit_session(pdb, report, out)

    text = out.read_text()
    assert "color forest" in text
    assert "color red, chirality_wrong" in text
    assert "chain A and resi 1" in text


def test_write_mirror_comparison_session(tmp_path):
    l_pdb = FIXTURE_DIR / "synthetic_l_ala_correct.pdb"
    d_pdb = FIXTURE_DIR / "synthetic_d_ala_inverted.pdb"
    out = tmp_path / "mirror.pml"

    write_mirror_comparison_session(l_pdb, d_pdb, out)

    text = out.read_text()
    assert "translate [35, 0, 0]" in text
    assert "color marine" in text
    assert "color magenta" in text


def test_write_af3_correction_session(tmp_path):
    before = FIXTURE_DIR / "synthetic_l_ala_inverted.pdb"
    after = tmp_path / "after.pdb"
    correct_af3_output(str(before), str(after))
    out = tmp_path / "correction.pml"

    write_af3_correction_session(before, after, out)

    text = out.read_text()
    assert "align" in text
    assert "corrected_before" in text
    assert "corrected_after" in text
    assert "chain A and resi 1" in text
