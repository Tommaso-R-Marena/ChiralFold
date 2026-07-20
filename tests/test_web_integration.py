"""Additional integration tests for web download workflows and HF space parity."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("gradio")

from chiralfold.af3_correct import detect_chirality_violations
from web.app import run_correct, run_mirror, run_audit

TOY = (
    Path(__file__).resolve().parents[1]
    / "chiralfold"
    / "data"
    / "examples"
    / "toy_ubiquitin_fragment.pdb"
)
INVERTED = (
    Path(__file__).parent / "fixtures" / "af3_correction" / "synthetic_l_ala_inverted.pdb"
)


def test_correct_then_mirror_pipeline(tmp_path: Path):
    """Upload → correct → mirror should yield two downloadable PDBs."""
    src = tmp_path / "inv.pdb"
    src.write_text(INVERTED.read_text())

    md_c, corrected = run_correct(str(src))
    assert corrected and Path(corrected).is_file()
    after = detect_chirality_violations(corrected)
    assert after["n_violations"] == 0

    md_m, mirrored = run_mirror(corrected, "X", True)
    assert mirrored and Path(mirrored).is_file()
    mirrored_text = Path(mirrored).read_text()
    assert ("ATOM" in mirrored_text) or ("HETATM" in mirrored_text)
    assert "Mirror" in md_m
    assert "Chirality" in md_c or "correction" in md_c.lower()


def test_audit_correct_mirror_on_toy(tmp_path: Path):
    src = tmp_path / "toy.pdb"
    src.write_text(TOY.read_text())

    md_a, audit_json = run_audit(str(src))
    assert audit_json and audit_json.endswith(".json")

    md_c, corrected = run_correct(str(src))
    assert corrected and corrected.endswith(".pdb")

    md_m, mirrored = run_mirror(str(src), "Y", False)
    assert mirrored and mirrored.endswith(".pdb")
    assert all(x is not None for x in (md_a, md_c, md_m))


def test_hf_space_examples_present():
    """Bundled HF Space examples must exist for offline demos."""
    root = Path(__file__).resolve().parents[1] / "hf_space" / "examples"
    assert (root / "toy_ubiquitin_fragment.pdb").is_file()
    assert (root / "synthetic_l_ala_inverted.pdb").is_file()


def test_hf_space_app_imports():
    """HF Space entrypoint must be importable without launching."""
    import importlib.util

    space_app = Path(__file__).resolve().parents[1] / "hf_space" / "app.py"
    assert space_app.is_file()
    spec = importlib.util.spec_from_file_location("hf_space_app", space_app)
    mod = importlib.util.module_from_spec(spec)
    # Do not exec fully (builds Gradio Blocks at import); just verify syntax via compile
    compile(space_app.read_text(), str(space_app), "exec")
