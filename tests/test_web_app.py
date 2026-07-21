"""Comprehensive tests for the ChiralFold web application."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

gradio = pytest.importorskip("gradio", exc_type=ImportError)

from web import app as web_app  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "af3_correction"
TOY_PDB = (
    Path(__file__).resolve().parents[1]
    / "chiralfold"
    / "data"
    / "examples"
    / "toy_ubiquitin_fragment.pdb"
)
INVERTED = FIXTURE_DIR / "synthetic_l_ala_inverted.pdb"
CORRECT = FIXTURE_DIR / "synthetic_l_ala_correct.pdb"


@pytest.fixture
def toy_pdb(tmp_path: Path) -> str:
    dest = tmp_path / "toy.pdb"
    dest.write_text(TOY_PDB.read_text())
    return str(dest)


@pytest.fixture
def inverted_pdb(tmp_path: Path) -> str:
    dest = tmp_path / "inverted.pdb"
    dest.write_text(INVERTED.read_text())
    return str(dest)


@pytest.fixture
def correct_pdb(tmp_path: Path) -> str:
    dest = tmp_path / "correct.pdb"
    dest.write_text(CORRECT.read_text())
    return str(dest)


def test_build_app():
    demo = web_app.build_app()
    assert demo is not None
    assert hasattr(demo, "launch")


def test_as_path_rejects_none():
    with pytest.raises(Exception):
        web_app._as_path(None)


def test_as_path_rejects_non_pdb(tmp_path: Path):
    bad = tmp_path / "notes.txt"
    bad.write_text("not a pdb")
    with pytest.raises(Exception):
        web_app._as_path(str(bad))


def test_as_path_rejects_empty(tmp_path: Path):
    empty = tmp_path / "empty.pdb"
    empty.write_text("")
    with pytest.raises(Exception):
        web_app._as_path(str(empty))


def test_as_path_accepts_pdb(toy_pdb: str):
    assert web_app._as_path(toy_pdb) == toy_pdb


def test_run_audit_toy(toy_pdb: str):
    md, path = web_app.run_audit(toy_pdb)
    assert "Stereochemistry audit" in md or "Overall score" in md
    assert path is not None
    assert path.endswith("_audit.json")
    data = json.loads(Path(path).read_text())
    assert "overall_score" in data
    assert "chirality" in data
    assert data["chirality"]["n_wrong"] == 0


def test_run_correct_inverted(inverted_pdb: str):
    md, path = web_app.run_correct(inverted_pdb)
    assert "Chirality correction" in md
    assert path is not None
    assert path.endswith("_corrected.pdb")
    assert Path(path).stat().st_size > 0
    assert "Violations" in md or "Fixed" in md


def test_run_correct_already_ok(correct_pdb: str):
    md, path = web_app.run_correct(correct_pdb)
    assert path is not None
    assert Path(path).is_file()
    # Already-correct structure should report zero residual violations
    assert "0" in md


def test_run_mirror_produces_download(toy_pdb: str):
    md, path = web_app.run_mirror(toy_pdb, "X", True)
    assert "Mirror" in md
    assert path is not None
    assert path.endswith("_mirror.pdb")
    text = Path(path).read_text()
    assert "ATOM" in text
    assert Path(path).stat().st_size > 0


@pytest.mark.parametrize("axis", ["X", "Y", "Z"])
def test_run_mirror_axes(toy_pdb: str, axis: str):
    md, path = web_app.run_mirror(toy_pdb, axis, False)
    assert path is not None
    assert axis in md or axis.lower() in md.lower()


def test_run_audit_no_file():
    with pytest.raises(Exception):
        web_app.run_audit(None)


def test_load_example_toy():
    path, msg = web_app.load_example("toy")
    assert path is not None
    assert Path(path).is_file()
    assert "ubiquitin" in msg.lower() or "toy" in msg.lower()


def test_load_example_inverted():
    path, msg = web_app.load_example("inverted")
    assert path is not None
    assert Path(path).is_file()
    assert Path(path).name == "synthetic_l_ala_inverted.pdb"


def test_load_example_unknown():
    path, msg = web_app.load_example("nope")
    assert path is None
    assert "not found" in msg.lower() or "unknown" in msg.lower()


def test_fetch_pdb_id_invalid():
    path, msg = web_app.fetch_pdb_id("BAD")
    assert path is None
    assert len(msg) > 0


def test_fetch_pdb_id_empty():
    path, msg = web_app.fetch_pdb_id("")
    assert path is None


def test_kpi_html_structure():
    html = web_app._kpi_html([("Score", "90", "ok"), ("Violations", "0", "ok")])
    assert "cf-kpi" in html
    assert "Score" in html
    assert "90" in html


def test_stem_and_copy(tmp_path: Path, toy_pdb: str):
    stem = web_app._stem_from_path(toy_pdb)
    assert stem
    named = web_app._copy_with_name(toy_pdb, stem, "_mirror.pdb")
    assert named.endswith("_mirror.pdb")
    assert Path(named).is_file()


@pytest.mark.slow
def test_fetch_pdb_id_1crn():
    path, msg = web_app.fetch_pdb_id("1CRN")
    assert path is not None, msg
    assert Path(path).is_file()
    assert "1CRN" in msg or "1crn" in msg.lower()
