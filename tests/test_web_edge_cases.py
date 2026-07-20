"""Edge-case tests for web upload validation and helpers."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("gradio")

from web import app as web_app


def test_max_upload_rejects_oversized(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(web_app, "MAX_UPLOAD_BYTES", 64)
    big = tmp_path / "big.pdb"
    big.write_bytes(b"ATOM" + b"x" * 200)
    with pytest.raises(Exception):
        web_app._as_path(str(big))


def test_ent_extension_accepted(tmp_path: Path):
    ent = tmp_path / "struct.ent"
    ent.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
    assert web_app._as_path(str(ent)) == str(ent)


def test_stem_strips_suffixes():
    assert "foo" in web_app._stem_from_path("/tmp/foo_corrected.pdb")
    assert web_app._stem_from_path("/tmp/bar_mirror.pdb")


def test_run_correct_rejects_missing():
    with pytest.raises(Exception):
        web_app.run_correct(None)


def test_run_mirror_rejects_missing():
    with pytest.raises(Exception):
        web_app.run_mirror(None, "X", True)
