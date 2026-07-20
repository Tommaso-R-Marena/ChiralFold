"""CLI and packaging smoke tests for chiralfold-web entry point."""
from __future__ import annotations

import importlib
from pathlib import Path

import pytest


def test_web_package_exports_main():
    pytest.importorskip("gradio", exc_type=ImportError)
    mod = importlib.import_module("web.app")
    assert callable(mod.main)
    assert callable(mod.build_app)


def test_pyproject_declares_web_extra():
    text = Path("pyproject.toml").read_text()
    assert 'chiralfold-web = "web.app:main"' in text
    assert "gradio" in text
    assert "[project.optional-dependencies]" in text
