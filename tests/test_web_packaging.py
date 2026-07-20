"""CLI and packaging smoke tests for chiralfold-web entry point."""
from __future__ import annotations

import importlib

import pytest


def test_web_package_exports_main():
    mod = importlib.import_module("web.app")
    assert callable(mod.main)
    assert callable(mod.build_app)


def test_pyproject_declares_web_extra():
    # setuptools may not expose extras easily; check pyproject text
    from pathlib import Path

    text = Path("pyproject.toml").read_text()
    assert 'chiralfold-web = "web.app:main"' in text
    assert "gradio" in text
    assert "[project.optional-dependencies]" in text
