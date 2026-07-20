"""Smoke test for web app (optional gradio dependency)."""
import pytest

gradio = pytest.importorskip("gradio")


def test_build_app():
    from web.app import build_app

    demo = build_app()
    assert demo is not None
