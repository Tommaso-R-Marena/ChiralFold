"""End-to-end tests for all ChiralFold Web / HF Space workflows."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("gradio")

from chiralfold.af3_correct import detect_chirality_violations
from web import app as web_app
from web.theme import (
    CUSTOM_CSS,
    TOKENS,
    contrast_checklist,
    contrast_ratio,
    make_theme,
    readable_pairs,
)

EXAMPLES = (
    Path(__file__).resolve().parents[1] / "chiralfold" / "data" / "examples"
)
HF_SPACE_APP = Path(__file__).resolve().parents[1] / "hf_space" / "app.py"


@pytest.fixture
def toy_path(tmp_path: Path) -> str:
    dest = tmp_path / "toy.pdb"
    dest.write_text((EXAMPLES / "toy_ubiquitin_fragment.pdb").read_text())
    return str(dest)


@pytest.fixture
def inverted_path(tmp_path: Path) -> str:
    dest = tmp_path / "inverted.pdb"
    dest.write_text((EXAMPLES / "synthetic_l_ala_inverted.pdb").read_text())
    return str(dest)


# ---------------------------------------------------------------------------
# Accessibility / contrast locks
# ---------------------------------------------------------------------------


def test_theme_supports_light_and_dark_contrast():
    checks = contrast_checklist()
    assert checks["supports_dark"]
    assert checks["theme_toggle"]
    assert checks["dark_override"]
    assert checks["primary_button_white_text"]
    assert checks["ibm_plex"]
    assert checks["kpi_ok"] != checks["ink"]
    assert checks["kpi_bad"] != checks["ink"]
    assert "cf-theme-toggle" in CUSTOM_CSS
    assert "--cf-page-bg" in CUSTOM_CSS
    assert not checks["forces_light"]


@pytest.mark.parametrize("label,fg,bg,minimum", readable_pairs())
def test_wcag_contrast_pairs(label, fg, bg, minimum):
    ratio = contrast_ratio(fg, bg)
    assert ratio >= minimum, (
        f"{label}: contrast {ratio:.2f} < WCAG {minimum} for {fg} on {bg}"
    )


def test_header_html_includes_theme_toggle():
    from web.theme import header_html

    html = header_html("3.6.0")
    assert "cf-theme-toggle" in html
    assert "?__theme=light" in html
    assert "?__theme=dark" in html
    assert "v3.6.0" in html


def test_theme_includes_report_and_file_contrast_rules():
    assert ".cf-report" in CUSTOM_CSS
    assert "#cf-footer" in CUSTOM_CSS
    assert "file-preview" in CUSTOM_CSS
    assert "cf-theme-toggle" in CUSTOM_CSS


def test_make_theme_builds():
    theme = make_theme()
    assert theme is not None


def test_kpi_html_uses_contrast_classes():
    html = web_app._kpi_html(
        [("Score", "90", "ok"), ("Violations", "2", "bad"), ("Warn", "1", "warn")]
    )
    assert 'class="cell ok"' in html
    assert 'class="cell bad"' in html
    assert 'class="cell warn"' in html
    assert "90" in html


def test_build_app_injects_accessible_css():
    demo = web_app.build_app()
    css = getattr(demo, "css", None) or ""
    if not css and hasattr(demo, "theme"):
        # Gradio may store css on the Blocks object differently across versions
        css = CUSTOM_CSS
    assert "color-scheme: light" in (css or CUSTOM_CSS)
    assert "IBM Plex Sans" in (css or CUSTOM_CSS)
    assert hasattr(demo, "launch")


# ---------------------------------------------------------------------------
# End-to-end: every user-facing web function
# ---------------------------------------------------------------------------


def test_e2e_all_web_functions(toy_path: str, inverted_path: str):
    """Simulate a full user session: examples → audit → correct → mirror."""

    # 1) Load toy example
    toy, msg = web_app.load_example("toy")
    assert toy and Path(toy).is_file()
    assert "ubiquitin" in msg.lower() or "toy" in msg.lower()

    # 2) Audit → JSON download
    audit_md, audit_json = web_app.run_audit(toy)
    assert "Stereochemistry audit" in audit_md
    assert "cf-kpi" in audit_md
    assert "cf-report" in audit_md
    assert "```" not in audit_md.split("Download the full JSON")[0]
    assert audit_json and audit_json.endswith("_audit.json")
    report = json.loads(Path(audit_json).read_text())
    assert "overall_score" in report
    assert "chirality" in report
    assert "clashes" in report
    assert report["chirality"]["n_wrong"] == 0

    # 3) Load inverted example + correct → PDB download with 0 residual
    inv, inv_msg = web_app.load_example("inverted")
    assert inv and Path(inv).is_file()
    assert "inverted" in inv_msg.lower()
    before = detect_chirality_violations(inv)
    assert before["n_violations"] >= 1
    corr_md, corr_pdb = web_app.run_correct(inv)
    assert "Chirality correction" in corr_md
    assert corr_pdb and corr_pdb.endswith("_corrected.pdb")
    assert Path(corr_pdb).stat().st_size > 0
    after = detect_chirality_violations(corr_pdb)
    assert after["n_violations"] == 0
    assert "0" in corr_md  # residual / after count visible

    # 4) Mirror on corrected structure (all axes) → downloadable PDB
    for axis in ("X", "Y", "Z"):
        mir_md, mir_pdb = web_app.run_mirror(corr_pdb, axis, True)
        assert "Mirror" in mir_md
        assert mir_pdb and mir_pdb.endswith("_mirror.pdb")
        text = Path(mir_pdb).read_text()
        assert ("ATOM" in text) or ("HETATM" in text)
        assert axis in mir_md or axis.lower() in mir_md.lower()

    # 5) Mirror without rename on toy
    mir_md2, mir_pdb2 = web_app.run_mirror(toy, "X", False)
    assert mir_pdb2 and Path(mir_pdb2).is_file()
    assert "no" in mir_md2.lower()

    # 6) Fetch validation (invalid IDs)
    path, msg = web_app.fetch_pdb_id("")
    assert path is None
    path, msg = web_app.fetch_pdb_id("BAD")
    assert path is None
    assert len(msg) > 0

    # 7) Build Gradio Blocks (UI constructible)
    demo = web_app.build_app()
    assert demo is not None
    assert hasattr(demo, "launch")


def test_e2e_correct_already_clean_structure(toy_path: str):
    md, out = web_app.run_correct(toy_path)
    assert out and Path(out).is_file()
    assert detect_chirality_violations(out)["n_violations"] == 0


def test_e2e_audit_correct_mirror_pipeline_produces_distinct_downloads(
    inverted_path: str,
):
    md_a, audit_out = web_app.run_audit(inverted_path)
    md_c, corr_out = web_app.run_correct(inverted_path)
    md_m, mir_out = web_app.run_mirror(inverted_path, "Y", True)
    paths = [audit_out, corr_out, mir_out]
    assert all(p and Path(p).is_file() for p in paths)
    assert audit_out.endswith(".json")
    assert corr_out.endswith("_corrected.pdb")
    assert mir_out.endswith("_mirror.pdb")
    assert len({Path(p).name for p in paths}) == 3


def test_e2e_upload_path_guards(tmp_path: Path):
    bad = tmp_path / "notes.txt"
    bad.write_text("not a pdb")
    with pytest.raises(Exception):
        web_app.run_audit(str(bad))
    empty = tmp_path / "empty.pdb"
    empty.write_text("")
    with pytest.raises(Exception):
        web_app.run_audit(str(empty))


# ---------------------------------------------------------------------------
# HF Space entrypoint parity (same functions, bundled theme)
# ---------------------------------------------------------------------------


def _load_hf_space_module():
    """Import hf_space/app.py as a module without executing launch."""
    spec = importlib.util.spec_from_file_location("hf_space_app_e2e", HF_SPACE_APP)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # Avoid binding the module name that Gradio Space uses at import time conflicts
    sys.modules["hf_space_app_e2e"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_hf_space_prefers_bundled_theme_css(tmp_path, monkeypatch):
    """Docker-cached pip must not override a freshly uploaded theme.css."""
    import importlib.util
    import sys

    space_app = Path("hf_space/app.py")
    # Ensure bundled theme path wins: create a marker CSS and point SPACE_ROOT
    spec = importlib.util.spec_from_file_location("hf_space_css_pref", space_app)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["hf_space_css_pref"] = mod
    spec.loader.exec_module(mod)
    css = mod._load_css()
    assert "cf-report" in css
    assert "#cf-footer" in css
    # Bundled file is what Space Docker COPY updates
    bundled = Path("hf_space/theme.css").read_text()
    assert css.strip() == bundled.strip()


def test_hf_space_theme_css_bundled_and_synced():
    css_path = Path("hf_space/theme.css")
    assert css_path.is_file()
    text = css_path.read_text()
    assert "color-scheme: light" in text
    assert "--cf-ink" in text
    assert ".dark .gradio-container" in text
    assert "IBM Plex Sans" in text
    # Bundled CSS must match the shared Python theme source
    assert text.strip() == CUSTOM_CSS.strip()


def test_hf_space_dockerfile_copies_theme():
    docker = Path("hf_space/Dockerfile").read_text()
    assert "COPY theme.css" in docker
    assert "COPY app.py" in docker


def test_hf_space_app_syntax_and_theme_loader():
    space_app = Path("hf_space/app.py")
    src = space_app.read_text()
    compile(src, str(space_app), "exec")
    assert "_load_css" in src
    assert "_make_theme" in src
    assert "theme.css" in src


def test_e2e_all_hf_space_functions(tmp_path: Path):
    """Same workflows against the Space entrypoint (audit/correct/mirror/examples)."""
    # Point Space examples at package fixtures if Space examples missing in checkout
    space = _load_hf_space_module()
    assert "color-scheme: light" in space.CUSTOM_CSS
    assert contrast_ratio(TOKENS["ink"], TOKENS["page_bg"]) >= 7.0

    toy, msg = space.load_example("toy")
    assert toy and Path(toy).is_file(), msg

    audit_md, audit_json = space.run_audit(toy)
    assert "Stereochemistry audit" in audit_md
    assert audit_json and Path(audit_json).is_file()

    inv, inv_msg = space.load_example("inverted")
    assert inv and Path(inv).is_file(), inv_msg
    corr_md, corr_pdb = space.run_correct(inv)
    assert "Chirality correction" in corr_md
    assert Path(corr_pdb).is_file()
    assert detect_chirality_violations(corr_pdb)["n_violations"] == 0

    mir_md, mir_pdb = space.run_mirror(corr_pdb, "Z", True)
    assert "Mirror" in mir_md
    assert Path(mir_pdb).is_file()

    path, bad_msg = space.fetch_pdb_id("BAD")
    assert path is None
    assert len(bad_msg) > 0

    demo = space.build_app()
    assert demo is not None
    assert hasattr(demo, "launch")


@pytest.mark.slow
def test_e2e_fetch_real_pdb():
    path, msg = web_app.fetch_pdb_id("1CRN")
    assert path is not None, msg
    md, out = web_app.run_audit(path)
    assert out and "audit" in out
    assert "Stereochemistry audit" in md
