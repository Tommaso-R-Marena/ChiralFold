"""Tests for shared web UI formatting helpers."""
from web.ui_format import audit_result_markdown, footer_html, report_pre_html


def test_report_pre_html_escapes_and_uses_class():
    html = report_pre_html("Score: 80\nWrong: 0")
    assert 'class="cf-report"' in html
    assert "Score: 80" in html
    assert "<script>" not in report_pre_html("<script>")


def test_audit_result_markdown_uses_pre_not_fence():
    md = audit_result_markdown("<div></div>", "1UBQ.pdb", "Residues: 76")
    assert "cf-report" in md
    assert "```" not in md
    assert "1UBQ.pdb" in md


def test_footer_html_shows_full_pip_command():
    html = footer_html("3.5.1")
    assert "pip install chiralfold" in html
    assert "cf-install-cmd" in html
    assert "cf-footer-note" in html
