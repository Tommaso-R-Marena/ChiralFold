"""Shared HTML/Markdown fragments for ChiralFold web + HF Space."""
from __future__ import annotations

import html
import traceback


def report_pre_html(text: str) -> str:
    """Monospace audit report without Gradio syntax-highlight spans."""
    return f'<pre class="cf-report">{html.escape(text)}</pre>'


def audit_result_markdown(kpis_html: str, filename: str, summary: str) -> str:
    return (
        f"## Stereochemistry audit\n\n{kpis_html}\n\n"
        f"**File:** `{filename}`\n\n"
        f"{report_pre_html(summary)}\n\n"
        "Download the full JSON report below."
    )


def error_markdown(exc: BaseException) -> str:
    tb = html.escape(traceback.format_exc())
    return (
        f"**Error:** {html.escape(str(exc))}\n\n"
        f'<pre class="cf-error">{tb}</pre>'
    )


def footer_html(version: str, *, include_hf_link: bool = False) -> str:
    hf = (
        ' · <a href="https://huggingface.co/spaces/The-Philosopher/ChiralFold-App">'
        "Hugging Face Space</a>"
        if include_hf_link
        else ""
    )
    return f"""
<div id="cf-footer">
  <p class="cf-footer-line">
    <strong>ChiralFold v{html.escape(version)}</strong>
    · <a href="https://github.com/Tommaso-R-Marena/ChiralFold">GitHub</a>{hf}
    · <code class="cf-install-cmd">pip install chiralfold</code>
  </p>
  <p class="cf-footer-note">Processed in-session; nothing stored permanently.</p>
</div>
"""
