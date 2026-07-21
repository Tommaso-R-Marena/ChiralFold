#!/usr/bin/env python3
"""
ChiralFold Hugging Face Space entrypoint.

Self-contained Gradio UI that depends on the `chiralfold` package
(installed via requirements.txt from GitHub). Theme/CSS is loaded from
``web.theme`` when available, otherwise from bundled ``theme.css``.
"""
from __future__ import annotations

import json
import os
import shutil
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import gradio as gr

from chiralfold import __version__
from chiralfold.af3_correct import correct_af3_output, detect_chirality_violations
from chiralfold.auditor import audit_pdb, format_report
from chiralfold.pdb_pipeline import mirror_pdb

MAX_UPLOAD_BYTES = int(os.environ.get("CHIRALFOLD_MAX_UPLOAD_MB", "25")) * 1024 * 1024
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
SPACE_ROOT = Path(__file__).resolve().parent
EXAMPLE_TOY = SPACE_ROOT / "examples" / "toy_ubiquitin_fragment.pdb"
EXAMPLE_INV = SPACE_ROOT / "examples" / "synthetic_l_ala_inverted.pdb"


def _load_css() -> str:
    """Prefer packaged theme (post-install); fall back to bundled theme.css."""
    try:
        from web.theme import CUSTOM_CSS as css  # type: ignore

        return css
    except Exception:
        css_path = SPACE_ROOT / "theme.css"
        if css_path.is_file():
            return css_path.read_text()
        return "/* missing theme */"


def _make_theme():
    try:
        from web.theme import make_theme

        return make_theme()
    except Exception:
        return gr.themes.Soft(
            primary_hue="teal",
            secondary_hue="slate",
            neutral_hue="slate",
            font=["IBM Plex Sans", "Segoe UI", "Helvetica", "Arial", "sans-serif"],
        )


try:
    from web.ui_format import audit_result_markdown, error_markdown, footer_html
except ImportError:
    import html as _html
    import traceback as _traceback

    def report_pre_html(text: str) -> str:
        return f'<pre class="cf-report">{_html.escape(text)}</pre>'

    def audit_result_markdown(kpis_html: str, filename: str, summary: str) -> str:
        return (
            f"## Stereochemistry audit\n\n{kpis_html}\n\n"
            f"**File:** `{filename}`\n\n"
            f"{report_pre_html(summary)}\n\n"
            "Download the full JSON report below."
        )

    def error_markdown(exc: BaseException) -> str:
        tb = _html.escape(_traceback.format_exc())
        return f"**Error:** {_html.escape(str(exc))}\n\n<pre class=\"cf-error\">{tb}</pre>"

    def footer_html(version: str, *, include_hf_link: bool = False) -> str:
        return (
            f'<div id="cf-footer"><p class="cf-footer-line">'
            f"<strong>ChiralFold v{_html.escape(version)}</strong> · "
            f'<a href="https://github.com/Tommaso-R-Marena/ChiralFold">GitHub</a> · '
            f'<code class="cf-install-cmd">pip install chiralfold</code></p>'
            f'<p class="cf-footer-note">Processed in-session; nothing stored permanently.</p></div>'
        )


CUSTOM_CSS = _load_css()


def _as_path(uploaded) -> str:
    if uploaded is None:
        raise gr.Error("Upload a PDB file, fetch a PDB ID, or load an example first.")
    path = uploaded if isinstance(uploaded, str) else getattr(uploaded, "name", None)
    if not path:
        raise gr.Error("Could not read the uploaded file.")
    path = str(path)
    if not path.lower().endswith((".pdb", ".ent")):
        raise gr.Error("Only PDB files are supported (.pdb or .ent).")
    size = Path(path).stat().st_size
    if size > MAX_UPLOAD_BYTES:
        raise gr.Error(
            f"File too large ({size / 1e6:.1f} MB). "
            f"Limit is {MAX_UPLOAD_BYTES / 1e6:.0f} MB."
        )
    if size == 0:
        raise gr.Error("Uploaded file is empty.")
    return path


def _copy_with_name(src: str, stem: str, suffix: str) -> str:
    out_dir = Path(tempfile.mkdtemp(prefix="chiralfold_dl_"))
    dest = out_dir / f"{stem}{suffix}"
    shutil.copy2(src, dest)
    return str(dest)


def _stem_from_path(path: str) -> str:
    name = Path(path).stem
    for junk in ("_corrected", "_mirror", "_audit"):
        name = name.replace(junk, "")
    return name or "structure"


def _kpi_html(cells):
    parts = ['<div class="cf-kpi">']
    for label, value, cls in cells:
        parts.append(
            f'<div class="cell {cls}"><div class="label">{label}</div>'
            f'<div class="value">{value}</div></div>'
        )
    parts.append("</div>")
    return "".join(parts)


def fetch_pdb_id(pdb_id: str) -> Tuple[Optional[str], str]:
    pid = (pdb_id or "").strip().upper()
    if len(pid) != 4 or not pid.isalnum():
        return None, "Enter a valid 4-character PDB ID (e.g. 1UBQ)."
    url = RCSB_DOWNLOAD.format(pdb_id=pid)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
    except urllib.error.HTTPError as exc:
        return None, f"RCSB returned HTTP {exc.code} for {pid}."
    except Exception as exc:  # noqa: BLE001
        return None, f"Could not fetch {pid}: {exc}"
    if not data or (b"ATOM" not in data and b"HETATM" not in data):
        return None, f"Downloaded file for {pid} does not look like a PDB."
    out = Path(tempfile.mkdtemp(prefix="chiralfold_rcsb_")) / f"{pid}.pdb"
    out.write_bytes(data)
    return str(out), f"Loaded **{pid}** from RCSB ({len(data) / 1024:.0f} KB)."


def load_example(which: str) -> Tuple[Optional[str], str]:
    mapping = {
        "toy": (EXAMPLE_TOY, "Toy ubiquitin fragment (10 residues) — bundled, offline."),
        "inverted": (
            EXAMPLE_INV,
            "Synthetic inverted L-Ala — demonstrates chirality correction.",
        ),
    }
    path, msg = mapping.get(which, (None, "Unknown example."))
    if path is None or not Path(path).is_file():
        return None, f"Example not found: {which}"
    dest = Path(tempfile.mkdtemp(prefix="chiralfold_ex_")) / path.name
    shutil.copy2(path, dest)
    return str(dest), msg


def run_audit(uploaded) -> Tuple[str, Optional[str]]:
    try:
        src = _as_path(uploaded)
        report = audit_pdb(src)
        chir = report["chirality"]
        rama = report["ramachandran"]
        n_wrong = chir.get("n_wrong", 0)
        score = report["overall_score"]
        score_cls = "ok" if score >= 80 else ("warn" if score >= 60 else "bad")
        chir_cls = "ok" if n_wrong == 0 else "bad"
        kpis = _kpi_html(
            [
                ("Overall score", f"{score:.1f}", score_cls),
                ("Cα correct %", f"{chir['pct_correct']:.1f}", chir_cls),
                ("Chirality violations", str(n_wrong), chir_cls),
                ("Ramachandran favored %", f"{rama['pct_favored']:.1f}", "ok"),
            ]
        )
        summary = format_report(report)
        md = audit_result_markdown(kpis, Path(src).name, summary)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False, prefix="chiralfold_audit_"
        ) as fp:
            json.dump(report, fp, indent=2, default=float)
            raw_json = fp.name
        return md, _copy_with_name(raw_json, _stem_from_path(src), "_audit.json")
    except gr.Error:
        raise
    except Exception as exc:  # noqa: BLE001
        return error_markdown(exc), None


def run_correct(uploaded) -> Tuple[str, Optional[str]]:
    try:
        src = _as_path(uploaded)
        before = detect_chirality_violations(src)
        out_raw = tempfile.NamedTemporaryFile(
            suffix="_corrected.pdb", delete=False, prefix="chiralfold_"
        )
        out_raw.close()
        result = correct_af3_output(src, out_raw.name)
        corr = result.get("correction", {})
        after = result.get("after", {})
        n_before = before.get(
            "n_violations", result.get("before", {}).get("n_violations", 0)
        )
        n_after = after.get("n_violations", 0)
        n_fixed = corr.get("n_corrected", max(0, n_before - n_after))
        before_cls = "ok" if n_before == 0 else "bad"
        after_cls = "ok" if n_after == 0 else "bad"
        kpis = _kpi_html(
            [
                ("Before", str(n_before), before_cls),
                ("Fixed", str(n_fixed), "ok"),
                ("After (residual)", str(n_after), after_cls),
            ]
        )
        md = (
            f"## Chirality correction\n\n{kpis}\n\n"
            f"**File:** `{Path(src).name}`\n\n"
            "**Download the corrected PDB below.**"
        )
        return md, _copy_with_name(out_raw.name, _stem_from_path(src), "_corrected.pdb")
    except gr.Error:
        raise
    except Exception as exc:  # noqa: BLE001
        return error_markdown(exc), None


def run_mirror(uploaded, axis: str, rename: bool) -> Tuple[str, Optional[str]]:
    try:
        src = _as_path(uploaded)
        out_raw = tempfile.NamedTemporaryFile(
            suffix="_mirror.pdb", delete=False, prefix="chiralfold_"
        )
        out_raw.close()
        info: Dict[str, Any] = mirror_pdb(
            src, out_raw.name, axis=axis.lower(), rename_residues=rename
        )
        kpis = _kpi_html(
            [
                ("Atoms mirrored", str(info.get("n_atoms", "?")), "ok"),
                ("Residues", str(info.get("n_residues", "?")), "ok"),
                ("Axis", axis.upper(), "ok"),
                ("Rename D↔L", "yes" if rename else "no", "ok"),
            ]
        )
        md = (
            f"## Mirror-image (L↔D) transform\n\n{kpis}\n\n"
            "Exact coordinate reflection for mirror-image binder design.\n\n"
            "**Download the mirrored PDB below.**"
        )
        return md, _copy_with_name(out_raw.name, _stem_from_path(src), "_mirror.pdb")
    except gr.Error:
        raise
    except Exception as exc:  # noqa: BLE001
        return error_markdown(exc), None


def build_app() -> gr.Blocks:
    with gr.Blocks(
        title=f"ChiralFold Web v{__version__}",
        css=CUSTOM_CSS,
        theme=_make_theme(),
    ) as demo:
        gr.HTML(
            f"""
            <div id="cf-header">
              <h1>ChiralFold</h1>
              <p class="cf-tagline">
                Upload any PDB — audit stereochemistry, fix chirality errors,
                or generate an exact mirror-image structure. No command line.
              </p>
              <div class="cf-meta">
                <span class="cf-badge">v{__version__}</span>
                <span class="cf-badge">0% residual chirality violations</span>
                <span class="cf-badge">open source</span>
              </div>
            </div>
            """
        )
        with gr.Row():
            with gr.Column(scale=5, elem_classes=["cf-panel"]):
                gr.Markdown("### 1 · Load a structure")
                pdb_in = gr.File(
                    label="Upload PDB (.pdb / .ent)",
                    file_types=[".pdb", ".ent"],
                    type="filepath",
                )
                with gr.Row():
                    pdb_id = gr.Textbox(
                        label="Or fetch from RCSB",
                        placeholder="1UBQ",
                        max_lines=1,
                        scale=3,
                    )
                    fetch_btn = gr.Button("Fetch PDB", scale=1)
                with gr.Row():
                    ex_toy = gr.Button("Example: ubiquitin fragment", size="sm")
                    ex_inv = gr.Button("Example: inverted Ala", size="sm")
                status = gr.Markdown(
                    "Upload a file, fetch a PDB ID, or load an example."
                )
                gr.Markdown("### 2 · Mirror options")
                axis = gr.Radio(
                    choices=["X", "Y", "Z"], value="X", label="Reflection axis"
                )
                rename = gr.Checkbox(
                    value=True, label="Rename residue codes L↔D after mirroring"
                )
            with gr.Column(scale=7, elem_classes=["cf-panel"]):
                gr.Markdown("### 3 · Choose an action")
                with gr.Tabs():
                    with gr.Tab("Audit"):
                        gr.Markdown(
                            "Full stereochemistry report with downloadable JSON."
                        )
                        audit_btn = gr.Button(
                            "Run audit → download JSON", variant="primary", size="lg"
                        )
                    with gr.Tab("Correct chirality"):
                        gr.Markdown(
                            "Fix inverted stereocenters (AF3 / diffusion outputs)."
                        )
                        correct_btn = gr.Button(
                            "Correct → download PDB", variant="primary", size="lg"
                        )
                    with gr.Tab("Mirror image"):
                        gr.Markdown(
                            "Exact L↔D coordinate reflection for binder design."
                        )
                        mirror_btn = gr.Button(
                            "Mirror → download PDB", variant="primary", size="lg"
                        )
                report_out = gr.Markdown(
                    value="Results appear here after you run an action."
                )
                download_out = gr.File(label="Download result")

        fetch_btn.click(fetch_pdb_id, inputs=pdb_id, outputs=[pdb_in, status])
        ex_toy.click(lambda: load_example("toy"), outputs=[pdb_in, status])
        ex_inv.click(lambda: load_example("inverted"), outputs=[pdb_in, status])
        audit_btn.click(run_audit, inputs=pdb_in, outputs=[report_out, download_out])
        correct_btn.click(
            run_correct, inputs=pdb_in, outputs=[report_out, download_out]
        )
        mirror_btn.click(
            run_mirror,
            inputs=[pdb_in, axis, rename],
            outputs=[report_out, download_out],
        )

        gr.HTML(footer_html(__version__))
    return demo


demo = build_app()

if __name__ == "__main__":
    demo.queue(default_concurrency_limit=4)
    demo.launch(server_name="0.0.0.0", server_port=7860, show_error=True)
