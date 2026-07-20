#!/usr/bin/env python3
"""
ChiralFold Web — upload a PDB, audit, correct chirality, or mirror L↔D.

Launch:
    pip install -e ".[web]"
    python web/app.py

Docker:
    docker compose up web
"""
from __future__ import annotations

import json
import os
import tempfile
import traceback
from pathlib import Path
from typing import Optional, Tuple

import gradio as gr

from chiralfold import __version__
from chiralfold.af3_correct import correct_af3_output
from chiralfold.auditor import audit_pdb, format_report
from chiralfold.pdb_pipeline import mirror_pdb

CUSTOM_CSS = """
:root {
  --cf-teal: #0d9488;
  --cf-navy: #0f172a;
  --cf-slate: #334155;
  --cf-bg: #f8fafc;
}
.gradio-container {
  font-family: 'Inter', 'Segoe UI', system-ui, sans-serif !important;
  background: linear-gradient(160deg, #f0fdfa 0%, #f8fafc 45%, #eff6ff 100%) !important;
}
#cf-header {
  text-align: center;
  padding: 1.5rem 1rem 0.5rem;
}
#cf-header h1 {
  color: var(--cf-navy);
  font-weight: 700;
  letter-spacing: -0.02em;
  margin-bottom: 0.25rem;
}
#cf-header p {
  color: var(--cf-slate);
  max-width: 42rem;
  margin: 0 auto;
}
.cf-card {
  border: 1px solid #e2e8f0 !important;
  border-radius: 16px !important;
  box-shadow: 0 4px 24px rgba(15, 23, 42, 0.06) !important;
  background: white !important;
}
footer { visibility: hidden; }
"""

EXAMPLES = {
    "Ubiquitin fragment (1LDF)": "https://files.rcsb.org/download/1LDF.pdb",
    "Hirulog–thrombin (1ABI)": "https://files.rcsb.org/download/1ABI.pdb",
}


def _save_upload(uploaded) -> str:
    if uploaded is None:
        raise gr.Error("Please upload a .pdb file first.")
    path = uploaded if isinstance(uploaded, str) else uploaded.name
    if not str(path).lower().endswith((".pdb", ".ent")):
        raise gr.Error("Only PDB files are supported (.pdb or .ent).")
    return path


def run_audit(uploaded) -> Tuple[str, Optional[str]]:
    try:
        src = _save_upload(uploaded)
        report = audit_pdb(src)
        summary = format_report(report)
        chir = report["chirality"]
        md = (
            f"## Audit complete\n\n"
            f"**Overall score:** {report['overall_score']:.1f}/100  \n"
            f"**Cα chirality:** {chir['pct_correct']:.1f}% correct "
            f"({chir['n_wrong']} violations)  \n"
            f"**Ramachandran favored:** {report['ramachandran']['pct_favored']:.1f}%  \n\n"
            f"```\n{summary}\n```"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix="_audit.json", delete=False, prefix="chiralfold_"
        ) as fp:
            json.dump(report, fp, indent=2, default=float)
            json_path = fp.name
        return md, json_path
    except Exception as exc:
        return f"**Error:** {exc}\n\n```\n{traceback.format_exc()}\n```", None


def run_correct(uploaded) -> Tuple[str, Optional[str]]:
    try:
        src = _save_upload(uploaded)
        out = tempfile.NamedTemporaryFile(suffix="_corrected.pdb", delete=False, prefix="chiralfold_")
        out.close()
        result = correct_af3_output(src, out.name)
        corr = result.get("correction", {})
        before = result.get("before", {})
        after = result.get("after", {})
        md = (
            "## Chirality correction\n\n"
            f"**Violations fixed:** {corr.get('n_corrected', 0)}  \n"
            f"**Before:** {before.get('n_violations', '?')} violations  \n"
            f"**After:** {after.get('n_violations', '?')} violations  \n\n"
            "Download the corrected PDB below."
        )
        return md, out.name
    except Exception as exc:
        return f"**Error:** {exc}", None


def run_mirror(uploaded, axis: str, rename: bool) -> Tuple[str, Optional[str]]:
    try:
        src = _save_upload(uploaded)
        out = tempfile.NamedTemporaryFile(suffix="_mirror.pdb", delete=False, prefix="chiralfold_")
        out.close()
        mirror_pdb(src, out.name, axis=axis.lower(), rename_residues=rename)
        md = (
            "## Mirror-image transform\n\n"
            f"**Axis:** reflect across {axis.upper()}  \n"
            f"**Rename D↔L codes:** {'yes' if rename else 'no'}  \n\n"
            "Mathematically exact L↔D coordinate reflection. Download below."
        )
        return md, out.name
    except Exception as exc:
        return f"**Error:** {exc}", None


def build_app() -> gr.Blocks:
    with gr.Blocks(
        title="ChiralFold Web",
        css=CUSTOM_CSS,
        theme=gr.themes.Soft(
            primary_hue="teal",
            neutral_hue="slate",
            font=[gr.themes.GoogleFont("Inter"), "system-ui", "sans-serif"],
        ),
    ) as demo:
        gr.HTML(
            f"""
            <div id="cf-header">
              <h1>ChiralFold</h1>
              <p>Upload any PDB structure — audit stereochemistry, fix AF3 chirality errors,
              or generate an exact mirror-image (L↔D) in one click. No command line required.
              <br><small>v{__version__} · open source · CUA</small></p>
            </div>
            """
        )

        with gr.Row(equal_height=True):
            with gr.Column(scale=1, elem_classes=["cf-card"]):
                pdb_in = gr.File(
                    label="Upload PDB",
                    file_types=[".pdb", ".ent"],
                    type="filepath",
                )
                gr.Markdown(
                    "**Workflow:** upload → choose an action → download results."
                )
                axis = gr.Radio(
                    choices=["X", "Y", "Z"],
                    value="X",
                    label="Mirror reflection axis",
                    info="Used only for Mirror tab",
                )
                rename = gr.Checkbox(
                    value=True,
                    label="Rename D-amino acid residue codes (mirror)",
                )

            with gr.Column(scale=2, elem_classes=["cf-card"]):
                with gr.Tabs():
                    with gr.Tab("Audit"):
                        gr.Markdown(
                            "Full stereochemistry report: Cα chirality, Ramachandran, "
                            "planarity, clashes, composite score."
                        )
                        audit_btn = gr.Button("Run audit", variant="primary")
                    with gr.Tab("Correct chirality"):
                        gr.Markdown(
                            "Detect and fix Cα chirality violations (ideal for AlphaFold 3 "
                            "or diffused structures)."
                        )
                        correct_btn = gr.Button("Correct & download", variant="primary")
                    with gr.Tab("Mirror image"):
                        gr.Markdown(
                            "Exact L↔D mirror transform for binder design and D-peptide workflows."
                        )
                        mirror_btn = gr.Button("Mirror & download", variant="primary")

                report_out = gr.Markdown(label="Results")
                download_out = gr.File(label="Download output")

        audit_btn.click(run_audit, pdb_in, [report_out, download_out])
        correct_btn.click(run_correct, pdb_in, [report_out, download_out])
        mirror_btn.click(run_mirror, [pdb_in, axis, rename], [report_out, download_out])

        gr.Markdown(
            "---\n"
            "**Citation:** Marena T.R. ChiralFold (2026). "
            "[GitHub](https://github.com/Tommaso-R-Marena/ChiralFold) · "
            "`pip install chiralfold` · "
            "Childs et al. (2025) bioRxiv for AF3 D-peptide context."
        )

    return demo


def main() -> None:
    host = os.environ.get("CHIRALFOLD_WEB_HOST", "0.0.0.0")
    port = int(os.environ.get("CHIRALFOLD_WEB_PORT", "7860"))
    share = os.environ.get("CHIRALFOLD_WEB_SHARE", "").lower() in ("1", "true", "yes")
    demo = build_app()
    demo.queue(default_concurrency_limit=4)
    demo.launch(server_name=host, server_port=port, share=share, show_error=True)


if __name__ == "__main__":
    main()
