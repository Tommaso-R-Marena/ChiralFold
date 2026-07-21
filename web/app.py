#!/usr/bin/env python3
"""
ChiralFold Web — upload a PDB, audit, correct chirality, or mirror L↔D.

Launch locally:
    pip install -e ".[web]"
    chiralfold-web
    # → http://localhost:7860

Hugging Face Space:
    https://huggingface.co/spaces/The-Philosopher/ChiralFold-App

Docker:
    docker compose up web
"""
from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import gradio as gr

from chiralfold import __version__
from chiralfold.af3_correct import correct_af3_output, detect_chirality_violations
from chiralfold.auditor import audit_pdb, format_report
from chiralfold.fetch import FetchError, fetch_structure
from chiralfold.pdb_pipeline import mirror_pdb
from chiralfold.structure_io import SUPPORTED_UPLOAD_SUFFIXES
from web.theme import CUSTOM_CSS, header_html, make_theme
from web.ui_format import audit_result_markdown, error_markdown, footer_html

# ---------------------------------------------------------------------------
# Limits & paths
# ---------------------------------------------------------------------------

MAX_UPLOAD_BYTES = int(os.environ.get("CHIRALFOLD_MAX_UPLOAD_MB", "25")) * 1024 * 1024

_PKG_ROOT = Path(__file__).resolve().parent.parent
_EXAMPLE_TOY = _PKG_ROOT / "chiralfold" / "data" / "examples" / "toy_ubiquitin_fragment.pdb"
_FIXTURE_INVERTED = (
    _PKG_ROOT / "chiralfold" / "data" / "examples" / "synthetic_l_ala_inverted.pdb"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _as_path(uploaded) -> str:
    if uploaded is None:
        raise gr.Error(
            "Upload a PDB/mmCIF/FASTA file, fetch an ID, or load an example first."
        )
    path = uploaded if isinstance(uploaded, str) else getattr(uploaded, "name", None)
    if not path:
        raise gr.Error("Could not read the uploaded file.")
    path = str(path)
    size = Path(path).stat().st_size
    if size > MAX_UPLOAD_BYTES:
        mb = MAX_UPLOAD_BYTES / (1024 * 1024)
        raise gr.Error(f"File too large ({size / 1e6:.1f} MB). Limit is {mb:.0f} MB.")
    if size == 0:
        raise gr.Error("Uploaded file is empty.")
    try:
        result = fetch_structure(path, max_bytes=MAX_UPLOAD_BYTES)
    except FetchError as exc:
        raise gr.Error(str(exc)) from exc
    return result.path


def _copy_with_name(src: str, stem: str, suffix: str) -> str:
    """Copy to a stable download filename in a temp dir."""
    out_dir = Path(tempfile.mkdtemp(prefix="chiralfold_dl_"))
    dest = out_dir / f"{stem}{suffix}"
    shutil.copy2(src, dest)
    return str(dest)


def _stem_from_path(path: str) -> str:
    name = Path(path).stem
    for junk in ("_corrected", "_mirror", "_audit"):
        name = name.replace(junk, "")
    return name or "structure"


def _kpi_html(cells: List[Tuple[str, str, str]]) -> str:
    """cells: list of (label, value, css_class)."""
    parts = ['<div class="cf-kpi">']
    for label, value, cls in cells:
        parts.append(
            f'<div class="cell {cls}"><div class="label">{label}</div>'
            f'<div class="value">{value}</div></div>'
        )
    parts.append("</div>")
    return "".join(parts)


def fetch_pdb_id(query: str) -> Tuple[Optional[str], str]:
    """Download from RCSB or AlphaFold DB (PDB ID / UniProt / AFDB)."""
    q = (query or "").strip()
    if not q:
        return (
            None,
            "Enter a PDB ID (1UBQ), UniProt accession (P04637), or AFDB id (AF-P04637-F1).",
        )
    try:
        result = fetch_structure(q, max_bytes=MAX_UPLOAD_BYTES)
    except FetchError as exc:
        return None, str(exc)
    except Exception as exc:  # noqa: BLE001
        return None, f"Could not fetch {q}: {exc}"
    kb = (result.bytes_downloaded or Path(result.path).stat().st_size) / 1024
    return (
        result.path,
        f"Loaded **{result.label}** from `{result.source}` ({kb:.0f} KB).",
    )


def load_example(which: str) -> Tuple[Optional[str], str]:
    mapping = {
        "toy": (_EXAMPLE_TOY, "Toy ubiquitin fragment (10 residues) — bundled, offline."),
        "inverted": (
            _FIXTURE_INVERTED,
            "Synthetic inverted L-Ala — demonstrates chirality correction.",
        ),
    }
    path, msg = mapping.get(which, (None, "Unknown example."))
    if path is None or not Path(path).is_file():
        return None, f"Example not found: {which}"
    dest = Path(tempfile.mkdtemp(prefix="chiralfold_ex_")) / path.name
    shutil.copy2(path, dest)
    return str(dest), msg


# ---------------------------------------------------------------------------
# Core actions
# ---------------------------------------------------------------------------


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
        named = _copy_with_name(raw_json, _stem_from_path(src), "_audit.json")
        return md, named
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
        n_before = before.get("n_violations", result.get("before", {}).get("n_violations", 0))
        n_after = after.get("n_violations", 0)
        n_fixed = corr.get("n_corrected", max(0, n_before - n_after))
        before_cls = "ok" if n_before == 0 else "bad"
        after_cls = "ok" if n_after == 0 else "warn"
        kpis = _kpi_html(
            [
                ("Violations before", str(n_before), before_cls),
                ("Violations after", str(n_after), after_cls),
                ("Fixed", str(n_fixed), "ok"),
                ("Residual rate", f"{after.get('pct_violations', 0):.1f}%", after_cls),
            ]
        )
        md = (
            f"## Chirality correction\n\n{kpis}\n\n"
            "Detects inverted Cα stereocenters (common in AlphaFold 3 D-peptide "
            "outputs) and rewrites coordinates to the expected chirality class.\n\n"
            "**Download the corrected PDB below.**"
        )
        named = _copy_with_name(out_raw.name, _stem_from_path(src), "_corrected.pdb")
        return md, named
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
        n_atoms = info.get("n_atoms", "?")
        n_res = info.get("n_residues", "?")
        kpis = _kpi_html(
            [
                ("Atoms mirrored", str(n_atoms), "ok"),
                ("Residues", str(n_res), "ok"),
                ("Axis", axis.upper(), "ok"),
                ("Rename D↔L", "yes" if rename else "no", "ok"),
            ]
        )
        md = (
            f"## Mirror-image (L↔D) transform\n\n{kpis}\n\n"
            "Exact coordinate reflection — mathematically inverts chirality at "
            "every stereocenter. Ideal for mirror-image binder / D-peptide design.\n\n"
            "**Download the mirrored PDB below.**"
        )
        named = _copy_with_name(out_raw.name, _stem_from_path(src), "_mirror.pdb")
        return md, named
    except gr.Error:
        raise
    except Exception as exc:  # noqa: BLE001
        return error_markdown(exc), None


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------


def build_app() -> gr.Blocks:
    with gr.Blocks(
        title=f"ChiralFold Web v{__version__}",
        css=CUSTOM_CSS,
        theme=make_theme(),
    ) as demo:
        gr.HTML(header_html(__version__))

        with gr.Row(equal_height=False):
            with gr.Column(scale=5, elem_classes=["cf-panel"]):
                gr.Markdown("### 1 · Load a structure")
                pdb_in = gr.File(
                    label="Upload PDB / mmCIF / FASTA",
                    file_types=list(SUPPORTED_UPLOAD_SUFFIXES),
                    type="filepath",
                    height=120,
                )
                with gr.Row():
                    pdb_id = gr.Textbox(
                        label="Or fetch (RCSB / AlphaFold DB)",
                        placeholder="1UBQ · P04637 · AF-P04637-F1",
                        max_lines=1,
                        scale=3,
                    )
                    fetch_btn = gr.Button("Fetch", scale=1)
                with gr.Row():
                    ex_toy = gr.Button("Example: ubiquitin fragment", size="sm")
                    ex_inv = gr.Button("Example: inverted Ala (fix demo)", size="sm")
                status = gr.Markdown(
                    "Upload PDB/mmCIF/FASTA, fetch an ID, or load an example."
                )

                gr.Markdown("### 2 · Mirror options")
                axis = gr.Radio(
                    choices=["X", "Y", "Z"],
                    value="X",
                    label="Reflection axis",
                    info="Used only for Mirror Image",
                )
                rename = gr.Checkbox(
                    value=True,
                    label="Rename residue codes L↔D after mirroring",
                )

            with gr.Column(scale=7, elem_classes=["cf-panel"]):
                gr.Markdown("### 3 · Choose an action")
                with gr.Tabs():
                    with gr.Tab("Audit"):
                        gr.Markdown(
                            "Full stereochemistry report: **Cα chirality**, Ramachandran, "
                            "planarity, clashes, and a composite quality score."
                        )
                        audit_btn = gr.Button(
                            "Run audit → download JSON", variant="primary", size="lg"
                        )
                    with gr.Tab("Correct chirality"):
                        gr.Markdown(
                            "Detect and **fix inverted stereocenters** — built for "
                            "AlphaFold 3 / diffusion outputs with D-peptide chirality errors "
                            "(Childs et al. 2025)."
                        )
                        correct_btn = gr.Button(
                            "Correct → download PDB", variant="primary", size="lg"
                        )
                    with gr.Tab("Mirror image"):
                        gr.Markdown(
                            "Exact **L↔D coordinate reflection** for mirror-image binder "
                            "and D-peptide design. RMSD to the reflected structure is 0.0 Å."
                        )
                        mirror_btn = gr.Button(
                            "Mirror → download PDB", variant="primary", size="lg"
                        )

                report_out = gr.Markdown(value="Results appear here after you run an action.")
                download_out = gr.File(label="Download result")

        fetch_btn.click(fetch_pdb_id, inputs=pdb_id, outputs=[pdb_in, status])
        ex_toy.click(lambda: load_example("toy"), outputs=[pdb_in, status])
        ex_inv.click(lambda: load_example("inverted"), outputs=[pdb_in, status])

        audit_btn.click(run_audit, inputs=pdb_in, outputs=[report_out, download_out])
        correct_btn.click(run_correct, inputs=pdb_in, outputs=[report_out, download_out])
        mirror_btn.click(
            run_mirror, inputs=[pdb_in, axis, rename], outputs=[report_out, download_out]
        )

        gr.HTML(footer_html(__version__, include_hf_link=True))

    return demo


def main() -> None:
    host = os.environ.get("CHIRALFOLD_WEB_HOST", "0.0.0.0")
    port = int(os.environ.get("CHIRALFOLD_WEB_PORT", "7860"))
    share = os.environ.get("CHIRALFOLD_WEB_SHARE", "").lower() in ("1", "true", "yes")
    demo = build_app()
    demo.queue(default_concurrency_limit=4)
    demo.launch(
        server_name=host,
        server_port=port,
        share=share,
        show_error=True,
        favicon_path=None,
    )


if __name__ == "__main__":
    main()
