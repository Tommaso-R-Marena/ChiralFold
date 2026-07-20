#!/usr/bin/env python3
"""Build results/colab_integrated_manifest.json from repo + Colab run artefacts."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def _load_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def build_manifest() -> dict:
    d_survey = _load_json(RESULTS / "d_residue_verification_summary.json")
    exp_val = _load_json(RESULTS / "experimental_validation_report.json")
    af3 = _load_json(RESULTS / "af3_resource_benchmark.json")
    chainfix = _load_json(RESULTS / "ramachandran_279struct_chainfix_summary.json")
    colab_100 = _load_json(RESULTS / "colab_runs/ramachandran_100struct/ramachandran_100struct_summary.json")
    colab_200 = _load_json(RESULTS / "colab_runs/ramachandran_200struct/ramachandran_200struct_summary.json")
    error_cls = _load_json(RESULTS / "error_classification.json")

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "chiralfold_version": "3.5.0",
        "dashboard_notebook": "demos/ChiralFold_Results_Dashboard.ipynb",
        "interpretation": {
            "paper_headline_ramachandran": (
                "Manuscript headline: n=362 clean (363 minus 5M2K glycopeptide), "
                "Spearman ρ=0.521, Pearson r=0.853 on Ramachandran outlier percentage."
            ),
            "colab_replications": (
                "Colab Expanded Ramachandran runs (n=155, n=285) are independent "
                "replications at different sample sizes/seeds — not interchangeable "
                "with the chain-ID-fixed n=362 benchmark."
            ),
            "d_residue_survey": (
                "PDB-wide survey: 12,573 checkable D-residues in 4,616 files; "
                "29 D-label/L-coordinate mismatches in 16 structures (0.23%)."
            ),
            "experimental_validation": (
                "14/14 non-borderline error structures pass automated RCSB+CCD criteria; "
                "2 borderline (1KO0, 2H9E) require manual review."
            ),
            "af3_correction": (
                "Synthetic AF3-mimetic benchmark: 100% detection and correction of "
                "inverted stereocenters; 0% residual violations; ~37 ms full pipeline."
            ),
            "aristotle_lean": (
                "Machine-checked Lean 4 chirality no-go theorems in "
                "formal/chirality_nogo/ (Aristotle / Harmonic): distance-only "
                "representations cannot recover signed orientation; lake build clean."
            ),
        },
        "headline_metrics": {
            "d_residue_survey": {
                "residues_checked": d_survey["checkable_residues"],
                "pdb_files": d_survey.get("pdb_files_scanned", 4616),
                "errors": d_survey["l_error"],
                "structures_with_errors": len(d_survey.get("errors_by_structure", {})),
                "error_rate_pct": d_survey["error_rate_pct"],
            },
            "ramachandran_paper": {
                "n_clean": chainfix["n_success_clean"],
                "spearman_rho": chainfix["spearman_rho_clean"],
                "pearson_r": chainfix["pearson_r_clean"],
                "excluded": chainfix["excluded_entries"],
            },
            "experimental_validation": {
                "n_structures": exp_val["n_structures"],
                "n_automated_pass": exp_val["n_automated_pass"],
                "n_borderline": exp_val["n_borderline"],
            },
            "af3_correction": {
                "detection_recall_pct": af3["aggregate_detection_recall_pct"],
                "residual_violation_rate_pct": af3["post_correction_residual_violation_rate_pct"],
                "full_pipeline_ms": round(
                    af3["operations"]["full_correct_af3_output_pipeline"]["total_wall_time_s"] * 1000,
                    1,
                ),
            },
        },
        "runs": [
            {
                "id": "d_residue_survey",
                "label": "PDB-wide D-residue survey",
                "source": "repo",
                "notebook": None,
                "artifacts": [
                    "results/d_residue_verification.csv",
                    "results/d_residue_verification_summary.json",
                    "results/error_classification.json",
                ],
                "metrics": {
                    "residues": d_survey["checkable_residues"],
                    "files": 4616,
                    "errors": d_survey["l_error"],
                    "structures": 16,
                },
            },
            {
                "id": "experimental_validation",
                "label": "D-residue experimental validation",
                "source": "colab+repo",
                "notebook": "demos/D_Residue_Experimental_Validation.ipynb",
                "artifacts": [
                    "results/experimental_validation_report.json",
                    "results/experimental_validation_summary.csv",
                ],
                "metrics": {
                    "n_automated_pass": exp_val["n_automated_pass"],
                    "n_borderline": exp_val["n_borderline"],
                    "n_structures": exp_val["n_structures"],
                },
            },
            {
                "id": "ramachandran_chainfix",
                "label": "Ramachandran era-representative (paper headline)",
                "source": "repo",
                "notebook": None,
                "role": "manuscript_headline",
                "artifacts": [
                    "results/ramachandran_279struct_chainfix_summary.json",
                    "results/ramachandran_279struct_chainfix_comparison.csv",
                    "results/ramachandran_279struct_chainfix_plot.png",
                ],
                "metrics": {
                    "n_success": chainfix["n_success_clean"],
                    "spearman_rho": chainfix["spearman_rho_clean"],
                    "pearson_r": chainfix["pearson_r_clean"],
                },
            },
            {
                "id": "ramachandran_colab_100",
                "label": "Ramachandran Colab replication (~100 target)",
                "source": "colab",
                "notebook": "demos/Expanded_Ramachandran_Benchmark.ipynb",
                "role": "independent_replication",
                "artifacts": [
                    "results/colab_runs/ramachandran_100struct/ramachandran_100struct_summary.json",
                    "results/colab_runs/ramachandran_100struct/ramachandran_100struct_comparison.csv",
                    "results/colab_runs/ramachandran_100struct/ramachandran_100struct_plot.png",
                ],
                "metrics": {
                    "n_success": colab_100["n_success"],
                    "spearman_rho": colab_100["spearman_rho"],
                    "pearson_r": colab_100["pearson_r"],
                },
            },
            {
                "id": "ramachandran_colab_200",
                "label": "Ramachandran Colab replication (~200 target)",
                "source": "colab",
                "notebook": "demos/Expanded_Ramachandran_Benchmark.ipynb",
                "role": "independent_replication",
                "artifacts": [
                    "results/colab_runs/ramachandran_200struct/ramachandran_200struct_summary.json",
                ],
                "metrics": {
                    "n_success": colab_200["n_success"],
                    "spearman_rho": colab_200["spearman_rho"],
                    "pearson_r": colab_200["pearson_r"],
                },
            },
            {
                "id": "af3_resource_benchmark",
                "label": "AF3 synthetic correction benchmark",
                "source": "repo",
                "notebook": "demos/ChiralFold_Toy_Dataset_Demo.ipynb",
                "artifacts": [
                    "results/af3_resource_benchmark.json",
                    "results/af3_accuracy_comparison.png",
                ],
                "metrics": {
                    "n_systems": af3["n_systems"],
                    "detection_recall_pct": af3["aggregate_detection_recall_pct"],
                    "residual_violation_rate_pct": af3["post_correction_residual_violation_rate_pct"],
                },
            },
            {
                "id": "quick_demo",
                "label": "ChiralFold Quick Demo (Colab)",
                "source": "colab",
                "notebook": "demos/ChiralFold_Quick_Demo.ipynb",
                "metrics": {
                    "audit_pdb": "1LDF",
                    "audit_score": 83.6,
                    "chirality_violations": 0,
                    "conformers_converged": "45/45",
                },
            },
            {
                "id": "toy_demo",
                "label": "ChiralFold Toy Dataset Demo (Colab)",
                "source": "colab",
                "notebook": "demos/ChiralFold_Toy_Dataset_Demo.ipynb",
                "metrics": {
                    "af3_correction": "1 violation → 0 after correct_af3_output",
                },
            },
        ],
        "error_taxonomy": error_cls.get("classification", {}),
        "ramachandran_progression": [
            {"label": "Original (31 struct)", "n": 31, "spearman_rho": 0.49, "source": "repo_historical"},
            {"label": "Expanded (104 struct)", "n": 104, "spearman_rho": 0.564, "source": "repo"},
            {"label": "Colab (~100 target)", "n": colab_100["n_success"], "spearman_rho": colab_100["spearman_rho"], "source": "colab"},
            {"label": "Expanded (279 struct)", "n": 279, "spearman_rho": 0.478, "source": "repo_historical"},
            {"label": "Colab (~200 target)", "n": colab_200["n_success"], "spearman_rho": colab_200["spearman_rho"], "source": "colab"},
            {"label": "Paper headline (362 clean)", "n": chainfix["n_success_clean"], "spearman_rho": chainfix["spearman_rho_clean"], "source": "repo_manuscript"},
        ],
    }


def main() -> None:
    manifest = build_manifest()
    out = RESULTS / "colab_integrated_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
