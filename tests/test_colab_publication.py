"""Schema / smoke tests that Colab dashboards and notebooks stay publication-ready."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _read_text(path: Path) -> str:
    """Read text as UTF-8 (Windows CI defaults to cp1252 without this)."""
    return path.read_text(encoding="utf-8")


def test_af3_benchmark_schema_matches_dashboard():
    """Dashboard cell 9 must use correction.violations_after (not s['after'])."""
    af3 = json.loads(_read_text(ROOT / "results/af3_resource_benchmark.json"))
    assert af3["aggregate_detection_recall_pct"] == 100.0
    for s in af3["systems"]:
        assert "before" in s and "n_violations" in s["before"]
        assert "correction" in s
        assert "violations_after" in s["correction"]
        assert "after" not in s  # guard against regressing to broken key


def test_mmcif_reverification_matches_legacy_survey():
    legacy = json.loads(_read_text(ROOT / "results/d_residue_verification_summary.json"))
    mmcif = json.loads(_read_text(ROOT / "results/mmcif_d_residue_expansion_summary.json"))
    known = mmcif.get("known_error_cohort") or {}
    # Back-compat: older summaries only had top-level n_structures == 16
    n_known = known.get("n_structures", mmcif.get("n_structures"))
    n_known_errors = known.get("n_errors", mmcif.get("n_errors"))
    error_pdbs = set(known.get("error_pdbs") or mmcif.get("error_pdbs") or [])
    assert n_known == 16
    assert n_known_errors == legacy["l_error"] == 29
    assert error_pdbs == set(legacy["errors_by_structure"])
    assert mmcif.get("matches_legacy_survey_errors") in (True, None)


def test_mmcif_universe_survey_artifacts():
    mmcif = json.loads(_read_text(ROOT / "results/mmcif_d_residue_expansion_summary.json"))
    assert mmcif["mode"] in ("universe", "both", "known-errors")
    universe_ids = ROOT / "results/mmcif_only_universe_ids.json"
    if mmcif["mode"] in ("universe", "both"):
        assert universe_ids.is_file()
        disc = json.loads(_read_text(universe_ids))
        assert disc["n_mmcif_only"] >= 1
        assert disc["n_universe_entries"] >= disc["n_mmcif_only"]
        uni = mmcif.get("universe_cohort") or {}
        assert uni.get("n_structures_requested", 0) >= 1


def test_colab_manifest_is_current():
    m = json.loads(_read_text(ROOT / "results/colab_integrated_manifest.json"))
    assert m["chiralfold_version"] == "3.5.1"
    assert "in progress" not in m["interpretation"]["aristotle_lean"].lower()
    assert "mmcif_reverification" in m["interpretation"]
    assert "9BC4" in m["interpretation"]["mmcif_reverification"] or "universe" in m["interpretation"]["mmcif_reverification"].lower()
    run_ids = {r["id"] for r in m["runs"]}
    assert "mmcif_reverification" in run_ids
    assert m["headline_metrics"]["d_residue_survey"]["pdb_files"] == 4623
    assert m["headline_metrics"]["mmcif_universe"]["known_errors_recovered"] == 29
    assert m["headline_metrics"]["mmcif_universe"]["n_mmcif_only"] >= 1


def test_dashboard_notebook_has_fixed_af3_key():
    nb = json.loads(_read_text(ROOT / "demos/ChiralFold_Results_Dashboard.ipynb"))
    src = "\n".join("".join(c.get("source", [])) for c in nb["cells"])
    assert "violations_after" in src
    assert "s['after']['n_violations']" not in src
    assert "chiralfold==3.5.1" in src or "chiralfold==3.5.1" in src.replace(" ", "")
    assert "Formal verification in progress" not in src
    # Cell 4 must use d_survey (not undefined alias d)
    assert "d_survey['checkable_residues']" in src
    assert "d['checkable_residues']" not in src
    assert "colab_integrated_manifest.json" in src
    assert "mmcif_d_residue_expansion_summary.json" in src
    assert "experimental_validation_report.json" in src
    assert "af3_resource_benchmark.json" in src


def test_mmcif_colab_notebook_exists():
    path = ROOT / "demos/Reproduce_mmCIF_D_Residue_Survey.ipynb"
    assert path.is_file()
    nb = json.loads(_read_text(path))
    src = "\n".join("".join(c.get("source", [])) for c in nb["cells"])
    assert "mmcif_d_residue_expansion.py" in src
    assert "gemmi" in src
    assert "--mode" in src or "universe" in src.lower()


def test_demo_notebooks_prefer_pypi_install():
    demos = list((ROOT / "demos").glob("*.ipynb"))
    assert demos
    stale = []
    for p in demos:
        text = _read_text(p)
        if "PyPI package pending" in text:
            stale.append(p.name)
        if "cursor/aristotle-formal-proofs-9901" in text and "in progress" in text.lower():
            stale.append(p.name)
    assert not stale, f"Stale notebook text: {stale}"


@pytest.mark.parametrize(
    "rel",
    [
        "demos/ChiralFold_Results_Dashboard.ipynb",
        "demos/Reproduce_mmCIF_D_Residue_Survey.ipynb",
        "demos/Reproduce_PDB_D_Residue_Errors_5min.ipynb",
        "demos/Demo_Unusual_Cases_Clash_Safety.ipynb",
        "demos/D_Residue_Experimental_Validation.ipynb",
        "demos/ChiralFold_Quick_Demo.ipynb",
        "demos/ChiralFold_Toy_Dataset_Demo.ipynb",
        "demos/Expanded_Ramachandran_Benchmark.ipynb",
    ],
)
def test_demo_notebooks_are_valid_json(rel):
    nb = json.loads(_read_text(ROOT / rel))
    assert nb["nbformat"] == 4
    assert nb["cells"]
