"""Schema / smoke tests that Colab dashboards and notebooks stay publication-ready."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_af3_benchmark_schema_matches_dashboard():
    """Dashboard cell 9 must use correction.violations_after (not s['after'])."""
    af3 = json.loads((ROOT / "results/af3_resource_benchmark.json").read_text())
    assert af3["aggregate_detection_recall_pct"] == 100.0
    for s in af3["systems"]:
        assert "before" in s and "n_violations" in s["before"]
        assert "correction" in s
        assert "violations_after" in s["correction"]
        assert "after" not in s  # guard against regressing to broken key


def test_mmcif_reverification_matches_legacy_survey():
    legacy = json.loads((ROOT / "results/d_residue_verification_summary.json").read_text())
    mmcif = json.loads((ROOT / "results/mmcif_d_residue_expansion_summary.json").read_text())
    assert mmcif["n_structures"] == 16
    assert mmcif["n_errors"] == legacy["l_error"] == 29
    assert set(mmcif["error_pdbs"]) == set(legacy["errors_by_structure"])


def test_colab_manifest_is_current():
    m = json.loads((ROOT / "results/colab_integrated_manifest.json").read_text())
    assert m["chiralfold_version"] == "3.5.1"
    assert "in progress" not in m["interpretation"]["aristotle_lean"].lower()
    assert "mmcif_reverification" in m["interpretation"]
    run_ids = {r["id"] for r in m["runs"]}
    assert "mmcif_reverification" in run_ids
    assert m["headline_metrics"]["d_residue_survey"]["pdb_files"] == 4623


def test_dashboard_notebook_has_fixed_af3_key():
    nb = json.loads((ROOT / "demos/ChiralFold_Results_Dashboard.ipynb").read_text())
    src = "\n".join("".join(c.get("source", [])) for c in nb["cells"])
    assert "violations_after" in src
    assert "s['after']['n_violations']" not in src
    assert "chiralfold==3.5.1" in src or "chiralfold==3.5.1" in src.replace(" ", "")
    assert "Formal verification in progress" not in src


def test_mmcif_colab_notebook_exists():
    path = ROOT / "demos/Reproduce_mmCIF_D_Residue_Survey.ipynb"
    assert path.is_file()
    nb = json.loads(path.read_text())
    src = "\n".join("".join(c.get("source", [])) for c in nb["cells"])
    assert "mmcif_d_residue_expansion.py" in src
    assert "gemmi" in src


def test_demo_notebooks_prefer_pypi_install():
    demos = list((ROOT / "demos").glob("*.ipynb"))
    assert demos
    stale = []
    for p in demos:
        text = p.read_text()
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
    nb = json.loads((ROOT / rel).read_text())
    assert nb["nbformat"] == 4
    assert nb["cells"]
