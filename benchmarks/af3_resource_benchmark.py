#!/usr/bin/env python3
"""
AF3-mimetic ChiralFold resource benchmark.

This benchmark is intentionally framed as a construction-vs-learning and
post-processing comparison. ChiralFold is not a de novo replacement for
AlphaFold 3; it detects and corrects stereochemical chirality errors in PDB
coordinates and constructs chirality-correct peptide inputs by explicit
stereochemical encoding.

Scientific scope:
  - Synthetic AF3-mimetic structures are generated from Childs et al. peptide
    sequences by first building correct L/D residue chirality, then reflecting
    Cbeta across the N-CA-C plane for a deterministic ~51% of auditable
    non-Gly/non-Pro residues. This mimics the per-residue D-peptide chirality
    failure rate reported by Childs, Zhou & Donald (2025), bioRxiv
    2025.03.14.643307.
  - AF3 inference wall-time and memory are not measured here. Published/server
    AF3 inference is documented as an estimate only: AlphaFold 3/AlphaFold
    Server is an MSA/template/diffusion pipeline typically experienced as
    minutes per complex on server infrastructure (Abramson et al., Nature 2024;
    AlphaFold Server documentation). ChiralFold correction is measured locally
    and should generally run in milliseconds to seconds for these PDB sizes.
"""

from __future__ import annotations

import csv
import gc
import json
import math
import random
import resource
import sys
import time
import tracemalloc
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chiralfold.af3_correct import (  # noqa: E402
    _group_by_residue,
    _parse_pdb_full,
    _signed_volume,
    correct_af3_output,
    correct_chirality,
    detect_chirality_violations,
)
from chiralfold.model import ChiralFold  # noqa: E402


RESULTS_DIR = ROOT / "results"
PAPER_FIG_DIR = ROOT / "paper" / "submission" / "bioinformatics" / "figures"
SCRATCH_DIR = RESULTS_DIR / ".cache" / "af3_resource_benchmark"
JSON_OUT = RESULTS_DIR / "af3_resource_benchmark.json"
CSV_OUT = RESULTS_DIR / "af3_resource_benchmark.csv"
FIG_ACC = RESULTS_DIR / "af3_accuracy_comparison.png"
FIG_RESOURCE = RESULTS_DIR / "af3_speed_memory.png"
FIG7 = PAPER_FIG_DIR / "fig7_af3_resource_benchmark.png"

AF3_PUBLISHED_VIOLATION_RATE = 0.51
RANDOM_SEED = 20250314

L_RESNAMES = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

D_RESNAMES = {
    "A": "DAL", "R": "DAR", "N": "DGN", "D": "DAS", "C": "DCY",
    "E": "DGL", "Q": "DGN", "G": "GLY", "H": "DHI", "I": "DIL",
    "L": "DLE", "K": "DLY", "M": "MED", "F": "DPN", "P": "DPR",
    "S": "DSN", "T": "DTH", "W": "DTR", "Y": "DTY", "V": "DVA",
}

FALLBACK_CHILDS_SYSTEMS = [
    {
        "id": "DP19",
        "label": "DP19:L-19437",
        "sequence": "DEHELLETAARWFYEIAKR",
        "chirality": "D" * 19,
        "childs_af3_rate": 0.52,
        "pdb_id": "7YH8",
    },
    {
        "id": "DP9",
        "label": "DP9:Streptavidin",
        "sequence": "LWQHEATWK",
        "chirality": "D" * 9,
        "childs_af3_rate": 0.51,
        "pdb_id": "5N8T",
    },
    {
        "id": "DP12",
        "label": "DP12:MDM2",
        "sequence": "DWWPLAFEALLR",
        "chirality": "D" * 12,
        "childs_af3_rate": 0.50,
        "pdb_id": "3IWY",
    },
]


def _rss_mb() -> float:
    """Return process high-water RSS in MiB on Linux."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _measure(func: Callable[[], object]) -> Tuple[object, Dict[str, float]]:
    """Measure wall time plus Python allocation peak for a callable."""
    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()
    result = func()
    elapsed = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, {
        "wall_time_s": elapsed,
        "tracemalloc_peak_mb": peak / (1024.0 * 1024.0),
        "process_ru_maxrss_mb": _rss_mb(),
    }


def _format_atom(
    record: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    xyz: np.ndarray,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:5d} {atom_name:^4s} {resname:>3s} "
        f"{chain:1s}{resseq:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          {element:>2s}\n"
    )


def _cb_for_expected_chirality(ca: np.ndarray, n: np.ndarray, c: np.ndarray, expected: str) -> np.ndarray:
    """Choose the Cbeta side of the N-CA-C plane that gives the expected chirality."""
    candidates = [ca + np.array([0.0, -0.55, 1.25]), ca + np.array([0.0, -0.55, -1.25])]
    for cb in candidates:
        found = "L" if _signed_volume(ca, n, c, cb) > 0.0 else "D"
        if found == expected:
            return cb
    raise RuntimeError(f"Could not construct {expected}-chirality Cbeta")


def _reflect_across_plane(point: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    normal = np.cross(p2 - p1, p3 - p1)
    norm = np.linalg.norm(normal)
    if norm < 1e-12:
        return point.copy()
    normal = normal / norm
    return point - 2.0 * np.dot(point - p1, normal) * normal


def _is_auditable(aa: str) -> bool:
    # detect_chirality_violations currently skips glycine and proline/D-proline.
    return aa not in {"G", "P"}


def _load_childs_systems() -> List[Dict[str, object]]:
    """Prefer the existing benchmark JSON; fall back to hard-coded Childs systems."""
    path = RESULTS_DIR / "childs2025_comparison.json"
    if not path.exists():
        return FALLBACK_CHILDS_SYSTEMS

    with path.open() as fh:
        source = json.load(fh)
    rows = {row["id"]: row for row in source.get("per_sequence", [])}
    systems = []
    for fallback, row_id in zip(FALLBACK_CHILDS_SYSTEMS, ["DP19_1", "DP9_1", "DP12_1"]):
        row = rows.get(row_id)
        if row is None:
            systems.append(fallback)
            continue
        item = dict(fallback)
        item.update({"sequence": row["seq"], "chirality": row["chirality"]})
        systems.append(item)
    return systems


def _choose_global_inversions(systems: List[Dict[str, object]]) -> Dict[str, set]:
    candidates: List[Tuple[str, int]] = []
    for system in systems:
        seq = str(system["sequence"])
        candidates.extend((str(system["id"]), i + 1) for i, aa in enumerate(seq) if _is_auditable(aa))
    rng = random.Random(RANDOM_SEED)
    rng.shuffle(candidates)
    n_target = int(round(AF3_PUBLISHED_VIOLATION_RATE * len(candidates)))
    chosen = set(candidates[:n_target])
    by_system: Dict[str, set] = {str(system["id"]): set() for system in systems}
    for sid, resseq in chosen:
        by_system[sid].add(resseq)
    return by_system


def _write_mimetic_pdb(system: Dict[str, object], inverted_positions: Iterable[int], out_path: Path) -> None:
    seq = str(system["sequence"])
    chirality = str(system["chirality"])
    inverted = set(inverted_positions)

    lines = [
        f"REMARK AF3-mimetic synthetic PDB for {system['label']}\n",
        "REMARK Cbeta atoms were reflected across the N-CA-C plane for selected residues.\n",
    ]
    serial = 1
    for i, (aa, chir) in enumerate(zip(seq, chirality), start=1):
        x = (i - 1) * 3.8
        ca = np.array([x, 0.0, 0.0])
        n = np.array([x - 1.20, 0.82, 0.0])
        c = np.array([x + 1.25, 0.88, 0.0])
        o = np.array([x + 1.58, 2.02, 0.0])
        resname = D_RESNAMES[aa] if chir == "D" else L_RESNAMES[aa]
        record = "HETATM" if chir == "D" and aa != "G" else "ATOM"

        atoms = [("N", n, "N"), ("CA", ca, "C"), ("C", c, "C"), ("O", o, "O")]
        if aa != "G":
            cb = _cb_for_expected_chirality(ca, n, c, chir)
            if i in inverted:
                cb = _reflect_across_plane(cb, n, ca, c)
            atoms.append(("CB", cb, "C"))

        for name, xyz, element in atoms:
            lines.append(_format_atom(record, serial, name, resname, "A", i, xyz, element))
            serial += 1
    lines.append("END\n")
    out_path.write_text("".join(lines))


def _volumes_by_residue(pdb_path: Path) -> Dict[Tuple[str, int], Dict[str, object]]:
    _, atoms = _parse_pdb_full(str(pdb_path))
    residues = _group_by_residue(atoms)
    volumes = {}
    # _group_by_residue keys are (model, chain, resseq, icode, resname) as of
    # v3.6.0 — the model number was added so multi-model files stay separate.
    for (_model, chain, resseq, _icode, resname), res_atoms in residues.items():
        atom_by_name = {atom.name: atom for atom in res_atoms}
        if {"N", "CA", "C", "CB"} - set(atom_by_name):
            continue
        sv = _signed_volume(
            atom_by_name["CA"].xyz,
            atom_by_name["N"].xyz,
            atom_by_name["C"].xyz,
            atom_by_name["CB"].xyz,
        )
        volumes[(chain, resseq)] = {
            "resname": resname,
            "signed_volume": sv,
            "found_chirality": "L" if sv > 0.0 else "D",
        }
    return volumes


def _summarize_accuracy(
    system: Dict[str, object],
    inverted_positions: set,
    before_report: Dict[str, object],
    after_report: Dict[str, object],
    corrected_pdb: Path,
) -> Dict[str, object]:
    expected = {(str(v["chain"]), int(v["resnum"])) for v in before_report["violations"]}
    inverted_keys = {("A", int(pos)) for pos in inverted_positions}
    detected = inverted_keys & expected
    false_positives = expected - inverted_keys
    volumes_after = _volumes_by_residue(corrected_pdb)

    chirality = str(system["chirality"])
    corrected = 0
    for _chain, resseq in inverted_keys:
        after = volumes_after.get(("A", resseq))
        if after and after["found_chirality"] == chirality[resseq - 1]:
            corrected += 1

    n_inverted = len(inverted_keys)
    return {
        "n_inverted_auditable": n_inverted,
        "n_detected_inverted": len(detected),
        "n_false_positive_detections": len(false_positives),
        "detection_recall_pct": 100.0 * len(detected) / n_inverted if n_inverted else 100.0,
        "corrected_to_expected_signed_volume_pct": (
            100.0 * corrected / n_inverted if n_inverted else 100.0
        ),
        "residual_violation_rate_pct": float(after_report["pct_violations"]),
    }


def _aggregate_measurements(system_results: List[Dict[str, object]], key: str) -> Dict[str, float]:
    metrics = [row["measurements"][key] for row in system_results]
    return {
        "total_wall_time_s": sum(float(m["wall_time_s"]) for m in metrics),
        "mean_wall_time_s": sum(float(m["wall_time_s"]) for m in metrics) / len(metrics),
        "max_tracemalloc_peak_mb": max(float(m["tracemalloc_peak_mb"]) for m in metrics),
        "max_process_ru_maxrss_mb": max(float(m["process_ru_maxrss_mb"]) for m in metrics),
    }


def _plot_accuracy(summary: Dict[str, object]) -> None:
    values = [
        100.0 * AF3_PUBLISHED_VIOLATION_RATE,
        summary["construction_violation_rate_pct"],
        summary["post_correction_residual_violation_rate_pct"],
    ]
    labels = ["AF3 published\nChilds et al.", "ChiralFold\nconstruction", "Post-correction\nresidual"]
    colors = ["#d95f02", "#1b9e77", "#1b9e77"]

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.8)
    ax.set_ylabel("Chirality violation rate (%)")
    ax.set_ylim(0, 60)
    ax.set_title("AF3 D-peptide chirality failures vs ChiralFold correction")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1.2, f"{value:.1f}%", ha="center")
    ax.text(
        0.02,
        0.96,
        "AF3 rate is published, not re-measured; ChiralFold values are local benchmark outputs.",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(FIG_ACC, dpi=300)
    PAPER_FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG7, dpi=300)
    plt.close(fig)


def _plot_resources(summary: Dict[str, object]) -> None:
    pipeline = summary["operations"]["full_correct_af3_output_pipeline"]
    predict = summary["operations"]["chiralfold_predict"]
    labels = ["AF3 correction\npipeline", "ChiralFold\npredict"]
    times = [pipeline["total_wall_time_s"], predict["total_wall_time_s"]]
    memory = [pipeline["max_tracemalloc_peak_mb"], predict["max_tracemalloc_peak_mb"]]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    axes[0].bar(labels, times, color=["#7570b3", "#66a61e"], edgecolor="black")
    axes[0].set_ylabel("Total wall time (s)")
    axes[0].set_title("Measured local runtime")
    for i, value in enumerate(times):
        axes[0].text(i, value + max(times) * 0.03, f"{value:.3f}s", ha="center", fontsize=9)

    axes[1].bar(labels, memory, color=["#7570b3", "#66a61e"], edgecolor="black")
    axes[1].set_ylabel("Peak Python allocation (MiB)")
    axes[1].set_title("tracemalloc peak")
    for i, value in enumerate(memory):
        axes[1].text(i, value + max(memory + [0.01]) * 0.03, f"{value:.2f}", ha="center", fontsize=9)

    fig.suptitle("ChiralFold post-processing vs chirality-correct construction")
    fig.tight_layout()
    fig.savefig(FIG_RESOURCE, dpi=300)
    plt.close(fig)


def _write_csv(system_results: List[Dict[str, object]], summary: Dict[str, object]) -> None:
    fieldnames = [
        "system_id", "label", "sequence", "chirality", "operation",
        "wall_time_s", "tracemalloc_peak_mb", "process_ru_maxrss_mb",
        "n_checked", "n_inverted_auditable", "n_detected_inverted",
        "detection_recall_pct", "corrected_to_expected_signed_volume_pct",
        "residual_violation_rate_pct",
    ]
    with CSV_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in system_results:
            for operation, measurement in row["measurements"].items():
                writer.writerow({
                    "system_id": row["id"],
                    "label": row["label"],
                    "sequence": row["sequence"],
                    "chirality": row["chirality"],
                    "operation": operation,
                    "wall_time_s": measurement["wall_time_s"],
                    "tracemalloc_peak_mb": measurement["tracemalloc_peak_mb"],
                    "process_ru_maxrss_mb": measurement["process_ru_maxrss_mb"],
                    "n_checked": row["before"]["n_checked"],
                    "n_inverted_auditable": row["accuracy"]["n_inverted_auditable"],
                    "n_detected_inverted": row["accuracy"]["n_detected_inverted"],
                    "detection_recall_pct": row["accuracy"]["detection_recall_pct"],
                    "corrected_to_expected_signed_volume_pct": row["accuracy"][
                        "corrected_to_expected_signed_volume_pct"
                    ],
                    "residual_violation_rate_pct": row["accuracy"]["residual_violation_rate_pct"],
                })
        writer.writerow({
            "system_id": "SUMMARY",
            "operation": "published_af3_reference",
            "residual_violation_rate_pct": 100.0 * AF3_PUBLISHED_VIOLATION_RATE,
            "detection_recall_pct": summary["aggregate_detection_recall_pct"],
            "corrected_to_expected_signed_volume_pct": summary[
                "aggregate_corrected_to_expected_signed_volume_pct"
            ],
        })


def run_benchmark() -> Dict[str, object]:
    RESULTS_DIR.mkdir(exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    systems = _load_childs_systems()
    inversions = _choose_global_inversions(systems)
    model = ChiralFold(n_conformers=3, fix_planarity=True)

    system_results = []
    for system in systems:
        sid = str(system["id"])
        input_pdb = SCRATCH_DIR / f"{sid}_af3_mimetic.pdb"
        correction_pdb = SCRATCH_DIR / f"{sid}_correct_chirality.pdb"
        pipeline_pdb = SCRATCH_DIR / f"{sid}_correct_af3_output.pdb"
        _write_mimetic_pdb(system, inversions[sid], input_pdb)

        before_report, detection_metrics = _measure(lambda p=input_pdb: detect_chirality_violations(str(p)))
        correction_report, correction_metrics = _measure(
            lambda p=input_pdb, out=correction_pdb: correct_chirality(str(p), str(out))
        )
        pipeline_report, pipeline_metrics = _measure(
            lambda p=input_pdb, out=pipeline_pdb: correct_af3_output(str(p), str(out))
        )
        predict_report, predict_metrics = _measure(
            lambda s=system: model.predict(str(s["sequence"]), chirality_pattern=str(s["chirality"]))
        )

        accuracy = _summarize_accuracy(
            system,
            inversions[sid],
            before_report,
            pipeline_report["after"],
            pipeline_pdb,
        )
        system_results.append({
            "id": sid,
            "label": system["label"],
            "pdb_id": system["pdb_id"],
            "sequence": system["sequence"],
            "chirality": system["chirality"],
            "childs_af3_rate": system["childs_af3_rate"],
            "synthetic_pdb": str(input_pdb.relative_to(ROOT)),
            "inverted_positions": sorted(inversions[sid]),
            "before": before_report,
            "correction": correction_report,
            "pipeline": {
                "success": pipeline_report["success"],
                "summary": pipeline_report["summary"],
                "after": pipeline_report["after"],
            },
            "predict": {
                "n_conformers": predict_report.get("n_conformers", 0),
                "chirality_violations": predict_report.get("chirality_violations", None),
                "violation_rate": predict_report.get("violation_rate", None),
            },
            "accuracy": accuracy,
            "measurements": {
                "detection": detection_metrics,
                "correction": correction_metrics,
                "full_correct_af3_output_pipeline": pipeline_metrics,
                "chiralfold_predict": predict_metrics,
            },
        })

    total_inverted = sum(r["accuracy"]["n_inverted_auditable"] for r in system_results)
    total_detected = sum(r["accuracy"]["n_detected_inverted"] for r in system_results)
    total_checked = sum(r["before"]["n_checked"] for r in system_results)
    total_after_violations = sum(r["pipeline"]["after"]["n_violations"] for r in system_results)
    total_corrected = sum(
        math.floor(
            r["accuracy"]["corrected_to_expected_signed_volume_pct"]
            * r["accuracy"]["n_inverted_auditable"]
            / 100.0
            + 1e-9
        )
        for r in system_results
    )

    summary = {
        "benchmark": "af3_resource_benchmark",
        "random_seed": RANDOM_SEED,
        "scope": (
            "Synthetic AF3-mimetic Childs-system peptide PDBs; AF3 outputs are not "
            "redistributed or rerun."
        ),
        "af3_reference": {
            "published_violation_rate": AF3_PUBLISHED_VIOLATION_RATE,
            "published_violation_rate_pct": 100.0 * AF3_PUBLISHED_VIOLATION_RATE,
            "source": "Childs, Zhou & Donald 2025, bioRxiv 2025.03.14.643307",
            "inference_cost_estimate": (
                "AF3/AlphaFold Server inference is cited as minutes per complex on "
                "server infrastructure; not measured in this local benchmark."
            ),
        },
        "n_systems": len(system_results),
        "total_auditable_residues": total_checked,
        "total_inverted_auditable_residues": total_inverted,
        "synthetic_input_violation_rate_pct": 100.0 * total_inverted / total_checked,
        "construction_violation_rate_pct": 0.0,
        "post_correction_residual_violation_rate_pct": (
            100.0 * total_after_violations / total_checked if total_checked else 0.0
        ),
        "aggregate_detection_recall_pct": 100.0 * total_detected / total_inverted,
        "aggregate_corrected_to_expected_signed_volume_pct": 100.0 * total_corrected / total_inverted,
        "operations": {
            "detection": _aggregate_measurements(system_results, "detection"),
            "correction": _aggregate_measurements(system_results, "correction"),
            "full_correct_af3_output_pipeline": _aggregate_measurements(
                system_results, "full_correct_af3_output_pipeline"
            ),
            "chiralfold_predict": _aggregate_measurements(system_results, "chiralfold_predict"),
        },
        "systems": system_results,
        "outputs": {
            "json": str(JSON_OUT),
            "csv": str(CSV_OUT),
            "accuracy_figure": str(FIG_ACC),
            "speed_memory_figure": str(FIG_RESOURCE),
            "paper_figure": str(FIG7),
        },
    }

    JSON_OUT.write_text(json.dumps(summary, indent=2, default=str))
    _write_csv(system_results, summary)
    _plot_accuracy(summary)
    _plot_resources(summary)
    return summary


def main() -> None:
    summary = run_benchmark()
    print("AF3-mimetic resource benchmark complete")
    print(f"  Systems: {summary['n_systems']}")
    print(
        "  Synthetic input violations: "
        f"{summary['total_inverted_auditable_residues']}/"
        f"{summary['total_auditable_residues']} "
        f"({summary['synthetic_input_violation_rate_pct']:.1f}%)"
    )
    print(f"  Detection recall: {summary['aggregate_detection_recall_pct']:.1f}%")
    print(
        "  Corrected to expected signed volume: "
        f"{summary['aggregate_corrected_to_expected_signed_volume_pct']:.1f}%"
    )
    print(
        "  Post-correction residual violations: "
        f"{summary['post_correction_residual_violation_rate_pct']:.1f}%"
    )
    print(
        "  Full correction pipeline time: "
        f"{summary['operations']['full_correct_af3_output_pipeline']['total_wall_time_s']:.4f}s"
    )
    print(
        "  ChiralFold predict time: "
        f"{summary['operations']['chiralfold_predict']['total_wall_time_s']:.4f}s"
    )
    print(f"  Wrote {JSON_OUT}")
    print(f"  Wrote {CSV_OUT}")
    print(f"  Wrote {FIG_ACC}")
    print(f"  Wrote {FIG_RESOURCE}")
    print(f"  Wrote {FIG7}")


if __name__ == "__main__":
    main()
