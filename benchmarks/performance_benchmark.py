#!/usr/bin/env python3
"""ChiralFold performance benchmark — throughput and scaling of the core paths.

Measures wall time and peak memory for the operations that dominate a
PDB-scale survey:

* ``audit_pdb`` and each of its five checks
* clash detection alone (the auditor's hot path)
* ``score_interface`` on receptor:ligand complexes
* ``detect_chirality_violations`` / ``correct_af3_output``
* ``mirror_pdb``

Scaling inputs are built by tiling a real deposited structure on a lattice, so
the benchmark exercises realistic covalent topology, altloc handling and
neighbour density rather than random point clouds. Copies are separated far
enough that no new inter-copy contacts appear, which keeps the per-atom work
constant and makes the wall-time curve a clean read on algorithmic scaling.

Usage::

    python benchmarks/performance_benchmark.py
    python benchmarks/performance_benchmark.py --repeats 7 --out results/perf.json
    python benchmarks/performance_benchmark.py --quick     # small sizes only

Output: JSON report (default ``results/performance_benchmark.json``) plus a
human-readable table on stdout. No network access required.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402

from chiralfold import auditor as auditor_mod  # noqa: E402
from chiralfold.af3_correct import (  # noqa: E402
    correct_af3_output,
    detect_chirality_violations,
)
from chiralfold.auditor import audit_pdb  # noqa: E402
from chiralfold.interface_scorer import score_interface  # noqa: E402
from chiralfold.pdb_pipeline import mirror_pdb  # noqa: E402

#: Real structure tiled to build the scaling series (p53:MDM2 D-peptide complex).
BASE_STRUCTURE = REPO_ROOT / "results" / "3IWY.pdb"

#: Lattice spacing (Å) between tiled copies — well beyond any vdW contact.
TILE_SPACING = 80.0

#: Copy counts for the scaling series.
SCALING_COPIES = (1, 2, 4, 8, 16)
QUICK_COPIES = (1, 2, 4)


# ─────────────────────────────────────────────────────────────────────────────
# Timing helpers
# ─────────────────────────────────────────────────────────────────────────────

def _time_call(fn: Callable[[], object], repeats: int, warmup: int = 1) -> Dict:
    """Run *fn* and return timing statistics in milliseconds.

    The minimum is reported alongside the median: for CPU-bound deterministic
    work the minimum is the least noise-contaminated estimate, while the median
    reflects what a user actually experiences.
    """
    for _ in range(warmup):
        fn()
    samples: List[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return {
        "median_ms": round(statistics.median(samples), 4),
        "min_ms": round(min(samples), 4),
        "max_ms": round(max(samples), 4),
        "stdev_ms": round(statistics.stdev(samples), 4) if len(samples) > 1 else 0.0,
        "repeats": repeats,
    }


def _peak_rss_mb() -> float:
    """Peak resident set size of this process in MB."""
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports KiB, macOS reports bytes.
    return peak / 1024.0 if sys.platform != "darwin" else peak / (1024.0 * 1024.0)


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic scaling inputs built from a real structure
# ─────────────────────────────────────────────────────────────────────────────

def _tile_structure(src: Path, dest: Path, n_copies: int,
                    spacing: float = TILE_SPACING) -> int:
    """Write *n_copies* translated copies of *src* to *dest*; return atom count.

    Each copy gets its own chain-ID suffix so residues stay distinguishable, and
    copies are placed on a cubic lattice with *spacing* Å between neighbours.
    """
    chain_alphabet = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" "abcdefghijklmnopqrstuvwxyz" "0123456789"
    )
    src_lines = [
        ln for ln in src.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    ]
    side = int(np.ceil(n_copies ** (1.0 / 3.0)))

    out: List[str] = []
    serial = 0
    chain_map: Dict[str, str] = {}
    next_chain = 0

    for copy_idx in range(n_copies):
        cx = (copy_idx % side) * spacing
        cy = ((copy_idx // side) % side) * spacing
        cz = (copy_idx // (side * side)) * spacing
        for ln in src_lines:
            padded = ln.ljust(80)
            try:
                x = float(padded[30:38]) + cx
                y = float(padded[38:46]) + cy
                z = float(padded[46:54]) + cz
            except ValueError:
                continue
            src_chain = padded[21]
            key = f"{copy_idx}:{src_chain}"
            if key not in chain_map:
                chain_map[key] = chain_alphabet[next_chain % len(chain_alphabet)]
                next_chain += 1
            serial += 1
            out.append(
                f"{padded[:6]}{serial % 100000:>5}{padded[11:21]}"
                f"{chain_map[key]}{padded[22:30]}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}{padded[54:80].rstrip()}"
            )
    dest.write_text("\n".join(out) + "\nEND\n")
    return serial


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark sections
# ─────────────────────────────────────────────────────────────────────────────

def bench_audit_scaling(workdir: Path, repeats: int, copies) -> List[Dict]:
    """Full audit + per-check breakdown across the tiled scaling series."""
    rows: List[Dict] = []
    for n in copies:
        path = workdir / f"tiled_{n}.pdb"
        n_atoms = _tile_structure(BASE_STRUCTURE, path, n)
        p = str(path)

        report = audit_pdb(p)
        atoms = auditor_mod._parse_pdb(p)
        groups = auditor_mod._group_by_residue(atoms)
        keys = list(groups.keys())
        table = auditor_mod._ResidueTable(groups, keys)

        row = {
            "copies": n,
            "n_atoms": report["n_atoms"],
            "n_residues": report["n_residues"],
            "source_records": n_atoms,
            "audit_pdb": _time_call(lambda p=p: audit_pdb(p), repeats),
            "parse": _time_call(lambda p=p: auditor_mod._parse_pdb(p), repeats),
            "residue_table": _time_call(
                lambda g=groups, k=keys: auditor_mod._ResidueTable(g, k), repeats
            ),
            "chirality": _time_call(
                lambda t=table: auditor_mod._check_chirality(t), repeats
            ),
            "bond_geometry": _time_call(
                lambda t=table: auditor_mod._check_bond_geometry(t), repeats
            ),
            "ramachandran": _time_call(
                lambda t=table: auditor_mod._check_ramachandran(t), repeats
            ),
            "planarity": _time_call(
                lambda t=table: auditor_mod._check_planarity(t), repeats
            ),
            "clashes": _time_call(
                lambda a=atoms: auditor_mod._check_clashes(a), repeats
            ),
            "overall_score": report["overall_score"],
            "clash_score": report["clashes"]["clash_score"],
        }
        row["us_per_atom"] = round(
            row["audit_pdb"]["median_ms"] * 1000.0 / max(report["n_atoms"], 1), 3
        )
        rows.append(row)
    return rows


def bench_interface(workdir: Path, repeats: int, copies) -> List[Dict]:
    """Interface scoring between two tiled copies placed in contact."""
    rows: List[Dict] = []
    for n in copies:
        rec = workdir / f"iface_rec_{n}.pdb"
        lig = workdir / f"iface_lig_{n}.pdb"
        _tile_structure(BASE_STRUCTURE, rec, n)
        _tile_structure(BASE_STRUCTURE, lig, n)

        # Shift the ligand copy slightly so the two sets interpenetrate and the
        # contact search has real work to do.
        shifted = []
        for ln in lig.read_text().splitlines():
            if ln.startswith(("ATOM  ", "HETATM")):
                padded = ln.ljust(80)
                x = float(padded[30:38]) + 6.0
                shifted.append(
                    f"{padded[:21]}Z{padded[22:30]}{x:8.3f}{padded[38:80].rstrip()}"
                )
            else:
                shifted.append(ln)
        lig.write_text("\n".join(shifted) + "\n")

        result = score_interface(str(rec), str(lig))
        rows.append({
            "copies": n,
            "n_interface_pairs": result["n_interface_pairs"],
            "salt_bridges": result["salt_bridges"],
            "hbonds": result["hbonds"],
            "score_interface": _time_call(
                lambda r=str(rec), l=str(lig): score_interface(r, l), repeats
            ),
        })
    return rows


def bench_pipelines(workdir: Path, repeats: int) -> Dict:
    """Chirality detection, AF3 correction and mirroring on the base structure."""
    src = str(BASE_STRUCTURE)
    out_corrected = str(workdir / "corrected.pdb")
    out_mirror = str(workdir / "mirrored.pdb")

    detect = detect_chirality_violations(src)
    return {
        "structure": BASE_STRUCTURE.name,
        "n_checked": detect["n_checked"],
        "n_violations": detect["n_violations"],
        "n_unassignable": detect["n_unassignable"],
        "detect_chirality_violations": _time_call(
            lambda: detect_chirality_violations(src), repeats
        ),
        "correct_af3_output": _time_call(
            lambda: correct_af3_output(src, out_corrected), repeats
        ),
        "mirror_pdb": _time_call(lambda: mirror_pdb(src, out_mirror), repeats),
    }


def bench_ramachandran_batch(repeats: int) -> Dict:
    """Batched vs per-residue Ramachandran classification."""
    from chiralfold.ramachandran import classify_regions, score_ramachandran

    rng = np.random.default_rng(20260731)
    n = 20000
    phi = rng.uniform(-180, 180, size=n)
    psi = rng.uniform(-180, 180, size=n)
    rtypes = ["general"] * n

    batched = _time_call(lambda: classify_regions(phi, psi, rtypes), repeats)
    scalar = _time_call(
        lambda: [score_ramachandran(a, b, "general") for a, b in zip(phi, psi)],
        max(1, repeats // 3),
    )
    return {
        "n_points": n,
        "batched": batched,
        "per_residue": scalar,
        "speedup": round(scalar["median_ms"] / max(batched["median_ms"], 1e-9), 2),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────────────

def _git_commit() -> Optional[str]:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT, stderr=subprocess.DEVNULL, text=True,
        ).strip()
    except Exception:
        return None


def _environment() -> Dict:
    import scipy

    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "git_commit": _git_commit(),
    }


def _print_table(report: Dict) -> None:
    print("\n" + "=" * 78)
    print("ChiralFold performance benchmark")
    print("=" * 78)

    print("\n-- audit_pdb scaling (tiled 3IWY) " + "-" * 43)
    print(f"{'atoms':>8} {'residues':>9} {'audit ms':>10} {'us/atom':>9} "
          f"{'clash ms':>9} {'rama ms':>8} {'parse ms':>9}")
    for row in report["audit_scaling"]:
        print(
            f"{row['n_atoms']:>8} {row['n_residues']:>9} "
            f"{row['audit_pdb']['median_ms']:>10.2f} {row['us_per_atom']:>9.2f} "
            f"{row['clashes']['median_ms']:>9.2f} "
            f"{row['ramachandran']['median_ms']:>8.2f} "
            f"{row['parse']['median_ms']:>9.2f}"
        )

    print("\n-- score_interface scaling " + "-" * 50)
    print(f"{'copies':>7} {'contacts':>10} {'ms':>10}")
    for row in report["interface_scaling"]:
        print(f"{row['copies']:>7} {row['n_interface_pairs']:>10} "
              f"{row['score_interface']['median_ms']:>10.2f}")

    print("\n-- pipelines (3IWY) " + "-" * 57)
    p = report["pipelines"]
    for key in ("detect_chirality_violations", "correct_af3_output", "mirror_pdb"):
        print(f"{key:>30}: {p[key]['median_ms']:8.2f} ms")

    r = report["ramachandran_batch"]
    print(f"\n-- Ramachandran classification ({r['n_points']} points) " + "-" * 24)
    print(f"{'batched':>30}: {r['batched']['median_ms']:8.2f} ms")
    print(f"{'per-residue':>30}: {r['per_residue']['median_ms']:8.2f} ms")
    print(f"{'speedup':>30}: {r['speedup']:8.2f}x")
    print(f"\npeak RSS: {report['peak_rss_mb']:.0f} MB")
    print("=" * 78)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5,
                        help="timed repetitions per measurement (default 5)")
    parser.add_argument("--quick", action="store_true",
                        help="use the small scaling series only")
    parser.add_argument(
        "--out", default=str(REPO_ROOT / "results" / "performance_benchmark.json"),
        help="path for the JSON report",
    )
    args = parser.parse_args(argv)

    if not BASE_STRUCTURE.is_file():
        print(f"ERROR: base structure not found: {BASE_STRUCTURE}", file=sys.stderr)
        return 1

    copies = QUICK_COPIES if args.quick else SCALING_COPIES

    with tempfile.TemporaryDirectory(prefix="chiralfold_perf_") as tmp:
        workdir = Path(tmp)
        report = {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "environment": _environment(),
            "base_structure": BASE_STRUCTURE.name,
            "tile_spacing_angstrom": TILE_SPACING,
            "audit_scaling": bench_audit_scaling(workdir, args.repeats, copies),
            "interface_scaling": bench_interface(workdir, args.repeats, copies),
            "pipelines": bench_pipelines(workdir, args.repeats),
            "ramachandran_batch": bench_ramachandran_batch(args.repeats),
        }
        report["peak_rss_mb"] = round(_peak_rss_mb(), 1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n")

    _print_table(report)
    print(f"\nJSON report written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
