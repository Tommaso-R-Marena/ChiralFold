#!/usr/bin/env python3
"""Compare ChiralFold public-API performance against an earlier git revision.

Creates a throwaway ``git worktree`` at the requested baseline revision, copies
this script into it, and runs the *same* measurement code in both checkouts as
separate subprocesses. Only public API entry points are timed, so the script
runs unchanged against older revisions whose internals differ.

Usage::

    python benchmarks/compare_baseline_performance.py --baseline v3.5.1
    python benchmarks/compare_baseline_performance.py --baseline master --repeats 7

Output: JSON (default ``results/performance_comparison.json``) with per-operation
before/after medians and speedup factors, plus a table on stdout.

Both checkouts import the *same* installed third-party dependencies; only
ChiralFold's own code differs, which is what the comparison is about.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

#: Structures timed in both checkouts (repo-relative, present in every revision).
AUDIT_TARGETS = (
    "results/3IWY.pdb",
    "results/1UBQ_D_mirror.pdb",
    "chiralfold/data/examples/toy_ubiquitin_fragment.pdb",
)

#: Tiled-copy counts for the scaling comparison.
TILE_COPIES = (1, 4, 16)

TILE_SPACING = 80.0


# ─────────────────────────────────────────────────────────────────────────────
# Child mode — runs inside each checkout, public API only
# ─────────────────────────────────────────────────────────────────────────────

def _tile_structure(src: Path, dest: Path, n_copies: int) -> None:
    """Tile *src* on a cubic lattice (identical logic in both checkouts)."""
    import numpy as np

    alphabet = (
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
        cx = (copy_idx % side) * TILE_SPACING
        cy = ((copy_idx // side) % side) * TILE_SPACING
        cz = (copy_idx // (side * side)) * TILE_SPACING
        for ln in src_lines:
            padded = ln.ljust(80)
            try:
                x = float(padded[30:38]) + cx
                y = float(padded[38:46]) + cy
                z = float(padded[46:54]) + cz
            except ValueError:
                continue
            key = f"{copy_idx}:{padded[21]}"
            if key not in chain_map:
                chain_map[key] = alphabet[next_chain % len(alphabet)]
                next_chain += 1
            serial += 1
            out.append(
                f"{padded[:6]}{serial % 100000:>5}{padded[11:21]}"
                f"{chain_map[key]}{padded[22:30]}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}{padded[54:80].rstrip()}"
            )
    dest.write_text("\n".join(out) + "\nEND\n")


def _time(fn, repeats: int, warmup: int = 1) -> Dict:
    for _ in range(warmup):
        fn()
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return {
        "median_ms": round(statistics.median(samples), 4),
        "min_ms": round(min(samples), 4),
        "repeats": repeats,
    }


def run_child(repo_root: Path, repeats: int) -> Dict:
    """Time the public API inside whichever checkout is on ``sys.path``."""
    import resource

    from chiralfold.af3_correct import correct_af3_output, detect_chirality_violations
    from chiralfold.auditor import audit_pdb
    from chiralfold.interface_scorer import score_interface
    from chiralfold.pdb_pipeline import mirror_pdb

    results: Dict = {"audit": {}, "audit_scaling": [], "interface_scaling": []}

    for rel in AUDIT_TARGETS:
        path = repo_root / rel
        if not path.is_file():
            continue
        report = audit_pdb(str(path))
        results["audit"][rel] = {
            "n_atoms": report["n_atoms"],
            "timing": _time(lambda p=str(path): audit_pdb(p), repeats),
        }

    with tempfile.TemporaryDirectory(prefix="cf_cmp_") as tmp:
        work = Path(tmp)
        base = repo_root / "results" / "3IWY.pdb"

        for n in TILE_COPIES:
            tiled = work / f"tiled_{n}.pdb"
            _tile_structure(base, tiled, n)
            report = audit_pdb(str(tiled))
            results["audit_scaling"].append({
                "copies": n,
                "n_atoms": report["n_atoms"],
                "timing": _time(lambda p=str(tiled): audit_pdb(p), repeats),
            })

            rec = work / f"rec_{n}.pdb"
            lig = work / f"lig_{n}.pdb"
            _tile_structure(base, rec, n)
            shifted = []
            for ln in rec.read_text().splitlines():
                if ln.startswith(("ATOM  ", "HETATM")):
                    padded = ln.ljust(80)
                    x = float(padded[30:38]) + 6.0
                    shifted.append(
                        f"{padded[:21]}Z{padded[22:30]}{x:8.3f}"
                        f"{padded[38:80].rstrip()}"
                    )
            lig.write_text("\n".join(shifted) + "\n")
            results["interface_scaling"].append({
                "copies": n,
                "timing": _time(
                    lambda r=str(rec), l=str(lig): score_interface(r, l), repeats
                ),
            })

        src = str(base)
        results["detect_chirality_violations"] = _time(
            lambda: detect_chirality_violations(src), repeats
        )
        out_c = str(work / "corrected.pdb")
        results["correct_af3_output"] = _time(
            lambda: correct_af3_output(src, out_c), repeats
        )
        out_m = str(work / "mirrored.pdb")
        results["mirror_pdb"] = _time(lambda: mirror_pdb(src, out_m), repeats)

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    results["peak_rss_mb"] = round(
        peak / 1024.0 if sys.platform != "darwin" else peak / (1024.0 * 1024.0), 1
    )
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────────────────────

def _run_in(checkout: Path, script: Path, repeats: int) -> Dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(checkout) + os.pathsep + env.get("PYTHONPATH", "")
    out = subprocess.run(
        [sys.executable, str(script), "--child", "--repo-root", str(checkout),
         "--repeats", str(repeats)],
        cwd=str(checkout), env=env, capture_output=True, text=True, check=True,
    )
    return json.loads(out.stdout)


def _speedup(before: Optional[Dict], after: Optional[Dict]) -> Optional[float]:
    if not before or not after:
        return None
    b = before.get("median_ms")
    a = after.get("median_ms")
    if not b or not a:
        return None
    return round(b / a, 2)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--repo-root", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--baseline", default="master",
                        help="git revision to compare against (default: master)")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--out", default=None,
                        help="JSON output path (default results/performance_comparison.json)")
    args = parser.parse_args(argv)

    if args.child:
        root = Path(args.repo_root) if args.repo_root else Path.cwd()
        json.dump(run_child(root, args.repeats), sys.stdout)
        return 0

    here = Path(__file__).resolve()
    repo_root = here.parent.parent
    out_path = Path(args.out) if args.out else (
        repo_root / "results" / "performance_comparison.json"
    )

    after = _run_in(repo_root, here, args.repeats)

    with tempfile.TemporaryDirectory(prefix="cf_baseline_") as tmp:
        worktree = Path(tmp) / "baseline"
        subprocess.run(
            ["git", "worktree", "add", "--detach", str(worktree), args.baseline],
            cwd=str(repo_root), check=True, capture_output=True, text=True,
        )
        try:
            shutil.copy2(here, worktree / "benchmarks" / here.name)
            before = _run_in(worktree, worktree / "benchmarks" / here.name,
                             args.repeats)
            baseline_sha = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=str(worktree), text=True,
            ).strip()
        finally:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(worktree)],
                cwd=str(repo_root), check=False, capture_output=True,
            )

    current_sha = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], cwd=str(repo_root), text=True
    ).strip()

    comparison = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "baseline_ref": args.baseline,
        "baseline_commit": baseline_sha,
        "current_commit": current_sha,
        "repeats": args.repeats,
        "audit": {},
        "audit_scaling": [],
        "interface_scaling": [],
        "pipelines": {},
        "peak_rss_mb": {"before": before.get("peak_rss_mb"),
                        "after": after.get("peak_rss_mb")},
    }

    for rel in AUDIT_TARGETS:
        b = before["audit"].get(rel)
        a = after["audit"].get(rel)
        if not b or not a:
            continue
        comparison["audit"][rel] = {
            "n_atoms": a["n_atoms"],
            "before_ms": b["timing"]["median_ms"],
            "after_ms": a["timing"]["median_ms"],
            "speedup": _speedup(b["timing"], a["timing"]),
        }

    for b, a in zip(before["audit_scaling"], after["audit_scaling"]):
        comparison["audit_scaling"].append({
            "copies": a["copies"],
            "n_atoms": a["n_atoms"],
            "before_ms": b["timing"]["median_ms"],
            "after_ms": a["timing"]["median_ms"],
            "speedup": _speedup(b["timing"], a["timing"]),
        })

    for b, a in zip(before["interface_scaling"], after["interface_scaling"]):
        comparison["interface_scaling"].append({
            "copies": a["copies"],
            "before_ms": b["timing"]["median_ms"],
            "after_ms": a["timing"]["median_ms"],
            "speedup": _speedup(b["timing"], a["timing"]),
        })

    for key in ("detect_chirality_violations", "correct_af3_output", "mirror_pdb"):
        comparison["pipelines"][key] = {
            "before_ms": before[key]["median_ms"],
            "after_ms": after[key]["median_ms"],
            "speedup": _speedup(before[key], after[key]),
        }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(comparison, indent=2) + "\n")

    print("=" * 76)
    print(f"ChiralFold performance: {comparison['baseline_commit']} "
          f"→ {comparison['current_commit']}")
    print("=" * 76)
    print(f"\n{'structure':<52}{'before':>8}{'after':>8}{'x':>6}")
    for rel, row in comparison["audit"].items():
        print(f"{rel:<52}{row['before_ms']:>8.1f}{row['after_ms']:>8.1f}"
              f"{row['speedup']:>6.1f}")
    print(f"\n{'audit_pdb, tiled 3IWY':<52}{'before':>8}{'after':>8}{'x':>6}")
    for row in comparison["audit_scaling"]:
        print(f"{'  ' + str(row['n_atoms']) + ' atoms':<52}"
              f"{row['before_ms']:>8.1f}{row['after_ms']:>8.1f}{row['speedup']:>6.1f}")
    print(f"\n{'score_interface, tiled 3IWY':<52}{'before':>8}{'after':>8}{'x':>6}")
    for row in comparison["interface_scaling"]:
        print(f"{'  ' + str(row['copies']) + ' copies':<52}"
              f"{row['before_ms']:>8.1f}{row['after_ms']:>8.1f}{row['speedup']:>6.1f}")
    print(f"\n{'pipelines':<52}{'before':>8}{'after':>8}{'x':>6}")
    for key, row in comparison["pipelines"].items():
        print(f"{'  ' + key:<52}{row['before_ms']:>8.1f}{row['after_ms']:>8.1f}"
              f"{row['speedup']:>6.1f}")
    print(f"\nJSON written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
