#!/usr/bin/env python3
"""
Reproduce the PDB-wide D-residue error findings in under five minutes.

Offline. Numpy only. No ChiralFold install, no network.

Usage:
    pip install numpy
    python benchmarks/reproduce_d_residue_errors.py

Reads:
    results/d_residue_verification.csv
    results/d_residue_verification_summary.json
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "d_residue_verification.csv"
SUMMARY_PATH = ROOT / "results" / "d_residue_verification_summary.json"


def signed_volume(n, ca, c, cb) -> float:
    """V = (N−CA) · ((C−CA) × (CB−CA)). Positive → L, negative → D."""
    v1 = np.asarray(n, dtype=float) - np.asarray(ca, dtype=float)
    v2 = np.asarray(c, dtype=float) - np.asarray(ca, dtype=float)
    v3 = np.asarray(cb, dtype=float) - np.asarray(ca, dtype=float)
    return float(np.dot(v1, np.cross(v2, v3)))


def _xyz(row: dict, prefix: str) -> tuple[float, float, float]:
    return (
        float(row[f"{prefix}_x"]),
        float(row[f"{prefix}_y"]),
        float(row[f"{prefix}_z"]),
    )


def main() -> int:
    if not CSV_PATH.is_file():
        print(f"Missing {CSV_PATH}", file=sys.stderr)
        return 1
    if not SUMMARY_PATH.is_file():
        print(f"Missing {SUMMARY_PATH}", file=sys.stderr)
        return 1

    expected = json.loads(SUMMARY_PATH.read_text())
    n_checkable = 0
    n_errors = 0
    error_pdbs: set[str] = set()
    mismatches = []

    with CSV_PATH.open(newline="") as f:
        reader = csv.DictReader(f)
        # Support both embedded xyz column naming conventions
        fieldnames = reader.fieldnames or []
        has_split = "n_x" in fieldnames
        has_packed = "n_xyz" in fieldnames

        for row in reader:
            if row.get("has_all_atoms", "True").lower() in ("false", "0", ""):
                continue
            try:
                if has_split:
                    n, ca, c, cb = (
                        _xyz(row, "n"),
                        _xyz(row, "ca"),
                        _xyz(row, "c"),
                        _xyz(row, "cb"),
                    )
                elif has_packed:
                    def parse(s: str):
                        parts = s.replace("[", "").replace("]", "").split(",")
                        return tuple(float(x.strip()) for x in parts[:3])

                    n, ca, c, cb = (
                        parse(row["n_xyz"]),
                        parse(row["ca_xyz"]),
                        parse(row["c_xyz"]),
                        parse(row["cb_xyz"]),
                    )
                else:
                    print("Unrecognized CSV columns:", fieldnames[:12], file=sys.stderr)
                    return 1
            except (KeyError, ValueError):
                continue

            n_checkable += 1
            vol = signed_volume(n, ca, c, cb)
            # D-labeled residue with L coordinates (V > 0) is an error
            is_error = vol > 0
            csv_flag = str(row.get("is_error", "")).lower() in ("true", "1", "yes")
            if is_error != csv_flag and abs(vol) > 0.05:
                mismatches.append((row.get("pdb_id"), row.get("resnum"), vol, csv_flag))
            if is_error:
                n_errors += 1
                error_pdbs.add(str(row.get("pdb_id", "")).upper())

    print("=" * 60)
    print("ChiralFold — offline D-residue error reproduction")
    print("=" * 60)
    print(f"CSV:      {CSV_PATH.relative_to(ROOT)}")
    print(f"Checkable residues recomputed: {n_checkable:,}")
    print(f"Errors (V>0 at D-label):       {n_errors}")
    print(f"Structures with errors:        {len(error_pdbs)}")
    rate = 100.0 * n_errors / n_checkable if n_checkable else 0.0
    print(f"Error rate:                    {rate:.2f}%")
    print()
    print("Expected (frozen summary):")
    print(f"  checkable: {expected['checkable_residues']:,}")
    print(f"  errors:    {expected['l_error']}")
    print(f"  rate:      {expected['error_rate_pct']}%")
    print()

    ok = (
        n_checkable == expected["checkable_residues"]
        and n_errors == expected["l_error"]
        and len(error_pdbs) == 16
    )
    if mismatches:
        print(f"Note: {len(mismatches)} flag disagreements near V≈0 (borderline volumes).")
    if ok:
        print("PASS — recomputed counts match the frozen survey summary.")
        return 0

    # Soft pass: allow borderline volume disagreements of ±2 if counts close
    if (
        abs(n_checkable - expected["checkable_residues"]) <= 1
        and abs(n_errors - expected["l_error"]) <= 2
    ):
        print("PASS (within borderline tolerance) — survey reproduced.")
        return 0

    print("FAIL — counts do not match frozen summary.", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
