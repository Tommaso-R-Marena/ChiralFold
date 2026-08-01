#!/usr/bin/env python3
"""Quantify how the v3.6.0 alternate-location policy changes audit results.

Up to v3.5.1 ChiralFold selected alternate conformations **per atom**, accepting
only altloc labels ``{' ', 'A', '1'}``. That had two consequences:

1. A residue modelled only as altloc ``B``/``C`` lost *every* atom.
2. Where alternates existed, the ``A`` conformer was kept even when it was the
   minority occupancy, and — because the filter ran per atom — a residue could
   end up with its backbone from one alternate and its side chain from another.

v3.6.0 selects **per residue**: blank-altloc atoms are always kept, and among
lettered alternates the label with the highest summed occupancy wins.

This script audits every structure it is given under both policies and reports
the differences, so the scope of the change is measured rather than asserted.
Residues without alternate conformations are provably unaffected, which bounds
the effect on any previously published aggregate.

Usage::

    python benchmarks/altloc_policy_sensitivity.py
    python benchmarks/altloc_policy_sensitivity.py path/to/*.pdb --out report.json

No network access required.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chiralfold import _pdbio  # noqa: E402
from chiralfold import auditor as auditor_mod  # noqa: E402

#: Report fields compared between the two policies.
_COMPARED_FIELDS = (
    ("n_atoms", lambda r: r["n_atoms"]),
    ("n_residues", lambda r: r["n_residues"]),
    ("chirality_pct_correct", lambda r: r["chirality"]["pct_correct"]),
    ("chirality_n_wrong", lambda r: r["chirality"]["n_wrong"]),
    ("rama_pct_outlier", lambda r: r["ramachandran"]["pct_outlier"]),
    ("rama_pct_favored", lambda r: r["ramachandran"]["pct_favored"]),
    ("rama_n_evaluated", lambda r: r["ramachandran"]["n_evaluated"]),
    ("planarity_pct_within_6deg", lambda r: r["planarity"]["pct_within_6deg"]),
    ("clash_score", lambda r: r["clashes"]["clash_score"]),
    ("bl_rmsd", lambda r: r["bond_geometry"]["bl_rmsd"]),
    ("overall_score", lambda r: r["overall_score"]),
)


def _audit_with_policy(path: str, policy: str) -> dict:
    """Audit *path* forcing a specific alternate-location policy."""
    original = auditor_mod.read_atom_records

    def _patched(p, **kwargs):
        kwargs.setdefault("altloc_policy", policy)
        return original(p, **kwargs)

    auditor_mod.read_atom_records = _patched
    try:
        return auditor_mod.audit_pdb(path)
    finally:
        auditor_mod.read_atom_records = original


def _altloc_census(path: str) -> Dict[str, int]:
    """Count alternate-location usage in *path*."""
    with open(path) as fh:
        records = list(_pdbio.iter_atom_records(fh))
    lettered = [r for r in records if r.altloc not in _pdbio.BLANK_ALTLOCS]
    residues_with_alts = {r.residue_key for r in lettered}
    winners = _pdbio.winning_altlocs(records)
    # Residues where the majority conformer is not the legacy 'A'/'1' choice.
    changed = {
        rkey for rkey, label in winners.items()
        if label not in _pdbio.LEGACY_ALTLOC_ACCEPT
    }
    # Residues the legacy filter would have deleted entirely.
    legacy_kept = {
        r.residue_key for r in records
        if r.altloc in _pdbio.LEGACY_ALTLOC_ACCEPT
    }
    dropped = {r.residue_key for r in records} - legacy_kept
    return {
        "n_atom_records": len(records),
        "n_lettered_altloc_atoms": len(lettered),
        "n_residues": len({r.residue_key for r in records}),
        "n_residues_with_alternates": len(residues_with_alts),
        "n_residues_majority_not_A": len(changed),
        "n_residues_lost_under_legacy": len(dropped),
    }


def analyse(paths: List[Path]) -> dict:
    rows: List[dict] = []
    n_identical = 0
    for path in paths:
        p = str(path)
        try:
            legacy = _audit_with_policy(p, "legacy")
            residue = _audit_with_policy(p, "residue")
        except Exception as exc:  # pragma: no cover - diagnostic path
            rows.append({"structure": path.name, "error": str(exc)})
            continue

        deltas = {}
        for name, getter in _COMPARED_FIELDS:
            before, after = getter(legacy), getter(residue)
            if before != after:
                deltas[name] = {"legacy": before, "residue": after}

        if not deltas:
            n_identical += 1
        rows.append({
            "structure": path.name,
            "census": _altloc_census(p),
            "identical": not deltas,
            "deltas": deltas,
        })

    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_structures": len(rows),
        "n_identical_under_both_policies": n_identical,
        "policies": list(_pdbio.ALTLOC_POLICIES),
        "structures": rows,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="PDB files (default: repo corpus)")
    parser.add_argument(
        "--out",
        default=str(REPO_ROOT / "results" / "altloc_policy_sensitivity.json"),
    )
    args = parser.parse_args(argv)

    if args.paths:
        paths = [Path(p) for p in args.paths]
    else:
        paths = sorted(
            list((REPO_ROOT / "results").glob("*.pdb"))
            + list((REPO_ROOT / "chiralfold" / "data" / "examples").glob("*.pdb"))
            + list((REPO_ROOT / "tests" / "fixtures").rglob("*.pdb"))
        )

    report = analyse(paths)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n")

    print("=" * 72)
    print("Alternate-location policy sensitivity (legacy per-atom → per-residue)")
    print("=" * 72)
    print(f"{'structure':<34}{'alt res':>8}{'maj≠A':>7}{'identical':>11}")
    for row in report["structures"]:
        if "error" in row:
            print(f"{row['structure']:<34}  ERROR: {row['error']}")
            continue
        c = row["census"]
        print(
            f"{row['structure']:<34}{c['n_residues_with_alternates']:>8}"
            f"{c['n_residues_majority_not_A']:>7}"
            f"{'yes' if row['identical'] else 'NO':>11}"
        )
    for row in report["structures"]:
        if row.get("deltas"):
            print(f"\nchanged fields for {row['structure']}:")
            for name, d in row["deltas"].items():
                print(f"  {name:<28} {d['legacy']} → {d['residue']}")
    print(
        f"\n{report['n_identical_under_both_policies']}/{report['n_structures']} "
        f"structures identical under both policies"
    )
    print(f"JSON written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
