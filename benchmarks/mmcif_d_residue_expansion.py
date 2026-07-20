#!/usr/bin/env python3
"""
Expand D-residue verification to mmCIF-only PDB entries (Colab / Rockfish friendly).

Reads PDB IDs flagged as mmCIF-only in results/ccd_code_coverage_summary.csv,
downloads mmCIF from RCSB, and computes signed Cα volume for D-labeled residues.

Usage:
    pip install gemmi numpy
    python benchmarks/mmcif_d_residue_expansion.py
    python benchmarks/mmcif_d_residue_expansion.py --limit 20   # smoke test

Output:
    results/mmcif_d_residue_expansion.csv
    results/mmcif_d_residue_expansion_summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
import urllib.request
from pathlib import Path

import numpy as np

try:
    import gemmi
except ImportError as e:
    raise SystemExit("Install gemmi: pip install gemmi") from e

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
CACHE = RESULTS / "mmcif_cache"
CACHE.mkdir(parents=True, exist_ok=True)

D_AA = {
    "DAL", "DAR", "DSG", "DAS", "DCY", "DGL", "DGN", "DHI", "DIL",
    "DLE", "DLY", "MED", "DPN", "DPR", "DSN", "DTH", "DTR", "DTY", "DVA",
}


def signed_volume(n, ca, c, cb):
    v1 = n - ca
    v2 = c - ca
    v3 = cb - ca
    return float(np.dot(v1, np.cross(v2, v3)))


def pos(atom):
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z])


def download_mmcif(pdb_id: str) -> Path:
    path = CACHE / f"{pdb_id.upper()}.cif"
    if path.exists():
        return path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    with urllib.request.urlopen(url, timeout=60) as r:
        path.write_bytes(r.read())
    return path


def mmcif_only_pdb_ids(limit: int | None) -> list[str]:
    """PDB IDs from error structures + high mmcif-only counts in coverage table."""
    ids = set()
    cov = RESULTS / "ccd_code_coverage_summary.csv"
    if cov.exists():
        with cov.open() as fp:
            for row in csv.DictReader(fp):
                n = int(row.get("mmcif_only_unavailable", 0) or 0)
                if n > 0:
                    for pdb in (row.get("error_structures") or "").split(";"):
                        pdb = pdb.strip()
                        if pdb:
                            ids.add(pdb)
    # Also scan verification CSV for structures we already know
    ver = RESULTS / "d_residue_verification.csv"
    if ver.exists():
        with ver.open() as fp:
            for row in csv.DictReader(fp):
                try:
                    if float(row.get("signed_volume") or 0) > 0:
                        ids.add(row["pdb_id"])
                except ValueError:
                    pass
    out = sorted(ids)
    return out[:limit] if limit else out


def scan_mmcif(pdb_id: str) -> list[dict]:
    path = download_mmcif(pdb_id)
    st = gemmi.read_structure(str(path))
    rows = []
    model = st[0]  # first model only (NMR ensembles otherwise duplicate residues)
    for chain in model:
        for res in chain:
                if res.name not in D_AA:
                    continue
                n = res.find_atom("N", "*")
                ca = res.find_atom("CA", "*")
                c = res.find_atom("C", "*")
                cb = res.find_atom("CB", "*")
                if not (n and ca and c and cb):
                    continue
                vol = signed_volume(pos(n), pos(ca), pos(c), pos(cb))
                rows.append({
                    "pdb_id": pdb_id.upper(),
                    "chain": chain.name,
                    "resnum": res.seqid.num,
                    "resname": res.name,
                    "signed_volume": vol,
                    "chirality": "L" if vol > 0 else "D",
                    "is_error": vol > 0,
                    "source_format": "mmCIF",
                })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="Max structures (smoke test)")
    args = ap.parse_args()

    pdb_ids = mmcif_only_pdb_ids(args.limit)
    print(f"Scanning {len(pdb_ids)} structures via mmCIF...")
    all_rows: list[dict] = []
    for i, pdb_id in enumerate(pdb_ids, 1):
        try:
            rows = scan_mmcif(pdb_id)
            all_rows.extend(rows)
            print(f"  [{i}/{len(pdb_ids)}] {pdb_id}: {len(rows)} D-residues")
        except Exception as exc:
            print(f"  [{i}/{len(pdb_ids)}] {pdb_id}: FAILED ({exc})")

    out_csv = RESULTS / "mmcif_d_residue_expansion.csv"
    if all_rows:
        with out_csv.open("w", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)

    errors = [r for r in all_rows if r["is_error"]]
    summary = {
        "n_structures": len(pdb_ids),
        "n_d_residues": len(all_rows),
        "n_errors": len(errors),
        "error_pdbs": sorted({r["pdb_id"] for r in errors}),
        "note": "mmCIF-native re-verification of known error structures and mmcif-flagged entries",
    }
    out_json = RESULTS / "mmcif_d_residue_expansion_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {out_csv} ({len(all_rows)} rows, {len(errors)} errors)")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
