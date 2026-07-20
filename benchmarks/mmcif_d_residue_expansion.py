#!/usr/bin/env python3
"""
mmCIF D-residue survey (Colab / CI friendly).

Modes
-----
  known-errors  Re-verify the 16 known legacy-PDB error structures in native mmCIF.
  universe      Discover and scan the full live mmCIF-only universe of entries that
                contain any of the 18 standard D-amino acid CCD codes (RCSB search
                + legacy-.pdb HEAD filter).
  both          Run known-errors then universe (default for publication notebooks).

Usage:
    pip install gemmi numpy
    python benchmarks/mmcif_d_residue_expansion.py
    python benchmarks/mmcif_d_residue_expansion.py --mode universe
    python benchmarks/mmcif_d_residue_expansion.py --mode universe --limit 20
    python benchmarks/mmcif_d_residue_expansion.py --mode known-errors

Output:
    results/mmcif_d_residue_expansion.csv
    results/mmcif_d_residue_expansion_summary.json
    results/mmcif_only_universe_ids.json   (universe / both)

Notes:
    The historically cited "~245 mmCIF-only" figure counted failed legacy-.pdb
    downloads against a broad DPN/DAS/DAL CCD enumeration. Live RCSB enumeration
    of entries that actually contain the 18 D-AA CCD codes yields a smaller
    mmCIF-only set (tens of structures); this script surveys that full live set.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
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

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_HEAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

D_AA = {
    "DAL", "DAR", "DSG", "DAS", "DCY", "DGL", "DGN", "DHI", "DIL",
    "DLE", "DLY", "MED", "DPN", "DPR", "DSN", "DTH", "DTR", "DTY", "DVA",
}

# Default max CIF download size (~25 MB). Larger cryo-EM assemblies are skipped
# but recorded so Colab stays responsive; re-run with --max-cif-mb 0 to disable.
DEFAULT_MAX_CIF_MB = 25.0


def signed_volume(n, ca, c, cb):
    v1 = n - ca
    v2 = c - ca
    v3 = cb - ca
    return float(np.dot(v1, np.cross(v2, v3)))


def pos(atom):
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z])


def _http_json(payload: dict, timeout: int = 120) -> dict:
    req = urllib.request.Request(
        SEARCH_URL,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def search_entries_for_ccd(code: str) -> set[str]:
    """Return experimental PDB entry IDs containing CCD ``code`` as polymer or ligand."""
    ids: set[str] = set()
    for attr in (
        "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
        "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
    ):
        payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": attr,
                    "operator": "exact_match",
                    "value": code,
                },
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 10000},
                "results_content_type": ["experimental"],
            },
        }
        try:
            data = _http_json(payload)
            ids.update(x["identifier"].upper() for x in data.get("result_set", []))
        except Exception as exc:
            print(f"  [search warn] {code} via {attr.split('.')[-1]}: {exc}")
        time.sleep(0.05)
    return ids


def legacy_pdb_available(pdb_id: str, timeout: int = 20) -> bool:
    url = PDB_HEAD_URL.format(pdb_id=pdb_id.upper())
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return 200 <= getattr(resp, "status", 200) < 300
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return False
        return False
    except Exception:
        return False


def discover_mmcif_only_universe(
    codes: set[str] | None = None,
    workers: int = 16,
) -> dict:
    """Enumerate live mmCIF-only entries containing D-AA CCD codes."""
    codes = codes or (D_AA - {"MED"})  # MED optional; keep parity with coverage CSV
    per_code: dict[str, list[str]] = {}
    all_ids: set[str] = set()
    print(f"Discovering RCSB entries for {len(codes)} CCD codes...")
    for code in sorted(codes):
        found = search_entries_for_ccd(code)
        per_code[code] = sorted(found)
        all_ids |= found
        print(f"  {code}: {len(found)}")
    print(f"Unique entries: {len(all_ids)}")
    print("Filtering legacy .pdb availability (HEAD)...")
    mmcif_only: list[str] = []
    legacy: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(legacy_pdb_available, pid): pid for pid in sorted(all_ids)}
        for i, fut in enumerate(as_completed(futs), 1):
            pid = futs[fut]
            if fut.result():
                legacy.append(pid)
            else:
                mmcif_only.append(pid)
            if i % 200 == 0 or i == len(futs):
                print(f"  checked {i}/{len(futs)}  mmcif_only={len(mmcif_only)}")
    mmcif_only = sorted(mmcif_only)
    legacy = sorted(legacy)
    report = {
        "n_ccd_codes": len(codes),
        "n_universe_entries": len(all_ids),
        "n_legacy_pdb": len(legacy),
        "n_mmcif_only": len(mmcif_only),
        "mmcif_only": mmcif_only,
        "per_code_counts": {k: len(v) for k, v in per_code.items()},
        "note": (
            "Live RCSB polymer+ligand CCD search, then HEAD filter for missing "
            "legacy .pdb. Historical '~245' counted failed downloads against a "
            "broader DPN/DAS/DAL enumeration; this is the checkable live set."
        ),
    }
    out = RESULTS / "mmcif_only_universe_ids.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out} ({len(mmcif_only)} mmCIF-only IDs)")
    return report


def known_error_pdb_ids() -> list[str]:
    """PDB IDs of the 16 known D-label/L-coordinate mismatches."""
    ids: set[str] = set()
    legacy = RESULTS / "d_residue_verification_summary.json"
    if legacy.exists():
        data = json.loads(legacy.read_text(encoding="utf-8"))
        ids.update(str(x).upper() for x in data.get("errors_by_structure", {}))
    ver = RESULTS / "d_residue_verification.csv"
    if ver.exists():
        with ver.open(encoding="utf-8") as fp:
            for row in csv.DictReader(fp):
                try:
                    if float(row.get("signed_volume") or 0) > 0:
                        ids.add(row["pdb_id"].upper())
                except ValueError:
                    pass
    # Fallback frozen set if summaries missing
    if not ids:
        ids = {
            "1ABI", "1BG0", "1D7T", "1HHZ", "1KO0", "1MCB", "1OF6", "1P52",
            "1UHG", "1XT7", "2AOU", "2ATS", "2H9E", "2RMI", "2W76", "3RIT",
        }
    return sorted(ids)


def download_mmcif(pdb_id: str, max_bytes: int | None) -> Path | None:
    path = CACHE / f"{pdb_id.upper()}.cif"
    if path.exists():
        if max_bytes and path.stat().st_size > max_bytes:
            return None
        return path
    url = CIF_URL.format(pdb_id=pdb_id.upper())
    if max_bytes:
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=30) as resp:
                cl = resp.headers.get("Content-Length")
                if cl and int(cl) > max_bytes:
                    return None
        except Exception:
            pass
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    if max_bytes and len(data) > max_bytes:
        return None
    path.write_bytes(data)
    return path


def scan_mmcif(pdb_id: str, max_bytes: int | None) -> tuple[list[dict], str | None]:
    """Return (rows, skip_reason)."""
    path = download_mmcif(pdb_id, max_bytes=max_bytes)
    if path is None:
        return [], "too_large_or_unavailable"
    try:
        st = gemmi.read_structure(str(path))
    except Exception as exc:
        return [], f"parse_error:{exc}"
    rows: list[dict] = []
    if len(st) == 0:
        return [], "empty_structure"
    model = st[0]
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
    return rows, None


def scan_ids(
    pdb_ids: list[str],
    *,
    cohort: str,
    max_cif_mb: float,
) -> tuple[list[dict], dict]:
    max_bytes = None if max_cif_mb <= 0 else int(max_cif_mb * 1024 * 1024)
    print(f"Scanning {len(pdb_ids)} structures via mmCIF [{cohort}]...")
    all_rows: list[dict] = []
    skipped: dict[str, str] = {}
    failed: dict[str, str] = {}
    for i, pdb_id in enumerate(pdb_ids, 1):
        try:
            rows, reason = scan_mmcif(pdb_id, max_bytes=max_bytes)
            if reason:
                skipped[pdb_id] = reason
                print(f"  [{i}/{len(pdb_ids)}] {pdb_id}: SKIP ({reason})")
            else:
                all_rows.extend(rows)
                n_err = sum(1 for r in rows if r["is_error"])
                print(
                    f"  [{i}/{len(pdb_ids)}] {pdb_id}: "
                    f"{len(rows)} D-residues ({n_err} errors)"
                )
        except Exception as exc:
            failed[pdb_id] = str(exc)
            print(f"  [{i}/{len(pdb_ids)}] {pdb_id}: FAILED ({exc})")
    meta = {
        "cohort": cohort,
        "n_requested": len(pdb_ids),
        "n_scanned": len(pdb_ids) - len(skipped) - len(failed),
        "n_skipped": len(skipped),
        "n_failed": len(failed),
        "skipped": skipped,
        "failed": failed,
    }
    return all_rows, meta


def write_outputs(all_rows: list[dict], summary: dict) -> None:
    out_csv = RESULTS / "mmcif_d_residue_expansion.csv"
    if all_rows:
        with out_csv.open("w", newline="", encoding="utf-8") as fp:
            w = csv.DictWriter(fp, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    else:
        out_csv.write_text("", encoding="utf-8")
    out_json = RESULTS / "mmcif_d_residue_expansion_summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {out_csv} ({len(all_rows)} rows)")
    print(f"Wrote {out_json}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--mode",
        choices=("known-errors", "universe", "both"),
        default="both",
        help="Survey mode (default: both)",
    )
    ap.add_argument("--limit", type=int, default=None, help="Max structures per cohort (smoke test)")
    ap.add_argument(
        "--max-cif-mb",
        type=float,
        default=DEFAULT_MAX_CIF_MB,
        help="Skip CIF downloads larger than this many MB (0 = no limit)",
    )
    ap.add_argument("--workers", type=int, default=16, help="Parallel HEAD workers for discovery")
    args = ap.parse_args()

    all_rows: list[dict] = []
    cohorts: dict[str, dict] = {}
    known_ids: list[str] = []
    universe_ids: list[str] = []

    if args.mode in ("known-errors", "both"):
        known_ids = known_error_pdb_ids()
        if args.limit is not None:
            known_ids = known_ids[: args.limit]
        rows, meta = scan_ids(known_ids, cohort="known-errors", max_cif_mb=args.max_cif_mb)
        for r in rows:
            r = dict(r)
            r["cohort"] = "known-errors"
            all_rows.append(r)
        # attach cohort without mutating scan output schema inconsistently
        cohorts["known-errors"] = meta

    if args.mode in ("universe", "both"):
        # Reuse cached discovery when present and fresh enough for smoke tests
        disc = discover_mmcif_only_universe(workers=args.workers)
        universe_ids = list(disc["mmcif_only"])
        if args.limit is not None:
            universe_ids = universe_ids[: args.limit]
        rows, meta = scan_ids(universe_ids, cohort="universe", max_cif_mb=args.max_cif_mb)
        for r in rows:
            r = dict(r)
            r["cohort"] = "universe"
            all_rows.append(r)
        cohorts["universe"] = meta
        cohorts["universe"]["discovery"] = {
            "n_universe_entries": disc["n_universe_entries"],
            "n_legacy_pdb": disc["n_legacy_pdb"],
            "n_mmcif_only": disc["n_mmcif_only"],
        }

    # Normalize CSV columns (cohort optional for back-compat readers)
    if all_rows and "cohort" not in all_rows[0]:
        pass

    errors = [r for r in all_rows if r["is_error"]]
    error_pdbs = sorted({r["pdb_id"] for r in errors})
    known_errors = [r for r in all_rows if r.get("cohort") == "known-errors" and r["is_error"]]
    universe_errors = [r for r in all_rows if r.get("cohort") == "universe" and r["is_error"]]

    # Deduplicate residue rows for combined CSV display counts
    scanned_structures = sorted({r["pdb_id"] for r in all_rows})
    if args.mode == "known-errors":
        n_structures = len(known_ids)
    elif args.mode == "universe":
        n_structures = len(universe_ids)
    else:
        n_structures = len(set(known_ids) | set(universe_ids))

    summary = {
        "mode": args.mode,
        "n_structures": n_structures,
        "n_structures_with_d_residues": len(scanned_structures),
        "n_d_residues": len(all_rows),
        "n_errors": len(errors),
        "error_pdbs": error_pdbs,
        "known_error_cohort": {
            "n_structures": len(known_ids),
            "n_errors": len(known_errors),
            "error_pdbs": sorted({r["pdb_id"] for r in known_errors}),
        },
        "universe_cohort": {
            "n_structures_requested": len(universe_ids),
            "n_errors": len(universe_errors),
            "error_pdbs": sorted({r["pdb_id"] for r in universe_errors}),
        },
        "cohorts": cohorts,
        "max_cif_mb": args.max_cif_mb,
        "matches_legacy_survey_errors": None,
        "legacy_error_count": 29,
        "note": (
            "Native mmCIF survey. known-errors re-verifies the 16 legacy mismatches; "
            "universe scans all live mmCIF-only entries containing the 18 D-AA CCD codes."
        ),
    }

    # Compare known-error cohort to frozen legacy survey when available
    if known_ids and args.mode in ("known-errors", "both"):
        legacy_path = RESULTS / "d_residue_verification_summary.json"
        if legacy_path.exists():
            legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
            legacy_pdbs = set(legacy.get("errors_by_structure", {}))
            got = set(summary["known_error_cohort"]["error_pdbs"])
            summary["matches_legacy_survey_errors"] = (
                got == legacy_pdbs
                and summary["known_error_cohort"]["n_errors"] == legacy.get("l_error", 29)
            )

    write_outputs(all_rows, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
