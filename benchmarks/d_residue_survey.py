"""
D-residue PDB Chirality Survey
================================
Downloads up to 200 PDB structures containing D-amino acid residues,
audits each with ChiralFold, and records chirality violations.

Supports resume: skips PDB IDs already recorded in d_survey_results.json.
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.request
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────
WORKSPACE     = Path("/home/user/workspace/chiralfold")
BENCHMARKS    = WORKSPACE / "benchmarks"
RESULTS_DIR   = WORKSPACE / "results" / "d_survey"
RESULTS_JSON  = WORKSPACE / "results" / "d_survey_results.json"
PDB_IDS_FILE  = BENCHMARKS / "d_residue_pdb_ids.json"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(WORKSPACE))

# ── parameters ─────────────────────────────────────────────────────────────
MAX_STRUCTURES   = 200
MAX_FILE_SIZE_B  = 2 * 1024 * 1024   # 2 MB
TIMEOUT_S        = 30                 # per-structure audit timeout (seconds)
SAVE_EVERY       = 5                  # checkpoint frequency (very frequent for resume)
DOWNLOAD_TIMEOUT = 20                 # urllib timeout (seconds)

# D-amino acid residue names (mirrors auditor.py)
D_RESNAMES = {
    "DAL", "DAR", "DSG", "DAS", "DCY", "DGL", "DGN",
    "DHI", "DIL", "DLE", "DLY", "MED", "DPN", "DPR",
    "DSN", "DTH", "DTR", "DTY", "DVA",
}


# ── helpers ────────────────────────────────────────────────────────────────

def load_pdb_ids() -> list[str]:
    with open(PDB_IDS_FILE) as f:
        data = json.load(f)
    return data["all_ids"]


def load_existing_results() -> dict:
    """Load existing results for resume support."""
    if RESULTS_JSON.exists():
        try:
            with open(RESULTS_JSON) as f:
                data = json.load(f)
            print(f"Resuming from checkpoint: {data.get('total_structures_audited', 0)} already audited")
            # Ensure skipped_ids list exists
            if "skipped_ids" not in data:
                data["skipped_ids"] = []
            return data
        except Exception as e:
            print(f"Could not load checkpoint ({e}), starting fresh")
    return {
        "total_structures_attempted": 0,
        "total_structures_audited": 0,
        "total_with_perfect_chirality": 0,
        "total_with_errors": 0,
        "total_skipped_too_large": 0,
        "total_skipped_download_fail": 0,
        "total_skipped_timeout": 0,
        "total_skipped_audit_error": 0,
        "error_details": [],
        "all_audited": [],
        "skipped_ids": [],  # IDs that failed/timed out/too large — skip on resume
    }


def get_already_audited_ids(results: dict) -> set[str]:
    """Return set of PDB IDs already recorded in results (audited or skipped)."""
    audited = {r["pdb_id"] for r in results.get("all_audited", [])}
    skipped = set(results.get("skipped_ids", []))
    return audited | skipped


def get_file_size_bytes(pdb_id: str) -> int | None:
    """HEAD request to get Content-Length before downloading."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT) as resp:
            cl = resp.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def download_pdb(pdb_id: str, path: Path) -> bool:
    """Download a PDB file; return True on success."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, str(path))
        return True
    except Exception as exc:
        print(f"  [DOWNLOAD FAIL] {pdb_id}: {exc}")
        if path.exists():
            path.unlink()
        return False


def audit_with_timeout(pdb_path: str, timeout: int) -> dict | None:
    """
    Run audit_pdb in a separate process with a hard timeout.
    Returns the report dict or None if timed out / errored.
    """
    import multiprocessing

    def _worker(path, result_queue):
        try:
            sys.path.insert(0, str(WORKSPACE))
            from chiralfold.auditor import audit_pdb
            report = audit_pdb(path)
            result_queue.put(("ok", report))
        except Exception as exc:
            result_queue.put(("error", str(exc)))

    q = multiprocessing.Queue()
    p = multiprocessing.Process(target=_worker, args=(pdb_path, q), daemon=True)
    p.start()
    p.join(timeout)

    if p.is_alive():
        p.terminate()
        p.join(2)
        return None   # timed out

    if not q.empty():
        status, payload = q.get()
        if status == "ok":
            return payload
        else:
            print(f"  [AUDIT ERROR] {payload}")
            return None
    return None


def extract_violation_info(violations: list[dict]) -> list[dict]:
    """Enrich violation records with D-residue flag."""
    enriched = []
    for v in violations:
        enriched.append({
            "chain":    v.get("chain", "?"),
            "resseq":   v.get("resseq", 0),
            "resname":  v.get("resname", "???"),
            "expected": v.get("expected", "?"),
            "observed": v.get("observed", "?"),
            "is_d_residue": v.get("resname", "") in D_RESNAMES,
        })
    return enriched


def save_results(results: dict):
    with open(RESULTS_JSON, "w") as f:
        json.dump(results, f, indent=2)


# ── main survey ────────────────────────────────────────────────────────────

def main():
    print("=== D-residue PDB Chirality Survey ===\n")

    all_ids = load_pdb_ids()
    print(f"Loaded {len(all_ids)} PDB IDs from {PDB_IDS_FILE}\n")

    # ── Resume support ────────────────────────────────────────────────────
    results = load_existing_results()
    already_audited = get_already_audited_ids(results)

    audited_count   = results["total_structures_audited"]
    attempted_count = results["total_structures_attempted"]

    print(f"Already audited: {audited_count}, need {MAX_STRUCTURES - audited_count} more\n")

    new_audited_this_run = 0

    for pdb_id in all_ids:
        if audited_count >= MAX_STRUCTURES:
            break

        # Skip already processed
        if pdb_id in already_audited:
            continue

        attempted_count += 1
        results["total_structures_attempted"] = attempted_count
        pdb_path = RESULTS_DIR / f"{pdb_id}.pdb"

        # ── already downloaded? ──────────────────────────────────────────
        if pdb_path.exists():
            file_size = pdb_path.stat().st_size
            if file_size > MAX_FILE_SIZE_B:
                print(f"  [SKIP large cached] {pdb_id} ({file_size/1024:.0f} KB)")
                results["total_skipped_too_large"] += 1
                results["skipped_ids"].append(pdb_id)
                already_audited.add(pdb_id)
                continue
        else:
            # HEAD check before downloading
            file_size_estimate = get_file_size_bytes(pdb_id)
            if file_size_estimate is not None and file_size_estimate > MAX_FILE_SIZE_B:
                print(f"  [SKIP too large] {pdb_id} ({file_size_estimate/1024:.0f} KB)")
                results["total_skipped_too_large"] += 1
                results["skipped_ids"].append(pdb_id)
                already_audited.add(pdb_id)
                continue

            # Download
            ok = download_pdb(pdb_id, pdb_path)
            if not ok:
                results["total_skipped_download_fail"] += 1
                results["skipped_ids"].append(pdb_id)
                already_audited.add(pdb_id)
                continue

            # Re-check size after download
            if pdb_path.exists():
                file_size = pdb_path.stat().st_size
                if file_size > MAX_FILE_SIZE_B:
                    print(f"  [SKIP large after dl] {pdb_id} ({file_size/1024:.0f} KB)")
                    results["total_skipped_too_large"] += 1
                    results["skipped_ids"].append(pdb_id)
                    already_audited.add(pdb_id)
                    pdb_path.unlink()
                    continue
            else:
                results["total_skipped_download_fail"] += 1
                results["skipped_ids"].append(pdb_id)
                already_audited.add(pdb_id)
                continue

        # ── audit ─────────────────────────────────────────────────────────
        t0 = time.time()
        report = audit_with_timeout(str(pdb_path), TIMEOUT_S)
        elapsed = time.time() - t0

        if report is None:
            # Distinguish timeout vs other failure
            if elapsed >= TIMEOUT_S - 1:
                print(f"  [TIMEOUT] {pdb_id} after {elapsed:.1f}s")
                results["total_skipped_timeout"] += 1
            else:
                print(f"  [AUDIT FAIL] {pdb_id}")
                results["total_skipped_audit_error"] += 1
            results["skipped_ids"].append(pdb_id)
            already_audited.add(pdb_id)
            # Save checkpoint after every skip so we don't re-process them
            _add_summary(results)
            save_results(results)
            continue

        audited_count += 1
        new_audited_this_run += 1
        already_audited.add(pdb_id)
        results["total_structures_audited"] = audited_count

        chirality   = report["chirality"]
        n_wrong     = chirality.get("n_wrong", 0)
        n_correct   = chirality.get("n_correct", 0)
        pct_correct = chirality.get("pct_correct", 100.0)
        violations  = chirality.get("violations", [])

        record = {
            "pdb_id":       pdb_id,
            "n_residues":   report.get("n_residues", 0),
            "n_atoms":      report.get("n_atoms", 0),
            "overall_score": report.get("overall_score", 0),
            "n_correct":    n_correct,
            "n_wrong":      n_wrong,
            "pct_correct":  pct_correct,
            "n_violations": len(violations),
            "has_errors":   n_wrong > 0,
        }
        results["all_audited"].append(record)

        if n_wrong > 0:
            results["total_with_errors"] += 1
            enriched_violations = extract_violation_info(violations)
            n_d_residue_errors = sum(1 for v in enriched_violations if v["is_d_residue"])
            results["error_details"].append({
                "pdb_id":            pdb_id,
                "n_errors":          n_wrong,
                "n_d_residue_errors": n_d_residue_errors,
                "n_l_residue_errors": n_wrong - n_d_residue_errors,
                "pct_correct":       pct_correct,
                "error_residues":    enriched_violations,
            })
        else:
            results["total_with_perfect_chirality"] += 1

        # ── progress ──────────────────────────────────────────────────────
        if audited_count % 10 == 0:
            print(
                f"[{audited_count}/{MAX_STRUCTURES}] audited | "
                f"errors: {results['total_with_errors']} | "
                f"last: {pdb_id} ({n_wrong} wrong, {elapsed:.1f}s)"
            )

        # ── checkpoint ────────────────────────────────────────────────────
        if new_audited_this_run % SAVE_EVERY == 0:
            _add_summary(results)
            save_results(results)
            print(f"  >>> Checkpoint saved at {audited_count} structures <<<")

    # ── final summary ──────────────────────────────────────────────────────
    _add_summary(results)
    save_results(results)

    print("\n=== Survey Complete ===")
    print(f"Attempted      : {results['total_structures_attempted']}")
    print(f"Audited        : {results['total_structures_audited']}")
    print(f"Perfect chiral : {results['total_with_perfect_chirality']}")
    print(f"With errors    : {results['total_with_errors']}")
    print(f"Too large      : {results['total_skipped_too_large']}")
    print(f"Download fails : {results['total_skipped_download_fail']}")
    print(f"Timeouts       : {results['total_skipped_timeout']}")
    print(f"Audit errors   : {results['total_skipped_audit_error']}")
    print(f"\nResults saved to: {RESULTS_JSON}")


def _add_summary(results: dict):
    """Compute and attach summary statistics."""
    audited = results["all_audited"]
    n = len(audited)
    if n == 0:
        results["summary"] = {"n_audited": 0}
        return

    wrongs     = [r["n_wrong"]     for r in audited]
    pcts       = [r["pct_correct"] for r in audited]
    n_residues = [r["n_residues"]  for r in audited]
    scores     = [r["overall_score"] for r in audited]

    # D-residue error breakdown
    d_err_counts = [e["n_d_residue_errors"] for e in results["error_details"]]
    l_err_counts = [e["n_l_residue_errors"] for e in results["error_details"]]

    results["summary"] = {
        "n_audited":               n,
        "n_with_errors":           results["total_with_errors"],
        "error_rate_pct":          round(100 * results["total_with_errors"] / n, 2) if n else 0,
        "mean_pct_correct_chiral": round(sum(pcts) / n, 3) if n else 0,
        "min_pct_correct_chiral":  round(min(pcts), 3) if pcts else 0,
        "max_n_wrong":             max(wrongs) if wrongs else 0,
        "total_chirality_errors":  sum(wrongs),
        "total_d_residue_errors":  sum(d_err_counts),
        "total_l_residue_errors":  sum(l_err_counts),
        "mean_overall_score":      round(sum(scores) / n, 2) if n else 0,
        "mean_residues":           round(sum(n_residues) / n, 1) if n else 0,
        "top_error_structures":    sorted(
            [{"pdb_id": e["pdb_id"], "n_errors": e["n_errors"]} for e in results["error_details"]],
            key=lambda x: x["n_errors"], reverse=True
        )[:10],
    }


if __name__ == "__main__":
    main()
