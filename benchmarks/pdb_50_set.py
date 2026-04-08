"""
pdb_50_set.py
Download 50 diverse PDB structures and fetch wwPDB validation metrics.
Results saved to pdb50_metadata.json.
"""

import json
import os
import sys
import time
import urllib.request
import urllib.error

# ── 50 PDB IDs ────────────────────────────────────────────────────────────────
PDB_IDS = {
    "ultra_high_res_xray": ["1EJG", "2VB1", "1US0", "2OL9", "3NIR"],
    "high_res_xray": ["1CRN", "2LZM", "1L2Y", "3WBM", "4LZT", "2QMT", "1A6M", "4CI0", "3AYX", "5B3L"],
    "medium_res_xray": ["1UBQ", "1SHG", "3IWY", "1YCR", "2RGF", "4HHB", "1IGT", "5HHD", "2AKE", "3BLM"],
    "low_res_xray": ["3J3Q", "1AON", "4V5D", "2HHB", "3PQR", "1DKG", "3OGO", "5T89", "6VXX", "7BV2"],
    "nmr": ["1D3Z", "2KOD", "1G6J", "5PO5", "2N3G", "1TV0", "2JXR", "1J8B", "6I2T", "2RVD"],
    "cryo_em": ["7K3G", "7A4M", "8FZK", "6W41", "5VOS"],
}

ALL_PDB_IDS = []
for category, ids in PDB_IDS.items():
    ALL_PDB_IDS.extend(ids)

assert len(ALL_PDB_IDS) == 50, f"Expected 50 IDs, got {len(ALL_PDB_IDS)}"

# ── Paths ─────────────────────────────────────────────────────────────────────
PDB_DIR = "/home/user/workspace/chiralfold/results/pdb50"
OUTPUT_JSON = "/home/user/workspace/chiralfold/benchmarks/pdb50_metadata.json"

os.makedirs(PDB_DIR, exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────

def fetch_json(url, timeout=30):
    """Fetch a URL and return parsed JSON, or None on failure."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ChiralFold-Benchmark/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        print(f"    [WARN] JSON fetch failed for {url}: {e}")
        return None


def download_pdb(pdb_id, out_dir, timeout=60):
    """Download a PDB file. Returns True if successful."""
    pdb_lower = pdb_id.lower()
    out_path = os.path.join(out_dir, f"{pdb_lower}.pdb")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print(f"  [SKIP] {pdb_id}: already downloaded")
        return True

    # Try RCSB PDB file download
    url = f"https://files.rcsb.org/download/{pdb_lower}.pdb"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ChiralFold-Benchmark/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
        if len(data) == 0:
            raise ValueError("Empty response")
        with open(out_path, "wb") as f:
            f.write(data)
        print(f"  [OK]   {pdb_id}: downloaded {len(data):,} bytes")
        return True
    except Exception as e:
        print(f"  [FAIL] {pdb_id}: download failed — {e}")
        # Remove partial file if any
        if os.path.exists(out_path):
            os.remove(out_path)
        return False


def safe_first(lst):
    """Return first element of a list, or None."""
    if isinstance(lst, list) and len(lst) > 0:
        return lst[0]
    return lst  # might already be a scalar


def fetch_validation_metrics(pdb_id):
    """
    Fetch resolution, method, Ramachandran outliers %, and clashscore
    from the RCSB Data API.
    Returns a dict with keys: resolution, method, rama_outliers_wwpdb, clashscore_wwpdb
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    data = fetch_json(url)

    resolution = None
    method = None
    rama_outliers = None
    clashscore = None

    if data is None:
        return {
            "resolution": resolution,
            "method": method,
            "rama_outliers_wwpdb": rama_outliers,
            "clashscore_wwpdb": clashscore,
        }

    # Resolution
    try:
        res_list = data.get("rcsb_entry_info", {}).get("resolution_combined", None)
        resolution = safe_first(res_list)
    except Exception:
        pass

    # Experimental method
    try:
        exptl = data.get("exptl", [])
        if exptl:
            method = exptl[0].get("method", None)
    except Exception:
        pass

    # Validation summary — geometry metrics live in pdbx_vrpt_summary_geometry
    try:
        geom = data.get("pdbx_vrpt_summary_geometry", None)
        if isinstance(geom, list) and geom:
            geom = geom[0]
        elif isinstance(geom, dict):
            pass
        else:
            geom = None
        if geom:
            rama_outliers = geom.get("percent_ramachandran_outliers", None)
            clashscore = geom.get("clashscore", None)
    except Exception:
        pass

    return {
        "resolution": resolution,
        "method": method,
        "rama_outliers_wwpdb": rama_outliers,
        "clashscore_wwpdb": clashscore,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    results = []
    total = len(ALL_PDB_IDS)

    for i, pdb_id in enumerate(ALL_PDB_IDS, 1):
        print(f"\n[{i:02d}/{total}] Processing {pdb_id} ...")

        # Step 1: Download PDB file
        downloaded = download_pdb(pdb_id, PDB_DIR)

        # Step 2: Fetch validation metrics
        print(f"  Fetching validation metrics ...")
        metrics = fetch_validation_metrics(pdb_id)

        record = {
            "pdb_id": pdb_id.upper(),
            "resolution": metrics["resolution"],
            "method": metrics["method"],
            "rama_outliers_wwpdb": metrics["rama_outliers_wwpdb"],
            "clashscore_wwpdb": metrics["clashscore_wwpdb"],
            "downloaded": downloaded,
        }
        results.append(record)
        print(f"  => resolution={record['resolution']}, method={record['method']}, "
              f"rama={record['rama_outliers_wwpdb']}, clash={record['clashscore_wwpdb']}, "
              f"downloaded={record['downloaded']}")

        # Be polite to the API
        time.sleep(0.3)

    # Save results
    with open(OUTPUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✓ Saved metadata for {len(results)} entries to {OUTPUT_JSON}")

    # Summary
    downloaded_count = sum(1 for r in results if r["downloaded"])
    with_rama = sum(1 for r in results if r["rama_outliers_wwpdb"] is not None)
    with_clash = sum(1 for r in results if r["clashscore_wwpdb"] is not None)
    print(f"  Downloaded PDB files : {downloaded_count}/{total}")
    print(f"  With Ramachandran %  : {with_rama}/{total}")
    print(f"  With clashscore      : {with_clash}/{total}")


if __name__ == "__main__":
    main()
