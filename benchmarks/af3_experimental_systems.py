#!/usr/bin/env python3
"""
Audit experimental Childs et al. D-peptide systems from RCSB.

This script does not redistribute AlphaFold 3 predictions and does not claim to
measure AF3 on these structures. It audits deposited experimental PDB
coordinates for the systems discussed by Childs, Zhou & Donald (2025), while the
published AF3 D-peptide chirality failure rates are cited from that paper.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chiralfold.af3_correct import detect_chirality_violations  # noqa: E402
from chiralfold.auditor import audit_pdb  # noqa: E402


RESULTS_DIR = ROOT / "results"
CACHE_DIR = RESULTS_DIR / ".cache" / "childs_experimental_pdbs"
OUT_JSON = RESULTS_DIR / "af3_experimental_systems.json"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

SYSTEMS = [
    {
        "pdb_id": "3IWY",
        "system": "DP12:MDM2",
        "childs_sequence": "DWWPLAFEALLR",
        "published_af3_violation_rate": 0.50,
        "note": "Experimental dPMI-gamma/MDM2 structure used as the DP12 system.",
    },
    {
        "pdb_id": "5N8T",
        "system": "DP9:Streptavidin",
        "childs_sequence": "LWQHEATWK",
        "published_af3_violation_rate": 0.51,
        "note": "Experimental DP9 streptavidin-binding D-peptide system.",
    },
    {
        "pdb_id": "7YH8",
        "system": "DP19:L-19437",
        "childs_sequence": "DEHELLETAARWFYEIAKR",
        "published_af3_violation_rate": 0.52,
        "note": "Experimental DP19 system if available from RCSB.",
    },
]


def _download_pdb(pdb_id: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / f"{pdb_id}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path

    url = RCSB_URL.format(pdb_id=pdb_id)
    with urllib.request.urlopen(url, timeout=30) as response:
        content = response.read()
    if not content.startswith((b"HEADER", b"TITLE", b"REMARK", b"ATOM", b"HETATM")):
        raise ValueError(f"RCSB response for {pdb_id} did not look like a PDB file")
    path.write_bytes(content)
    return path


def _compact_audit(report: Dict[str, object]) -> Dict[str, object]:
    """Keep result JSON compact and stable by removing atom-level internals."""
    chirality = report.get("chirality", {})
    ramachandran = report.get("ramachandran", {})
    planarity = report.get("planarity", {})
    clashes = report.get("clashes", {})
    bond_geometry = report.get("bond_geometry", {})
    return {
        "n_residues": report.get("n_residues"),
        "overall_score": report.get("overall_score"),
        "chirality": {
            "n_correct": chirality.get("n_correct"),
            "n_wrong": chirality.get("n_wrong"),
            "n_glycine": chirality.get("n_glycine"),
            "pct_correct": chirality.get("pct_correct"),
            "violations": chirality.get("violations", []),
        },
        "ramachandran": {
            "n_evaluated": ramachandran.get("n_evaluated"),
            "pct_favored": ramachandran.get("pct_favored"),
            "pct_allowed": ramachandran.get("pct_allowed"),
            "pct_outlier": ramachandran.get("pct_outlier"),
        },
        "planarity": {
            "n_peptide_bonds": planarity.get("n_peptide_bonds"),
            "pct_within_6deg": planarity.get("pct_within_6deg"),
        },
        "clashes": {
            "n_clashes": clashes.get("n_clashes"),
            "clash_score": clashes.get("clash_score"),
        },
        "bond_geometry": {
            "n_bonds_checked": bond_geometry.get("n_bonds_checked"),
            "n_angles_checked": bond_geometry.get("n_angles_checked"),
            "n_outliers": len(bond_geometry.get("outliers", [])),
        },
    }


def audit_systems() -> Dict[str, object]:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows: List[Dict[str, object]] = []
    t0 = time.perf_counter()

    for system in SYSTEMS:
        pdb_id = system["pdb_id"]
        row = dict(system)
        row["rcsb_url"] = RCSB_URL.format(pdb_id=pdb_id)
        try:
            pdb_path = _download_pdb(pdb_id)
            audit = audit_pdb(str(pdb_path))
            af3_detector = detect_chirality_violations(str(pdb_path))
            row.update({
                "status": "ok",
                "cached_pdb": str(pdb_path.relative_to(ROOT)),
                "audit": _compact_audit(audit),
                "af3_detector": af3_detector,
            })
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError, OSError) as exc:
            row.update({
                "status": "unavailable",
                "error": f"{type(exc).__name__}: {exc}",
            })
        rows.append(row)

    output = {
        "benchmark": "af3_experimental_systems",
        "scope": "Experimental RCSB structures for Childs et al. systems; no AF3 outputs redistributed.",
        "af3_reference_note": (
            "Published AF3 failure rates are from Childs, Zhou & Donald 2025 "
            "bioRxiv 2025.03.14.643307, not from local AF3 predictions."
        ),
        "elapsed_s": time.perf_counter() - t0,
        "systems": rows,
    }
    OUT_JSON.write_text(json.dumps(output, indent=2, default=str))
    return output


def main() -> None:
    output = audit_systems()
    ok = sum(1 for row in output["systems"] if row["status"] == "ok")
    unavailable = sum(1 for row in output["systems"] if row["status"] != "ok")
    print("Experimental Childs-system audit complete")
    print(f"  Audited: {ok}")
    print(f"  Unavailable/errors: {unavailable}")
    for row in output["systems"]:
        if row["status"] == "ok":
            detector = row["af3_detector"]
            audit = row["audit"]
            print(
                f"  {row['pdb_id']} {row['system']}: "
                f"detector {detector['n_violations']}/{detector['n_checked']} violations; "
                f"audit chirality {audit['chirality']['pct_correct']}%"
            )
        else:
            print(f"  {row['pdb_id']} {row['system']}: {row['status']} ({row['error']})")
    print(f"  Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
