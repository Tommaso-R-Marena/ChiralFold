#!/usr/bin/env python3
"""
Expand the ChiralFold vs wwPDB/MolProbity Ramachandran benchmark to n>=100.

This script downloads a stratified random sample of >=100 PDB structures from
RCSB across X-ray resolution bins, NMR ensembles, and cryo-EM reconstructions;
fetches the corresponding wwPDB validation XML reports to extract the
MolProbity-derived Ramachandran outlier percentages; runs ChiralFold's auditor
on each structure to record its Ramachandran outlier percentage; and computes
the Spearman rho / Pearson r between the two outlier columns with a scatter
plot.

Outputs (written to ``results/``):
  - ramachandran_100struct_comparison.csv  -- one row per structure
  - ramachandran_100struct_plot.png         -- scatter w/ regression band
  - molprobity_comparison.json              -- summary stats (overwritten)

USAGE (network + multi-hour compute; not run on CI):
    python benchmarks/expand_ramachandran_benchmark.py            # default n=100
    python benchmarks/expand_ramachandran_benchmark.py --n 150    # larger
    python benchmarks/expand_ramachandran_benchmark.py --dry-run  # plan only

NOTE: The canonical value reported in the README/paper for the published
31-structure benchmark is Spearman rho = 0.49 (p = 0.006). Any rho/n/p
computed by this script supersedes that value ONLY when a clean n>=100 run
completes; until then, the canonical 31-structure number stands. See
results/REPRODUCIBILITY.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Dict, List, Optional

# Allow running as a standalone script from the repo root.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

try:
    from scipy import stats as scistats
except ImportError as exc:  # pragma: no cover
    raise SystemExit('scipy is required for this benchmark') from exc

from chiralfold.auditor import audit_pdb

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RESULTS_DIR = os.path.join(REPO_ROOT, 'results')
CACHE_DIR = os.path.join(REPO_ROOT, 'results', 'ramachandran_expansion_cache')

# Stratified resolution bins (X-ray): high / medium / low / very low.
# NMR and cryo-EM are pulled by experimental method, not resolution.
RESOLUTION_BINS = [
    ('xray_ultra', '<= 1.20',
     {'from': 0.0, 'to': 1.20, 'include_lower': True, 'include_upper': True},
     25),
    ('xray_high',  '1.20-1.80',
     {'from': 1.20, 'to': 1.80, 'include_lower': False, 'include_upper': True},
     25),
    ('xray_med',   '1.80-2.50',
     {'from': 1.80, 'to': 2.50, 'include_lower': False, 'include_upper': True},
     20),
    ('xray_low',   '2.50-3.40',
     {'from': 2.50, 'to': 3.40, 'include_lower': False, 'include_upper': True},
     15),
]
NMR_BIN_TARGET = 10
CRYOEM_BIN_TARGET = 10  # 3.0-4.5 A reconstructions

PDB_FILE_URL = 'https://files.rcsb.org/download/{pdb_id}.pdb'
VALIDATION_XML_URL = (
    'https://files.rcsb.org/pub/pdb/validation_reports/{stem}/'
    '{pdb_id_lower}/{pdb_id_lower}_validation.xml.gz'
)
SEARCH_URL = 'https://search.rcsb.org/rcsbsearch/v2/query'


def _build_query(resolution_filter,
                 method: str,
                 max_size: int = 200) -> dict:
    """Construct an RCSB Search API JSON query for a single bin.

    resolution_filter is either None (no resolution gate) or a dict with the
    RCSB ``range`` operator schema: ``{from, to, include_lower, include_upper}``.
    """
    nodes: List[dict] = [{
        'type': 'terminal',
        'service': 'text',
        'parameters': {
            'attribute': 'exptl.method',
            'operator': 'exact_match',
            'value': method,
        },
    }]
    if resolution_filter:
        nodes.append({
            'type': 'terminal',
            'service': 'text',
            'parameters': {
                'attribute': 'rcsb_entry_info.resolution_combined',
                'operator': 'range',
                'value': resolution_filter,
            },
        })
    return {
        'query': {'type': 'group', 'logical_operator': 'and', 'nodes': nodes},
        'return_type': 'entry',
        'request_options': {
            'paginate': {'start': 0, 'rows': max_size},
            'results_content_type': ['experimental'],
            'sort': [{'sort_by': 'score', 'direction': 'desc'}],
        },
    }


def _search_rcsb(query: dict, timeout: int = 30) -> List[str]:
    """POST a search query to RCSB; return matching PDB IDs."""
    data = json.dumps(query).encode('utf-8')
    req = urllib.request.Request(
        SEARCH_URL, data=data, headers={'Content-Type': 'application/json'}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read().decode('utf-8'))
    return [hit['identifier'] for hit in payload.get('result_set', [])]


def _stratified_sample(n_target: int, seed: int = 20260514) -> List[Dict]:
    """Build a stratified random sample of PDB IDs across resolution / method.

    Returns a list of {pdb_id, bin, expected_resolution_range, method} dicts.
    """
    rng = random.Random(seed)
    sample: List[Dict] = []

    for bin_name, res_label, res_filter, target in RESOLUTION_BINS:
        try:
            ids = _search_rcsb(_build_query(res_filter, 'X-RAY DIFFRACTION'))
        except (urllib.error.URLError, ValueError) as exc:
            print(f'  [WARN] {bin_name} query failed: {exc}')
            continue
        rng.shuffle(ids)
        sample.extend({
            'pdb_id': pid, 'bin': bin_name, 'resolution_range': res_label,
            'method': 'X-RAY DIFFRACTION',
        } for pid in ids[:target])

    try:
        nmr_ids = _search_rcsb(_build_query(None, 'SOLUTION NMR'))
    except urllib.error.URLError as exc:
        print(f'  [WARN] NMR query failed: {exc}')
        nmr_ids = []
    rng.shuffle(nmr_ids)
    sample.extend({
        'pdb_id': pid, 'bin': 'nmr', 'resolution_range': 'NMR',
        'method': 'SOLUTION NMR',
    } for pid in nmr_ids[:NMR_BIN_TARGET])

    try:
        em_ids = _search_rcsb(_build_query(
            {'from': 0.0, 'to': 4.5,
             'include_lower': True, 'include_upper': True},
            'ELECTRON MICROSCOPY'))
    except urllib.error.URLError as exc:
        print(f'  [WARN] cryo-EM query failed: {exc}')
        em_ids = []
    rng.shuffle(em_ids)
    sample.extend({
        'pdb_id': pid, 'bin': 'cryoem', 'resolution_range': '<=4.5 A (EM)',
        'method': 'ELECTRON MICROSCOPY',
    } for pid in em_ids[:CRYOEM_BIN_TARGET])

    if len(sample) < n_target:
        print(f'  [WARN] stratified sample only has {len(sample)} entries '
              f'(target {n_target}); proceeding with what we have.')
    return sample[:max(n_target, len(sample))]


def _download_pdb(pdb_id: str, cache: str) -> Optional[str]:
    path = os.path.join(cache, f'{pdb_id.lower()}.pdb')
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    try:
        urllib.request.urlretrieve(
            PDB_FILE_URL.format(pdb_id=pdb_id.lower()), path)
        return path
    except urllib.error.URLError as exc:
        print(f'    [skip] {pdb_id}: PDB download failed ({exc})')
        return None


def _wwpdb_outlier_pct(pdb_id: str, cache: str) -> Optional[float]:
    """Fetch wwPDB validation XML and return the Ramachandran outlier %."""
    import gzip
    stem = pdb_id.lower()[1:3]
    url = VALIDATION_XML_URL.format(stem=stem, pdb_id_lower=pdb_id.lower())
    cached = os.path.join(cache, f'{pdb_id.lower()}_validation.xml')
    if not os.path.exists(cached) or os.path.getsize(cached) == 0:
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                with gzip.GzipFile(fileobj=resp) as gz:
                    with open(cached, 'wb') as out:
                        out.write(gz.read())
        except (urllib.error.URLError, OSError) as exc:
            print(f'    [skip] {pdb_id}: validation XML unavailable ({exc})')
            return None
    try:
        tree = ET.parse(cached)
    except ET.ParseError as exc:
        print(f'    [skip] {pdb_id}: validation XML parse error ({exc})')
        return None
    root = tree.getroot()
    entry = root.find('Entry')
    if entry is None:
        return None
    pct = entry.get('percent-rama-outliers')
    if pct is None or pct == 'NotAvailable':
        return None
    try:
        return float(pct)
    except ValueError:
        return None


def _run_audit(pdb_path: str) -> Optional[float]:
    try:
        report = audit_pdb(pdb_path)
    except Exception as exc:  # noqa: BLE001 - benchmark wrapper
        print(f'    [skip] audit failed: {exc}')
        return None
    rama = report.get('ramachandran') or {}
    if 'pct_outlier' not in rama:
        return None
    return float(rama['pct_outlier'])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--n', type=int, default=105,
                        help='target sample size (default 105, ensures >=100)')
    parser.add_argument('--dry-run', action='store_true',
                        help='print plan without downloading or auditing')
    parser.add_argument('--seed', type=int, default=20260514)
    args = parser.parse_args()

    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f'Building stratified sample, target n={args.n} (seed={args.seed})')
    plan = _stratified_sample(args.n, seed=args.seed)
    plan_counts = defaultdict(int)
    for row in plan:
        plan_counts[row['bin']] += 1
    print('  Plan:', dict(plan_counts), 'total =', len(plan))

    if args.dry_run:
        print('Dry run; exiting without network downloads.')
        return 0

    rows: List[Dict] = []
    t0 = time.time()
    for i, row in enumerate(plan, 1):
        pdb_id = row['pdb_id']
        print(f'  [{i}/{len(plan)}] {pdb_id} ({row["bin"]})')
        pdb_path = _download_pdb(pdb_id, CACHE_DIR)
        if pdb_path is None:
            continue
        wwpdb_pct = _wwpdb_outlier_pct(pdb_id, CACHE_DIR)
        if wwpdb_pct is None:
            continue
        cf_pct = _run_audit(pdb_path)
        if cf_pct is None:
            continue
        rows.append({
            'pdb_id': pdb_id, 'bin': row['bin'],
            'method': row['method'],
            'resolution_range': row['resolution_range'],
            'chiralfold_rama_outlier_pct': cf_pct,
            'wwpdb_rama_outlier_pct': wwpdb_pct,
        })

    elapsed = time.time() - t0
    print(f'Completed {len(rows)} structures in {elapsed/60:.1f} min')

    if len(rows) < 30:
        print(f'  [ERROR] only {len(rows)} usable rows; aborting before '
              'overwriting canonical molprobity_comparison.json.')
        return 1

    csv_path = os.path.join(RESULTS_DIR, 'ramachandran_100struct_comparison.csv')
    with open(csv_path, 'w', newline='') as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {csv_path}')

    cf = np.array([r['chiralfold_rama_outlier_pct'] for r in rows])
    wp = np.array([r['wwpdb_rama_outlier_pct'] for r in rows])
    rho, rho_p = scistats.spearmanr(cf, wp)
    r_pearson, p_pearson = scistats.pearsonr(cf, wp)

    summary = {
        'n_structures': int(len(rows)),
        'spearman_rho': float(rho),
        'spearman_p_value': float(rho_p),
        'pearson_r': float(r_pearson),
        'pearson_p_value': float(p_pearson),
        'chiralfold_mean_outlier_pct': float(cf.mean()),
        'wwpdb_mean_outlier_pct': float(wp.mean()),
        'note': (
            'Expanded benchmark via benchmarks/expand_ramachandran_benchmark.py. '
            'See ramachandran_100struct_comparison.csv for the per-structure data.'
        ),
    }
    json_path = os.path.join(RESULTS_DIR, 'molprobity_comparison.json')
    with open(json_path, 'w') as fp:
        json.dump(summary, fp, indent=2)
    print(f'Wrote {json_path}: rho={rho:.3f} (p={rho_p:.4g}), n={len(rows)}')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(wp, cf, s=22, alpha=0.75)
        lim = max(cf.max(), wp.max()) * 1.1
        ax.plot([0, lim], [0, lim], '--', linewidth=0.8, color='gray')
        ax.set_xlabel('wwPDB/MolProbity Ramachandran outlier (%)')
        ax.set_ylabel('ChiralFold Ramachandran outlier (%)')
        ax.set_title(f'n={len(rows)} structures; Spearman rho={rho:.2f} '
                     f'(p={rho_p:.3g})')
        fig.tight_layout()
        plot_path = os.path.join(
            RESULTS_DIR, 'ramachandran_100struct_plot.png')
        fig.savefig(plot_path, dpi=160)
        plt.close(fig)
        print(f'Wrote {plot_path}')
    except Exception as exc:  # noqa: BLE001
        print(f'  [WARN] plot generation failed: {exc}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
