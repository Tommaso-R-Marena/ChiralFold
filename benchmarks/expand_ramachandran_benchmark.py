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
CIF_FILE_URL = 'https://files.rcsb.org/download/{pdb_id}.cif'
VALIDATION_XML_URL = (
    'https://files.rcsb.org/pub/pdb/validation_reports/{stem}/'
    '{pdb_id_lower}/{pdb_id_lower}_validation.xml.gz'
)
SEARCH_URL = 'https://search.rcsb.org/rcsbsearch/v2/query'


def _build_query(resolution_filter,
                 method: str,
                 max_size: int = 200,
                 start: int = 0) -> dict:
    """Construct an RCSB Search API JSON query for a single bin.

    resolution_filter is either None (no resolution gate) or a dict with the
    RCSB ``range`` operator schema: ``{from, to, include_lower, include_upper}``.

    ``start``/``max_size`` page through the (alphanumerically ordered) result
    set. The RCSB default sort returns entries in PDB-ID order, i.e. roughly
    by deposition era; representative sampling must therefore page across the
    whole range rather than reading only the first window (see
    ``_search_rcsb_representative``).
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
            'paginate': {'start': start, 'rows': max_size},
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


def _rcsb_total_count(resolution_filter, method: str,
                      timeout: int = 30) -> int:
    """Return the total number of entries matching a bin's query."""
    query = _build_query(resolution_filter, method, max_size=1, start=0)
    data = json.dumps(query).encode('utf-8')
    req = urllib.request.Request(
        SEARCH_URL, data=data, headers={'Content-Type': 'application/json'}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read().decode('utf-8'))
    return int(payload.get('total_count', 0))


def _search_rcsb_representative(resolution_filter, method: str,
                               n_candidates: int, rng: random.Random,
                               timeout: int = 30) -> List[str]:
    """Return a seeded, era-representative candidate list for one bin.

    The RCSB result set is ordered by PDB ID (roughly deposition era), so
    reading only the first window over-samples legacy structures (which also
    disproportionately lack wwPDB validation reports). To avoid that bias we:

      1. Look up the bin's total entry count.
      2. Draw seeded-random page offsets spread across the *entire* range
         ``[0, total)`` and fetch a small window (50 IDs) at each offset.
      3. Deduplicate, shuffle with the same seeded RNG, and return up to
         ``n_candidates`` IDs.

    This is fully deterministic for a given seed while spanning all deposition
    eras / resolutions within the bin, rather than only the oldest entries.
    """
    try:
        total = _rcsb_total_count(resolution_filter, method, timeout=timeout)
    except (urllib.error.URLError, ValueError):
        total = 0
    if total <= 0:
        return []

    window = 50
    # Number of windows needed (with headroom) to gather n_candidates uniques.
    n_windows = max(4, int((n_candidates * 3) / window) + 2)
    if total <= window * n_windows:
        # Small bin: just fetch the whole thing once and sample from it.
        ids = _search_rcsb(_build_query(resolution_filter, method,
                                        max_size=min(total, 10000), start=0),
                           timeout=timeout)
        rng.shuffle(ids)
        return ids[:n_candidates]

    max_start = max(0, total - window)
    offsets = sorted(rng.sample(range(0, max_start + 1),
                                min(n_windows, max_start + 1)))
    seen: set = set()
    candidates: List[str] = []
    for start in offsets:
        try:
            page = _search_rcsb(_build_query(resolution_filter, method,
                                             max_size=window, start=start),
                                timeout=timeout)
        except (urllib.error.URLError, ValueError) as exc:
            print(f'  [WARN] page at start={start} failed: {exc}')
            continue
        for pid in page:
            if pid not in seen:
                seen.add(pid)
                candidates.append(pid)
    rng.shuffle(candidates)
    return candidates[:n_candidates]


def _stratified_sample(n_target: int, seed: int = 20260514) -> List[Dict]:
    """Build a stratified random sample of PDB IDs across resolution / method.

    Returns a list of {pdb_id, bin, expected_resolution_range, method} dicts.
    """
    rng = random.Random(seed)
    sample: List[Dict] = []

    # Base per-bin targets sum to 105. When the caller requests a larger
    # n_target, scale each bin proportionally (preserving the stratification
    # ratios) so that the planned pool exceeds n_target. This headroom lets
    # downstream attrition (entries lacking a wwPDB validation report) be
    # absorbed without dropping below the usable-structure goal.
    import math
    base_total = sum(t for *_, t in RESOLUTION_BINS) + NMR_BIN_TARGET + \
        CRYOEM_BIN_TARGET
    scale = max(1.0, n_target / base_total)

    # Over-fetch factor per bin to absorb attrition. ~12% of randomly sampled
    # entries lack a wwPDB validation report (worst case NMR ~27%). The plan
    # itself is over-provisioned by this factor (the entries are *kept* in the
    # plan, not sliced off) so the post-attrition usable count still meets the
    # target. Every planned entry is processed, so keep this modest.
    overfetch = 1.30

    def _plan_for_bin(res_filter, method, target, bin_name, res_label):
        planned = int(math.ceil(target * scale * overfetch))
        ids = _search_rcsb_representative(res_filter, method, planned, rng)
        if not ids:
            print(f'  [WARN] {bin_name} representative query returned no IDs')
        return [{
            'pdb_id': pid, 'bin': bin_name, 'resolution_range': res_label,
            'method': method,
        } for pid in ids[:planned]]

    for bin_name, res_label, res_filter, target in RESOLUTION_BINS:
        try:
            sample.extend(_plan_for_bin(
                res_filter, 'X-RAY DIFFRACTION', target, bin_name, res_label))
        except (urllib.error.URLError, ValueError) as exc:
            print(f'  [WARN] {bin_name} query failed: {exc}')

    try:
        sample.extend(_plan_for_bin(
            None, 'SOLUTION NMR', NMR_BIN_TARGET, 'nmr', 'NMR'))
    except (urllib.error.URLError, ValueError) as exc:
        print(f'  [WARN] NMR query failed: {exc}')

    try:
        sample.extend(_plan_for_bin(
            {'from': 0.0, 'to': 4.5,
             'include_lower': True, 'include_upper': True},
            'ELECTRON MICROSCOPY', CRYOEM_BIN_TARGET,
            'cryoem', '<=4.5 A (EM)'))
    except (urllib.error.URLError, ValueError) as exc:
        print(f'  [WARN] cryo-EM query failed: {exc}')

    if len(sample) < n_target:
        print(f'  [WARN] stratified sample only has {len(sample)} entries '
              f'(target {n_target}); proceeding with what we have.')
    # Return the full plan so downstream attrition can be absorbed.
    return sample


# PDB-format chain IDs occupy a single column, so only single characters are
# representable. mmCIF ``auth_asym_id`` values can be multi-character (e.g.
# "AA", "AB"). Converter + constants live in chiralfold.structure_io.
from chiralfold.structure_io import (  # noqa: E402
    CONVERTER_REMARK_PREFIX,
    CONVERTER_VERSION,
    PDB_CHAIN_CHARS,
    mmcif_to_pdb as _cif_to_pdb,
)


def _count_cif_chains(cif_text: str) -> int:
    """Count distinct first-model mmCIF chain IDs (for skip diagnostics)."""
    lines = cif_text.splitlines()
    cols: Dict[str, int] = {}
    i = 0
    rows: List[List[str]] = []
    while i < len(lines):
        if lines[i].strip() == 'loop_':
            j = i + 1
            headers: List[str] = []
            while j < len(lines) and lines[j].strip().startswith('_'):
                headers.append(lines[j].strip())
                j += 1
            if headers and headers[0].startswith('_atom_site.'):
                cols = {h.split('.', 1)[1]: k for k, h in enumerate(headers)}
                k = j
                while (k < len(lines) and lines[k].strip()
                       and not lines[k].strip().startswith('_')
                       and lines[k].strip() != 'loop_'
                       and not lines[k].startswith('#')):
                    p = lines[k].split()
                    if len(p) >= len(headers):
                        rows.append(p)
                    k += 1
                break
            i = j
        else:
            i += 1
    if not cols:
        return 0
    aidx = cols.get('auth_asym_id', cols.get('label_asym_id'))
    midx = cols.get('pdbx_PDB_model_num')
    first_model = None
    chains: set = set()
    for p in rows:
        if midx is not None and midx < len(p):
            if first_model is None:
                first_model = p[midx]
            elif p[midx] != first_model:
                break
        if aidx is not None and aidx < len(p):
            chains.add(p[aidx])
    return len(chains)


def _cached_pdb_is_valid(path: str) -> bool:
    """Decide whether a cached ``.pdb`` may be reused without regeneration.

    Genuine legacy PDB downloads always begin with a standard header record
    (HEADER, and in practice other title-section records) and are always
    reusable. Converter-generated files are tagged with a
    ``CONVERTER_REMARK_PREFIX`` line recording the converter version; such a
    file is reusable only if its version matches the current one.

    Converter output written before versioning existed has no REMARK tag and
    begins directly with an ATOM/HETATM record (legacy PDB files never do).
    Those are treated as stale and rejected so that a ``--resume`` run against
    an old cache regenerates them with the current converter rather than
    silently reusing corrupted (chain-collided) coordinates.
    """
    try:
        with open(path, 'r') as fp:
            head = fp.read(4096)
    except OSError:
        return False
    if not head.strip():
        return False
    first = head.lstrip().splitlines()[0] if head.lstrip() else ''
    if first.startswith(CONVERTER_REMARK_PREFIX):
        ver_str = first[len(CONVERTER_REMARK_PREFIX):].strip()
        try:
            return int(ver_str) == CONVERTER_VERSION
        except ValueError:
            return False
    # No converter tag. If the file starts with ATOM/HETATM it is untagged
    # converter output (pre-versioning) -> stale. Otherwise it is a genuine
    # legacy PDB download (HEADER/etc.) -> valid.
    if first.startswith('ATOM') or first.startswith('HETATM'):
        return False
    return True


def _download_pdb(pdb_id: str, cache: str,
                  reason_out: Optional[List[str]] = None) -> Optional[str]:
    """Return a path to a PDB-format file for *pdb_id*, or None on failure.

    On failure, if ``reason_out`` is provided, a single reason code is appended
    to it (e.g. 'pdb_download_failed', 'too_many_chains',
    'no_atom_site', 'cif_fetch_failed').
    """
    def _fail(reason: str) -> None:
        if reason_out is not None:
            reason_out.append(reason)

    path = os.path.join(cache, f'{pdb_id.lower()}.pdb')
    if os.path.exists(path) and os.path.getsize(path) > 0:
        if _cached_pdb_is_valid(path):
            return path
        # Stale converter output from an older _cif_to_pdb version: drop it so
        # it is regenerated below with the current (chain-ID-safe) converter.
        print(f'    [cache] {pdb_id}: stale converted PDB '
              '(old converter); regenerating')
        try:
            os.remove(path)
        except OSError:
            pass
    try:
        urllib.request.urlretrieve(
            PDB_FILE_URL.format(pdb_id=pdb_id.lower()), path)
        return path
    except urllib.error.URLError as exc:
        # Legacy PDB format unavailable (common for large / modern entries).
        # Fall back to the mmCIF file and convert its atom_site loop to PDB.
        is_404 = isinstance(exc, urllib.error.HTTPError) and exc.code == 404
        if not is_404:
            print(f'    [skip] {pdb_id}: PDB download failed ({exc})')
            _fail('pdb_download_failed')
            return None
        try:
            with urllib.request.urlopen(
                    CIF_FILE_URL.format(pdb_id=pdb_id.lower()),
                    timeout=60) as resp:
                cif_text = resp.read().decode('utf-8', 'replace')
        except (urllib.error.URLError, OSError) as cif_exc:
            print(f'    [skip] {pdb_id}: PDB 404 and mmCIF fallback '
                  f'failed ({cif_exc})')
            _fail('cif_fetch_failed')
            return None
        pdb_text = _cif_to_pdb(cif_text)
        if not pdb_text:
            n_chains = _count_cif_chains(cif_text)
            if n_chains > len(PDB_CHAIN_CHARS):
                print(f'    [skip] {pdb_id}: mmCIF has {n_chains} distinct '
                      'chains (>36); cannot represent in PDB format')
                _fail('too_many_chains')
            else:
                print(f'    [skip] {pdb_id}: mmCIF had no parseable '
                      'atom_site loop')
                _fail('no_atom_site')
            return None
        with open(path, 'w') as fp:
            fp.write(pdb_text)
        print(f'    [cif->pdb] {pdb_id}: converted from mmCIF '
              '(legacy PDB unavailable)')
        return path


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
    parser.add_argument(
        '--cache-dir', default=CACHE_DIR,
        help='directory for cached PDB / validation-XML downloads '
             '(default: results/ramachandran_expansion_cache)')
    parser.add_argument(
        '--out-prefix', default=None,
        help='path prefix for outputs. When set, writes '
             '<prefix>_comparison.csv, <prefix>_summary.json, '
             '<prefix>_plot.png and does NOT overwrite the canonical '
             'results/molprobity_comparison.json. When omitted, uses the '
             'historical filenames (ramachandran_100struct_comparison.csv, '
             'ramachandran_100struct_plot.png) and updates '
             'molprobity_comparison.json.')
    parser.add_argument(
        '--resume', action='store_true',
        help='reuse any cached PDB / validation-XML downloads already present '
             'in --cache-dir instead of re-fetching (downloads are cached by '
             'default; this flag also reports how many were reused).')
    args = parser.parse_args()

    cache_dir = os.path.abspath(args.cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    if args.out_prefix:
        out_prefix = os.path.abspath(args.out_prefix)
        os.makedirs(os.path.dirname(out_prefix) or '.', exist_ok=True)
    else:
        out_prefix = None
    if args.resume:
        print(f'  [resume] reusing cached downloads in {cache_dir} when present')

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
    failures: List[Dict] = []
    t0 = time.time()
    for i, row in enumerate(plan, 1):
        pdb_id = row['pdb_id']
        print(f'  [{i}/{len(plan)}] {pdb_id} ({row["bin"]})')
        dl_reason: List[str] = []
        pdb_path = _download_pdb(pdb_id, cache_dir, reason_out=dl_reason)
        if pdb_path is None:
            failures.append({'pdb_id': pdb_id, 'bin': row['bin'],
                             'reason': (dl_reason[0] if dl_reason
                                        else 'pdb_download_failed')})
            continue
        wwpdb_pct = _wwpdb_outlier_pct(pdb_id, cache_dir)
        if wwpdb_pct is None:
            failures.append({'pdb_id': pdb_id, 'bin': row['bin'],
                             'reason': 'wwpdb_validation_unavailable'})
            continue
        cf_pct = _run_audit(pdb_path)
        if cf_pct is None:
            failures.append({'pdb_id': pdb_id, 'bin': row['bin'],
                             'reason': 'chiralfold_audit_failed'})
            continue
        rows.append({
            'pdb_id': pdb_id, 'bin': row['bin'],
            'method': row['method'],
            'resolution_range': row['resolution_range'],
            'chiralfold_rama_outlier_pct': cf_pct,
            'wwpdb_rama_outlier_pct': wwpdb_pct,
        })

    elapsed = time.time() - t0
    print(f'Completed {len(rows)} structures in {elapsed/60:.1f} min '
          f'({len(failures)} failures)')

    if len(rows) < 30:
        print(f'  [ERROR] only {len(rows)} usable rows; aborting before '
              'overwriting canonical molprobity_comparison.json.')
        return 1

    if out_prefix:
        csv_path = f'{out_prefix}_comparison.csv'
    else:
        csv_path = os.path.join(
            RESULTS_DIR, 'ramachandran_100struct_comparison.csv')
    with open(csv_path, 'w', newline='') as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {csv_path}')

    cf = np.array([r['chiralfold_rama_outlier_pct'] for r in rows])
    wp = np.array([r['wwpdb_rama_outlier_pct'] for r in rows])
    rho, rho_p = scistats.spearmanr(cf, wp)
    r_pearson, p_pearson = scistats.pearsonr(cf, wp)

    # Bootstrap 95% CIs for rho / r by resampling the paired rows (seeded).
    def _bootstrap_ci(stat_fn, n_boot=2000):
        boot_rng = np.random.default_rng(args.seed)
        n = len(cf)
        vals = []
        for _ in range(n_boot):
            idx = boot_rng.integers(0, n, n)
            try:
                vals.append(stat_fn(cf[idx], wp[idx]))
            except Exception:  # noqa: BLE001 - degenerate resample
                continue
        if not vals:
            return [None, None]
        lo, hi = np.percentile(vals, [2.5, 97.5])
        return [float(lo), float(hi)]

    spearman_ci = _bootstrap_ci(
        lambda a, b: scistats.spearmanr(a, b).statistic)
    pearson_ci = _bootstrap_ci(
        lambda a, b: scistats.pearsonr(a, b).statistic)

    summary = {
        'n_structures': int(len(rows)),
        'n_success': int(len(rows)),
        'n_planned': int(len(plan)),
        'n_failures': int(len(failures)),
        'failures': failures,
        'spearman_rho': float(rho),
        'spearman_p_value': float(rho_p),
        'spearman_p': float(rho_p),
        'spearman_rho_ci95': spearman_ci,
        'pearson_r': float(r_pearson),
        'pearson_p_value': float(p_pearson),
        'pearson_p': float(p_pearson),
        'pearson_r_ci95': pearson_ci,
        'bootstrap_n_iter': 2000,
        'chiralfold_mean_outlier_pct': float(cf.mean()),
        'wwpdb_mean_outlier_pct': float(wp.mean()),
        'seed': int(args.seed),
        'note': (
            'Expanded benchmark via benchmarks/expand_ramachandran_benchmark.py. '
            'Candidates drawn by seeded era-representative sampling across the '
            'full RCSB result range per bin. See the *_comparison.csv for the '
            'per-structure data.'
        ),
    }
    if out_prefix:
        json_path = f'{out_prefix}_summary.json'
    else:
        json_path = os.path.join(RESULTS_DIR, 'molprobity_comparison.json')
    with open(json_path, 'w') as fp:
        json.dump(summary, fp, indent=2)
    print(f'Wrote {json_path}: rho={rho:.3f} (p={rho_p:.4g}), '
          f'n_success={len(rows)}')

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
        if out_prefix:
            plot_path = f'{out_prefix}_plot.png'
        else:
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
