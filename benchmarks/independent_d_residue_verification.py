#!/usr/bin/env python3
"""
Independent D-Residue Chirality Verification
==============================================

Scans ALL downloaded PDB files containing D-amino acid residues,
extracts raw N/CA/C/CB coordinates for every D-labeled residue,
computes the signed tetrahedron volume, and classifies each as
D-correct (negative) or L-error (positive).

THIS SCRIPT USES ONLY:
  - Standard library (os, csv, glob)
  - numpy (for cross product and dot product)
  - Raw PDB file parsing

NO ChiralFold code is imported or used. This is an independent
ground-truth verification.

Output:
  results/d_residue_verification.csv — one row per D-amino acid residue
  results/d_residue_verification_summary.json — aggregate statistics
"""

import os
import csv
import json
import glob
import numpy as np
from collections import Counter

# Known D-amino acid 3-letter codes in PDB
D_AA_CODES = {
    'DAL', 'DAR', 'DSG', 'DAS', 'DCY', 'DGL', 'DGN', 'DHI', 'DIL',
    'DLE', 'DLY', 'MED', 'DPN', 'DPR', 'DSN', 'DTH', 'DTR', 'DTY', 'DVA',
}

# Map D-AA code to one-letter code for readability
D_TO_ONE = {
    'DAL': 'A', 'DAR': 'R', 'DSG': 'N', 'DAS': 'D', 'DCY': 'C',
    'DGL': 'E', 'DGN': 'Q', 'DHI': 'H', 'DIL': 'I', 'DLE': 'L',
    'DLY': 'K', 'MED': 'M', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S',
    'DTH': 'T', 'DTR': 'W', 'DTY': 'Y', 'DVA': 'V',
}


def signed_volume(n_pos, ca_pos, c_pos, cb_pos):
    """
    Compute signed volume of the tetrahedron formed by N, CA, C, CB
    relative to CA.

    Returns:
        float: negative = D-chirality, positive = L-chirality
    """
    v1 = n_pos - ca_pos
    v2 = c_pos - ca_pos
    v3 = cb_pos - ca_pos
    return float(np.dot(v1, np.cross(v2, v3)))


def extract_d_residues(pdb_path):
    """
    Parse a PDB file and extract all D-amino acid residues with
    their backbone atom coordinates.

    Returns:
        list of dicts with keys: pdb_id, chain, resnum, resname,
        one_letter, n_xyz, ca_xyz, c_xyz, cb_xyz, has_all_atoms
    """
    pdb_id = os.path.basename(pdb_path).replace('.pdb', '').replace('.PDB', '').upper()
    
    # Collect atoms by (chain, resnum, resname)
    residues = {}
    
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            if len(line) < 54:
                continue
            
            resname = line[17:20].strip()
            if resname not in D_AA_CODES:
                continue
            
            chain = line[21]
            try:
                resnum = int(line[22:26])
            except ValueError:
                # Handle non-integer residue numbers (insertion codes etc.)
                resnum_str = line[22:27].strip()
                try:
                    resnum = int(''.join(c for c in resnum_str if c.isdigit()))
                except ValueError:
                    continue
            
            aname = line[12:16].strip()
            
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            
            key = (chain, resnum, resname)
            if key not in residues:
                residues[key] = {
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'resnum': resnum,
                    'resname': resname,
                    'one_letter': D_TO_ONE.get(resname, '?'),
                    'atoms': {},
                }
            residues[key]['atoms'][aname] = np.array([x, y, z])
    
    # Process each residue
    results = []
    for key, res in sorted(residues.items()):
        atoms = res['atoms']
        has_n = 'N' in atoms
        has_ca = 'CA' in atoms
        has_c = 'C' in atoms
        has_cb = 'CB' in atoms
        has_all = has_n and has_ca and has_c and has_cb
        
        entry = {
            'pdb_id': res['pdb_id'],
            'chain': res['chain'],
            'resnum': res['resnum'],
            'resname': res['resname'],
            'one_letter': res['one_letter'],
            'has_n': has_n,
            'has_ca': has_ca,
            'has_c': has_c,
            'has_cb': has_cb,
            'has_all_atoms': has_all,
            'signed_volume': None,
            'chirality': None,
            'is_error': None,
        }
        
        if has_all:
            vol = signed_volume(atoms['N'], atoms['CA'], atoms['C'], atoms['CB'])
            entry['signed_volume'] = round(vol, 4)
            
            if abs(vol) < 0.01:
                entry['chirality'] = 'flat'
                entry['is_error'] = 'ambiguous'
            elif vol < 0:
                entry['chirality'] = 'D'
                entry['is_error'] = False
            else:
                entry['chirality'] = 'L'
                entry['is_error'] = True
            
            # Store raw coordinates for verification
            entry['n_xyz'] = f"{atoms['N'][0]:.3f},{atoms['N'][1]:.3f},{atoms['N'][2]:.3f}"
            entry['ca_xyz'] = f"{atoms['CA'][0]:.3f},{atoms['CA'][1]:.3f},{atoms['CA'][2]:.3f}"
            entry['c_xyz'] = f"{atoms['C'][0]:.3f},{atoms['C'][1]:.3f},{atoms['C'][2]:.3f}"
            entry['cb_xyz'] = f"{atoms['CB'][0]:.3f},{atoms['CB'][1]:.3f},{atoms['CB'][2]:.3f}"
        else:
            missing = []
            if not has_n: missing.append('N')
            if not has_ca: missing.append('CA')
            if not has_c: missing.append('C')
            if not has_cb: missing.append('CB')
            entry['chirality'] = 'incomplete'
            entry['is_error'] = 'cannot_check'
            entry['n_xyz'] = ''
            entry['ca_xyz'] = ''
            entry['c_xyz'] = ''
            entry['cb_xyz'] = ''
        
        results.append(entry)
    
    return results


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import time
    t_start = time.time()
    
    results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    survey_dir = os.path.join(results_dir, 'd_survey')
    
    print("=" * 70)
    print("  Independent D-Residue Chirality Verification")
    print("  (no ChiralFold code — numpy + raw PDB only)")
    print("=" * 70)
    
    # Find all PDB files
    pdb_files = sorted(glob.glob(os.path.join(survey_dir, '*.pdb')))
    
    # Also check the results/ directory for individually downloaded PDBs
    for f in glob.glob(os.path.join(results_dir, '*.pdb')):
        if f not in pdb_files:
            pdb_files.append(f)
    
    print(f"\n  PDB files to scan: {len(pdb_files)}")
    
    all_residues = []
    files_with_d_aa = 0
    
    for i, pdb_path in enumerate(pdb_files):
        try:
            residues = extract_d_residues(pdb_path)
            if residues:
                all_residues.extend(residues)
                files_with_d_aa += 1
            if (i + 1) % 50 == 0:
                print(f"  Scanned {i+1}/{len(pdb_files)} files, "
                      f"found {len(all_residues)} D-AA residues so far...")
        except Exception as e:
            pass
    
    print(f"\n  Total D-AA residues found: {len(all_residues)}")
    print(f"  Files containing D-AAs: {files_with_d_aa}")
    
    # ── Statistics ─────────────────────────────────────────────────────────
    checkable = [r for r in all_residues if r['has_all_atoms']]
    d_correct = [r for r in checkable if r['chirality'] == 'D']
    l_error = [r for r in checkable if r['chirality'] == 'L']
    flat = [r for r in checkable if r['chirality'] == 'flat']
    incomplete = [r for r in all_residues if not r['has_all_atoms']]
    
    print(f"\n  With all 4 atoms (N/CA/C/CB): {len(checkable)}")
    print(f"  Missing atoms (cannot check): {len(incomplete)}")
    print(f"  D-chirality (correct):        {len(d_correct)}")
    print(f"  L-chirality (ERROR):          {len(l_error)}")
    print(f"  Flat/ambiguous:               {len(flat)}")
    
    error_rate = len(l_error) / len(checkable) * 100 if checkable else 0
    print(f"\n  Error rate: {len(l_error)}/{len(checkable)} = {error_rate:.1f}%")
    
    # ── Error details ─────────────────────────────────────────────────────
    if l_error:
        print(f"\n  {'PDB':<6}{'Chain':<6}{'Res#':<7}{'Name':<5}{'AA':<4}{'Volume':>9}")
        print(f"  {'-'*40}")
        for r in sorted(l_error, key=lambda x: (x['pdb_id'], x['chain'], x['resnum'])):
            print(f"  {r['pdb_id']:<6}{r['chain']:<6}{r['resnum']:<7}"
                  f"{r['resname']:<5}{r['one_letter']:<4}{r['signed_volume']:>+9.4f}")
    
    # ── By residue type ───────────────────────────────────────────────────
    print(f"\n  Errors by D-AA type:")
    type_counts = Counter(r['resname'] for r in l_error)
    for resname, count in type_counts.most_common():
        total_of_type = sum(1 for r in checkable if r['resname'] == resname)
        print(f"    {resname}: {count}/{total_of_type} errors")
    
    # ── By PDB structure ──────────────────────────────────────────────────
    error_pdbs = Counter(r['pdb_id'] for r in l_error)
    print(f"\n  Errors by PDB structure:")
    for pdb_id, count in error_pdbs.most_common():
        print(f"    {pdb_id}: {count} error(s)")
    
    # ── Save CSV ──────────────────────────────────────────────────────────
    csv_path = os.path.join(results_dir, 'd_residue_verification.csv')
    fieldnames = [
        'pdb_id', 'chain', 'resnum', 'resname', 'one_letter',
        'has_all_atoms', 'signed_volume', 'chirality', 'is_error',
        'n_xyz', 'ca_xyz', 'c_xyz', 'cb_xyz',
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in sorted(all_residues, key=lambda x: (x['pdb_id'], x['chain'], x['resnum'])):
            row = {k: r.get(k, '') for k in fieldnames}
            writer.writerow(row)
    print(f"\n  Saved: {csv_path} ({len(all_residues)} rows)")
    
    # ── Save summary JSON ─────────────────────────────────────────────────
    summary = {
        'description': 'Independent D-amino acid chirality verification using signed tetrahedron volume',
        'method': 'signed_volume = dot(N-CA, cross(C-CA, CB-CA)). Negative = D, Positive = L.',
        'code_used': 'numpy only, no ChiralFold',
        'pdb_files_scanned': len(pdb_files),
        'files_with_d_aa': files_with_d_aa,
        'total_d_residues': len(all_residues),
        'checkable_residues': len(checkable),
        'incomplete_residues': len(incomplete),
        'd_correct': len(d_correct),
        'l_error': len(l_error),
        'flat_ambiguous': len(flat),
        'error_rate_pct': round(error_rate, 2),
        'errors': [
            {
                'pdb_id': r['pdb_id'],
                'chain': r['chain'],
                'resnum': r['resnum'],
                'resname': r['resname'],
                'one_letter': r['one_letter'],
                'signed_volume': r['signed_volume'],
            }
            for r in sorted(l_error, key=lambda x: (x['pdb_id'], x['chain'], x['resnum']))
        ],
        'errors_by_type': dict(type_counts.most_common()),
        'errors_by_structure': dict(error_pdbs.most_common()),
    }
    
    json_path = os.path.join(results_dir, 'd_residue_verification_summary.json')
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {json_path}")
    
    elapsed = time.time() - t_start
    print(f"\n  Done in {elapsed:.1f}s")
    print("=" * 70)
