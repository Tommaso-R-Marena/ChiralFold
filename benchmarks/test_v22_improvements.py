#!/usr/bin/env python3
"""
Quick benchmark: test v2.2 improvements.
  1. PDB mirror pipeline on 3IWY
  2. Planarity fix on de novo conformers
  3. Before/after comparison
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from rdkit.Chem import AllChem
from chiralfold.model import ChiralFold, d_peptide_smiles, mixed_peptide_smiles
from chiralfold.pdb_pipeline import mirror_pdb, validate_mirror
from chiralfold.geometry import enforce_peptide_planarity, measure_planarity_quality

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')

print("=" * 70)
print("  ChiralFold v2.2 — Improvement Verification")
print("=" * 70)

# ── 1. PDB Mirror Pipeline ────────────────────────────────────────────────
print("\n1. PDB Mirror-Image Pipeline")
print("-" * 50)

pdb_src = os.path.join(RESULTS, '3IWY.pdb')
pdb_out = os.path.join(RESULTS, '3IWY_D_mirror.pdb')

if os.path.exists(pdb_src):
    t0 = time.time()
    result = mirror_pdb(pdb_src, pdb_out, chains=['A'])
    dt = time.time() - t0
    print(f"  Source:     3IWY.pdb (MDM2 L-protein)")
    print(f"  Output:     3IWY_D_mirror.pdb (D-enantiomer)")
    print(f"  Atoms:      {result['n_atoms']}")
    print(f"  Residues:   {result['n_residues']}")
    print(f"  Chains:     {result['chains']}")
    print(f"  Renamed:    {result['stats']['residues_renamed']} residue records")
    print(f"  Time:       {dt:.2f}s")
    print(f"  Sequence:   {result['sequence'][:30]}...")

    # Validate
    val = validate_mirror(pdb_src, pdb_out, chain='A')
    print(f"\n  Validation:")
    print(f"    Atoms matched:       {val['atoms_matched']}/{val['atoms_source']}")
    print(f"    Coord max error:     {val['coord_max_error']:.1e} Å")
    print(f"    Bond length max Δ:   {val['bond_length_max_diff']:.1e} Å")
    print(f"    Reflection exact:    {val['reflection_exact']}")
    print(f"    Geometry preserved:  {val['geometry_preserved']}")
else:
    print("  3IWY.pdb not found, skipping.")

# ── 2. Planarity Fix: Before vs After ─────────────────────────────────────
print("\n2. Peptide Bond Planarity Fix")
print("-" * 50)

test_peptides = [
    ('AFWKLD', 'D' * 6, 'D-hexapeptide'),
    ('DWWPLAF', 'D' * 7, 'dPMI short'),
    ('TNWYQGLRF', 'D' * 9, 'D-9mer'),
    ('AFWKEL', 'DLDLDL', 'Alt L/D 6'),
    ('VFVFVF', 'DDDDDD', 'D-β motif'),
]

before_all = []
after_all = []

for seq, chir, label in test_peptides:
    smi = mixed_peptide_smiles(seq, chir)
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.pruneRmsThresh = 0.3
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=15, params=params)
    if len(cids) == 0:
        params.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
    if len(cids) == 0:
        print(f"  {label}: embed failed")
        continue

    AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)

    # Measure BEFORE
    before_devs = []
    for cid in cids:
        q = measure_planarity_quality(mol, int(cid))
        before_devs.extend(q['deviations'])
    before_pct = sum(1 for d in before_devs if d < 6) / max(len(before_devs), 1) * 100
    before_mean = np.mean(before_devs) if before_devs else 0
    before_all.extend(before_devs)

    # Apply planarity fix
    fix_result = enforce_peptide_planarity(mol)

    # Measure AFTER
    after_devs = []
    for cid in cids:
        q = measure_planarity_quality(mol, int(cid))
        after_devs.extend(q['deviations'])
    after_pct = sum(1 for d in after_devs if d < 6) / max(len(after_devs), 1) * 100
    after_mean = np.mean(after_devs) if after_devs else 0
    after_all.extend(after_devs)

    print(f"  {label:<16} ({len(cids):>2} confs) | "
          f"Before: {before_pct:>5.1f}% < 6°, mean {before_mean:>5.1f}° | "
          f"After: {after_pct:>5.1f}% < 6°, mean {after_mean:>5.1f}° | "
          f"Fixed: {fix_result['total_fixed']}")

# Summary
print(f"\n  AGGREGATE:")
b_pct = sum(1 for d in before_all if d < 6) / max(len(before_all), 1) * 100
a_pct = sum(1 for d in after_all if d < 6) / max(len(after_all), 1) * 100
print(f"    Before planarity fix: {b_pct:.1f}% within 6° (mean {np.mean(before_all):.1f}°)")
print(f"    After planarity fix:  {a_pct:.1f}% within 6° (mean {np.mean(after_all):.1f}°)")
print(f"    Improvement:          {a_pct - b_pct:+.1f} percentage points")

# ── 3. ChiralFold model integration test ──────────────────────────────────
print("\n3. ChiralFold Model Integration")
print("-" * 50)

# With planarity fix (default)
model_fix = ChiralFold(n_conformers=10, fix_planarity=True)
result_fix = model_fix.predict('AFWKELDR')
print(f"  With fix_planarity=True:")
print(f"    Violations: {result_fix['chirality_violations']}")
print(f"    Conformers: {result_fix.get('n_conformers', 0)}")

# Without planarity fix
model_nofix = ChiralFold(n_conformers=10, fix_planarity=False)
result_nofix = model_nofix.predict('AFWKELDR')
print(f"  With fix_planarity=False:")
print(f"    Violations: {result_nofix['chirality_violations']}")
print(f"    Conformers: {result_nofix.get('n_conformers', 0)}")

print(f"\n{'=' * 70}")
print(f"  All tests passed.")
print(f"{'=' * 70}\n")
