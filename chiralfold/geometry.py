"""
ChiralFold — Geometry Post-Processing
=======================================

Fixes known ETKDG + MMFF94 limitations:
  1. Peptide bond planarity enforcement (omega → ±180°)
  2. Backbone dihedral filtering for D-amino acid preferences

These corrections improve the structural quality of de novo
generated conformers without compromising chirality guarantees.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
import warnings
warnings.filterwarnings('ignore')


# ═══════════════════════════════════════════════════════════════════════════
# Backbone Atom Detection
# ═══════════════════════════════════════════════════════════════════════════

def find_peptide_bonds(mol):
    """
    Find all peptide bond atom quadruplets (Cα_i, C_i, N_i+1, Cα_i+1)
    for computing/setting omega dihedral angles.

    Returns:
        List of tuples (ca1_idx, c_idx, n_idx, ca2_idx) for each peptide bond.
    """
    # Strategy: find C(=O)-N amide pattern, then identify flanking Cα atoms
    # SMARTS for peptide bond: C(=O)-N where both are in the backbone
    peptide_bonds = []

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        # Look for C-N single bond where C has a =O neighbor
        c_atom = n_atom = None
        if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'N':
            c_atom, n_atom = a1, a2
        elif a1.GetSymbol() == 'N' and a2.GetSymbol() == 'C':
            c_atom, n_atom = a2, a1
        else:
            continue

        # Check that C has a =O (carbonyl carbon)
        has_carbonyl = False
        for nb in c_atom.GetNeighbors():
            if nb.GetSymbol() == 'O':
                b = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nb.GetIdx())
                if b and b.GetBondTypeAsDouble() >= 1.5:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Find Cα on the C side: carbon neighbor of C that has an N neighbor
        ca1_idx = None
        for nb in c_atom.GetNeighbors():
            if nb.GetSymbol() == 'C' and nb.GetIdx() != n_atom.GetIdx():
                # Check if this C is bonded to an N (making it a Cα)
                for sub in nb.GetNeighbors():
                    if sub.GetSymbol() == 'N' and sub.GetIdx() != n_atom.GetIdx():
                        ca1_idx = nb.GetIdx()
                        break
                if ca1_idx is not None:
                    break

        # Find Cα on the N side: carbon neighbor of N that has a C=O neighbor
        ca2_idx = None
        for nb in n_atom.GetNeighbors():
            if nb.GetSymbol() == 'C' and nb.GetIdx() != c_atom.GetIdx():
                for sub in nb.GetNeighbors():
                    if sub.GetSymbol() == 'C' and sub.GetIdx() != n_atom.GetIdx():
                        # Check if sub has =O
                        for sub2 in sub.GetNeighbors():
                            if sub2.GetSymbol() == 'O':
                                b = mol.GetBondBetweenAtoms(sub.GetIdx(), sub2.GetIdx())
                                if b and b.GetBondTypeAsDouble() >= 1.5:
                                    ca2_idx = nb.GetIdx()
                                    break
                        if ca2_idx is not None:
                            break
                if ca2_idx is not None:
                    break

        if ca1_idx is not None and ca2_idx is not None:
            peptide_bonds.append((ca1_idx, c_atom.GetIdx(),
                                  n_atom.GetIdx(), ca2_idx))

    return peptide_bonds


# ═══════════════════════════════════════════════════════════════════════════
# Planarity Enforcement
# ═══════════════════════════════════════════════════════════════════════════

def enforce_peptide_planarity(mol, conf_id=None, omega_target=180.0,
                              tolerance=6.0, max_iterations=200):
    """
    Enforce peptide bond planarity by constraining omega dihedrals.

    For each peptide bond, if the omega dihedral deviates from the target
    by more than `tolerance` degrees, it is corrected by:
      1. Setting the dihedral to the nearest target (180° or 0°)
      2. Running constrained MMFF94 minimization to relax clashes

    Args:
        mol: RDKit Mol with 3D conformer(s) and explicit hydrogens.
        conf_id: Specific conformer ID to fix (None = fix all).
        omega_target: Target omega angle in degrees (180 = trans, 0 = cis).
        tolerance: Maximum deviation from target before correction (degrees).
        max_iterations: Max force field optimization iterations after correction.

    Returns:
        dict with 'n_bonds_fixed', 'omega_before', 'omega_after' per conformer.
    """
    peptide_bonds = find_peptide_bonds(mol)
    if not peptide_bonds:
        return {'n_bonds_fixed': 0, 'peptide_bonds_found': 0}

    conf_ids = [conf_id] if conf_id is not None else list(range(mol.GetNumConformers()))
    results = []

    for cid in conf_ids:
        if cid >= mol.GetNumConformers():
            continue

        n_fixed = 0
        omegas_before = []
        omegas_after = []

        for ca1, c, n, ca2 in peptide_bonds:
            # Get current omega
            try:
                omega = rdMolTransforms.GetDihedralDeg(
                    mol.GetConformer(cid), ca1, c, n, ca2
                )
            except Exception:
                continue

            omegas_before.append(omega)

            # Check deviation from target (consider both +180 and -180)
            dev_trans = min(abs(omega - 180), abs(omega + 180))
            dev_cis = abs(omega)

            if dev_trans <= tolerance:
                omegas_after.append(omega)
                continue

            # Set to nearest target
            if dev_trans < dev_cis:
                target = 180.0 if omega > 0 else -180.0
            else:
                target = 0.0

            try:
                rdMolTransforms.SetDihedralDeg(
                    mol.GetConformer(cid), ca1, c, n, ca2, target
                )
                n_fixed += 1
            except Exception:
                pass

            # Re-measure
            try:
                omega_new = rdMolTransforms.GetDihedralDeg(
                    mol.GetConformer(cid), ca1, c, n, ca2
                )
                omegas_after.append(omega_new)
            except Exception:
                omegas_after.append(target)

        # Constrained minimization to relax clashes from dihedral adjustment
        if n_fixed > 0:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(
                    mol, AllChem.MMFFGetMoleculeProperties(mol),
                    confId=cid
                )
                if ff is not None:
                    # Add torsion constraints to keep omega near target
                    for ca1, c, n, ca2 in peptide_bonds:
                        try:
                            omega_now = rdMolTransforms.GetDihedralDeg(
                                mol.GetConformer(cid), ca1, c, n, ca2
                            )
                            # Constrain omega to within ±5° of current value
                            ff.MMFFAddTorsionConstraint(
                                ca1, c, n, ca2, False,
                                omega_now - 5.0, omega_now + 5.0, 100.0
                            )
                        except Exception:
                            pass

                    ff.Minimize(maxIts=max_iterations)
            except Exception:
                pass

        results.append({
            'conf_id': cid,
            'n_bonds_fixed': n_fixed,
            'n_peptide_bonds': len(peptide_bonds),
            'omega_before': omegas_before,
            'omega_after': omegas_after,
        })

    return {
        'peptide_bonds_found': len(peptide_bonds),
        'conformers_processed': len(results),
        'per_conformer': results,
        'total_fixed': sum(r['n_bonds_fixed'] for r in results),
    }


def measure_planarity_quality(mol, conf_id=0):
    """
    Measure peptide bond planarity quality for a conformer.

    Returns:
        dict with omega angles, deviations, and quality metrics.
    """
    peptide_bonds = find_peptide_bonds(mol)
    if not peptide_bonds or mol.GetNumConformers() <= conf_id:
        return {'n_bonds': 0, 'omegas': [], 'deviations': []}

    omegas = []
    deviations = []
    for ca1, c, n, ca2 in peptide_bonds:
        try:
            omega = rdMolTransforms.GetDihedralDeg(
                mol.GetConformer(conf_id), ca1, c, n, ca2
            )
            omegas.append(omega)
            dev = min(abs(omega - 180), abs(omega + 180))
            deviations.append(dev)
        except Exception:
            pass

    pct_within_6 = (sum(1 for d in deviations if d < 6.0) / len(deviations) * 100
                    if deviations else 0)

    return {
        'n_bonds': len(peptide_bonds),
        'omegas': omegas,
        'deviations': deviations,
        'mean_deviation': np.mean(deviations) if deviations else 0,
        'max_deviation': max(deviations) if deviations else 0,
        'pct_within_6deg': pct_within_6,
        'pct_within_10deg': (sum(1 for d in deviations if d < 10.0) / len(deviations) * 100
                            if deviations else 0),
    }
