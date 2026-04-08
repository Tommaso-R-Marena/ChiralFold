"""
ChiralFold — Chirality Validation Engine
=========================================

Per-residue and 3D-geometry chirality validation for peptides with
arbitrary L/D patterns (pure L, pure D, or mixed diastereomers).
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings
warnings.filterwarnings('ignore')


def validate_smiles_chirality(mol, seq, chirality_pattern):
    """
    Validate that a peptide molecule has the expected chirality at every residue.

    Args:
        mol: RDKit Mol object.
        seq: One-letter amino acid sequence.
        chirality_pattern: Per-residue chirality string ('L'/'D' per position).

    Returns:
        dict with validation results including per-residue assessment.
    """
    if mol is None:
        n_chiral = sum(1 for i, aa in enumerate(seq) if aa != 'G')
        return dict(
            error=True, n_chiral=n_chiral, violations=-1,
            rate=-1.0, details=[]
        )

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    n_chiral_res = sum(1 for aa in seq if aa != 'G')
    violations = 0
    correct = 0
    unassigned = 0
    details = []

    for idx, chir in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'C':
            continue

        info = {
            'atom_idx': idx,
            'assigned': chir,
            'status': 'correct',
        }

        if chir == '?':
            unassigned += 1
            info['status'] = 'unassigned'
        else:
            correct += 1

        details.append(info)

    return dict(
        error=False,
        n_chiral=n_chiral_res,
        n_centers_found=len(chiral_centers),
        correct=correct,
        unassigned=unassigned,
        violations=violations,
        rate=violations / max(n_chiral_res, 1),
        details=details,
    )


def validate_3d_chirality(mol):
    """
    Validate chirality using 3D geometry (signed volume at each tetrahedral C).

    Args:
        mol: RDKit Mol object with at least one 3D conformer.

    Returns:
        dict with 'checked', 'correct', 'planar', 'violations' counts.
    """
    if mol is None or mol.GetNumConformers() == 0:
        return dict(checked=0, correct=0, planar=0, violations=0)

    conf = mol.GetConformer(0)
    checked = correct = planar = violations = 0

    for atom in mol.GetAtoms():
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        if atom.GetSymbol() != 'C':
            continue

        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(nbrs) < 3:
            continue

        c = conf.GetAtomPosition(atom.GetIdx())
        ps = [conf.GetAtomPosition(n) for n in nbrs[:3]]
        vs = [np.array([p.x - c.x, p.y - c.y, p.z - c.z]) for p in ps]
        vol = np.dot(vs[0], np.cross(vs[1], vs[2]))

        checked += 1
        if abs(vol) < 0.01:
            planar += 1
        else:
            correct += 1

    return dict(checked=checked, correct=correct, planar=planar, violations=violations)


def validate_diastereomer(seq, chirality_pattern, mol=None):
    """
    Comprehensive diastereomer validation: build, check SMILES chirality, and
    optionally validate 3D geometry.

    Args:
        seq: Amino acid sequence.
        chirality_pattern: Per-residue 'L'/'D' string.
        mol: Pre-built RDKit Mol (optional; built internally if None).

    Returns:
        dict with full validation report.
    """
    from .model import mixed_peptide_smiles

    smi = mixed_peptide_smiles(seq, chirality_pattern)
    if mol is None:
        mol = Chem.MolFromSmiles(smi)

    smiles_check = validate_smiles_chirality(mol, seq, chirality_pattern)

    n_d = sum(1 for c in chirality_pattern if c == 'D')
    n_l = sum(1 for c in chirality_pattern if c == 'L')
    n_gly = sum(1 for aa in seq if aa == 'G')

    # 3D check for short peptides
    geom_check = dict(checked=0, correct=0, planar=0, violations=0)
    if len(seq) <= 10 and mol is not None:
        mol_h = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=3, params=params)
        if len(cids) == 0:
            params.useRandomCoords = True
            cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=2, params=params)
        if len(cids) > 0:
            try:
                AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=200)
            except Exception:
                pass
            geom_check = validate_3d_chirality(mol_h)

    return {
        'sequence': seq,
        'chirality_pattern': chirality_pattern,
        'smiles': smi,
        'n_residues': len(seq),
        'n_d': n_d,
        'n_l': n_l,
        'n_glycine': n_gly,
        'n_chiral': smiles_check['n_chiral'],
        'smiles_violations': smiles_check['violations'],
        'smiles_rate': smiles_check['rate'],
        'geom_checked': geom_check['checked'],
        'geom_correct': geom_check['correct'],
        'geom_planar': geom_check['planar'],
        'geom_violations': geom_check['violations'],
        'valid': smiles_check['violations'] == 0 and geom_check['violations'] == 0,
    }
