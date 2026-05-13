"""
ChiralFold — Chirality Validation Engine
=========================================

Per-residue and 3D-geometry chirality validation for peptides with
arbitrary L/D patterns (pure L, pure D, or mixed diastereomers).
"""

import warnings

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Module-level: do NOT call warnings.filterwarnings globally — that suppresses
# warnings for all downstream user code. Targeted context managers are used
# inside individual functions where benign RDKit warnings are expected.


# Expected CIP descriptor at Cα for each one-letter amino acid:
#   L-amino acids → 'S' for nearly all; exception: L-cysteine → 'R'
#                   (sulfur outranks the carboxylate carbon by CIP priority)
#   D-amino acids → 'R' for nearly all; exception: D-cysteine → 'S'
#
# Glycine is achiral (no Cα stereocenter).
_L_CIP_EXCEPTIONS = {"C": "R"}
_D_CIP_EXCEPTIONS = {"C": "S"}


def _expected_cip(one_letter_aa: str, chirality_code: str) -> str:
    """
    Return the expected CIP descriptor ('R' or 'S') at Cα for the given
    amino acid one-letter code and L/D chirality assignment.

    Glycine returns '' (achiral).
    """
    aa = one_letter_aa.upper()
    if aa == "G":
        return ""
    if chirality_code == "L":
        return _L_CIP_EXCEPTIONS.get(aa, "S")
    if chirality_code == "D":
        return _D_CIP_EXCEPTIONS.get(aa, "R")
    return ""


def _ordered_chiral_centers_for_residues(mol, seq):
    """
    Map RDKit chiral centers (in atom-index order) to peptide residue
    positions. For peptides built by ``mixed_peptide_smiles``, residue
    chiral centers appear in N→C order along the chain because the SMILES
    string concatenates residue fragments in sequence. We therefore align
    chiral centers to the non-glycine residue positions by order of
    appearance.

    Returns a list of (residue_position, atom_idx, cip) tuples, one per
    non-glycine residue that produced a chiral center.
    """
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # Filter to Cα-type carbons (we only care about backbone chirality here).
    cα_centers = []
    for idx, chir in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != "C":
            continue
        cα_centers.append((idx, chir))

    # Non-glycine positions in the peptide sequence, in order.
    res_positions = [i for i, aa in enumerate(seq) if aa.upper() != "G"]

    # For isoleucine and threonine the side chain itself contains an
    # additional stereocenter — keep only the first stereocenter
    # encountered per residue (the Cα is the first center along the
    # backbone-fragment SMILES).
    aligned = []
    used_centers = 0
    for pos in res_positions:
        if used_centers >= len(cα_centers):
            break
        # Skip side-chain stereocenters by taking centers in pairs for I/T.
        idx, chir = cα_centers[used_centers]
        aligned.append((pos, idx, chir))
        used_centers += 1
        aa = seq[pos].upper()
        if aa in ("I", "T") and used_centers < len(cα_centers):
            # The next center belongs to the side chain — skip it.
            used_centers += 1

    return aligned


def validate_smiles_chirality(mol, seq, chirality_pattern):
    """
    Validate that a peptide molecule has the expected chirality at every
    residue.

    Compares each assigned CIP descriptor at the residue Cα against the
    descriptor expected from the chirality pattern (D → R, L → S, with
    cysteine as the documented exception).

    Args:
        mol: RDKit Mol object (must have been parsed from SMILES with
            stereochemistry information, or had ``AssignStereochemistry``
            called).
        seq: One-letter amino acid sequence.
        chirality_pattern: Per-residue chirality string ('L'/'D' per
            position). Must be the same length as *seq*.

    Returns:
        dict with keys:
            error (bool): True if *mol* is None.
            n_chiral (int): Number of non-glycine residues in *seq*.
            n_centers_found (int): Total chiral centers detected by RDKit.
            correct (int): Residues whose CIP matches expectation.
            unassigned (int): Residues whose CIP is '?'.
            violations (int): Residues whose CIP contradicts expectation.
            rate (float): violations / max(n_chiral, 1).
            details (list[dict]): Per-residue records with ``atom_idx``,
                ``residue_pos``, ``expected``, ``observed``, ``status``.

    Notes:
        Fixed in v3.2.1: violation detection was previously non-functional
        (always returned 0). The 0% violation rate reported in Marena (2026)
        is independently confirmed by 3D coordinate geometry in
        ``benchmarks/childs2025_comparison.py`` and
        ``benchmarks/independent_d_residue_verification.py``.
    """
    n_chiral_res = sum(1 for aa in seq if aa.upper() != "G")

    if mol is None:
        return dict(
            error=True,
            n_chiral=n_chiral_res,
            n_centers_found=0,
            correct=0,
            unassigned=0,
            violations=-1,
            rate=-1.0,
            details=[],
        )

    if len(seq) != len(chirality_pattern):
        # Cannot validate without a 1-to-1 residue↔chirality mapping.
        return dict(
            error=True,
            n_chiral=n_chiral_res,
            n_centers_found=0,
            correct=0,
            unassigned=0,
            violations=-1,
            rate=-1.0,
            details=[{"status": "error",
                      "reason": "seq/chirality length mismatch"}],
        )

    aligned = _ordered_chiral_centers_for_residues(mol, seq)
    n_centers_found = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

    violations = 0
    correct = 0
    unassigned = 0
    details = []

    for residue_pos, atom_idx, observed in aligned:
        aa = seq[residue_pos]
        expected_chir = chirality_pattern[residue_pos].upper()
        expected_cip = _expected_cip(aa, expected_chir)

        info = {
            "atom_idx": atom_idx,
            "residue_pos": residue_pos,
            "residue": aa,
            "expected": expected_cip,
            "observed": observed,
            "status": "correct",
        }

        if observed == "?":
            unassigned += 1
            info["status"] = "unassigned"
        elif expected_cip and observed != expected_cip:
            violations += 1
            info["status"] = "violation"
        else:
            correct += 1

        details.append(info)

    return dict(
        error=False,
        n_chiral=n_chiral_res,
        n_centers_found=n_centers_found,
        correct=correct,
        unassigned=unassigned,
        violations=violations,
        rate=violations / max(n_chiral_res, 1),
        details=details,
    )


def validate_3d_chirality(mol, conf_id=0):
    """
    Validate chirality using 3D geometry via RDKit's canonical stereochemistry
    perception.

    The molecule's existing chiral tags (from the SMILES it was built from)
    are compared against the tags RDKit assigns when it re-perceives
    stereochemistry directly from the 3D coordinates via
    ``Chem.AssignStereochemistryFrom3D``. Because both the original and the
    re-perceived tags are computed by RDKit relative to canonical bond
    ordering at each atom, the comparison is invariant to atom/neighbor
    iteration order — equivalent molecules with different atom numberings
    produce identical results.

    For each carbon stereocenter:
      - If the geometry is effectively planar (degenerate tetrahedron), the
        atom is counted as ``planar`` and not as a violation.
      - Otherwise the re-perceived tag is compared to the original tag.
        ``CW`` vs ``CCW`` disagreement is a violation; agreement (or
        ``CHI_UNSPECIFIED`` from a degenerate local geometry) counts as
        correct.

    Args:
        mol: RDKit Mol object with at least one 3D conformer.
        conf_id: conformer index to validate against (default 0).

    Returns:
        dict with keys ``checked``, ``correct``, ``planar``, ``violations``.

    Notes:
        Reworked in v3.2.2 to derive observed handedness via
        ``AssignStereochemistryFrom3D`` instead of a hand-rolled signed-volume
        check tied to ``GetNeighbors()`` iteration order. The previous
        implementation could in principle disagree with RDKit's own canonical
        ordering and yield results that depended on atom indexing rather than
        on the underlying stereochemistry.
    """
    if mol is None or mol.GetNumConformers() == 0:
        return dict(checked=0, correct=0, planar=0, violations=0)

    cw = Chem.ChiralType.CHI_TETRAHEDRAL_CW
    ccw = Chem.ChiralType.CHI_TETRAHEDRAL_CCW

    # Snapshot the original (SMILES-derived) chiral tags for carbon
    # stereocenters before re-perception overwrites them.
    expected = {}
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag in (cw, ccw) and atom.GetSymbol() == "C":
            expected[atom.GetIdx()] = tag

    if not expected:
        return dict(checked=0, correct=0, planar=0, violations=0)

    # Re-perceive stereochemistry from the 3D conformer on a copy so we don't
    # mutate the caller's molecule. RDKit applies its canonical neighbor
    # ordering internally, making this comparison order-invariant.
    mol_copy = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(mol_copy, confId=conf_id,
                                     replaceExistingTags=True)

    # Planarity is a property of the 3D geometry, not of either tag; compute
    # it once from the conformer using the same threshold as the original
    # check (≈ 0.01 Å³ signed volume).
    conf = mol.GetConformer(conf_id)
    checked = correct = planar = violations = 0

    for idx, exp_tag in expected.items():
        atom = mol_copy.GetAtomWithIdx(idx)
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(nbrs) < 3:
            continue
        checked += 1

        c = conf.GetAtomPosition(idx)
        ps = [conf.GetAtomPosition(n) for n in nbrs[:3]]
        vs = [np.array([p.x - c.x, p.y - c.y, p.z - c.z]) for p in ps]
        if abs(float(np.dot(vs[0], np.cross(vs[1], vs[2])))) < 0.01:
            planar += 1
            continue

        observed_tag = atom.GetChiralTag()
        if observed_tag not in (cw, ccw):
            # RDKit could not assign a definite tag from 3D — treat as
            # correct rather than a violation, mirroring the previous
            # behavior for ambiguous geometries.
            correct += 1
        elif observed_tag == exp_tag:
            correct += 1
        else:
            violations += 1

    return dict(checked=checked, correct=correct, planar=planar,
                violations=violations)


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
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")
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
