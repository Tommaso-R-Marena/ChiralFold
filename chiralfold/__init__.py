"""
ChiralFold — Chirality-Preserving Peptide Structure Prediction
===============================================================

A model that guarantees 0% chirality violations for D-peptides and
mixed L/D diastereomers, compared to AlphaFold 3's 51% violation rate.

Quick start::

    from chiralfold import ChiralFold, MirrorImagePredictor

    model = ChiralFold()

    # Pure D-peptide (with planarity fix)
    result = model.predict('AFWKELDR')

    # Mixed L/D diastereomer
    result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')

    # Mirror-image PDB transformation (highest quality)
    MirrorImagePredictor.from_pdb_id('1SHG', 'D_SH3.pdb')

    # Mirror-image from coordinates
    import numpy as np
    l_coords = np.random.randn(100, 3)
    result = model.predict_from_mirror(l_coords, 'AEAAAKEAAA')

Reference:
    Childs, Zhou & Donald (2025). "Has AlphaFold 3 Solved the Protein
    Folding Problem for D-Peptides?" bioRxiv 2025.03.14.643307
"""

__version__ = "3.0.0"
__author__ = "Tommaso R. Marena"

from .model import (
    ChiralFold,
    MirrorImagePredictor,
    d_peptide_smiles,
    l_peptide_smiles,
    mixed_peptide_smiles,
    D_AMINO_ACID_SMILES,
    L_AMINO_ACID_SMILES,
)
from .validator import (
    validate_smiles_chirality,
    validate_3d_chirality,
    validate_diastereomer,
)
from .pdb_pipeline import (
    mirror_pdb,
    mirror_pdb_string,
    fetch_and_mirror,
    validate_mirror,
)
from .geometry import (
    enforce_peptide_planarity,
    measure_planarity_quality,
)
from .ramachandran import (
    score_ramachandran,
    filter_conformers_by_ramachandran,
)
from .auditor import (
    audit_pdb,
    format_report,
)
from .rotamers import validate_rotamers
from .threading import thread_sequence, thread_and_mirror, find_template
from .fragments import assemble_protein, predict_secondary_structure
from .af3_correct import correct_af3_output, detect_chirality_violations, correct_chirality
from .interface_scorer import score_interface, compare_interfaces
from .enumerate import enumerate_diastereomers

__all__ = [
    'ChiralFold',
    'MirrorImagePredictor',
    'd_peptide_smiles',
    'l_peptide_smiles',
    'mixed_peptide_smiles',
    'D_AMINO_ACID_SMILES',
    'L_AMINO_ACID_SMILES',
    'validate_smiles_chirality',
    'validate_3d_chirality',
    'validate_diastereomer',
    'mirror_pdb',
    'mirror_pdb_string',
    'fetch_and_mirror',
    'validate_mirror',
    'enforce_peptide_planarity',
    'measure_planarity_quality',
    'score_ramachandran',
    'filter_conformers_by_ramachandran',
    'audit_pdb',
    'format_report',
    'validate_rotamers',
    'thread_sequence',
    'thread_and_mirror',
    'find_template',
    'assemble_protein',
    'predict_secondary_structure',
    'correct_af3_output',
    'detect_chirality_violations',
    'correct_chirality',
    'score_interface',
    'compare_interfaces',
    'enumerate_diastereomers',
]
