"""
ChiralFold Unit Tests
======================

Tests for chirality correctness across pure D, pure L, and mixed L/D peptides.
"""

import pytest
import numpy as np
from rdkit import Chem

from chiralfold.model import (
    ChiralFold,
    d_peptide_smiles,
    l_peptide_smiles,
    mixed_peptide_smiles,
    D_AMINO_ACID_SMILES,
    L_AMINO_ACID_SMILES,
    MirrorImagePredictor,
)
from chiralfold.validator import (
    validate_smiles_chirality,
    validate_3d_chirality,
    validate_diastereomer,
)
from chiralfold.data.test_sequences import PURE_D_SEQS, DIASTEREOMER_SEQS


# ── Amino acid library tests ──────────────────────────────────────────────

class TestAminoAcidLibrary:
    def test_d_library_has_20_amino_acids(self):
        assert len(D_AMINO_ACID_SMILES) == 20

    def test_l_library_has_20_amino_acids(self):
        assert len(L_AMINO_ACID_SMILES) == 20

    def test_d_amino_acids_are_valid_smiles(self):
        for aa, smi in D_AMINO_ACID_SMILES.items():
            mol = Chem.MolFromSmiles(smi)
            assert mol is not None, f"Invalid SMILES for D-{aa}: {smi}"

    def test_l_amino_acids_are_valid_smiles(self):
        for aa, smi in L_AMINO_ACID_SMILES.items():
            mol = Chem.MolFromSmiles(smi)
            assert mol is not None, f"Invalid SMILES for L-{aa}: {smi}"

    def test_glycine_is_achiral(self):
        mol = Chem.MolFromSmiles(D_AMINO_ACID_SMILES['G'])
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        cc = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        assert len(cc) == 0

    def test_d_amino_acids_have_stereocenters(self):
        for aa, smi in D_AMINO_ACID_SMILES.items():
            if aa == 'G':
                continue
            mol = Chem.MolFromSmiles(smi)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            cc = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            assert len(cc) >= 1, f"D-{aa} has no stereocenters"


# ── SMILES builder tests ─────────────────────────────────────────────────

class TestSMILESBuilders:
    def test_d_peptide_smiles_short(self):
        smi = d_peptide_smiles('AFK')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_l_peptide_smiles_short(self):
        smi = l_peptide_smiles('AFK')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_mixed_peptide_smiles(self):
        smi = mixed_peptide_smiles('AFK', 'DLD')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_mixed_raises_on_length_mismatch(self):
        with pytest.raises(ValueError, match="length"):
            mixed_peptide_smiles('AFK', 'DL')

    def test_mixed_raises_on_invalid_chirality(self):
        with pytest.raises(ValueError, match="Invalid chirality"):
            mixed_peptide_smiles('AFK', 'DLX')

    def test_mixed_raises_on_unknown_amino_acid(self):
        with pytest.raises(ValueError, match="Unknown amino acid"):
            mixed_peptide_smiles('AZK', 'DDD')

    def test_d_equals_all_d_mixed(self):
        seq = 'AFWKELDR'
        smi_d = d_peptide_smiles(seq)
        smi_m = mixed_peptide_smiles(seq, 'D' * len(seq))
        assert smi_d == smi_m

    def test_l_equals_all_l_mixed(self):
        seq = 'AFWKELDR'
        smi_l = l_peptide_smiles(seq)
        smi_m = mixed_peptide_smiles(seq, 'L' * len(seq))
        assert smi_l == smi_m

    def test_glycine_in_mixed_peptide(self):
        smi = mixed_peptide_smiles('GAG', 'DDD')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_proline_in_mixed_peptide(self):
        smi = mixed_peptide_smiles('APK', 'DDD')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_all_20_amino_acids_in_d_peptide(self):
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        smi = d_peptide_smiles(seq)
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None


# ── Chirality validation tests ────────────────────────────────────────────

class TestChiralityValidation:
    def test_pure_d_zero_violations(self):
        for sid, seq in PURE_D_SEQS.items():
            smi = d_peptide_smiles(seq)
            mol = Chem.MolFromSmiles(smi)
            sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
            assert sv['violations'] == 0, f"{sid}: {sv['violations']} violations"

    def test_diastereomer_zero_violations(self):
        for sid, data in DIASTEREOMER_SEQS.items():
            report = validate_diastereomer(data['seq'], data['chirality'])
            assert report['smiles_violations'] == 0, (
                f"{sid}: {report['smiles_violations']} violations"
            )
            assert report['valid'], f"{sid}: validation failed"

    def test_none_mol_returns_error(self):
        sv = validate_smiles_chirality(None, 'AFK', 'DDD')
        assert sv['error'] is True


# ── ChiralFold model tests ───────────────────────────────────────────────

class TestChiralFoldModel:
    def test_predict_pure_d(self):
        model = ChiralFold(n_conformers=3)
        result = model.predict('AFK')
        assert result['chirality_violations'] == 0
        assert result['violation_rate'] == 0.0
        assert result['n_d_residues'] == 3

    def test_predict_mixed(self):
        model = ChiralFold(n_conformers=3)
        result = model.predict('AFK', chirality_pattern='DLD')
        assert result['chirality_violations'] == 0
        assert result['n_d_residues'] == 2
        assert result['n_l_residues'] == 1

    def test_predict_all_l(self):
        model = ChiralFold(n_conformers=3)
        result = model.predict('AFK', chirality_pattern='LLL')
        assert result['chirality_violations'] == 0
        assert result['n_l_residues'] == 3

    def test_predict_defaults_to_all_d(self):
        model = ChiralFold(n_conformers=3)
        result = model.predict('AFK')
        assert result['chirality_pattern'] == 'DDD'

    def test_predict_mirror(self):
        model = ChiralFold()
        coords = np.random.randn(50, 3)
        result = model.predict_from_mirror(coords, 'AEAAA')
        assert result['chirality_preserved'] is True
        assert result['rmsd_to_ideal_mirror'] == 0.0


# ── Mirror-image predictor tests ─────────────────────────────────────────

class TestMirrorImagePredictor:
    def test_reflect_inverts_x(self):
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        reflected = MirrorImagePredictor.reflect_structure(coords, axis='x')
        np.testing.assert_array_equal(reflected[:, 0], -coords[:, 0])
        np.testing.assert_array_equal(reflected[:, 1:], coords[:, 1:])

    def test_reflect_inverts_y(self):
        coords = np.array([[1.0, 2.0, 3.0]])
        reflected = MirrorImagePredictor.reflect_structure(coords, axis='y')
        assert reflected[0, 1] == -2.0

    def test_verify_mirror_rmsd_zero(self):
        coords = np.random.randn(100, 3)
        reflected = MirrorImagePredictor.reflect_structure(coords, axis='x')
        v = MirrorImagePredictor.verify_mirror_chirality(coords, reflected)
        assert v['rmsd_to_expected'] == pytest.approx(0.0, abs=1e-10)
        assert v['chirality_inverted'] is True

    def test_predict_d_structure(self):
        coords = np.random.randn(50, 3)
        result = MirrorImagePredictor.predict_d_structure(coords)
        assert result['n_atoms'] == 50
        assert result['method'] == 'mirror_image_reflection'


# ── Integration test ──────────────────────────────────────────────────────

class TestIntegration:
    def test_full_30_sequence_suite_zero_violations(self):
        """The key benchmark: all 30 pure D sequences must have 0 violations."""
        total_chiral = 0
        total_viol = 0
        for sid, seq in PURE_D_SEQS.items():
            smi = d_peptide_smiles(seq)
            mol = Chem.MolFromSmiles(smi)
            sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
            total_chiral += sv['n_chiral']
            total_viol += sv['violations']
        assert total_viol == 0
        assert total_chiral > 250  # Should be 302

    def test_full_15_diastereomer_suite_zero_violations(self):
        """All 15 diastereomer sequences must have 0 violations."""
        total_chiral = 0
        total_viol = 0
        for sid, data in DIASTEREOMER_SEQS.items():
            report = validate_diastereomer(data['seq'], data['chirality'])
            total_chiral += report['n_chiral']
            total_viol += report['smiles_violations']
        assert total_viol == 0
        assert total_chiral > 100


# ── External PDB validation test ──────────────────────────────────────

class TestExternalPDB:
    """Validate chirality on a real PDB structure fetched from RCSB."""

    def test_ubiquitin_1ubq_chirality(self):
        """PDB 1UBQ (ubiquitin, 76 residues) must have 100% Cα chirality."""
        import os
        import urllib.request
        import tempfile

        # Download 1UBQ if not cached
        cache_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
        pdb_path = os.path.join(cache_dir, '1UBQ.pdb')

        if not os.path.exists(pdb_path):
            os.makedirs(cache_dir, exist_ok=True)
            url = 'https://files.rcsb.org/download/1UBQ.pdb'
            try:
                urllib.request.urlretrieve(url, pdb_path)
            except Exception:
                pytest.skip('Could not download 1UBQ.pdb from RCSB')

        if not os.path.exists(pdb_path):
            pytest.skip('1UBQ.pdb not available')

        from chiralfold.auditor import audit_pdb
        report = audit_pdb(pdb_path)

        assert report['chirality']['pct_correct'] == 100.0, (
            f"1UBQ chirality: {report['chirality']['pct_correct']}% "
            f"(expected 100%)"
        )
        assert report['n_residues'] >= 70, (
            f"1UBQ has {report['n_residues']} residues (expected >= 70)"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
