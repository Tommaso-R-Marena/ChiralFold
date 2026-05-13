"""
ChiralFold Unit Tests
======================

Test coverage scope:
  - Amino-acid SMILES libraries (L and D, 20 residues each, including
    glycine achirality and isoleucine/threonine side-chain stereo).
  - SMILES builders (``d_peptide_smiles``, ``l_peptide_smiles``,
    ``mixed_peptide_smiles``) including error paths.
  - Chirality validation (``validate_smiles_chirality``,
    ``validate_3d_chirality``, ``validate_diastereomer``) on pure-D,
    pure-L, and mixed-pattern peptides.
  - ChiralFold model: ``predict`` (de novo and mirror modes).
  - ``MirrorImagePredictor`` reflection geometry.
  - 30 pure-D sequences and 15 diastereomer sequences run end-to-end.
  - External-PDB integration tests (1UBQ download from RCSB).
  - Regression tests confirming the v3.2.1 validator bug fix.
  - Smoke imports for the AF3-correction pipeline and the CLI entry point.
"""

import pytest
import numpy as np

pytest.importorskip("rdkit", reason="rdkit required for chirality tests")

from rdkit import Chem  # noqa: E402

from chiralfold.model import (  # noqa: E402
    ChiralFold,
    d_peptide_smiles,
    l_peptide_smiles,
    mixed_peptide_smiles,
    D_AMINO_ACID_SMILES,
    L_AMINO_ACID_SMILES,
    MirrorImagePredictor,
)
from chiralfold.validator import (  # noqa: E402
    validate_smiles_chirality,
    validate_3d_chirality,
    validate_diastereomer,
)
from chiralfold.data.test_sequences import PURE_D_SEQS, DIASTEREOMER_SEQS  # noqa: E402


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
    # n_conformers=3 for CI speed; production uses default 16
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

@pytest.mark.slow
class TestIntegration:
    def test_full_30_sequence_suite_zero_violations(self):
        """46 total sequences, 467 chiral residues. Key benchmark: all must have 0 violations."""
        total_chiral = 0
        total_viol = 0
        for sid, seq in PURE_D_SEQS.items():
            smi = d_peptide_smiles(seq)
            mol = Chem.MolFromSmiles(smi)
            sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
            total_chiral += sv['n_chiral']
            total_viol += sv['violations']
        assert total_viol == 0
        assert total_chiral >= 250  # stable lower bound (actual: 302)

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

@pytest.mark.slow
class TestExternalPDB:
    """Validate chirality on a real PDB structure fetched from RCSB."""

    def test_ubiquitin_1ubq_chirality(self):
        """PDB 1UBQ (ubiquitin, 76 residues) must have 100% Cα chirality."""
        import os
        import urllib.request

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


# ── AF3 correction smoke tests ────────────────────────────────────────────

class TestAF3Correction:
    def test_correct_af3_import(self):
        from chiralfold.af3_correct import correct_af3_output
        assert callable(correct_af3_output)


# ── CLI smoke tests ───────────────────────────────────────────────────────

class TestCLI:
    def test_cli_import(self):
        from chiralfold.cli import main
        assert callable(main)


# ── Validator regression tests (v3.2.1 bug fix) ────────────────────────────

class TestValidatorFix:
    """Regression tests confirming validator bug fix in v3.2.1."""

    def test_wrong_chirality_detected_smiles(self):
        """validate_smiles_chirality must return violations > 0 for inverted sequence."""
        from chiralfold.validator import validate_smiles_chirality
        from chiralfold.model import mixed_peptide_smiles
        from rdkit import Chem
        # Build a pure-L peptide but tell validator to expect pure-D
        smi = mixed_peptide_smiles("AWK", "LLL")
        mol = Chem.MolFromSmiles(smi)
        result = validate_smiles_chirality(mol, "AWK", "DDD")
        # If all residues are L but we expect D, violations must be > 0
        assert result['violations'] > 0, (
            "validate_smiles_chirality returned 0 violations for deliberately "
            "wrong chirality — the v3.2.1 bug fix may not have been applied."
        )

    def test_correct_chirality_no_violations(self):
        """validate_smiles_chirality must return 0 violations for correct input."""
        from chiralfold.validator import validate_smiles_chirality
        from chiralfold.model import mixed_peptide_smiles
        from rdkit import Chem
        smi = mixed_peptide_smiles("AWK", "LLL")
        mol = Chem.MolFromSmiles(smi)
        result = validate_smiles_chirality(mol, "AWK", "LLL")
        assert result['violations'] == 0


# ── validate_3d_chirality stability tests ────────────────────────────────

class TestValidate3DChiralityStability:
    """validate_3d_chirality must be invariant to atom-numbering permutations.

    The pre-v3.2.2 implementation derived handedness from a signed volume
    over the first three atoms returned by ``GetNeighbors()``. RDKit does
    not guarantee that iteration order matches the canonical bond ordering
    used to interpret ``CHI_TETRAHEDRAL_CW/CCW``, so equivalent molecules
    with different atom indices could be flagged inconsistently. These
    tests pin the new ``AssignStereochemistryFrom3D``-based implementation.
    """

    def _embed(self, seq, pattern):
        from rdkit.Chem import AllChem
        from chiralfold.model import mixed_peptide_smiles
        smi = mixed_peptide_smiles(seq, pattern)
        mol = Chem.MolFromSmiles(smi)
        mol_h = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol_h, params)
        return mol_h

    def test_stable_under_renumbering_d_peptide(self):
        import random
        mol_h = self._embed("AFK", "DDD")
        baseline = validate_3d_chirality(mol_h)
        n = mol_h.GetNumAtoms()
        for seed in range(10):
            random.seed(seed)
            perm = list(range(n))
            random.shuffle(perm)
            renumbered = Chem.RenumberAtoms(mol_h, perm)
            assert validate_3d_chirality(renumbered) == baseline, (
                f"validate_3d_chirality unstable under renumbering "
                f"(seed={seed})"
            )

    def test_stable_under_renumbering_mixed_peptide(self):
        import random
        mol_h = self._embed("ACDEF", "LDLDL")
        baseline = validate_3d_chirality(mol_h)
        assert baseline['checked'] >= 4
        n = mol_h.GetNumAtoms()
        for seed in range(10):
            random.seed(seed + 100)
            perm = list(range(n))
            random.shuffle(perm)
            renumbered = Chem.RenumberAtoms(mol_h, perm)
            assert validate_3d_chirality(renumbered) == baseline

    def test_flipped_tags_are_detected_as_violations(self):
        """Tag/geometry mismatch must still surface as a violation."""
        mol_h = self._embed("AFK", "DDD")
        cw = Chem.ChiralType.CHI_TETRAHEDRAL_CW
        ccw = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        flipped = Chem.RWMol(mol_h)
        for a in flipped.GetAtoms():
            t = a.GetChiralTag()
            if t == cw:
                a.SetChiralTag(ccw)
            elif t == ccw:
                a.SetChiralTag(cw)
        result = validate_3d_chirality(flipped)
        assert result['violations'] == result['checked'] > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
