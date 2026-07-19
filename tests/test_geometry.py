import pytest

pytest.importorskip("rdkit", reason="rdkit required for geometry tests")

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, rdMolTransforms  # noqa: E402

from chiralfold.geometry import (  # noqa: E402
    enforce_peptide_planarity,
    find_peptide_bonds,
    measure_planarity_quality,
)
from chiralfold.model import mixed_peptide_smiles  # noqa: E402


def _embedded_peptide(sequence="AFK", chirality="DLD"):
    mol = Chem.MolFromSmiles(mixed_peptide_smiles(sequence, chirality))
    assert mol is not None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 7
    assert AllChem.EmbedMolecule(mol, params) == 0
    return mol


def test_find_peptide_bonds_detects_short_mixed_peptide_backbone():
    mol = _embedded_peptide()

    bonds = find_peptide_bonds(mol)

    assert len(bonds) == 2
    for quad in bonds:
        assert len(quad) == 4
        assert all(isinstance(idx, int) for idx in quad)


def test_measure_planarity_quality_reports_expected_keys():
    mol = _embedded_peptide()

    metrics = measure_planarity_quality(mol)

    assert set(metrics) >= {
        "n_bonds",
        "omegas",
        "deviations",
        "mean_deviation",
        "max_deviation",
        "pct_within_6deg",
        "pct_within_10deg",
    }
    assert metrics["n_bonds"] == 2
    assert len(metrics["omegas"]) == metrics["n_bonds"]
    assert len(metrics["deviations"]) == metrics["n_bonds"]


def test_enforce_peptide_planarity_improves_intentionally_distorted_omega():
    mol = _embedded_peptide()
    ca1, c_atom, n_atom, ca2 = find_peptide_bonds(mol)[0]
    rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), ca1, c_atom, n_atom, ca2, 120.0)
    before = measure_planarity_quality(mol)

    result = enforce_peptide_planarity(mol, conf_id=0, tolerance=6.0, max_iterations=0)
    after = measure_planarity_quality(mol)

    assert result["peptide_bonds_found"] == 2
    assert result["conformers_processed"] == 1
    assert result["total_fixed"] >= 1
    assert after["max_deviation"] <= before["max_deviation"]
    assert after["mean_deviation"] <= 10.0
