import numpy as np
import pytest

pytest.importorskip("rdkit", reason="rdkit required for conformer filtering smoke test")

from rdkit import Chem  # noqa: E402

from chiralfold.ramachandran import (  # noqa: E402
    filter_conformers_by_ramachandran,
    score_ramachandran,
)


def test_score_ramachandran_classifies_known_l_regions():
    assert score_ramachandran(-60.0, -45.0, "general") == "favored"
    assert score_ramachandran(-120.0, 130.0, "general") in {"favored", "allowed"}
    assert score_ramachandran(0.0, -150.0, "general") == "outlier"


def test_score_ramachandran_classifies_d_mirror_region():
    assert score_ramachandran(60.0, 45.0, "D-general") == "favored"
    assert score_ramachandran(-60.0, -45.0, "D-general") == "outlier"


def test_score_ramachandran_handles_special_residue_types_and_nan():
    assert score_ramachandran(60.0, -120.0, "glycine") == "favored"
    assert score_ramachandran(-65.0, -35.0, "proline") == "favored"
    assert score_ramachandran(float("nan"), -45.0, "general") == "outlier"


def test_filter_conformers_by_ramachandran_smoke_with_minimal_conformer():
    rw = Chem.RWMol()
    for _ in range(6):
        rw.AddAtom(Chem.Atom("C"))
    for idx in range(5):
        rw.AddBond(idx, idx + 1, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(6)
    for idx, xyz in enumerate(
        [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.5, 1.0, 0.0),
            (2.5, 1.0, 0.5),
            (3.0, 2.0, 0.3),
            (4.0, 2.2, 1.0),
        ]
    ):
        conf.SetAtomPosition(idx, xyz)
    mol.AddConformer(conf)
    residues = [
        {"N": 0, "CA": 1, "C": 2},
        {"N": 3, "CA": 4, "C": 5},
    ]

    passing = filter_conformers_by_ramachandran(mol, residues, [0, 99], min_pct_favored=0.0)

    assert passing == [(0, np.float64(0.0))] or passing == [(0, 0.0)]
