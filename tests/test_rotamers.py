import pytest
import numpy as np

from chiralfold import rotamers as rotamer_mod
from chiralfold.rotamers import (
    _angle_deviation,
    _classify_chi1,
    _group_by_residue,
    _parse_pdb_atoms,
    compute_chi1,
    validate_rotamers,
)


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element):
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"ATOM  {serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _ala_leu_pdb():
    return "".join(
        [
            _atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
            _atom_line(2, "CA", "ALA", "A", 1, 1.458, 0.000, 0.000, "C"),
            _atom_line(3, "CB", "ALA", "A", 1, 1.936, 1.440, 0.000, "C"),
            _atom_line(4, "N", "LEU", "A", 2, 3.000, 0.000, 0.000, "N"),
            _atom_line(5, "CA", "LEU", "A", 2, 4.458, 0.000, 0.000, "C"),
            _atom_line(6, "CB", "LEU", "A", 2, 4.936, 1.440, 0.000, "C"),
            _atom_line(7, "CG", "LEU", "A", 2, 4.200, 2.500, 0.900, "C"),
            "END\n",
        ]
    )


def test_validate_rotamers_minimal_ala_leu_pdb(tmp_path):
    pdb_path = tmp_path / "ala_leu.pdb"
    pdb_path.write_text(_ala_leu_pdb())

    report = validate_rotamers(str(pdb_path), chain="A")

    assert set(report) == {
        "n_residues_checked",
        "n_favored",
        "n_allowed",
        "n_outlier",
        "pct_favored",
        "pct_outlier",
        "outliers",
    }
    assert report["n_residues_checked"] == 1
    assert report["n_favored"] + report["n_allowed"] + report["n_outlier"] == 1


def test_validate_rotamers_raises_for_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_rotamers(str(tmp_path / "missing.pdb"))


def test_compute_chi1_and_classification_helpers():
    n = np.array([0.000, 0.000, 0.000])
    ca = np.array([1.458, 0.000, 0.000])
    cb = np.array([1.936, 1.440, 0.000])
    cg = np.array([1.200, 2.500, 0.900])

    chi = compute_chi1(n, ca, cb, cg)

    assert -180.0 <= chi <= 180.0
    assert _angle_deviation(175.0, -180.0) == 5.0
    assert _classify_chi1(-55.0, [-60.0, 180.0]) == "favored"
    assert _classify_chi1(-105.0, [-60.0, 180.0]) == "allowed"
    assert _classify_chi1(20.0, [-60.0, 180.0]) == "outlier"


def test_compute_chi1_rejects_degenerate_geometry():
    with pytest.raises(ValueError, match="Degenerate dihedral"):
        compute_chi1(
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([3.0, 0.0, 0.0]),
        )


def test_parse_and_group_rotamer_atoms_skip_bad_records(tmp_path):
    pdb = tmp_path / "bad_records.pdb"
    pdb.write_text(
        "".join(
            [
                _atom_line(1, "N", "LEU", "A", 1, 0.0, 0.0, 0.0, "N"),
                "ATOM  BADLY FORMATTED\n",
                "REMARK ignored\n",
            ]
        )
    )

    atoms = _parse_pdb_atoms(str(pdb))
    grouped = _group_by_residue(atoms)

    assert len(atoms) == 1
    assert list(grouped) == [("A", 1, "", "LEU")]


def test_validate_rotamers_chain_filter_and_skip_paths(tmp_path, monkeypatch):
    pdb = tmp_path / "skip_paths.pdb"
    pdb.write_text(
        "".join(
            [
                _atom_line(1, "N", "LEU", "B", 1, 0.0, 0.0, 0.0, "N"),
                _atom_line(2, "CA", "LEU", "B", 1, 1.0, 0.0, 0.0, "C"),
                _atom_line(3, "CB", "LEU", "B", 1, 1.0, 1.0, 0.0, "C"),
                _atom_line(4, "CG", "LEU", "B", 1, 1.0, 1.0, 1.0, "C"),
                _atom_line(5, "N", "GLY", "A", 2, 0.0, 0.0, 0.0, "N"),
                _atom_line(6, "N", "LEU", "A", 3, 0.0, 0.0, 0.0, "N"),
                _atom_line(7, "CA", "LEU", "A", 3, 1.0, 0.0, 0.0, "C"),
                _atom_line(8, "CB", "LEU", "A", 3, 1.0, 1.0, 0.0, "C"),
                "END\n",
            ]
        )
    )

    assert validate_rotamers(str(pdb), chain="C")["n_residues_checked"] == 0

    monkeypatch.setitem(rotamer_mod._CG_ATOM_NAME, "LEU", None)
    assert validate_rotamers(str(pdb), chain="A")["n_residues_checked"] == 0


def test_validate_rotamers_handles_compute_error_and_outlier(tmp_path, monkeypatch):
    pdb = tmp_path / "leu.pdb"
    pdb.write_text(_ala_leu_pdb())

    monkeypatch.setattr(rotamer_mod, "compute_chi1", lambda *args: (_ for _ in ()).throw(ValueError))
    assert validate_rotamers(str(pdb))["n_residues_checked"] == 0

    monkeypatch.setattr(rotamer_mod, "compute_chi1", lambda *args: 20.0)
    report = validate_rotamers(str(pdb))

    assert report["n_outlier"] == 1
    assert report["pct_outlier"] == 100.0
    assert report["outliers"][0]["resname"] == "LEU"
