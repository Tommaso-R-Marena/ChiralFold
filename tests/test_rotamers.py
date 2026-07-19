import pytest

from chiralfold.rotamers import validate_rotamers


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
