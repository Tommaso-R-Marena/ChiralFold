import pytest

from chiralfold.interface_scorer import score_interface


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element, record="ATOM"):
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"{record:<6s}{serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _receptor_pdb():
    return "".join(
        [
            _atom_line(1, "N", "ALA", "A", 1, 0.000, 0.000, 0.000, "N"),
            _atom_line(2, "CA", "ALA", "A", 1, 1.458, 0.000, 0.000, "C"),
            _atom_line(3, "C", "ALA", "A", 1, 2.100, 1.200, 0.000, "C"),
            _atom_line(4, "O", "ALA", "A", 1, 1.800, 2.300, 0.000, "O"),
            _atom_line(5, "CB", "ALA", "A", 1, 1.900, -0.800, 1.200, "C"),
            "END\n",
        ]
    )


def _ligand_pdb():
    return "".join(
        [
            _atom_line(1, "N", "DAL", "B", 1, 3.000, 0.000, 0.000, "N", "HETATM"),
            _atom_line(2, "CA", "DAL", "B", 1, 4.458, 0.000, 0.000, "C", "HETATM"),
            _atom_line(3, "C", "DAL", "B", 1, 5.100, 1.200, 0.000, "C", "HETATM"),
            _atom_line(4, "O", "DAL", "B", 1, 4.800, 2.300, 0.000, "O", "HETATM"),
            _atom_line(5, "CB", "DAL", "B", 1, 4.900, -0.800, -1.200, "C", "HETATM"),
            "END\n",
        ]
    )


def test_score_interface_tiny_nearby_files(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.pdb"
    receptor.write_text(_receptor_pdb())
    ligand.write_text(_ligand_pdb())

    metrics = score_interface(str(receptor), str(ligand), receptor_chain="A", ligand_chain="B")

    assert set(metrics) >= {"bsa", "hbonds", "interface_score"}
    assert metrics["bsa"] > 0
    assert metrics["hbonds"] >= 1
    assert metrics["interface_score"] > 0
    if "buried_surface_area" in metrics:
        assert metrics["buried_surface_area"] == metrics["bsa"]
    if "n_hbonds" in metrics:
        assert metrics["n_hbonds"] == metrics["hbonds"]


def test_score_interface_missing_file_raises(tmp_path):
    existing = tmp_path / "receptor.pdb"
    existing.write_text(_receptor_pdb())

    with pytest.raises(FileNotFoundError):
        score_interface(str(existing), str(tmp_path / "missing.pdb"))
