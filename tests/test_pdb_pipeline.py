import numpy as np
import pytest

from chiralfold.pdb_pipeline import mirror_pdb, mirror_pdb_string, validate_mirror


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element, record="ATOM"):
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"{record:<6s}{serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _tiny_pdb():
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 1.000, 0.000, 0.000, "N"),
        _atom_line(2, "CA", "ALA", "A", 1, 2.000, 0.500, 0.000, "C"),
        _atom_line(3, "C", "ALA", "A", 1, 2.900, 1.600, 0.000, "C"),
        _atom_line(4, "O", "ALA", "A", 1, 2.600, 2.700, 0.000, "O"),
        _atom_line(5, "CB", "ALA", "A", 1, 2.100, -0.300, 1.200, "C"),
        _atom_line(6, "N", "PHE", "A", 2, 4.200, 1.300, 0.100, "N"),
        _atom_line(7, "CA", "PHE", "A", 2, 5.100, 2.300, 0.200, "C"),
        _atom_line(8, "C", "PHE", "A", 2, 6.400, 1.800, 0.100, "C"),
        _atom_line(9, "O", "PHE", "A", 2, 7.200, 2.500, 0.000, "O"),
        _atom_line(10, "CB", "PHE", "A", 2, 4.600, 3.100, 1.300, "C"),
        "END\n",
    ]
    return "".join(lines)


def _coords_by_serial(pdb_text):
    coords = {}
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords[int(line[6:11])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return coords


def test_mirror_pdb_string_reflects_x_and_renames_residues():
    source = _tiny_pdb()

    mirrored = mirror_pdb_string(source, axis="x")

    src = _coords_by_serial(source)
    dst = _coords_by_serial(mirrored)
    assert dst.keys() == src.keys()
    for serial, xyz in src.items():
        np.testing.assert_allclose(dst[serial], [-xyz[0], xyz[1], xyz[2]], atol=1e-6)
    assert "DAL" in mirrored
    assert "DPN" in mirrored


def test_mirror_pdb_file_validate_and_round_trip(tmp_path):
    source = tmp_path / "source.pdb"
    mirrored = tmp_path / "mirrored.pdb"
    roundtrip = tmp_path / "roundtrip.pdb"
    source.write_text(_tiny_pdb())

    result = mirror_pdb(str(source), str(mirrored), axis="x")
    validation = validate_mirror(str(source), str(mirrored), chain="A")
    mirror_pdb(str(mirrored), str(roundtrip), axis="x")

    assert result["n_atoms"] == 10
    assert result["n_residues"] == 2
    assert bool(validation["reflection_exact"]) is True
    assert validation["coord_max_error"] == pytest.approx(0.0, abs=1e-6)

    src = _coords_by_serial(source.read_text())
    rt = _coords_by_serial(roundtrip.read_text())
    for serial, xyz in src.items():
        np.testing.assert_allclose(rt[serial], xyz, atol=1e-6)
