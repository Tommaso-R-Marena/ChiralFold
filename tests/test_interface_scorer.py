import pytest

from chiralfold import interface_scorer as iface_mod
from chiralfold.interface_scorer import compare_interfaces, format_comparison_table, score_interface


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


def test_score_interface_chain_filters_and_empty_chain_errors(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.pdb"
    receptor.write_text(_receptor_pdb())
    ligand.write_text(_ligand_pdb())

    with pytest.raises(ValueError, match="No atoms found in receptor"):
        score_interface(str(receptor), str(ligand), receptor_chain="Z", ligand_chain="B")

    with pytest.raises(ValueError, match="No atoms found in ligand"):
        score_interface(str(receptor), str(ligand), receptor_chain="A", ligand_chain="Z")


def test_compare_interfaces_sorts_successes_and_formats_errors(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.pdb"
    missing = tmp_path / "missing.pdb"
    receptor.write_text(_receptor_pdb())
    ligand.write_text(_ligand_pdb())

    results = compare_interfaces(
        [
            (receptor, ligand, "near"),
            (receptor, missing, "missing"),
        ]
    )
    table = format_comparison_table(results)

    assert results[0]["label"] == "near"
    assert results[0]["interface_score"] > 0
    assert results[1]["label"] == "missing"
    assert "error" in results[1]
    assert "near" in table
    assert "ERROR:" in table
    assert format_comparison_table([]) == "(no results)"


def test_interface_parser_filters_and_infers_chains(tmp_path):
    pdb = tmp_path / "parse_edges.pdb"
    pdb.write_text(
        """\
ATOM      1  N   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  CA AALA A   2       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA BALA A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      4  CA AALA A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM  broken
END
"""
    )

    atoms = iface_mod._parse_atoms(str(pdb), chains=["A"])

    assert len(atoms) == 1
    assert atoms[0].element_sym == "C"
    assert iface_mod._infer_chains(str(pdb)) == ["A"]
    assert iface_mod._parse_atoms(str(pdb), chains=["Z"]) == []


def test_interface_metric_helpers_cover_empty_and_contact_categories():
    rec = iface_mod._Atom("ATOM", 1, "N", "LYS", "A", 1, 0.0, 0.0, 0.0, "N")
    lig_clash = iface_mod._Atom("ATOM", 2, "O", "ASP", "B", 1, 3.0, 0.0, 0.0, "O")
    lig_comp = iface_mod._Atom("ATOM", 3, "CA", "ALA", "B", 2, 4.0, 0.0, 0.0, "C")
    lig_far = iface_mod._Atom("ATOM", 4, "CB", "ALA", "B", 3, 6.0, 0.0, 0.0, "C")

    assert iface_mod._find_interface_pairs([], [lig_clash]) == []
    assert iface_mod._compute_shape_complementarity([]) == {
        "n_clashing": 0,
        "n_complementary": 0,
        "n_no_contact": 0,
        "sc_fraction": 0.0,
    }

    pairs = [(rec, lig_clash, 3.0), (rec, lig_comp, 4.0), (rec, lig_far, 6.0)]
    sc = iface_mod._compute_shape_complementarity(pairs)

    assert iface_mod._dist(rec.xyz, lig_comp.xyz) == 4.0
    assert iface_mod._compute_bsa(pairs) == 40.0
    assert sc["n_clashing"] == 1
    assert sc["n_complementary"] == 1
    assert sc["n_no_contact"] == 1
    assert iface_mod._compute_hbonds(pairs) == 1


def test_interface_residue_specific_metrics(tmp_path):
    receptor = tmp_path / "charged_hydrophobic_receptor.pdb"
    ligand = tmp_path / "charged_hydrophobic_ligand.pdb"
    receptor.write_text(
        "".join(
            [
                _atom_line(1, "CA", "LYS", "A", 1, 0.0, 0.0, 0.0, "C"),
                _atom_line(2, "CA", "ALA", "A", 2, 8.0, 0.0, 0.0, "C"),
                "END\n",
            ]
        )
    )
    ligand.write_text(
        "".join(
            [
                _atom_line(1, "CA", "ASP", "B", 1, 3.5, 0.0, 0.0, "C"),
                _atom_line(2, "CA", "VAL", "B", 2, 12.0, 0.0, 0.0, "C"),
                "END\n",
            ]
        )
    )
    rec_atoms = iface_mod._parse_atoms(str(receptor))
    lig_atoms = iface_mod._parse_atoms(str(ligand))

    assert iface_mod._compute_salt_bridges(rec_atoms, lig_atoms) == 1
    assert iface_mod._compute_hydrophobic(rec_atoms, lig_atoms) == 1
    assert iface_mod._ca_pairs([], lig_atoms, 4.0) == []
