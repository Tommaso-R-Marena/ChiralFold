from pathlib import Path

import numpy as np
import pytest

from chiralfold import fragments as fragment_mod
from chiralfold.fragments import (
    _place_atom,
    assemble_protein,
    build_backbone_from_fragments,
    predict_secondary_structure,
)
from chiralfold.threading import find_template, thread_and_mirror, thread_sequence


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element):
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"ATOM  {serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _template_backbone_pdb():
    coords = [
        (1, "N", 0.000, 0.000, 0.000),
        (1, "CA", 1.458, 0.000, 0.000),
        (1, "C", 2.130, 1.200, 0.000),
        (1, "O", 1.600, 2.200, 0.000),
        (2, "N", 3.460, 1.200, 0.000),
        (2, "CA", 4.200, 2.400, 0.000),
        (2, "C", 5.700, 2.400, 0.000),
        (2, "O", 6.300, 3.400, 0.000),
        (3, "N", 6.400, 1.300, 0.000),
        (3, "CA", 7.850, 1.100, 0.000),
        (3, "C", 8.500, 2.300, 0.000),
        (3, "O", 8.000, 3.300, 0.000),
    ]
    lines = [
        _atom_line(serial, atom, "ALA", "A", resseq, x, y, z, atom[0])
        for serial, (resseq, atom, x, y, z) in enumerate(coords, start=1)
    ]
    return "".join(lines + ["END\n"])


def test_predict_secondary_structure_short_sequence():
    ss = predict_secondary_structure("AVGP")

    assert ss == "HECC"


def test_assemble_protein_short_sequence_writes_backbone(tmp_path):
    output = tmp_path / "assembled.pdb"

    result = assemble_protein("AFK", output_pdb=str(output), seed=42)

    assert result["n_residues"] == 3
    assert result["n_atoms"] == 12
    assert result["ss_prediction"] == predict_secondary_structure("AFK")
    assert output.exists()
    assert output.read_text().count("\nATOM") >= 11


def test_thread_sequence_and_find_template_with_synthetic_template(tmp_path):
    template = tmp_path / "template.pdb"
    output = tmp_path / "threaded.pdb"
    template.write_text(_template_backbone_pdb())

    threaded = thread_sequence("AFK", str(template), "A", str(output), chirality="L")
    found = find_template("AAA", str(tmp_path), chain="A")

    assert threaded["n_residues"] == 3
    assert threaded["n_atoms_written"] == 12
    assert output.exists()
    text = output.read_text()
    assert "ALA" in text
    assert "PHE" in text
    assert "LYS" in text
    assert Path(found["template_pdb"]) == template
    assert found["template_seq"] == "AAA"


def test_thread_sequence_rejects_missing_template(tmp_path):
    with pytest.raises(FileNotFoundError):
        thread_sequence("AFK", str(tmp_path / "missing.pdb"), "A", str(tmp_path / "out.pdb"))


def test_thread_sequence_truncates_and_writes_d_residues(tmp_path):
    template = tmp_path / "template.pdb"
    output = tmp_path / "threaded_d.pdb"
    template.write_text(_template_backbone_pdb())

    with pytest.warns(UserWarning, match="Truncating to 3 residues"):
        result = thread_sequence("AFKW", str(template), "A", str(output), chirality="D")

    text = output.read_text()
    assert result["target_seq"] == "AFK"
    assert result["n_residues"] == 3
    assert "DAL" in text
    assert "DPN" in text
    assert "DLY" in text


def test_thread_sequence_rejects_bad_chirality_and_missing_chain(tmp_path):
    template = tmp_path / "template.pdb"
    template.write_text(_template_backbone_pdb())

    with pytest.raises(ValueError, match="chirality must be"):
        thread_sequence("AFK", str(template), "A", str(tmp_path / "bad.pdb"), chirality="X")

    with pytest.raises(ValueError, match="No backbone atoms"):
        thread_sequence("AFK", str(template), "Z", str(tmp_path / "missing_chain.pdb"))


def test_thread_and_mirror_writes_d_mirror_and_removes_intermediate(tmp_path):
    template = tmp_path / "template.pdb"
    output = tmp_path / "mirrored.pdb"
    template.write_text(_template_backbone_pdb())

    returned = thread_and_mirror("AFK", str(template), "A", str(output))

    text = output.read_text()
    assert returned == str(output)
    assert "DAL" in text
    assert "DPN" in text
    assert not Path(str(output) + ".L_intermediate.pdb").exists()


def test_find_template_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_template("AFK", str(tmp_path / "missing"), chain="A")

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="No .pdb files"):
        find_template("AFK", str(empty), chain="A")

    bad = tmp_path / "bad"
    bad.mkdir()
    (bad / "only_chain_b.pdb").write_text(_template_backbone_pdb().replace(" A", " B"))
    with pytest.raises(ValueError, match="Could not read any valid chain"):
        find_template("AFK", str(bad), chain="A")


def test_place_atom_and_build_backbone_edge_cases():
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([1.458, 0.0, 0.0])
    p3 = np.array([1.936, 1.440, 0.0])

    placed = _place_atom(p1, p2, p3, bond_length=1.5, bond_angle=109.5, dihedral=180.0)
    assert placed.shape == (3,)
    assert np.linalg.norm(placed - p3) == pytest.approx(1.5, abs=0.01)

    residues = build_backbone_from_fragments("AV", seed=1)
    assert len(residues) == 2
    assert set(residues[0]) == {"N", "CA", "C", "O"}

    with pytest.raises(ValueError, match="ss_prediction length"):
        build_backbone_from_fragments("AV", ss_prediction="H")


def test_assemble_protein_d_default_output_and_invalid_inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = assemble_protein("AF", chirality="D", seed=2)

    output = tmp_path / "fragment_2aa.pdb"
    text = output.read_text()
    assert result["output_pdb"] == "fragment_2aa.pdb"
    assert "DAL" in text
    assert "DPN" in text

    with pytest.raises(ValueError, match="non-empty"):
        assemble_protein("")

    with pytest.raises(ValueError, match="chirality must be"):
        assemble_protein("A", chirality="X")


def test_assemble_protein_skips_missing_backbone_atom(tmp_path, monkeypatch):
    output = tmp_path / "partial.pdb"

    def fake_backbone(sequence, ss_prediction, seed=None):
        return [
            {
                "N": np.array([0.0, 0.0, 0.0]),
                "CA": None,
                "C": np.array([1.0, 1.0, 0.0]),
                "O": np.array([1.0, 2.0, 0.0]),
            }
        ]

    monkeypatch.setattr(fragment_mod, "build_backbone_from_fragments", fake_backbone)

    result = assemble_protein("A", output_pdb=str(output))

    assert result["n_atoms"] == 3
    assert " CA " not in output.read_text()
