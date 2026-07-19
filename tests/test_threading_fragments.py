from pathlib import Path

import pytest

from chiralfold.fragments import assemble_protein, predict_secondary_structure
from chiralfold.threading import find_template, thread_sequence


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
