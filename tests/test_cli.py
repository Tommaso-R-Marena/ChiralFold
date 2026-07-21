import json
import sys
from pathlib import Path

import pytest

from chiralfold import __version__
import chiralfold.cli as cli
from chiralfold.cli import main


def _run_cli(monkeypatch, capsys, argv):
    monkeypatch.setattr(sys, "argv", ["chiralfold", *argv])
    main()
    return capsys.readouterr()


def test_cli_version(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["chiralfold", "--version"])

    with pytest.raises(SystemExit) as exc:
        main()

    assert exc.value.code == 0
    assert f"chiralfold {__version__}" in capsys.readouterr().out


def test_cli_help(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["chiralfold", "--help"])

    with pytest.raises(SystemExit) as exc:
        main()

    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "ChiralFold" in out
    assert "predict" in out
    assert "validate" in out


def test_cli_predict_afk_json(monkeypatch, capsys, tmp_path):
    _ = tmp_path

    captured = _run_cli(monkeypatch, capsys, ["predict", "AFK", "--json"])
    payload = json.loads(captured.out)

    assert payload["sequence"] == "AFK"
    assert payload["chirality_pattern"] == "DDD"
    assert payload["chirality_violations"] == 0


def test_cli_validate_afk_json(monkeypatch, capsys, tmp_path):
    _ = tmp_path

    captured = _run_cli(monkeypatch, capsys, ["validate", "AFK", "--chirality", "DLD", "--json"])
    payload = json.loads(captured.out)

    assert payload["valid"] is True
    assert payload["n_residues"] == 3
    assert payload["n_d"] == 2
    assert payload["n_l"] == 1


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element, record="ATOM"):
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"{record:<6s}{serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _tiny_ala_pdb(record="ATOM", resname="ALA", chain="A", x_offset=0.0):
    return "".join(
        [
            _atom_line(1, "N", resname, chain, 1, x_offset + 1.201, 0.847, 0.000, "N", record),
            _atom_line(2, "CA", resname, chain, 1, x_offset + 0.000, 0.000, 0.000, "C", record),
            _atom_line(3, "C", resname, chain, 1, x_offset - 1.250, 0.881, 0.000, "C", record),
            _atom_line(4, "O", resname, chain, 1, x_offset - 1.200, 2.095, 0.000, "O", record),
            _atom_line(5, "CB", resname, chain, 1, x_offset + 0.000, -0.500, 1.200, "C", record),
            "END\n",
        ]
    )


def test_cli_without_command_prints_help(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["chiralfold"])

    with pytest.raises(SystemExit) as exc:
        main()

    assert exc.value.code == 0
    assert "Available commands" in capsys.readouterr().out


def test_cli_predict_text_and_error(monkeypatch, capsys):
    class DummyModel:
        def __init__(self, n_conformers):
            assert n_conformers == 10

        def predict(self, sequence, chirality_pattern):
            assert (sequence, chirality_pattern) == ("AFK", "DLD")
            return {
                "sequence": sequence,
                "chirality_pattern": chirality_pattern,
                "n_d_residues": 2,
                "n_l_residues": 1,
                "smiles": "N[C@@H](C)C(=O)O",
                "chirality_violations": 0,
                "violation_rate": 0.0,
                "conformers": [object(), object()],
            }

    monkeypatch.setattr(cli, "ChiralFold", DummyModel)
    captured = _run_cli(monkeypatch, capsys, ["predict", "AFK", "--chirality", "DLD"])
    assert "Sequence:          AFK" in captured.out
    assert "3D conformers:     2" in captured.out

    class ErrorModel:
        def __init__(self, n_conformers):
            pass

        def predict(self, sequence, chirality_pattern):
            return {"error": "bad input"}

    monkeypatch.setattr(cli, "ChiralFold", ErrorModel)
    monkeypatch.setattr(sys, "argv", ["chiralfold", "predict", "AFK"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 1
    assert "Error: bad input" in capsys.readouterr().err


def test_cli_validate_text(monkeypatch, capsys):
    captured = _run_cli(monkeypatch, capsys, ["validate", "AFK", "--chirality", "DLD"])

    assert "[PASS] AFK (DLD)" in captured.out
    assert "Residues: 3" in captured.out


def test_cli_audit_json_text_and_missing_arg(monkeypatch, capsys, tmp_path):
    pdb = tmp_path / "audit.pdb"
    pdb.write_text(_tiny_ala_pdb())

    captured = _run_cli(monkeypatch, capsys, ["audit", str(pdb), "--json"])
    payload = json.loads(captured.out)
    assert payload["n_residues"] == 1
    assert payload["chirality"]["n_correct"] == 1

    captured = _run_cli(monkeypatch, capsys, ["audit", str(pdb)])
    assert "ChiralFold PDB Auditor" in captured.out

    monkeypatch.setattr(sys, "argv", ["chiralfold", "audit"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 1
    assert "provide a structure file" in capsys.readouterr().err


def test_cli_audit_rcsb_batch_writes_success_and_error_rows(monkeypatch, capsys, tmp_path):
    import chiralfold.auditor as auditor

    batch = tmp_path / "ids.txt"
    output = tmp_path / "audit.csv"
    batch.write_text("1ABC\nBADID\n2XYZ\n")

    def fake_urlretrieve(url, filename):
        if "2XYZ" in url:
            raise OSError("network unavailable")
        Path(filename).write_text(_tiny_ala_pdb())

    def fake_audit_pdb(path):
        return {
            "chirality": {"pct_correct": 100.0},
            "ramachandran": {"pct_favored": 0.0, "pct_outlier": 0.0},
            "planarity": {"pct_within_6deg": 100.0},
            "clashes": {"clash_score": 0.0},
            "overall_score": 42.0,
        }

    monkeypatch.setattr("urllib.request.urlretrieve", fake_urlretrieve)
    monkeypatch.setattr(auditor, "audit_pdb", fake_audit_pdb)

    cli._audit_rcsb_batch(str(batch), str(output))

    assert "Auditing 3 structures" in capsys.readouterr().out
    csv_text = output.read_text()
    assert "1ABC,100.0" in csv_text
    assert "2XYZ,ERROR" in csv_text
    assert "BADID" not in csv_text


def test_cli_correct_af3_mirror_and_mirror_id(monkeypatch, capsys, tmp_path):
    import chiralfold.pdb_pipeline as pdb_pipeline

    fixture = Path(__file__).parent / "fixtures" / "af3_correction" / "synthetic_l_ala_inverted.pdb"
    corrected = tmp_path / "corrected.pdb"
    captured = _run_cli(monkeypatch, capsys, ["correct-af3", str(fixture), "-o", str(corrected)])
    assert "Output:" in captured.out
    assert corrected.exists()

    pdb = tmp_path / "input.pdb"
    mirrored = tmp_path / "mirrored.pdb"
    pdb.write_text(_tiny_ala_pdb())
    captured = _run_cli(monkeypatch, capsys, ["mirror", str(pdb), "-o", str(mirrored), "--chains", "A"])
    assert "Mirrored" in captured.out
    assert mirrored.exists()

    def fake_fetch_and_mirror(pdb_id, output, chains=None):
        Path(output).write_text("END\n")
        return {"n_atoms": 3, "output_path": output}

    monkeypatch.setattr(pdb_pipeline, "fetch_and_mirror", fake_fetch_and_mirror)
    captured = _run_cli(monkeypatch, capsys, ["mirror-id", "1abc", "-o", str(tmp_path / "id.pdb")])
    assert "Downloaded 1ABC" in captured.out


def test_cli_enumerate_json_and_text(monkeypatch, capsys):
    captured = _run_cli(monkeypatch, capsys, ["enumerate", "AFK", "--json", "--top", "2"])
    payload = json.loads(captured.out)
    assert len(payload) == 2
    assert payload[0]["rank"] == 1

    import chiralfold.enumerate as enum_mod

    called = {}

    def fake_enumerate(sequence, top_n):
        called["args"] = (sequence, top_n)
        return [{"rank": 1, "chirality_pattern": "LLL", "n_d": 0, "n_l": 3, "score": 1.0}]

    def fake_format(results):
        called["formatted"] = results
        return "formatted"

    monkeypatch.setattr(enum_mod, "enumerate_diastereomers", fake_enumerate)
    monkeypatch.setattr(enum_mod, "format_enumeration_results", fake_format)

    _run_cli(monkeypatch, capsys, ["enumerate", "afk", "--top", "2"])
    assert called["args"] == ("AFK", 2)
    assert called["formatted"][0]["chirality_pattern"] == "LLL"


def test_cli_score_interface_text_and_json(monkeypatch, capsys, tmp_path):
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.pdb"
    receptor.write_text(_tiny_ala_pdb(chain="A"))
    ligand.write_text(_tiny_ala_pdb(record="HETATM", resname="DAL", chain="B", x_offset=3.0))

    captured = _run_cli(
        monkeypatch,
        capsys,
        ["score-interface", str(receptor), str(ligand), "--receptor-chain", "A", "--ligand-chain", "B"],
    )
    assert "Buried Surface Area" in captured.out
    assert "Interface score" in captured.out

    captured = _run_cli(
        monkeypatch,
        capsys,
        [
            "score-interface",
            str(receptor),
            str(ligand),
            "--receptor-chain",
            "A",
            "--ligand-chain",
            "B",
            "--json",
        ],
    )
    assert json.loads(captured.out)["bsa"] > 0
