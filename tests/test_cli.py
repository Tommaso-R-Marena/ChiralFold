import json
import sys

import pytest

from chiralfold import __version__
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
