"""Tests for structure_io + fetch helpers (offline; network marked slow)."""

from __future__ import annotations

import os
from pathlib import Path
from unittest import mock

import pytest

from chiralfold.fetch import (
    FetchError,
    classify_query,
    fetch_structure,
    resolve_local,
)
from chiralfold.structure_io import (
    detect_format,
    ensure_pdb_path,
    mmcif_to_pdb,
    parse_fasta,
    uniprot_from_text,
    write_mmcif_as_pdb,
)


def _mini_cif() -> str:
    return """data_TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1 A ALA CA 1.000 2.000 3.000 1.00 10.00 ? 1
ATOM 2 N N ALA A 1 1 A ALA N 1.100 2.100 3.100 1.00 10.00 ? 1
#
"""


def test_classify_query_kinds():
    assert classify_query("1ubq") == ("pdb", "1UBQ")
    assert classify_query("P04637") == ("uniprot", "P04637")
    assert classify_query("AF-P04637-F1") == ("afdb", "AF-P04637-F1")
    assert classify_query("AF-P04637-F1-model_v4") == ("afdb", "AF-P04637-F1")
    assert classify_query("")[0] == "unknown"


def test_uniprot_from_fasta_headers():
    assert uniprot_from_text("sp|P04637|P53_HUMAN") == "P04637"
    assert uniprot_from_text("tr|Q9Y6K1|FOO") == "Q9Y6K1"
    assert uniprot_from_text("random header") is None


def test_parse_fasta_and_mmcif_roundtrip(tmp_path: Path):
    records = parse_fasta(">sp|P04637|P53_HUMAN\nMEEPQSD\n>seq2\nAAA\n")
    assert records[0][1] == "MEEPQSD"
    assert records[1][1] == "AAA"

    cif = tmp_path / "mini.cif"
    cif.write_text(_mini_cif())
    assert detect_format(str(cif)) == "mmcif"
    pdb_path = write_mmcif_as_pdb(str(cif))
    text = Path(pdb_path).read_text()
    assert "ATOM" in text
    assert "ALA" in text
    assert ensure_pdb_path(pdb_path) == pdb_path


def test_ensure_pdb_rejects_fasta_without_fetch(tmp_path: Path):
    fa = tmp_path / "seq.fasta"
    fa.write_text(">onlyseq\nACDEFGHIK\n")
    with pytest.raises(ValueError, match="FASTA"):
        ensure_pdb_path(str(fa))


def test_resolve_local_fasta_needs_uniprot(tmp_path: Path):
    fa = tmp_path / "bare.fasta"
    fa.write_text(">peptide\nACDE\n")
    with pytest.raises(FetchError, match="UniProt"):
        resolve_local(str(fa))


def test_fetch_structure_mock_rcsb(tmp_path: Path):
    toy = Path("chiralfold/data/examples/toy_ubiquitin_fragment.pdb").read_bytes()

    def fake_get(url, timeout=45):
        assert "1UBQ" in url.upper() or "1ubq" in url
        return toy

    with mock.patch("chiralfold.fetch._http_get", side_effect=fake_get):
        result = fetch_structure("1UBQ")
    assert result.source == "rcsb"
    assert Path(result.path).is_file()
    assert "ATOM" in Path(result.path).read_text()


def test_fetch_structure_mock_afdb():
    toy = Path("chiralfold/data/examples/toy_ubiquitin_fragment.pdb").read_bytes()

    def fake_get(url, timeout=45):
        if "api/prediction" in url:
            raise OSError("skip api")
        if "AF-P04637" in url and url.endswith(".pdb"):
            return toy
        raise urllib_error(url)

    def urllib_error(url):
        raise Exception(f"unexpected {url}")

    with mock.patch("chiralfold.fetch._http_get", side_effect=fake_get):
        result = fetch_structure("P04637")
    assert result.source == "afdb"
    assert "P04637" in result.label


def test_resolve_local_fasta_fetches_afdb(tmp_path: Path):
    toy = Path("chiralfold/data/examples/toy_ubiquitin_fragment.pdb").read_bytes()
    fa = tmp_path / "p53.fasta"
    fa.write_text(">sp|P04637|P53_HUMAN Cellular tumor antigen p53\nMEEPQSD\n")

    def fake_get(url, timeout=45):
        if "api/prediction" in url:
            raise OSError("skip")
        if "AF-P04637" in url and url.endswith(".pdb"):
            return toy
        raise OSError(url)

    with mock.patch("chiralfold.fetch._http_get", side_effect=fake_get):
        result = resolve_local(str(fa))
    assert result.source == "fasta-afdb"
    assert Path(result.path).is_file()


@pytest.mark.slow
def test_live_fetch_rcsb_1crn():
    result = fetch_structure("1CRN")
    assert result.source == "rcsb"
    assert os.path.getsize(result.path) > 500


@pytest.mark.slow
def test_live_fetch_afdb_p01308():
    # Insulin — small, stable AFDB entry
    result = fetch_structure("P01308")
    assert result.source == "afdb"
    assert os.path.getsize(result.path) > 500
