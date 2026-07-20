from chiralfold import enumerate as enum_mod
from chiralfold.enumerate import enumerate_diastereomers, format_enumeration_results


def test_enumerate_diastereomers_short_sequence_top_n():
    results = enumerate_diastereomers("AFK", top_n=3, n_conformers=1, seed=123)

    assert 0 < len(results) <= 3
    scores = [entry["score"] for entry in results]
    assert scores == sorted(scores, reverse=True)
    for rank, entry in enumerate(results, start=1):
        assert entry["rank"] == rank
        assert len(entry["chirality_pattern"]) == 3
        assert set(entry["chirality_pattern"]).issubset({"L", "D"})
        assert isinstance(entry["score"], float)
        assert entry["valid"] is True


def test_enumeration_pattern_helpers_cover_glycine_and_random_sampling():
    stereo_positions = enum_mod._stereo_positions("AGF")
    assert stereo_positions == [0, 2]
    assert enum_mod._expand_pattern("AGF", stereo_positions, ("D", "L")) == "DLL"
    assert enum_mod._score_diastereomer(0, 0, 2) == 101.0

    rng = enum_mod.random.Random(7)
    sampled = enum_mod._random_patterns(12, 5, rng)
    assert len(sampled) == 5
    assert len(set(sampled)) == 5
    assert all(len(pattern) == 12 for pattern in sampled)

    assert enum_mod._random_patterns(2, 10, enum_mod.random.Random(1)) == enum_mod._all_patterns(2)


def test_enumerate_invalid_sequence_marks_entries_invalid(monkeypatch):
    class DummyModel:
        def __init__(self, n_conformers, fix_planarity):
            assert (n_conformers, fix_planarity) == (3, True)

        def predict(self, sequence, chirality_pattern):
            raise AssertionError("predict should not run when SMILES generation fails")

    def fake_smiles(sequence, chirality):
        raise ValueError("unknown residue")

    monkeypatch.setattr("chiralfold.model.ChiralFold", DummyModel)
    monkeypatch.setattr("chiralfold.model.mixed_peptide_smiles", fake_smiles)

    results = enumerate_diastereomers("Z", top_n=1)

    assert results == [
        {
            "chirality_pattern": "L",
            "n_d": 0,
            "n_l": 1,
            "n_stereo": 1,
            "score": 100.0,
            "smiles": None,
            "n_conformers_generated": 0,
            "valid": False,
            "rank": 1,
        }
    ]


def test_enumerate_predict_exception_keeps_smiles_but_zero_conformers(monkeypatch):
    class DummyModel:
        def __init__(self, n_conformers, fix_planarity):
            pass

        def predict(self, sequence, chirality_pattern):
            raise RuntimeError("embedding failed")

    monkeypatch.setattr("chiralfold.model.ChiralFold", DummyModel)
    monkeypatch.setattr("chiralfold.model.mixed_peptide_smiles", lambda sequence, chirality: "C")

    results = enumerate_diastereomers("A", top_n=1)

    assert results[0]["smiles"] == "C"
    assert results[0]["n_conformers_generated"] == 0
    assert results[0]["valid"] is True


def test_format_enumeration_results_empty_and_truncated_pattern():
    assert format_enumeration_results([]) == "(no results)"

    table = format_enumeration_results(
        [
            {
                "rank": 1,
                "chirality_pattern": "L" * 40,
                "n_d": 0,
                "n_l": 40,
                "score": 101.5,
                "n_conformers_generated": 3,
                "valid": True,
            },
            {
                "chirality_pattern": "D",
                "n_d": 1,
                "n_l": 0,
                "score": 90.0,
                "valid": False,
            },
        ]
    )

    assert "ChiralFold" in table
    assert "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLL..." in table
    assert "yes" in table
    assert "no" in table
