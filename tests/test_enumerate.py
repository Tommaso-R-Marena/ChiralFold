from chiralfold.enumerate import enumerate_diastereomers


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
