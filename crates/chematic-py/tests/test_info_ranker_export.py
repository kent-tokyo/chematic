"""Snapshot ownership, table formatting and typed InfoTheory output paths."""

import io
import os

import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I


def populated(metric=I.InfoType.ENTROPY):
    ranker = I.InfoBitRanker(4, 2, metric)
    for label, ids in [(0, [0, 1]), (1, [1, 2])]:
        vector = ExplicitBitVect(4)
        for bit in ids:
            vector.SetBit(bit)
        ranker.AccumulateVotes(vector, label)
    return ranker


EXPECTED = (
    "         Bit InfoContent     class0     class1\n"
    "           0           1         1         0\n"
    "           2           1         0         1\n"
    "           1           0         1         1\n"
    "           3           0         0         0\n"
)


@pytest.mark.parametrize("kind", ["str", "bytes", "path", "custom_str", "custom_bytes"])
def test_exact_table_and_filesystem_paths(tmp_path, kind):
    ranker = populated()
    ranker.GetTopN(4)
    path = tmp_path / "résultat.txt"

    class Custom:
        def __fspath__(self):
            return os.fsencode(path) if kind == "custom_bytes" else str(path)

    value = {"str": str(path), "bytes": os.fsencode(path), "path": path}.get(
        kind, Custom()
    )
    path.write_text("previous content that must disappear\n" * 20)
    assert ranker.WriteTopBitsToFile(fileName=value) is None
    assert path.read_bytes() == EXPECTED.encode("ascii")


@pytest.mark.parametrize("with_votes", [False, True])
def test_before_ranking_refuses_without_touching_destination(tmp_path, with_votes):
    ranker = populated() if with_votes else I.InfoBitRanker(4, 2)
    path = tmp_path / "retained.txt"
    path.write_text("keep")
    with pytest.raises(ValueError, match="successful GetTopN"):
        ranker.WriteTopBitsToFile(path)
    assert path.read_text() == "keep"
    path.unlink()
    with pytest.raises(ValueError, match="successful GetTopN"):
        ranker.WriteTopBitsToFile(path)
    assert not path.exists()


def test_ranking_snapshot_survives_votes_mask_bias_and_return_array_edits(tmp_path):
    ranker = populated(I.InfoType.BIASENTROPY)
    rows = ranker.GetTopN(4)
    rows[:] = -123
    ranker.AccumulateVotes(ExplicitBitVect(4), 0)
    ranker.SetMaskBits([2])
    ranker.SetBiasList([1])
    path = tmp_path / "snapshot.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text() == EXPECTED
    reranked = ranker.GetTopN(1)
    ranker.WriteTopBitsToFile(path)
    assert len(path.read_text().splitlines()) == 2
    assert int(path.read_text().splitlines()[1].split()[0]) == int(reranked[0, 0]) == 2


def test_zero_ranking_replaces_snapshot_with_header_only(tmp_path):
    ranker = populated()
    ranker.GetTopN(4)
    ranker.GetTopN(0)
    path = tmp_path / "empty.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text() == EXPECTED.splitlines(keepends=True)[0]


@pytest.mark.parametrize("bad", [-1, 5, True, 1.5, "2"])
def test_failed_ranking_preserves_previous_snapshot(tmp_path, bad):
    ranker = populated()
    ranker.GetTopN(4)
    with pytest.raises((ValueError, TypeError)):
        ranker.GetTopN(bad)
    path = tmp_path / "unchanged.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text() == EXPECTED


@pytest.mark.parametrize("bad", [None, 42, True, io.StringIO(), bytearray(b"x")])
def test_wrong_output_path_type(tmp_path, bad, monkeypatch):
    ranker = populated()
    ranker.GetTopN(4)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(TypeError):
        ranker.WriteTopBitsToFile(bad)
    assert list(tmp_path.iterdir()) == []


def test_io_failure_preserves_ranked_data_and_snapshot(tmp_path):
    ranker = populated()
    expected = ranker.GetTopN(4)
    for path in [tmp_path, tmp_path / "missing" / "output.txt"]:
        with pytest.raises(OSError):
            ranker.WriteTopBitsToFile(path)
    with pytest.raises(ValueError):
        ranker.WriteTopBitsToFile("bad\x00name")
    path = tmp_path / "valid.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text() == EXPECTED
    np.testing.assert_array_equal(ranker.GetTopN(4), expected)


def test_multiclass_header_count_precision(tmp_path):
    ranker = I.InfoBitRanker(2, 12, I.InfoType.CHISQUARE)
    for label in range(12):
        fp = ExplicitBitVect(2)
        fp.SetBit(label % 2)
        ranker.AccumulateVotes(fp, label)
    rows = ranker.GetTopN(2)
    path = tmp_path / "classes.txt"
    ranker.WriteTopBitsToFile(path)
    lines = path.read_text().splitlines()
    assert lines[0].endswith("     class10     class11")
    for row, line in zip(rows, lines[1:]):
        tokens = line.split()
        assert len(tokens) == 14
        assert int(tokens[0]) == int(row[0])
        assert tokens[1] == format(float(row[1]), ".5g")
        assert [int(x) for x in tokens[2:]] == row[2:].astype(int).tolist()


def test_five_significant_digits_including_small_scores(tmp_path):
    ranker = I.InfoBitRanker(1, 2)
    for label in range(2):
        for row in range(128):
            fp = ExplicitBitVect(1)
            if row < 64 + label:
                fp.SetBit(0)
            ranker.AccumulateVotes(fp, label)
    score = ranker.GetTopN(1)[0, 1]
    assert 0 < score < 0.0001
    path = tmp_path / "small.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text().splitlines()[1].split()[1] == format(score, ".5g")
    assert "e-" in path.read_text()


def test_pathlike_exception_propagates_and_snapshot_is_retained(tmp_path):
    failure = RuntimeError("custom path failure")

    class Broken:
        def __fspath__(self):
            raise failure

    ranker = populated()
    ranker.GetTopN(4)
    with pytest.raises(RuntimeError) as caught:
        ranker.WriteTopBitsToFile(Broken())
    assert caught.value is failure
    path = tmp_path / "valid.txt"
    ranker.WriteTopBitsToFile(path)
    assert path.read_text() == EXPECTED
