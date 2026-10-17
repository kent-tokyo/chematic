"""Entropy ranker input/state contracts, motivated by RDKit #9601."""

import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I


class Sparse:
    def __init__(self, width, bits):
        self.width = width
        self.bits = bits

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return self.bits


def vector(width, bits):
    result = ExplicitBitVect(width)
    for bit in bits:
        result.SetBit(bit)
    return result


def populated():
    ranker = I.InfoBitRanker(4, 2)
    ranker.AccumulateVotes(vector(4, [0, 1]), 0)
    ranker.AccumulateVotes(Sparse(4, [1, 2]), 1)
    return ranker


def test_perfect_independent_and_constant_features():
    ranker = populated()
    rows = ranker.GetTopN(4)
    np.testing.assert_array_equal(rows[:, 0], [0.0, 2.0, 1.0, 3.0])
    np.testing.assert_allclose(rows[:, 1], [1.0, 1.0, 0.0, 0.0], atol=1e-15)
    np.testing.assert_array_equal(
        rows[:, 2:], [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]
    )
    assert rows.dtype == np.float64


def test_vote_order_and_vector_storage_do_not_change_scores():
    ranker = I.InfoBitRanker(4, 2)
    ranker.AccumulateVotes(Sparse(4, iter([2, 1])), 1)
    ranker.AccumulateVotes(Sparse(4, [1, 0]), 0)
    np.testing.assert_array_equal(ranker.GetTopN(4), populated().GetTopN(4))


def test_multiclass_counts_and_unobserved_classes():
    ranker = I.InfoBitRanker(3, 4)
    for label in range(3):
        ranker.AccumulateVotes(Sparse(3, [label]), label)
    rows = ranker.GetTopN(3)
    expected_gain = np.log2(3) - 2.0 / 3.0
    np.testing.assert_allclose(rows[:, 1], [expected_gain] * 3, atol=1e-15)
    np.testing.assert_array_equal(rows[:, 2:5], np.eye(3))
    assert np.all(rows[:, -1] == 0)


def test_empty_vectors_are_real_observations():
    ranker = I.InfoBitRanker(2, 2)
    ranker.AccumulateVotes(Sparse(2, []), 0)
    ranker.AccumulateVotes(Sparse(2, []), 1)
    np.testing.assert_array_equal(
        ranker.GetTopN(2), [[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]]
    )


def test_one_observed_class_has_zero_information_gain():
    ranker = I.InfoBitRanker(2, 1)
    ranker.AccumulateVotes(Sparse(2, [0]), 0)
    np.testing.assert_array_equal(ranker.GetTopN(2), [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])


def test_top_zero_and_result_isolation():
    ranker = populated()
    assert ranker.GetTopN(0).shape == (0, 4)
    original = ranker.GetTopN(4)
    altered = ranker.GetTopN(4)
    altered[:] = -99.0
    np.testing.assert_array_equal(ranker.GetTopN(4), original)


@pytest.mark.parametrize(
    "width,classes,exc",
    [
        (0, 2, ValueError),
        (-1, 2, ValueError),
        (65_537, 2, ValueError),
        (4, 0, ValueError),
        (4, 257, ValueError),
        (65_536, 17, ValueError),
        (True, 2, TypeError),
        (4, False, TypeError),
        (4.0, 2, TypeError),
        (4, "2", TypeError),
    ],
)
def test_constructor_bounds(width, classes, exc):
    with pytest.raises(exc):
        I.InfoBitRanker(width, classes)


@pytest.mark.parametrize(
    "metric,exc",
    [
        (I.InfoType.BIASENTROPY, ValueError),
        (I.InfoType.BIASCHISQUARE, ValueError),
        (99, ValueError),
        (True, TypeError),
    ],
)
def test_unsupported_metrics_are_explicit(metric, exc):
    with pytest.raises(exc):
        I.InfoBitRanker(4, 2, metric)


@pytest.mark.parametrize(
    "bad,exc",
    [
        (None, TypeError),
        ([], TypeError),
        (Sparse(3, []), ValueError),
        (Sparse(True, []), TypeError),
        (Sparse(4.0, []), TypeError),
        (Sparse(4, None), TypeError),
        (Sparse(4, [0, 0]), ValueError),
        (Sparse(4, [-1]), ValueError),
        (Sparse(4, [4]), ValueError),
        (Sparse(4, [True]), TypeError),
        (Sparse(4, [1.0]), TypeError),
        (Sparse(4, ["1"]), TypeError),
    ],
)
def test_invalid_votes_do_not_modify_counts(bad, exc):
    ranker = populated()
    before = ranker.GetTopN(4)
    with pytest.raises(exc):
        ranker.AccumulateVotes(bad, 0)
    np.testing.assert_array_equal(ranker.GetTopN(4), before)


def test_failing_iterator_and_vector_mutation_leave_no_borrowed_state():
    def broken():
        yield 0
        raise RuntimeError("source iteration failed")

    ranker = populated()
    before = ranker.GetTopN(4)
    with pytest.raises(RuntimeError, match="source iteration"):
        ranker.AccumulateVotes(Sparse(4, broken()), 0)
    np.testing.assert_array_equal(ranker.GetTopN(4), before)
    source = vector(4, [0])
    ranker.AccumulateVotes(source, 1)
    expected = ranker.GetTopN(4)
    source.SetBit(3)
    np.testing.assert_array_equal(ranker.GetTopN(4), expected)


@pytest.mark.parametrize(
    "label,exc",
    [(-1, ValueError), (2, ValueError), (True, TypeError), (1.0, TypeError)],
)
def test_invalid_labels_are_checked_before_vector_access(label, exc):
    class Unreadable:
        def GetNumBits(self):
            raise AssertionError("vector should not be accessed")

        def GetOnBits(self):
            raise AssertionError("vector should not be accessed")

    ranker = populated()
    before = ranker.GetTopN(4)
    with pytest.raises(exc):
        ranker.AccumulateVotes(Unreadable(), label)
    np.testing.assert_array_equal(ranker.GetTopN(4), before)


@pytest.mark.parametrize(
    "num,exc", [(-1, ValueError), (5, ValueError), (True, TypeError), (2.0, TypeError)]
)
def test_invalid_top_count(num, exc):
    with pytest.raises(exc):
        populated().GetTopN(num)


def test_empty_ranker_bias_and_masks_are_explicit_errors():
    ranker = I.InfoBitRanker(4, 2)
    with pytest.raises(ValueError, match="one record"):
        ranker.GetTopN(1)
    for method in (ranker.SetBiasList, ranker.SetMaskBits):
        with pytest.raises(NotImplementedError):
            method([1])


def test_numeric_index_protocol_and_dimension_boundaries():
    ranker = I.InfoBitRanker(np.int64(65_536), np.int64(1), np.int64(1))
    ranker.AccumulateVotes(Sparse(np.int64(65_536), [np.int64(65_535)]), np.int64(0))
    assert ranker.GetTopN(65_536)[-1, 2] == 1
    ranker = I.InfoBitRanker(1, 256)
    ranker.AccumulateVotes(Sparse(1, [0]), 255)
    assert ranker.GetTopN(1).shape == (1, 258)
    assert ranker.GetTopN(1)[0, -1] == 1


def test_record_limit_preserves_state(monkeypatch):
    import chematic._info_theory as implementation

    monkeypatch.setattr(implementation, "MAX_RANKER_RECORDS", 2)
    ranker = populated()
    before = ranker.GetTopN(4)
    with pytest.raises(ValueError, match="record limit"):
        ranker.AccumulateVotes(Sparse(4, [0]), 1)
    np.testing.assert_array_equal(ranker.GetTopN(4), before)
