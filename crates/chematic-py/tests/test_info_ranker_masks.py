"""Typed bit-mask selection and historical vote semantics for RDKit #9601."""

import itertools
import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I


class Sparse:
    def __init__(self, bits, width=4):
        self.bits, self.width = bits, width

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return self.bits


def rows(ranker, count):
    return {int(row[0]): row[1:] for row in ranker.GetTopN(count)}


def populated(metric=I.InfoType.ENTROPY):
    ranker = I.InfoBitRanker(4, 2, metric)
    ranker.AccumulateVotes(Sparse([0, 1]), 0)
    ranker.AccumulateVotes(Sparse([2, 3]), 1)
    return ranker


@pytest.mark.parametrize("metric", [I.InfoType.ENTROPY, I.InfoType.CHISQUARE])
def test_mask_before_votes_limits_candidates_and_retained_counts(metric):
    r = I.InfoBitRanker(4, 2, metric)
    r.SetMaskBits([0, 2])
    r.AccumulateVotes(Sparse([0, 1]), 0)
    r.AccumulateVotes(Sparse([2, 3]), 1)
    result = rows(r, 2)
    assert sorted(result) == [0, 2]
    np.testing.assert_array_equal(result[0][1:], [1, 0])
    np.testing.assert_array_equal(result[2][1:], [0, 1])
    r.SetMaskBits(range(4))
    result = rows(r, 4)
    np.testing.assert_array_equal(result[1][1:], [0, 0])
    np.testing.assert_array_equal(result[3][1:], [0, 0])


def test_mask_after_votes_does_not_erase_historical_counts():
    r = populated()
    r.SetMaskBits([0])
    r.AccumulateVotes(Sparse([0, 1, 2, 3]), 0)
    r.SetMaskBits(range(4))
    result = rows(r, 4)
    np.testing.assert_array_equal(result[0][1:], [2, 0])
    np.testing.assert_array_equal(result[1][1:], [1, 0])
    np.testing.assert_array_equal(result[2][1:], [0, 1])
    np.testing.assert_array_equal(result[3][1:], [0, 1])
    # All three observations still contribute to class totals.
    assert result[1][0] == pytest.approx(0.2516291673878229)


def test_replacement_is_not_union_and_restoring_all_ids_works():
    r = populated()
    r.SetMaskBits([0, 1])
    r.SetMaskBits([2])
    assert list(rows(r, 1)) == [2]
    r.SetMaskBits(range(4))
    assert sorted(rows(r, 4)) == [0, 1, 2, 3]


def test_empty_mask_counts_records_but_no_on_votes():
    r = I.InfoBitRanker(4, 2)
    r.SetMaskBits([])
    r.AccumulateVotes(Sparse([0, 1, 2, 3]), 0)
    assert r.GetTopN(0).shape == (0, 4)
    with pytest.raises(ValueError):
        r.GetTopN(1)
    r.SetMaskBits(range(4))
    np.testing.assert_array_equal(r.GetTopN(4)[:, 1:], np.zeros((4, 3)))


def test_duplicates_normalize_and_input_order_does_not_set_ties():
    r = populated()
    r.SetMaskBits([3, 0, 3, 0])
    np.testing.assert_array_equal(r.GetTopN(2)[:, 0], [0, 3])
    with pytest.raises(ValueError):
        r.GetTopN(3)


def test_numpy_integer_generator_and_source_ownership():
    source = [np.int64(3), np.int64(0)]
    r = populated()
    r.SetMaskBits(x for x in source)
    source[:] = [1, 2]
    assert sorted(rows(r, 2)) == [0, 3]
    result = r.GetTopN(2)
    result[:] = -1
    assert sorted(rows(r, 2)) == [0, 3]


@pytest.mark.parametrize(
    "value,exception",
    [
        (None, TypeError),
        (1, TypeError),
        ("01", TypeError),
        (b"\0\1", TypeError),
        (bytearray([0, 1]), TypeError),
        ([True], TypeError),
        ([np.bool_(False)], TypeError),
        ([1.0], TypeError),
        (["1"], TypeError),
        ([None], TypeError),
        ([-1], ValueError),
        ([4], ValueError),
        ([2**80], ValueError),
        ([0, 0, 0, 0, 0], ValueError),
    ],
)
def test_bad_mask_keeps_previous_mask_and_counts(value, exception):
    r = populated()
    r.SetMaskBits([0, 2])
    before = r.GetTopN(2)
    with pytest.raises(exception):
        r.SetMaskBits(value)
    np.testing.assert_array_equal(r.GetTopN(2), before)


def test_iteration_failure_is_atomic_and_preserves_original_exception():
    r = populated()
    r.SetMaskBits([0, 2])
    before = r.GetTopN(2)
    failure = RuntimeError("source failure")

    def failing():
        yield 1
        raise failure

    with pytest.raises(RuntimeError) as error:
        r.SetMaskBits(failing())
    assert error.value is failure
    np.testing.assert_array_equal(r.GetTopN(2), before)


def test_infinite_duplicate_mask_has_bounded_consumption():
    r = populated()
    r.SetMaskBits([0])
    with pytest.raises(ValueError):
        r.SetMaskBits(itertools.repeat(0))
    assert list(rows(r, 1)) == [0]


def test_mask_does_not_hide_invalid_vector_bits_or_partial_votes():
    r = populated()
    r.SetMaskBits([0])
    before = r.GetTopN(1)
    for bits in [[0, 4], [1, 1], [0, True]]:
        with pytest.raises((ValueError, TypeError)):
            r.AccumulateVotes(Sparse(bits), 0)
        np.testing.assert_array_equal(r.GetTopN(1), before)


def test_explicit_vector_and_index_objects():
    class Index:
        def __index__(self):
            return 2

    r = I.InfoBitRanker(4, 2)
    r.SetMaskBits([Index()])
    fp = ExplicitBitVect(4)
    fp.SetBit(2)
    r.AccumulateVotes(fp, 1)
    np.testing.assert_array_equal(r.GetTopN(1), [[2, 0, 0, 1]])


def test_no_records_rule_applies_even_to_empty_mask():
    r = I.InfoBitRanker(4, 2)
    r.SetMaskBits([])
    with pytest.raises(ValueError, match="one record"):
        r.GetTopN(0)
