"""Class preference semantics, typed input and transactional replacement."""

import itertools
from fractions import Fraction

import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I


class Sparse:
    def __init__(self, width, bits):
        self.width, self.bits = width, bits

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return self.bits


def populated(metric=I.InfoType.BIASENTROPY, bias=(0,)):
    # Bit 0 favours class 0, bit 1 class 1; bit 2 is constant.
    result = I.InfoBitRanker(3, 2, metric)
    result.SetBiasList(bias)
    for bits, label in [([0, 2], 0), ([0, 2], 0), ([1, 2], 1), ([1, 2], 1)]:
        result.AccumulateVotes(Sparse(3, bits), label)
    return result


@pytest.mark.parametrize("metric,score", [(2, 1.0), (4, 4.0)])
def test_directional_preference_retains_counts_and_bit_ids(metric, score):
    result = populated(metric)
    np.testing.assert_array_equal(
        result.GetTopN(3), [[0, score, 2, 0], [1, 0, 0, 2], [2, 0, 2, 2]]
    )
    result.SetBiasList([1])
    np.testing.assert_array_equal(
        result.GetTopN(3), [[1, score, 0, 2], [0, 0, 2, 0], [2, 0, 2, 2]]
    )
    result.SetBiasList([])
    np.testing.assert_array_equal(
        result.GetTopN(3), [[0, score, 2, 0], [1, score, 0, 2], [2, 0, 2, 2]]
    )


@pytest.mark.parametrize("metric", [2, 4])
def test_default_and_empty_bias_accept_all_features(metric):
    result = populated(metric, [])
    expected = result.GetTopN(3)
    default = I.InfoBitRanker(3, 2, metric)
    for bits, label in [([0, 2], 0), ([0, 2], 0), ([1, 2], 1), ([1, 2], 1)]:
        default.AccumulateVotes(Sparse(3, bits), label)
    np.testing.assert_array_equal(default.GetTopN(3), expected)


@pytest.mark.parametrize("metric", [1, 3])
def test_unbiased_modes_ignore_valid_bias(metric):
    result = populated(metric, [])
    before = result.GetTopN(3)
    result.SetBiasList([0])
    np.testing.assert_array_equal(result.GetTopN(3), before)


def test_fraction_rule_not_raw_count_and_equality_is_inclusive():
    result = I.InfoBitRanker(2, 3, I.InfoType.BIASCHISQUARE)
    # Totals 2, 4, 4. Bit 0: fractions 1/2, 1/2, 0 (tie accepted).
    # Bit 1: 1/2, 1/4, 0 (accepted despite matching raw counts).
    for label, total in enumerate([2, 4, 4]):
        for row in range(total):
            bits = []
            if row < [1, 2, 0][label]:
                bits.append(0)
            if row < [1, 1, 0][label]:
                bits.append(1)
            result.AccumulateVotes(Sparse(2, bits), label)
    result.SetBiasList([0])
    assert (result.GetTopN(2)[:, 1] > 0).all()
    result.SetBiasList([1])
    scores = {int(row[0]): row[1] for row in result.GetTopN(2)}
    assert scores[0] > 0
    assert scores[1] == 0
    result.SetBiasList([2, 1])  # any preferred class can pass
    assert result.GetTopN(2)[0, 1] > 0


@pytest.mark.parametrize("metric", [2, 4])
def test_empty_classes_and_single_class_are_finite(metric):
    result = I.InfoBitRanker(2, 3, metric)
    result.SetBiasList([2])  # unobserved preference; no divisions by zero
    result.AccumulateVotes(Sparse(2, [0]), 0)
    result.AccumulateVotes(Sparse(2, [1]), 1)
    with np.errstate(all="raise"):
        rows = result.GetTopN(2)
    np.testing.assert_array_equal(rows[:, 1], [0, 0])
    result.SetBiasList([0, 2])
    assert result.GetTopN(2)[0, 1] > 0
    single = I.InfoBitRanker(1, 1, metric)
    single.SetBiasList([])
    with pytest.raises(ValueError):
        single.SetBiasList([0])
    single.AccumulateVotes(Sparse(1, []), 0)
    np.testing.assert_array_equal(single.GetTopN(1), [[0, 0, 0]])


@pytest.mark.parametrize(
    "bad,exception",
    [
        (None, TypeError),
        (0, TypeError),
        ("0", TypeError),
        ([True], TypeError),
        ([np.bool_(False)], TypeError),
        ([0.0], TypeError),
        (["0"], TypeError),
        ([-1], ValueError),
        ([2], ValueError),
        ([0, 0], ValueError),
        ([0, 1], ValueError),
        ([0, 1, 0], ValueError),
    ],
)
def test_invalid_replacement_keeps_bias_and_votes(bad, exception):
    result = populated()
    before = result.GetTopN(3)
    with pytest.raises(exception):
        result.SetBiasList(bad)
    np.testing.assert_array_equal(result.GetTopN(3), before)


def test_broken_generator_and_infinite_input_preserve_state():
    result = populated()
    before = result.GetTopN(3)

    def broken():
        yield 1
        raise RuntimeError("broken preference")

    with pytest.raises(RuntimeError, match="broken preference"):
        result.SetBiasList(broken())
    with pytest.raises(ValueError):
        result.SetBiasList(itertools.cycle([0, 1]))
    np.testing.assert_array_equal(result.GetTopN(3), before)


def test_owns_bias_input_and_accepts_numpy_index_protocol():
    result = populated()
    preference = [np.int64(1)]
    result.SetBiasList(iter(preference))
    expected = result.GetTopN(3)
    preference[0] = np.int64(0)
    np.testing.assert_array_equal(result.GetTopN(3), expected)
    result.GetTopN(3)[:] = 99
    np.testing.assert_array_equal(result.GetTopN(3), expected)
    assert result.GetTopN(np.int64(0)).shape == (0, 4)


@pytest.mark.parametrize("metric", [2, 4])
def test_explicit_vectors_and_failed_votes_preserve_bias(metric):
    result = I.InfoBitRanker(3, 2, np.int64(metric))
    result.SetBiasList([1])
    vector = ExplicitBitVect(3)
    vector.SetBit(0)
    result.AccumulateVotes(vector, 0)
    vector = ExplicitBitVect(3)
    vector.SetBit(1)
    result.AccumulateVotes(vector, 1)
    before = result.GetTopN(3)
    for value in [Sparse(2, []), Sparse(3, [0, 0]), Sparse(3, [3])]:
        with pytest.raises(ValueError):
            result.AccumulateVotes(value, 0)
    np.testing.assert_array_equal(result.GetTopN(3), before)
    assert before[0, 0] == 1


def test_exhaustive_small_tables_against_exact_direction_rule():
    # All 2x3 on/off tables with 1..3 observations in each class.
    for totals in itertools.product(range(1, 4), repeat=3):
        for hits in itertools.product(*(range(n + 1) for n in totals)):
            for preferred in [(0,), (1,), (2,), (0, 2)]:
                result = I.InfoBitRanker(1, 3, I.InfoType.BIASCHISQUARE)
                result.SetBiasList(preferred)
                for label, count in enumerate(totals):
                    for row in range(count):
                        result.AccumulateVotes(
                            Sparse(1, [0] if row < hits[label] else []), label
                        )
                fractions = [Fraction(a, b) for a, b in zip(hits, totals)]
                accepts = max(fractions[c] for c in preferred) >= max(
                    fractions[c] for c in range(3) if c not in preferred
                )
                statistic = Fraction(0)
                table = [hits, tuple(n - k for n, k in zip(totals, hits))]
                for row in table:
                    for label, observed in enumerate(row):
                        expected = Fraction(sum(row) * totals[label], sum(totals))
                        if expected:
                            statistic += (observed - expected) ** 2 / expected
                expected_score = float(statistic) if accepts else 0.0
                assert result.GetTopN(1)[0, 1] == pytest.approx(
                    expected_score, abs=1e-14
                )


def test_maximum_class_preference_and_failed_boundary_replacement():
    result = I.InfoBitRanker(1, 256, I.InfoType.BIASCHISQUARE)
    result.SetBiasList(range(255))
    result.AccumulateVotes(Sparse(1, [0]), 254)
    result.AccumulateVotes(Sparse(1, []), 255)
    expected = result.GetTopN(1)
    assert expected[0, 1] == 2
    with pytest.raises(ValueError, match="fewer than"):
        result.SetBiasList(range(256))
    np.testing.assert_array_equal(result.GetTopN(1), expected)
