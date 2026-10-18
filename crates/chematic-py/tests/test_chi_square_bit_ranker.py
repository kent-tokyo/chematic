"""Pearson ranking, independently specified tables, and shared state contracts."""

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


def ranker(observations, width=4, classes=2):
    result = I.InfoBitRanker(width, classes, I.InfoType.CHISQUARE)
    for bits, label in observations:
        result.AccumulateVotes(Sparse(width, bits), label)
    return result


def test_perfect_independent_and_constant_features():
    result = ranker([([0, 1], 0), ([0], 0), ([1, 2], 1), ([2], 1)])
    np.testing.assert_array_equal(
        result.GetTopN(4),
        [[0, 4, 2, 0], [2, 4, 0, 2], [1, 0, 1, 1], [3, 0, 0, 0]],
    )


def test_multiclass_fraction_reference_and_empty_class():
    # Rows of the on/off contingency table, fourth class unobserved.
    table = [[4, 1, 2, 0], [1, 5, 3, 0]]
    observations = []
    for state, row in enumerate(table):
        for label, count in enumerate(row):
            observations.extend([([0] if state == 0 else [], label)] * count)
    result = ranker(observations, width=1, classes=4).GetTopN(1)
    total = sum(map(sum, table))
    expected = Fraction(0)
    for row in table:
        for label, value in enumerate(row):
            column = sum(r[label] for r in table)
            if column:
                e = Fraction(sum(row) * column, total)
                expected += (value - e) ** 2 / e
    assert result[0, 1] == pytest.approx(float(expected), abs=1e-14)
    np.testing.assert_array_equal(result[0, 2:], table[0])


@pytest.mark.parametrize("bits", [[], [0], [1], [0, 1]])
@pytest.mark.parametrize("classes,labels", [(1, [0, 0]), (3, [0, 2])])
def test_degenerate_tables_are_finite_zero(bits, classes, labels):
    with np.errstate(all="raise"):
        rows = ranker([(bits, label) for label in labels], 2, classes).GetTopN(2)
    np.testing.assert_array_equal(rows[:, 1], [0, 0])
    assert np.isfinite(rows).all()


def test_permutation_complement_and_duplication_laws():
    observations = [([0], 0), ([0, 1], 0), ([1], 1), ([1, 2], 1), ([], 1)]
    original = ranker(observations).GetTopN(4)
    np.testing.assert_array_equal(ranker(observations[::-1]).GetTopN(4), original)
    complement = [(sorted(set(range(4)) - set(bits)), label) for bits, label in observations]
    np.testing.assert_allclose(ranker(complement).GetTopN(4)[:, :2], original[:, :2])
    duplicated = ranker(observations * 2).GetTopN(4)
    np.testing.assert_array_equal(duplicated[:, 0], original[:, 0])
    np.testing.assert_allclose(duplicated[:, 1:], original[:, 1:] * 2)


@pytest.mark.parametrize("bad", [Sparse(4, [0, 0]), Sparse(4, [4]), Sparse(3, [0])])
def test_invalid_votes_preserve_chi_square_state(bad):
    result = ranker([([0], 0), ([1], 1)])
    before = result.GetTopN(4)
    with pytest.raises(ValueError):
        result.AccumulateVotes(bad, 0)
    np.testing.assert_array_equal(result.GetTopN(4), before)


def test_broken_generator_and_explicit_vector_ownership():
    def broken():
        yield 0
        raise RuntimeError("broken vector")

    result = ranker([([0], 0), ([1], 1)])
    before = result.GetTopN(4)
    with pytest.raises(RuntimeError, match="broken vector"):
        result.AccumulateVotes(Sparse(4, broken()), 0)
    np.testing.assert_array_equal(result.GetTopN(4), before)
    vector = ExplicitBitVect(4)
    vector.SetBit(2)
    result.AccumulateVotes(vector, 1)
    before = result.GetTopN(4)
    vector.SetBit(0)
    np.testing.assert_array_equal(result.GetTopN(4), before)


def test_empty_top_zero_and_result_copy():
    result = I.InfoBitRanker(4, 2, np.int64(3))
    with pytest.raises(ValueError, match="one record"):
        result.GetTopN(0)
    result.AccumulateVotes(Sparse(4, []), 0)
    assert result.GetTopN(0).shape == (0, 4)
    before = result.GetTopN(4)
    assert before.dtype == np.float64
    result.GetTopN(4)[:] = -1
    np.testing.assert_array_equal(result.GetTopN(4), before)


def test_counts_beyond_uint16_remain_exact():
    result = I.InfoBitRanker(1, 2, I.InfoType.CHISQUARE)
    for _ in range(65_536):
        result.AccumulateVotes(Sparse(1, [0]), 0)
    result.AccumulateVotes(Sparse(1, []), 1)
    np.testing.assert_allclose(result.GetTopN(1), [[0, 65_537, 65_536, 0]], rtol=1e-14)


def test_zero_signal_bit_keeps_its_id_and_full_ranking():
    # RDKit 2026.09.1 returns duplicate bit 0 and omits bit 2 for this table.
    observations = []
    for label, on_counts in enumerate(([8, 12, 10], [11, 7, 10])):
        for row in range(23):
            observations.append(([bit for bit, n in enumerate(on_counts) if row < n], label))
    rows = ranker(observations, width=3).GetTopN(3)
    np.testing.assert_array_equal(rows[:, 0], [1, 0, 2])
    np.testing.assert_array_equal(rows[:, 2:], [[12, 7], [8, 11], [10, 10]])
    assert rows[-1, 1] == 0
    assert np.all(rows[:-1, 1] > 0)


def test_record_limit_does_not_modify_chi_square_state(monkeypatch):
    import chematic._info_theory as implementation

    monkeypatch.setattr(implementation, "MAX_RANKER_RECORDS", 2)
    result = ranker([([0], 0), ([1], 1)])
    before = result.GetTopN(4)
    with pytest.raises(ValueError, match="record limit"):
        result.AccumulateVotes(Sparse(4, [0]), 1)
    np.testing.assert_array_equal(result.GetTopN(4), before)
