"""Independent binary pair counts and transactional InfoTheory state checks."""

import gc
import itertools
import weakref

import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I


class Bits:
    def __init__(self, width, bits):
        self.width, self.bits = width, bits

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return iter(self.bits)


def configured(bits=(5, 1, 3)):
    generator = I.BitCorrMatGenerator()
    generator.SetBitList(bits)
    return generator


def test_asymmetric_pair_order_and_unnormalized_counts():
    generator = configured()
    for ids in [[1, 3], [1, 5], [1, 3, 5], [3]]:
        generator.CollectVotes(Bits(8, ids))
    assert generator.GetCorrMatrix().tolist() == [2, 1, 2]
    assert generator.GetCorrMatrix().dtype == np.float64
    generator.CollectVotes(Bits(8, [1, 3, 5]))
    assert generator.GetCorrMatrix().tolist() == [3, 2, 3]


def test_all_small_binary_patterns_bit_permutations_and_unselected_bits():
    votes = []
    for mask in range(32):
        bits = {i for i in range(5) if mask & (1 << i)}
        votes.extend([bits] * (mask % 4 + 1))
    for selected in itertools.permutations([4, 1, 3]):
        generator = configured(iter(selected))
        for vote in reversed(votes):
            generator.CollectVotes(Bits(5, (i for i in reversed(sorted(vote)))))
        expected = [
            sum(selected[i] in vote and selected[j] in vote for vote in votes)
            for i in range(3)
            for j in range(i)
        ]
        assert generator.GetCorrMatrix().tolist() == expected


def test_explicit_numpy_and_high_unsigned_sparse_protocols():
    generator = configured(np.array([5, 1, 3], dtype=np.uint32))
    vector = ExplicitBitVect(8)
    for bit in [1, 3, 5]:
        vector.SetBit(bit)
    generator.CollectVotes(vector)
    generator.CollectVotes(Bits(np.uint32(8), np.array([5, 3, 1], dtype=np.uint32)))
    assert generator.GetCorrMatrix().tolist() == [2, 2, 2]
    sparse = configured([0, 2**32 - 2])
    sparse.CollectVotes(Bits(2**32 - 1, [np.uint32(0), np.uint32(2**32 - 2)]))
    assert sparse.GetCorrMatrix().tolist() == [1]


def test_configuration_required_and_zero_votes_are_exact():
    generator = I.BitCorrMatGenerator()
    with pytest.raises(ValueError, match="SetBitList"):
        generator.GetCorrMatrix()
    with pytest.raises(ValueError, match="SetBitList"):
        generator.CollectVotes(Bits(8, []))
    generator.SetBitList([0, 1, 2, 3])
    assert generator.GetCorrMatrix().tolist() == [0] * 6
    generator.CollectVotes(Bits(8, []))
    generator.CollectVotes(Bits(8, [7]))
    assert generator.GetCorrMatrix().tolist() == [0] * 6


def test_successful_reset_changes_bit_order_width_and_clears_counts():
    generator = configured()
    generator.CollectVotes(Bits(8, [1, 3, 5]))
    generator.SetBitList([7, 3])
    assert generator.GetCorrMatrix().tolist() == [0]
    generator.CollectVotes(Bits(16, [7, 3]))
    assert generator.GetCorrMatrix().tolist() == [1]


@pytest.mark.parametrize(
    "bits,error",
    [
        (None, TypeError),
        ("bits", TypeError),
        ([], ValueError),
        ([1], ValueError),
        ([1, 1], ValueError),
        ([0, -1], ValueError),
        ([0, 2**32 - 1], ValueError),
        ([True, 1], TypeError),
        ([np.bool_(False), 1], TypeError),
        ([1.0, 2], TypeError),
        (itertools.count(), ValueError),
    ],
)
def test_invalid_setter_preserves_configuration_and_votes(bits, error):
    generator = configured()
    generator.CollectVotes(Bits(8, [1, 3, 5]))
    before = generator.GetCorrMatrix()
    with pytest.raises(error):
        generator.SetBitList(bits)
    np.testing.assert_array_equal(generator.GetCorrMatrix(), before)
    with pytest.raises(ValueError, match="consistent"):
        generator.CollectVotes(Bits(16, [1, 3, 5]))
    generator.CollectVotes(Bits(8, [1, 3, 5]))
    assert generator.GetCorrMatrix().tolist() == [2, 2, 2]


@pytest.mark.parametrize(
    "source,error",
    [
        (None, TypeError),
        ([1, 3], TypeError),
        (Bits(True, [1, 3]), TypeError),
        (Bits(np.bool_(True), [1, 3]), TypeError),
        (Bits(8.0, [1, 3]), TypeError),
        (Bits(0, []), ValueError),
        (Bits(2**32, []), ValueError),
        (Bits(4, [1, 3]), ValueError),
        (Bits(8, [1, 1]), ValueError),
        (Bits(8, [1, 9]), ValueError),
        (Bits(8, [1, -1]), ValueError),
        (Bits(8, [1, True]), TypeError),
        (Bits(8, [1, np.bool_(False)]), TypeError),
        (Bits(8, [1, 3.0]), TypeError),
        (Bits(8, "13"), TypeError),
    ],
)
def test_invalid_vote_preserves_counts_and_first_width_binding(source, error):
    generator = configured()
    with pytest.raises(error):
        generator.CollectVotes(source)
    assert generator.GetCorrMatrix().tolist() == [0, 0, 0]
    generator.CollectVotes(Bits(16, [1, 3, 5]))
    assert generator.GetCorrMatrix().tolist() == [1, 1, 1]


def test_source_and_configuration_errors_propagate_without_partial_updates():
    generator = configured()
    generator.CollectVotes(Bits(8, [1, 3, 5]))
    error = RuntimeError("interrupted")

    def interrupted():
        yield 1
        raise error

    class Broken(Bits):
        def GetNumBits(self):
            raise error

    for action in [
        lambda: generator.SetBitList(interrupted()),
        lambda: generator.CollectVotes(Bits(8, interrupted())),
        lambda: generator.CollectVotes(Broken(8, [])),
    ]:
        with pytest.raises(RuntimeError) as caught:
            action()
        assert caught.value is error
        assert generator.GetCorrMatrix().tolist() == [1, 1, 1]


def test_source_and_output_ownership():
    selected = [5, 1, 3]
    positions = [1, 3, 5]
    generator = configured(selected)
    vector = Bits(8, positions)
    ref = weakref.ref(vector)
    generator.CollectVotes(vector)
    selected[:] = [0, 2]
    positions[:] = []
    del vector
    gc.collect()
    assert ref() is None
    a, b = generator.GetCorrMatrix(), generator.GetCorrMatrix()
    assert a.flags.owndata and b.flags.owndata
    a[:] = 999
    assert b.tolist() == generator.GetCorrMatrix().tolist() == [1, 1, 1]


def test_lazy_view_keeps_owner_alive_during_iteration():
    class Borrowed:
        def GetNumBits(self):
            return 8

        def GetOnBits(self):
            ref = weakref.ref(self)

            def bits():
                gc.collect()
                assert ref() is not None
                yield 1
                yield 3

            return bits()

    generator = configured()
    generator.CollectVotes(Borrowed())
    assert generator.GetCorrMatrix().tolist() == [0, 0, 1]


def test_actual_selected_and_raw_id_limits():
    generator = configured(range(1024))
    assert generator.GetCorrMatrix().shape == (523776,)
    with pytest.raises(ValueError, match="1024"):
        generator.SetBitList(range(1025))
    small = configured([0, 1])
    small.CollectVotes(Bits(65536, range(65536)))
    assert small.GetCorrMatrix().tolist() == [1]
    with pytest.raises(ValueError, match="65536"):
        small.CollectVotes(Bits(65536, range(65537)))
    assert small.GetCorrMatrix().tolist() == [1]


def test_actual_pair_work_limit_and_failed_vote_atomicity():
    generator = configured(range(1024))
    for _ in range(95):
        generator.CollectVotes(Bits(1024, range(1024)))
    assert np.all(generator.GetCorrMatrix() == 95)
    with pytest.raises(ValueError, match="50000000"):
        generator.CollectVotes(Bits(1024, range(1024)))
    assert np.all(generator.GetCorrMatrix() == 95)
    generator.SetBitList([0, 1])
    generator.CollectVotes(Bits(2, [0, 1]))
    assert generator.GetCorrMatrix().tolist() == [1]


def test_record_limit_and_reset_with_small_counter_budget(monkeypatch):
    import chematic._bit_correlation as impl

    monkeypatch.setattr(impl, "_MAX_RECORDS", 3)
    generator = configured()
    for _ in range(3):
        generator.CollectVotes(Bits(8, []))
    with pytest.raises(ValueError, match="record limit"):
        generator.CollectVotes(Bits(8, [1, 3, 5]))
    assert generator.GetCorrMatrix().tolist() == [0, 0, 0]
    generator.SetBitList([1, 3])
    generator.CollectVotes(Bits(8, [1, 3]))
    assert generator.GetCorrMatrix().tolist() == [1]
