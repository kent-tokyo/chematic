"""Public metric facade: native reuse, ownership and binary source contracts."""

from dataclasses import dataclass
import gc
import itertools
import math
import weakref

import numpy as np
import pytest

from chematic.rdkit_compat import ExplicitBitVect, rdMetricMatrixCalc as M


@dataclass
class Bits:
    width: object
    positions: object

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return iter(self.positions)


def test_public_euclidean_keyword_pair_order_and_owned_float64_result():
    rows = [[0, 0], [3, 4], [3, 0], [-1, -2]]
    want = [math.dist(rows[i], rows[j]) for i in range(4) for j in range(i)]
    for source in [
        rows,
        tuple(map(tuple, rows)),
        np.array(rows),
        (iter(r) for r in rows),
    ]:
        actual = M.GetEuclideanDistMat(descripMat=source)
        np.testing.assert_allclose(actual, want, rtol=1e-15, atol=0)
        assert actual.dtype == np.float64 and actual.shape == (6,)
        assert actual.flags.owndata and actual.flags.c_contiguous
    actual[:] = 99
    assert rows == [[0, 0], [3, 4], [3, 0], [-1, -2]]


def test_numpy_strided_descriptor_views_and_fresh_output():
    source = np.arange(30, dtype=np.float64).reshape(5, 6)[::-1, ::2]
    before = source.copy()
    result = M.GetEuclideanDistMat(source)
    want = [math.dist(source[i], source[j]) for i in range(5) for j in range(i)]
    np.testing.assert_allclose(result, want)
    result[:] = -1
    np.testing.assert_array_equal(source, before)
    assert np.all(M.GetEuclideanDistMat(source) >= 0)


@pytest.mark.parametrize("rows", [[], [[3, 4]], [[], []]])
def test_native_empty_singleton_and_zero_features(rows):
    result = M.GetEuclideanDistMat(rows)
    assert result.shape == ((1,) if len(rows) == 2 else (0,))
    assert np.all(result == 0)


@pytest.mark.parametrize("scale", [1e200, 1e-200, math.ulp(0.0)])
def test_native_scaled_norms_preserve_representable_extremes(scale):
    assert M.GetEuclideanDistMat([[0], [scale]]).tolist() == [scale]


@pytest.mark.parametrize(
    "rows,error",
    [
        ("text", TypeError),
        (["12", "34"], TypeError),
        ([1, 2], TypeError),
        ([[None]], TypeError),
        ([[1], [1, 2]], ValueError),
        ([[math.nan]], ValueError),
        ([[math.inf]], ValueError),
        ([[-1e308], [1e308]], ValueError),
        ([[0] * 257], ValueError),
        (itertools.repeat([]), ValueError),
    ],
)
def test_euclidean_forwards_native_input_and_work_limits(rows, error):
    with pytest.raises(error):
        M.GetEuclideanDistMat(rows)


def test_binary_methods_match_exact_sets_with_all_small_binary_inputs():
    sets = [set(i for i in range(4) if mask & (1 << i)) for mask in range(16)]
    # Include a repeated empty source to verify the facade's RDKit convention.
    sets += [set(), {0, 1}]
    sources = [Bits(4, list(reversed(sorted(s))) + list(s)) for s in sets]
    want = []
    for i in range(len(sets)):
        for j in range(i):
            union = sets[i] | sets[j]
            want.append(len(sets[i] & sets[j]) / len(union) if union else 0.0)
    np.testing.assert_array_equal(M.GetTanimotoSimMat(bitVectList=iter(sources)), want)
    np.testing.assert_array_equal(
        M.GetTanimotoDistMat(tuple(sources)), 1 - np.array(want)
    )


def test_chematic_explicit_vectors_are_structural_binary_sources():
    first, second = ExplicitBitVect(17), ExplicitBitVect(17)
    first.SetBit(1)
    first.SetBit(16)
    second.SetBit(16)
    assert M.GetTanimotoSimMat([first, second]).tolist() == [0.5]
    assert M.GetTanimotoDistMat([first, second]).tolist() == [0.5]
    assert first.GetOnBits() == [1, 16]


def test_high_unsigned_sparse_ids_do_not_expand_dense_storage():
    rows = [Bits(2**32 - 1, [2**32 - 2, np.uint32(7)]), Bits(2**32 - 1, [2**32 - 2])]
    assert M.GetTanimotoSimMat(rows).tolist() == [0.5]
    assert M.GetTanimotoDistMat(rows).tolist() == [0.5]


@pytest.mark.parametrize("rows", [[], [Bits(0, [])], [Bits(0, []), Bits(0, [])]])
def test_empty_vectors_and_small_rankings(rows):
    sim, dist = M.GetTanimotoSimMat(rows), M.GetTanimotoDistMat(rows)
    assert sim.dtype == dist.dtype == np.float64
    assert sim.shape == dist.shape == ((1,) if len(rows) == 2 else (0,))
    assert np.all(sim == 0) and np.all(dist == 1)


@pytest.mark.parametrize(
    "source,error",
    [
        (None, TypeError),
        ("rows", TypeError),
        ([[1, 2]], TypeError),
        ([Bits(True, [])], TypeError),
        ([Bits(np.bool_(True), [])], TypeError),
        ([Bits(4.0, [])], TypeError),
        ([Bits(-1, [])], ValueError),
        ([Bits(2**32, [])], ValueError),
        ([Bits(4, []), Bits(5, [])], ValueError),
        ([Bits(4, [True])], TypeError),
        ([Bits(4, [np.bool_(False)])], TypeError),
        ([Bits(4, [1.0])], TypeError),
        ([Bits(4, [-1])], OverflowError),
        ([Bits(4, [2**32])], OverflowError),
        ([Bits(4, [4])], ValueError),
        ([Bits(0, [0])], ValueError),
        ([Bits(4, "12")], TypeError),
    ],
)
def test_binary_types_dimensions_and_positions_refuse(source, error):
    for function in [M.GetTanimotoDistMat, M.GetTanimotoSimMat]:
        with pytest.raises(error):
            function(source)


def test_source_callback_iteration_errors_precede_native_calculation(monkeypatch):
    import chematic._metric_matrix as impl

    def forbidden(*args, **kwargs):
        raise AssertionError("native calculator must not start")

    monkeypatch.setattr(impl, "_tanimoto", forbidden)
    error = RuntimeError("source interrupted")

    def interrupted():
        yield 1
        raise error

    class BrokenWidth(Bits):
        def GetNumBits(self):
            raise error

    class BrokenGetter(Bits):
        def GetOnBits(self):
            raise error

    def outer():
        yield Bits(4, [1])
        raise error

    for rows in [
        [Bits(4, [1]), BrokenWidth(4, [])],
        [BrokenGetter(4, [])],
        [Bits(4, interrupted())],
        outer(),
    ]:
        with pytest.raises(RuntimeError) as caught:
            M.GetTanimotoSimMat(rows)
        assert caught.value is error


def test_lazy_sources_remain_alive_then_are_released():
    references = []

    class Borrowed:
        def __init__(self):
            references.append(weakref.ref(self))

        def GetNumBits(self):
            return 8

        def GetOnBits(self):
            owner = weakref.ref(self)

            def positions():
                gc.collect()
                assert owner() is not None
                yield 1

            return positions()

    result = M.GetTanimotoSimMat(Borrowed() for _ in range(2))
    gc.collect()
    assert all(ref() is None for ref in references)
    assert result.tolist() == [1.0]


def test_binary_results_are_owned_and_input_positions_are_unchanged():
    positions = np.array([3, 1, 1], dtype=np.uint32)
    rows = [Bits(8, positions), Bits(8, [1])]
    before = positions.copy()
    a, b = M.GetTanimotoSimMat(rows), M.GetTanimotoSimMat(rows)
    a[:] = 9
    assert b.tolist() == [0.5]
    assert a.flags.owndata and b.flags.owndata
    np.testing.assert_array_equal(positions, before)


def test_actual_binary_row_raw_id_total_and_work_boundaries():
    assert M.GetTanimotoDistMat(itertools.repeat(Bits(4, []), 1024)).size == 523776
    with pytest.raises(ValueError, match="1024 rows"):
        M.GetTanimotoDistMat(itertools.repeat(Bits(4, [])))
    assert M.GetTanimotoSimMat(
        [Bits(4, itertools.repeat(1, 65536)), Bits(4, [1])]
    ).tolist() == [1.0]
    with pytest.raises(ValueError, match="65536 raw"):
        M.GetTanimotoSimMat([Bits(4, itertools.repeat(1))])
    rows = [Bits(4, itertools.repeat(1, 65536)) for _ in range(15)]
    rows.append(Bits(4, itertools.repeat(1, 16960)))
    assert np.all(M.GetTanimotoSimMat(rows) == 1)
    with pytest.raises(ValueError, match="1000000 raw"):
        M.GetTanimotoSimMat(Bits(4, itertools.repeat(1, 65536)) for _ in range(16))
    with pytest.raises(ValueError, match="work"):
        M.GetTanimotoSimMat(Bits(64, range(50)) for _ in range(1024))


def test_public_namespace_and_source_protocol_are_exported():
    namespace = {}
    exec("from chematic.rdkit_compat import *", namespace)
    assert namespace["rdMetricMatrixCalc"] is M
    assert "MetricBitVector" in namespace
