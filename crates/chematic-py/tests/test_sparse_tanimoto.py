"""Sparse binary-set matrices: positions, ownership and bounded input."""

import itertools
import pytest
import chematic


@pytest.mark.parametrize(
    "metric,expected",
    [("similarity", [0, 0, 0.5, 0, 0, 0]), ("distance", [1, 1, 0.5, 1, 1, 1])],
)
def test_lower_triangle(metric, expected):
    assert (
        chematic.sparse_tanimoto_matrix_checked([[], [1], [1, 2], [3]], metric=metric)
        == expected
    )


@pytest.mark.parametrize(
    "rows,similarity,distance",
    [
        ([], [], []),
        ([[7]], [], []),
        ([[], []], [1], [0]),
        ([[], [], [7]], [1, 0, 0], [0, 1, 1]),
    ],
)
def test_empty_policy(rows, similarity, distance):
    assert (
        chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity") == similarity
    )
    assert chematic.sparse_tanimoto_matrix_checked(rows) == distance


def test_duplicates_large_ids_and_source_ownership():
    rows = [[2**32 - 1, 0, 2**32 - 1], [2**32 - 1], [2**32 - 2]]
    before = [row[:] for row in rows]
    result = chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity")
    assert result == [0.5, 0, 0]
    result[0] = -1
    assert rows == before
    assert chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity") == [
        0.5,
        0,
        0,
    ]


@pytest.mark.parametrize("value", [None, 1.5, "1", True, 1j])
def test_wrong_index_types(value):
    with pytest.raises(TypeError):
        chematic.sparse_tanimoto_matrix_checked([[value]])


@pytest.mark.parametrize("value", [-1, 2**32, 2**200])
def test_u32_range(value):
    with pytest.raises(OverflowError):
        chematic.sparse_tanimoto_matrix_checked([[value]])


@pytest.mark.parametrize(
    "rows", [None, 1, "12", b"12", bytearray(b"12"), [1], ["12"], [b"12"]]
)
def test_wrong_containers(rows):
    with pytest.raises(TypeError):
        chematic.sparse_tanimoto_matrix_checked(rows)


@pytest.mark.parametrize("metric", ["", "Similarity", "other"])
def test_invalid_metric_does_not_consume_input(metric):
    def source():
        raise AssertionError("must not consume")
        yield []

    with pytest.raises(ValueError, match="metric"):
        chematic.sparse_tanimoto_matrix_checked(source(), metric=metric)


@pytest.mark.parametrize("metric", [None, 1])
def test_metric_type(metric):
    with pytest.raises(TypeError):
        chematic.sparse_tanimoto_matrix_checked([], metric=metric)


@pytest.mark.parametrize(
    "rows,message",
    [
        ([[]] * 1025, "1024 fingerprints"),
        ([[7] * 65537], "65536 input"),
        ([[7] * 65536] * 16, "1000000 total"),
        ([list(range(64))] * 1024, "pair-bit"),
    ],
)
def test_limits(rows, message):
    with pytest.raises(ValueError, match=message):
        chematic.sparse_tanimoto_matrix_checked(rows)


def test_generators_are_bounded():
    seen = []

    def outer():
        for i in itertools.count():
            seen.append(i)
            yield []

    with pytest.raises(ValueError, match="1024 fingerprints"):
        chematic.sparse_tanimoto_matrix_checked(outer())
    assert len(seen) == 1025
    seen.clear()

    def inner():
        for i in itertools.count():
            seen.append(i)
            yield 7

    with pytest.raises(ValueError, match="65536 input"):
        chematic.sparse_tanimoto_matrix_checked([inner()])
    assert len(seen) == 65537
    assert chematic.sparse_tanimoto_matrix_checked(
        (iter(row) for row in [[], [1], [1, 2]])
    ) == [1, 1, 0.5]


@pytest.mark.parametrize("nested", [False, True])
def test_original_iterator_errors(nested):
    def broken():
        yield 7 if nested else []
        raise RuntimeError("bit source failed")

    with pytest.raises(RuntimeError, match="bit source failed"):
        chematic.sparse_tanimoto_matrix_checked([broken()] if nested else broken())


def test_numpy_integer_positions():
    np = pytest.importorskip("numpy")
    rows = np.array([[1, 9, 2, 9], [2, 9, 3, 9]], dtype=np.uint32)[:, ::2]
    assert chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity") == [1 / 3]
    assert chematic.sparse_tanimoto_matrix_checked(rows.astype(np.int64)) == [1 - 1 / 3]


def test_getonbits_sources_preserve_owner_lifetime():
    import weakref

    refs = []

    class Source:
        def __init__(self, ids):
            self.ids = ids
            refs.append(weakref.ref(self))

        def GetOnBits(self):
            reference = weakref.ref(self)

            def values():
                owner = reference()
                if owner is None:
                    raise RuntimeError("source freed before reading")
                yield from owner.ids

            return values()

    rows = (Source(ids) for ids in [[1], [1, 2]])
    assert chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity") == [0.5]
    assert all(ref() is None for ref in refs)


@pytest.mark.parametrize("getter", [False, True])
def test_getonbits_exceptions_are_preserved(getter):
    class Source:
        @property
        def GetOnBits(self):
            if getter:
                raise RuntimeError("on-bit getter failed")

            def call():
                raise RuntimeError("on-bit call failed")

            return call

    with pytest.raises(RuntimeError, match="on-bit"):
        chematic.sparse_tanimoto_matrix_checked([Source()])


def test_noncallable_getonbits_is_refused():
    class Source:
        GetOnBits = 7

    with pytest.raises(TypeError):
        chematic.sparse_tanimoto_matrix_checked([Source()])
