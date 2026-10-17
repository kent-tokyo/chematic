"""Binary on-bit sources and NumPy output ownership/transactional failure."""

import itertools

import numpy as np
import pytest
from chematic.rdkit_compat import DataStructs, ExplicitBitVect


class Source:
    def __init__(self, width, ids):
        self.width, self.ids = width, ids

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return self.ids


def test_dirty_destination_is_fully_overwritten():
    fp = ExplicitBitVect(8)
    fp.SetBit(2)
    dest = np.full(8, 7, dtype=np.int8)
    result = DataStructs.ConvertToNumpyArray(fp, dest)
    assert result is dest
    assert dest.tolist() == [0, 0, 1, 0, 0, 0, 0, 0]
    DataStructs.ConvertToNumpyArray(ExplicitBitVect(8), dest)
    assert not dest.any()


def test_source_error_leaves_destination_unchanged():
    def broken():
        yield 1
        raise RuntimeError("source failed")

    dest = np.full(8, 7, dtype=np.int8)
    with pytest.raises(RuntimeError, match="source failed"):
        DataStructs.ConvertToNumpyArray(Source(8, broken()), dest)
    assert dest.tolist() == [7] * 8


@pytest.mark.parametrize(
    "dtype",
    [
        np.bool_,
        np.int8,
        np.uint8,
        np.int16,
        np.uint32,
        np.int64,
        np.float16,
        np.float32,
        np.float64,
        ">i4",
        ">f8",
    ],
)
def test_numeric_dtypes_and_noncontiguous_views(dtype):
    storage = np.full(17, 1, dtype=dtype)
    dest = storage[1::2][::-1]
    assert not dest.flags.c_contiguous
    out = DataStructs.ConvertToNumpyArray(
        Source(np.int64(8), [np.uint32(0), 3, 3, 7]), dest
    )
    assert out is dest
    assert dest.tolist() == [1, 0, 0, 1, 0, 0, 0, 1]
    assert np.all(storage[::2] == 1)
    assert dest.dtype == storage.dtype


@pytest.mark.parametrize("width", [0, 1, 7, 8, 9, 31])
def test_optional_destination_owns_its_result(width):
    ids = [] if not width else [width - 1]
    first = DataStructs.ConvertToNumpyArray(Source(width, iter(ids)))
    second = DataStructs.ConvertToNumpyArray(Source(width, ids))
    assert first.dtype == np.int8
    assert first.shape == (width,)
    assert not np.shares_memory(first, second)
    assert first.tolist() == [int(i in ids) for i in range(width)]


@pytest.mark.parametrize(
    "dest",
    [
        [],
        [0] * 8,
        bytearray(8),
        np.zeros(8, dtype=complex),
        np.zeros(8, dtype=object),
        np.zeros(8, dtype="S1"),
        np.zeros(8, dtype="datetime64[D]"),
        np.zeros(8).view(np.ma.MaskedArray),
    ],
)
def test_wrong_destination_types(dest):
    with pytest.raises(TypeError):
        DataStructs.ConvertToNumpyArray(Source(8, [2]), dest)


@pytest.mark.parametrize(
    "dest", [np.full(7, 7), np.full(9, 7), np.full((2, 4), 7), np.array(7)]
)
def test_wrong_shape_is_refused_before_source_consumption(dest):
    class Guard:
        def GetNumBits(self):
            return 8

        def GetOnBits(self):
            raise AssertionError("must not read source")

    before = dest.copy()
    with pytest.raises(ValueError):
        DataStructs.ConvertToNumpyArray(Guard(), dest)
    assert np.array_equal(dest, before)


def test_readonly_and_overlapping_views_are_refused():
    base = np.full(16, 7, dtype=np.uint8)
    readonly = base[:8]
    readonly.flags.writeable = False
    overlap = np.lib.stride_tricks.as_strided(
        base, shape=(8,), strides=(0,), writeable=True
    )
    subitem = np.lib.stride_tricks.as_strided(
        base.view(np.uint32), shape=(8,), strides=(1,), writeable=True
    )
    for dest in [readonly, overlap, subitem]:
        with pytest.raises(ValueError):
            DataStructs.ConvertToNumpyArray(Source(8, [1]), dest)
        assert base.tolist() == [7] * 16


@pytest.mark.parametrize("index", [True, np.bool_(False), 1.5, "1", None])
def test_wrong_bit_index_does_not_mutate_output(index):
    dest = np.full(8, 7)
    with pytest.raises(TypeError):
        DataStructs.ConvertToNumpyArray(Source(8, [1, index]), dest)
    assert dest.tolist() == [7] * 8


@pytest.mark.parametrize("index", [-1, 8, 2**200])
def test_out_of_range_bit_does_not_mutate_output(index):
    dest = np.full(8, 7)
    with pytest.raises(ValueError):
        DataStructs.ConvertToNumpyArray(Source(8, [1, index]), dest)
    assert dest.tolist() == [7] * 8


@pytest.mark.parametrize(
    "width,error",
    [
        (-1, ValueError),
        (1048577, ValueError),
        (True, TypeError),
        (np.bool_(False), TypeError),
        (1.5, TypeError),
        ("8", TypeError),
    ],
)
def test_width_validation(width, error):
    with pytest.raises(error):
        DataStructs.ConvertToNumpyArray(Source(width, []))


@pytest.mark.parametrize(
    "source",
    [
        None,
        1,
        object(),
        type("Missing", (), {"GetNumBits": lambda self: 8})(),
        type("NotCallable", (), {"GetNumBits": 8, "GetOnBits": []})(),
    ],
)
def test_missing_source_protocol(source):
    with pytest.raises(TypeError):
        DataStructs.ConvertToNumpyArray(source)


def test_raw_source_limit_and_exact_width_cap():
    maximum = 1048576
    result = DataStructs.ConvertToNumpyArray(Source(maximum, [0, maximum - 1]))
    assert result.size == maximum and result.sum() == 2
    assert DataStructs.ConvertToNumpyArray(
        Source(1, itertools.repeat(0, maximum))
    ).tolist() == [1]
    seen = 0

    def unbounded():
        nonlocal seen
        while True:
            seen += 1
            yield 0

    dest = np.full(1, 7)
    with pytest.raises(ValueError, match="input bit IDs"):
        DataStructs.ConvertToNumpyArray(Source(1, unbounded()), dest)
    assert seen == maximum + 1
    assert dest.tolist() == [7]


def test_source_lifetime_covers_its_lazy_view():
    import weakref

    class Ephemeral:
        def GetNumBits(self):
            return 8

        def GetOnBits(self):
            owner = weakref.ref(self)

            def bits():
                assert owner() is not None
                yield 2

            return bits()

    assert DataStructs.ConvertToNumpyArray(Ephemeral()).tolist() == [
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
    ]


@pytest.mark.parametrize("method", ["GetNumBits", "GetOnBits"])
def test_source_call_errors_are_preserved(method):
    class Broken(Source):
        pass

    def fail(self):
        raise RuntimeError("source call failed")

    setattr(Broken, method, fail)
    dest = np.full(8, 7)
    with pytest.raises(RuntimeError, match="source call failed"):
        DataStructs.ConvertToNumpyArray(Broken(8, [1]), dest)
    assert dest.tolist() == [7] * 8


def test_source_ids_can_alias_the_destination():
    dest = np.array([1, 0, 1, 0], dtype=np.int8)
    DataStructs.ConvertToNumpyArray(Source(4, dest[:2]), dest)
    assert dest.tolist() == [1, 1, 0, 0]
