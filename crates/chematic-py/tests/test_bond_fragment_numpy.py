"""Caller-owned NumPy cut counts: layout, isolation and failure boundaries."""

import numpy as np
import pytest
import chematic as ch


@pytest.mark.parametrize("layout", ["plain", "stride", "reverse", "big_endian"])
def test_numpy_counts_replace_in_logical_order(layout):
    backing = np.full(8, 77, dtype=np.uint32)
    counts = {
        "plain": backing[:4],
        "stride": backing[::2],
        "reverse": backing[6::-2],
        "big_endian": np.full(4, 77, dtype=">u4"),
    }[layout]
    source = ch.from_smiles("CC(C)O")
    before = source.smiles
    result = ch.fragment_on_bonds_checked(source, [0, 1], cuts_per_atom=counts)
    assert counts.tolist() == [1, 2, 1, 0]
    assert result.cuts_per_atom == (1, 2, 1, 0)
    assert source.smiles == before
    if layout == "stride":
        assert backing[1::2].tolist() == [77] * 4
    if layout == "reverse":
        assert backing[1::2].tolist() == [77] * 4
        assert backing[::2].tolist() == [0, 1, 2, 1]
    counts[:] = 42
    assert result.cuts_per_atom == (1, 2, 1, 0)


@pytest.mark.parametrize(
    "dtype", [np.int32, np.uint8, np.uint64, np.float32, bool, object]
)
def test_wrong_dtype_leaves_storage_unchanged(dtype):
    counts = np.full(3, 1, dtype=dtype)
    before = counts.copy()
    with pytest.raises(TypeError, match="uint32"):
        ch.fragment_on_bonds_checked(ch.from_smiles("CCO"), [0], cuts_per_atom=counts)
    np.testing.assert_array_equal(counts, before)


@pytest.mark.parametrize("shape", [(), (1, 3), (0,), (2,), (4,)])
def test_wrong_shape_leaves_storage_unchanged(shape):
    counts = np.full(shape, 77, dtype=np.uint32)
    with pytest.raises((TypeError, ValueError)):
        ch.fragment_on_bonds_checked(ch.from_smiles("CCO"), [0], cuts_per_atom=counts)
    assert np.all(counts == 77)


def test_read_only_and_overlapping_views_are_refused():
    source = ch.from_smiles("CCO")
    counts = np.full(3, 77, dtype=np.uint32)
    counts.flags.writeable = False
    with pytest.raises(ValueError, match="writable"):
        ch.fragment_on_bonds_checked(source, [0], cuts_per_atom=counts)
    storage = np.full(12, 77, dtype=np.uint8)
    for stride in (0, 1, 2, 3):
        overlap = np.ndarray((3,), dtype=np.uint32, buffer=storage, strides=(stride,))
        with pytest.raises(ValueError, match="overlap"):
            ch.fragment_on_bonds_checked(source, [0], cuts_per_atom=overlap)
    assert np.all(storage == 77)


def test_ndarray_subclass_cannot_intercept_writes():
    class Trap(np.ndarray):
        def __setitem__(self, key, value):
            pytest.fail("user-defined write must not run")

    counts = np.full(3, 77, dtype=np.uint32).view(Trap)
    with pytest.raises(TypeError, match="plain"):
        ch.fragment_on_bonds_checked(ch.from_smiles("CCO"), [0], cuts_per_atom=counts)


@pytest.mark.parametrize(
    "smiles,cuts",
    [
        ("CCO", [0, 99]),
        ("CCO", [0, 0]),
        ("C1CC1", [0]),
        ("C=C", [0]),
        ("N[C@H](C)O", [0]),
    ],
)
def test_native_refusal_preserves_array_and_source(smiles, cuts):
    source = ch.from_smiles(smiles)
    before = source.smiles
    counts = np.full(source.atom_count, 77, dtype=np.uint32)
    with pytest.raises(ValueError):
        ch.fragment_on_bonds_checked(source, cuts, cuts_per_atom=counts)
    assert np.all(counts == 77)
    assert source.smiles == before


def test_empty_cuts_zero_counts_and_empty_molecule():
    counts = np.full(3, 2**32 - 1, dtype=np.uint32)
    ch.fragment_on_bonds_checked(ch.from_smiles("CCO"), [], cuts_per_atom=counts)
    assert counts.tolist() == [0, 0, 0]
    result = ch.fragment_on_bonds_checked(
        ch.RWMol().GetMol(), [], cuts_per_atom=np.empty(0, dtype=np.uint32)
    )
    assert result.cuts_per_atom == ()


def test_validation_precedes_native_call(monkeypatch):
    def forbidden(*args):
        pytest.fail("native fragmentation must not start")

    monkeypatch.setattr(ch.Mol, "_fragment_single_bonds_checked", forbidden)
    with pytest.raises(ValueError):
        ch.fragment_on_bonds_checked(
            ch.from_smiles("CCO"), [0], cuts_per_atom=np.zeros(2, dtype=np.uint32)
        )


def test_conversion_failure_preserves_destination(monkeypatch):
    import chematic._bond_fragment as implementation

    source = ch.from_smiles("CCO")
    before = source.smiles
    counts = np.full(3, 77, dtype=np.uint32)
    error = MemoryError("staging failed")

    def fail_conversion(*args, **kwargs):
        raise error

    with monkeypatch.context() as patch:
        patch.setattr(implementation.np, "asarray", fail_conversion)
        with pytest.raises(MemoryError) as caught:
            ch.fragment_on_bonds_checked(source, [0], cuts_per_atom=counts)
        assert caught.value is error
    assert counts.tolist() == [77, 77, 77]
    assert source.smiles == before
