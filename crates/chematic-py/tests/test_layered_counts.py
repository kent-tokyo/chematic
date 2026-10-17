"""Typed in-place Layered counters; original inputs, native core tested separately."""

import itertools

import numpy as np
import pytest

from chematic import rdkit_compat as C


@pytest.fixture
def mol():
    return C.MolFromSmiles("CCO")


def test_accumulates_and_keeps_trailing_entries_and_list_identity(mol):
    counts = [10, 20, 30, 99]
    original_id = id(counts)
    fp = C.LayeredFingerprint(mol, atomCounts=counts)
    assert counts == [12, 23, 32, 99] and id(counts) == original_id
    assert fp.GetNumBits() == 2048 and fp.GetNumOnBits() == 11
    again = C.LayeredFingerprint(mol, atomCounts=counts)
    assert counts == [14, 26, 34, 99]
    assert again.ToBitString() == fp.ToBitString()


@pytest.mark.parametrize("dtype", ["uint32", ">u4", "<u4"])
@pytest.mark.parametrize("view", ["plain", "stride", "reverse"])
def test_numpy_views_update_only_destination(mol, dtype, view):
    storage = np.full(9, 99, dtype=dtype)
    counts = (
        storage[:4]
        if view == "plain"
        else storage[::2] if view == "stride" else storage[6:2:-1]
    )
    counts[:3] = [10, 20, 30]
    before = storage.copy()
    result = C.LayeredFingerprint(mol, atomCounts=counts)
    assert result.GetNumBits() == 2048
    assert counts.tolist() == [12, 23, 32] + [99] * (len(counts) - 3)
    expected = before.copy()
    target = (
        expected[:4]
        if view == "plain"
        else expected[::2] if view == "stride" else expected[6:2:-1]
    )
    target[:3] = [12, 23, 32]
    np.testing.assert_array_equal(storage, expected)


@pytest.mark.parametrize("array", [False, True])
def test_wrapping_counts(mol, array):
    counts = np.full(4, 2**32 - 1, dtype=np.uint32) if array else [2**32 - 1] * 4
    C.LayeredFingerprint(mol, atomCounts=counts)
    assert list(counts) == [1, 2, 1, 2**32 - 1]


@pytest.mark.parametrize("mask_on", [False, True])
def test_mask_changes_counts_not_just_output(mol, mask_on):
    mask = C.ExplicitBitVect(1)
    if mask_on:
        mask.SetBit(0)
    counts = [10, 20, 30]
    fp = C.LayeredFingerprint(mol, fpSize=1, atomCounts=counts, setOnlyBits=mask)
    assert counts == ([12, 23, 32] if mask_on else [10, 20, 30])
    assert fp.GetNumOnBits() == int(mask_on)
    assert mask.GetNumOnBits() == int(mask_on)


def test_no_layers_no_count_increase_and_numpy_integer_options(mol):
    counts = [1, 2, 3]
    fp = C.LayeredFingerprint(
        mol, layerFlags=np.uint32(0x80000000), fpSize=np.int32(9), atomCounts=counts
    )
    assert counts == [1, 2, 3] and fp.GetNumOnBits() == 0 and fp.GetNumBits() == 9


@pytest.mark.parametrize(
    "kwargs",
    [
        {"fpSize": 0},
        {"fpSize": 65537},
        {"minPath": 0},
        {"maxPath": 13},
        {"minPath": 3, "maxPath": 2},
        {"layerFlags": -1},
        {"layerFlags": 2**32},
        {"maxVisits": 0},
        {"maxPaths": 0},
        {"maxVisits": 1},
        {"maxPaths": 1},
        {"branchedPaths": False},
        {"fromAtoms": [3]},
        {"fromAtoms": [-1]},
        {"fromAtoms": [0] * 4097},
        {"fromAtoms": itertools.repeat(0)},
        {"fpSize": 9, "setOnlyBits": C.ExplicitBitVect(10)},
    ],
)
def test_error_paths_preserve_seed_counters(mol, kwargs):
    counts = [10, 20, 30]
    with pytest.raises((ValueError, TypeError)):
        C.LayeredFingerprint(mol, atomCounts=counts, **kwargs)
    assert counts == [10, 20, 30]


@pytest.mark.parametrize(
    "key", ["layerFlags", "minPath", "maxPath", "fpSize", "maxVisits", "maxPaths"]
)
@pytest.mark.parametrize("value", [True, np.bool_(True), 1.5, "2"])
def test_scalar_types_checked_before_mutation(mol, key, value):
    counts = [10, 20, 30]
    with pytest.raises(TypeError):
        C.LayeredFingerprint(mol, atomCounts=counts, **{key: value})
    assert counts == [10, 20, 30]


@pytest.mark.parametrize(
    "counts",
    [[], [0], [0, 0, True], [0, 0, -1], [0, 0, 2**32], [0, 0, 1.5], [0, 0, "2"]],
)
def test_invalid_list_entry_leaves_original_unchanged(mol, counts):
    before = counts.copy()
    with pytest.raises((ValueError, TypeError)):
        C.LayeredFingerprint(mol, atomCounts=counts)
    assert counts == before


@pytest.mark.parametrize(
    "counts",
    [
        (0, 0, 0),
        np.zeros(3, dtype=np.uint64),
        np.zeros(3, dtype=np.int32),
        np.zeros(3, dtype=np.float64),
        np.zeros(3, dtype=bool),
        np.zeros((3, 1), dtype=np.uint32),
        np.ma.array([0, 0, 0], dtype=np.uint32),
    ],
)
def test_bad_destination_type(mol, counts):
    with pytest.raises(TypeError):
        C.LayeredFingerprint(mol, atomCounts=counts)


@pytest.mark.parametrize("stride", [0, 1, 2])
def test_overlapping_uint32_array_rejected(mol, stride):
    original = np.zeros(4, dtype=np.uint32)
    view = np.lib.stride_tricks.as_strided(
        original, shape=(3,), strides=(stride,), writeable=True
    )
    with pytest.raises(ValueError, match="overlap"):
        C.LayeredFingerprint(mol, atomCounts=view)
    assert not np.any(original)


def test_readonly_and_oversized_buffers(mol):
    counts = np.zeros(3, dtype=np.uint32)
    counts.flags.writeable = False
    with pytest.raises(ValueError, match="writable"):
        C.LayeredFingerprint(mol, atomCounts=counts)
    counts = [0] * 65537
    with pytest.raises(ValueError, match="65536"):
        C.LayeredFingerprint(mol, atomCounts=counts)
    assert not any(counts)


def test_destination_subclass_callbacks_never_run(mol):
    class CallbackList(list):
        def __setitem__(self, key, value):
            raise AssertionError("unexpected mutation")

    with pytest.raises(TypeError, match="plain"):
        C.LayeredFingerprint(mol, atomCounts=CallbackList([0, 0, 0]))


def test_mask_protocol_exceptions_and_duplicates_do_not_write(mol):
    class Mask:
        def GetNumBits(self):
            return 2048

        def GetOnBits(self):
            yield 1
            raise RuntimeError("source failed")

    counts = [10, 20, 30]
    with pytest.raises(RuntimeError, match="source failed"):
        C.LayeredFingerprint(mol, atomCounts=counts, setOnlyBits=Mask())
    assert counts == [10, 20, 30]

    class Duplicate(Mask):
        def GetOnBits(self):
            return itertools.repeat(1)

    with pytest.raises(ValueError, match="duplicate"):
        C.LayeredFingerprint(mol, atomCounts=counts, setOnlyBits=Duplicate())
    assert counts == [10, 20, 30]


def test_fingerprint_construction_finishes_before_counter_write(mol, monkeypatch):
    counts = np.array([10, 20, 30], dtype=np.uint32)

    def fail(*args):
        raise MemoryError("allocation failed")

    monkeypatch.setattr(C.ExplicitBitVect, "_from_bytes", fail)
    with pytest.raises(MemoryError, match="allocation failed"):
        C.LayeredFingerprint(mol, atomCounts=counts)
    np.testing.assert_array_equal(counts, [10, 20, 30])


def test_roots_duplicates_match_existing_core(mol):
    for roots in ([0], [0, 0], (0, 1, 2)):
        for branched in (False, True):
            seeds = [2, 3, 4, 99]
            expected = mol._mol.rdkit_layered_config(
                from_atoms=roots, branched_paths=branched, atom_counts=seeds
            )
            result = C.LayeredFingerprint(
                mol, fromAtoms=roots, branchedPaths=branched, atomCounts=seeds
            )
            assert (
                result._to_bytes() == expected["bits"]
                and seeds == expected["atom_counts"]
            )


def test_unsupported_graph_and_graph_cap_do_not_mutate():
    for smiles in ("N->[Cu+2]", ".".join(["C"] * 513)):
        mol = C.MolFromSmiles(smiles)
        counts = [7] * mol.GetNumAtoms()
        with pytest.raises(ValueError):
            C.LayeredFingerprint(mol, atomCounts=counts)
        assert counts == [7] * len(counts)


def test_no_counts_and_empty_molecule():
    assert C.LayeredFingerprint(C.MolFromSmiles("CCO")).GetNumOnBits() == 11
    assert C.LayeredFingerprint(C.RWMol().GetMol(), atomCounts=[]).GetNumOnBits() == 0
    with pytest.raises(TypeError, match="mol"):
        C.LayeredFingerprint(None)


def test_boolean_and_root_type_checks(mol):
    for kwargs in ({"branchedPaths": 1}, {"fromAtoms": "0"}, {"fromAtoms": [True]}):
        with pytest.raises(TypeError):
            C.LayeredFingerprint(mol, **kwargs)
