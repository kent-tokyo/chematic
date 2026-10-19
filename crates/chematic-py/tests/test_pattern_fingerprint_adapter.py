"""Typed Pattern facade: exact folding, explicit no-ops and failure staging."""

import numpy as np
import pytest
import chematic as ch
from chematic import rdkit_compat as C


@pytest.mark.parametrize("width", [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048])
def test_power_of_two_widths_fold_native_bits(width):
    mol = C.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    full = C.PatternFingerprint(mol)
    result = C.PatternFingerprint(mol, np.int64(width))
    assert result.GetNumBits() == width
    assert set(result.GetOnBits()) == {bit % width for bit in full.GetOnBits()}
    result.SetBit(0)
    assert C.PatternFingerprint(mol).ToBitString() == full.ToBitString()


@pytest.mark.parametrize(
    "layout", ["list", "plain", "stride", "reverse", "empty", "empty_array"]
)
def test_counter_seeds_and_storage_are_unchanged(layout):
    mol = C.MolFromSmiles("CCO")
    storage = np.arange(10, dtype=np.uint32)
    counts = {
        "list": [77, 88, 99, 42],
        "plain": storage[:4],
        "stride": storage[::2],
        "reverse": storage[::-2],
        "empty": [],
        "empty_array": np.empty(0, dtype=np.uint32),
    }[layout]
    before = list(counts)
    C.PatternFingerprint(mol, atomCounts=counts)
    assert list(counts) == before
    np.testing.assert_array_equal(storage, np.arange(10, dtype=np.uint32))


@pytest.mark.parametrize("bits", [(), (0,), tuple(range(128))])
def test_validated_mask_does_not_filter_pinned_pattern_behavior(bits):
    class Mask:
        def GetNumBits(self):
            return 128

        def GetOnBits(self):
            return iter(bits)

    mol = C.MolFromSmiles("CCO")
    assert (
        C.PatternFingerprint(mol, 128, setOnlyBits=Mask()).ToBitString()
        == C.PatternFingerprint(mol, 128).ToBitString()
    )


@pytest.mark.parametrize(
    "options,error",
    [
        ({"fpSize": True}, TypeError),
        ({"fpSize": 1.0}, TypeError),
        ({"fpSize": 0}, ValueError),
        ({"fpSize": 3}, ValueError),
        ({"fpSize": 2049}, ValueError),
        ({"fpSize": 4096}, ValueError),
        ({"tautomerFingerprints": True}, ValueError),
        ({"tautomerFingerprints": 0}, TypeError),
        ({"atomCounts": [1]}, ValueError),
        ({"atomCounts": [1, 2, -1]}, ValueError),
        ({"atomCounts": [1, 2, 2**32]}, ValueError),
        ({"atomCounts": [1, 2, True]}, TypeError),
        ({"atomCounts": [1, 2, 1.0]}, TypeError),
        ({"atomCounts": (0, 0, 0)}, TypeError),
        ({"atomCounts": np.zeros(3, dtype=np.int64)}, TypeError),
        ({"atomCounts": np.zeros((1, 3), dtype=np.uint32)}, TypeError),
        ({"atomCounts": [0] * 65537}, ValueError),
        ({"setOnlyBits": object()}, TypeError),
    ],
)
def test_invalid_options_precede_native_calculation(options, error, monkeypatch):
    def forbidden(*args):
        pytest.fail("native Pattern calculation must not start")

    monkeypatch.setattr(ch.Mol, "rdkit_pattern_fp", forbidden)
    with pytest.raises(error):
        C.PatternFingerprint(C.MolFromSmiles("CCO"), **options)


def test_mask_failures_preserve_destination_and_do_not_start_native(monkeypatch):
    def forbidden(*args):
        pytest.fail("native Pattern calculation must not start")

    monkeypatch.setattr(ch.Mol, "rdkit_pattern_fp", forbidden)
    error = RuntimeError("source failed")

    class Mask:
        def GetNumBits(self):
            return 2048

        def GetOnBits(self):
            yield 2
            raise error

    counts = [77, 88, 99]
    with pytest.raises(RuntimeError) as caught:
        C.PatternFingerprint(
            C.MolFromSmiles("CCO"), atomCounts=counts, setOnlyBits=Mask()
        )
    assert caught.value is error
    assert counts == [77, 88, 99]


@pytest.mark.parametrize(
    "width,bits",
    [(1024, [0]), (2048, [0, 0]), (2048, [-1]), (2048, [2048]), (2048, [True])],
)
def test_mask_width_ids_and_duplicates(width, bits):
    class Mask:
        def GetNumBits(self):
            return width

        def GetOnBits(self):
            return bits

    with pytest.raises((TypeError, ValueError)):
        C.PatternFingerprint(C.MolFromSmiles("CCO"), setOnlyBits=Mask())


def test_native_error_preserves_counters(monkeypatch):
    error = ValueError("native failure")

    def fail(*args):
        raise error

    monkeypatch.setattr(ch.Mol, "rdkit_pattern_fp", fail)
    counts = np.array([1, 2, 3], dtype=np.uint32)
    with pytest.raises(ValueError) as caught:
        C.PatternFingerprint(C.MolFromSmiles("CCO"), atomCounts=counts)
    assert caught.value is error
    assert counts.tolist() == [1, 2, 3]


def test_readonly_overlap_and_subclass_counters_refuse():
    readonly = np.zeros(3, dtype=np.uint32)
    readonly.flags.writeable = False
    overlap = np.lib.stride_tricks.as_strided(np.zeros(3, dtype=np.uint32), (3,), (0,))

    class Array(np.ndarray):
        pass

    subclass = np.zeros(3, dtype=np.uint32).view(Array)
    for counts in (readonly, overlap, subclass):
        with pytest.raises((TypeError, ValueError)):
            C.PatternFingerprint(C.MolFromSmiles("CCO"), atomCounts=counts)


@pytest.mark.parametrize(
    "value", [None, "CCO", ch.from_smiles("CCO"), C.MolFromSmarts("[#6]")]
)
def test_only_concrete_compatibility_molecules(value):
    with pytest.raises(TypeError):
        C.PatternFingerprint(value)


def test_empty_molecule_has_empty_fingerprint():
    mol = C.Mol(ch.RWMol().GetMol())
    assert C.PatternFingerprint(mol).GetOnBits() == []


def test_graph_cap_precedes_native_calculation(monkeypatch):
    mol = C.MolFromSmiles("C")
    monkeypatch.setattr(C.Mol, "GetNumAtoms", lambda self: 1025)
    with pytest.raises(ValueError, match="1024"):
        C.PatternFingerprint(mol)


def test_public_namespace_exports():
    assert {"PatternFingerprint", "PatternAtomCounts", "PatternBitMask"} <= set(
        C.__all__
    )
