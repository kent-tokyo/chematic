import gc
import itertools
import weakref

import pytest
import chematic


def test_masks_follow_order_and_cover_every_atom_in_match():
    mol = chematic.from_smiles("CC(=O)O")
    patterns = chematic.FeaturePatternSet(["[#8]", "C=O", "[#6]"])
    assert patterns.invariants(mol) == [4, 6, 3, 1]
    assert patterns.patterns == ["[#8]", "C=O", "[#6]"]


@pytest.mark.parametrize("factory", [list, tuple, iter])
def test_iterable_patterns_are_owned_and_reusable(factory):
    definitions = ["[#8]", "[#6]"]
    patterns = chematic.FeaturePatternSet(factory(definitions))
    definitions.clear()
    gc.collect()
    assert patterns.invariants(chematic.from_smiles("CO")) == [2, 1]
    assert patterns.invariants(chematic.from_smiles("OC")) == [1, 2]
    returned = patterns.patterns
    returned[0] = "N"
    assert patterns.patterns == ["[#8]", "[#6]"]


def test_generator_owner_not_retained_and_no_raw_query_lifetimes():
    class Source:
        def __iter__(self):
            yield "O"
    source = Source()
    reference = weakref.ref(source)
    patterns = chematic.FeaturePatternSet(source)
    del source
    gc.collect()
    assert reference() is None
    assert patterns.invariants(chematic.from_smiles("CO")) == [0, 1]


def test_full_width_masks_and_overflow_cap():
    patterns = chematic.FeaturePatternSet(["O"] * 32)
    assert patterns.invariants(chematic.from_smiles("CO")) == [0, 0xFFFFFFFF]
    with pytest.raises(ValueError, match="32"):
        chematic.FeaturePatternSet(["O"] * 33)
    consumed = []
    def endless():
        for index in itertools.count():
            consumed.append(index)
            yield "O"
    with pytest.raises(ValueError, match="32"):
        chematic.FeaturePatternSet(endless())
    assert len(consumed) == 33


@pytest.mark.parametrize("count", [0, 32])
def test_original_generator_error_propagates(count):
    error = RuntimeError("source failed")
    def broken():
        yield from ["O"] * count
        raise error
    with pytest.raises(RuntimeError) as caught:
        chematic.FeaturePatternSet(broken())
    assert caught.value is error


@pytest.mark.parametrize("source", [None, 42, "O", b"O", [None], [b"O"], [1]])
def test_wrong_input_types_refused(source):
    with pytest.raises(TypeError):
        chematic.FeaturePatternSet(source)


@pytest.mark.parametrize("text", ["", "[bad", "C" * 33, "C" * 4097])
def test_invalid_or_oversized_patterns_refused(text):
    with pytest.raises(ValueError):
        chematic.FeaturePatternSet([text])


def test_empty_set_means_zero_features():
    patterns = chematic.FeaturePatternSet([])
    assert patterns.invariants(chematic.from_smiles("CCO")) == [0, 0, 0]
    empty = chematic.RWMol().GetMol()
    assert patterns.invariants(empty) == []
    assert not any(patterns.fingerprint(empty))


def test_pattern_changes_and_high_bits_affect_fingerprint():
    mol = chematic.from_smiles("CCO")
    oxygen = chematic.FeaturePatternSet(["O"])
    carbon = chematic.FeaturePatternSet(["C"])
    high = chematic.FeaturePatternSet(["N"] * 31 + ["O"])
    assert high.invariants(mol) == [0, 0, 1 << 31]
    assert oxygen.fingerprint(mol) != carbon.fingerprint(mol)
    assert oxygen.fingerprint(mol) != high.fingerprint(mol)


@pytest.mark.parametrize("radius", [0, 1, 2, 3])
@pytest.mark.parametrize("nbits", [128, 256, 512, 1024, 2048])
def test_supported_fingerprint_shapes_and_stability(radius, nbits):
    patterns = chematic.FeaturePatternSet(["O", "C"])
    mol = chematic.from_smiles("CCO")
    result = patterns.fingerprint(mol, radius=radius, nbits=nbits)
    assert isinstance(result, bytes)
    assert len(result) == nbits // 8
    assert result == patterns.fingerprint(mol, radius=radius, nbits=nbits)


@pytest.mark.parametrize("kwargs", [{"radius": 4}, {"nbits": 0}, {"nbits": 4096},
                                    {"max_visits": 0}, {"max_visits": 1_000_001}])
def test_invalid_configuration_refused(kwargs):
    with pytest.raises(ValueError):
        chematic.FeaturePatternSet(["O"]).fingerprint(chematic.from_smiles("CCO"), **kwargs)


def test_budget_exhaustion_refuses_partial_masks():
    patterns = chematic.FeaturePatternSet(["CC"])
    mol = chematic.from_smiles("CCCCCC")
    with pytest.raises(ValueError, match="budget"):
        patterns.invariants(mol, max_visits=1)
    with pytest.raises(ValueError, match="budget"):
        patterns.fingerprint(mol, max_visits=1)
    assert patterns.invariants(mol) == [1] * 6


def test_target_size_guard():
    with pytest.raises(ValueError, match="256"):
        chematic.FeaturePatternSet(["C"]).invariants(chematic.from_smiles("C" * 257))


def test_aromatic_and_kekule_forms_agree():
    patterns = chematic.FeaturePatternSet(["[a]", "[#7]", "[#6]"])
    aromatic = chematic.from_smiles("n1ccccc1")
    kekule = chematic.from_smiles("N1=CC=CC=C1")
    assert patterns.invariants(aromatic) == patterns.invariants(kekule)
    assert patterns.fingerprint(aromatic) == patterns.fingerprint(kekule)
