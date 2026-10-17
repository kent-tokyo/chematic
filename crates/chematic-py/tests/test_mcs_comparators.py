"""Typed native comparator choices and cross-element pruning regressions."""
import pytest
import chematic as ch


@pytest.mark.parametrize("function", [ch.find_mcs, ch.find_mcs_checked])
@pytest.mark.parametrize("mode", list(ch.McsAtomCompare))
def test_enum_and_legacy_string_agree(function, mode):
    mols = [ch.from_smiles("CCCO"), ch.from_smiles("NCCN")]
    enum_result = function(mols, atom_compare=mode)
    string_result = function(mols, atom_compare=mode.value)
    if function is ch.find_mcs_checked:
        enum_result, enum_timeout = enum_result
        string_result, string_timeout = string_result
        assert enum_timeout is string_timeout is False
    assert enum_result.atom_count == string_result.atom_count
    assert enum_result.smiles == string_result.smiles


@pytest.mark.parametrize("inputs,size", [
    (["CCC", "NNN"], 3),
    (["CCCO", "NCCN"], 4),
    (["COC", "NNN"], 3),
    (["CCCC", "NNN", "OOO"], 3),
    (["c1ccccc1", "n1cncnc1"], 6),
])
@pytest.mark.parametrize("mode", [ch.McsAtomCompare.ANY, ch.McsAtomCompare.ANY_HEAVY_ATOM])
def test_cross_element_topology(inputs, size, mode):
    for strings in [inputs, inputs[::-1]]:
        result, timeout = ch.find_mcs_checked(
            [ch.from_smiles(s) for s in strings], atom_compare=mode
        )
        assert timeout is False
        assert result.atom_count == size


@pytest.mark.parametrize("function", [ch.find_mcs, ch.find_mcs_checked])
@pytest.mark.parametrize("mode", list(ch.McsBondCompare))
def test_bond_enum_and_legacy_string_agree(function, mode):
    mols = [ch.from_smiles("CCC"), ch.from_smiles("CC=C")]
    a = function(mols, bond_compare=mode)
    b = function(mols, bond_compare=mode.value)
    if function is ch.find_mcs_checked:
        a, at = a
        b, bt = b
        assert at is bt is False
    assert a.atom_count == b.atom_count == (3 if mode is ch.McsBondCompare.ANY else 2)


def test_heavy_comparator_excludes_explicit_hydrogen():
    mols = [ch.from_smiles("[H][H]"), ch.from_smiles("NN")]
    assert ch.find_mcs(mols, atom_compare=ch.McsAtomCompare.ANY_HEAVY_ATOM) is None
    assert ch.find_mcs([mols[0], mols[0]], atom_compare=ch.McsAtomCompare.ANY_HEAVY_ATOM) is None
    result = ch.find_mcs(mols, atom_compare=ch.McsAtomCompare.ANY)
    assert result.atom_count == 2


@pytest.mark.parametrize("argument", ["atom_compare", "bond_compare"])
@pytest.mark.parametrize("value", [None, True, 0, object(), lambda: True])
def test_non_string_non_enum_comparator_is_a_type_error(argument, value):
    mols = [ch.from_smiles("CC"), ch.from_smiles("CC")]
    with pytest.raises(TypeError):
        ch.find_mcs_checked(mols, **{argument: value})


@pytest.mark.parametrize("argument", ["atom_compare", "bond_compare"])
def test_invalid_legacy_string_is_a_value_error(argument):
    with pytest.raises(ValueError):
        ch.find_mcs_checked([ch.from_smiles("CC")], **{argument: "wrong"})


def test_default_comparison_is_still_elements():
    mols = [ch.from_smiles("CCC"), ch.from_smiles("NNN")]
    assert ch.find_mcs(mols) is None
    assert ch.find_mcs(mols, atom_compare=ch.McsAtomCompare.ELEMENTS) is None
