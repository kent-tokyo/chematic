"""Explicit core maps, attachment labels and checked refusal boundaries."""
from dataclasses import FrozenInstanceError
import json

import pytest
import chematic as ch


def test_explicit_map_attachment_provenance_and_frozen_results():
    source, core = ch.from_smiles("CCOC"), ch.from_smiles("CO")
    before = source.smiles, core.smiles
    matching = [1, 2]
    result = ch.replace_core_checked(source, core, matching)
    assert result.original_to_output == (0, None, None, 1)
    assert [(a.core_atom, a.source_core_atom, a.source_anchor, a.source_bond, a.anchor, a.dummy)
            for a in result.attachments] == [(0, 1, 0, 0, 0, 2), (1, 2, 3, 2, 1, 3)]
    assert result.molecule.atom_count == 4
    assert "[1*]" in result.molecule.smiles
    assert matching == [1, 2]
    assert (source.smiles, core.smiles) == before
    with pytest.raises(FrozenInstanceError):
        result.original_to_output = ()
    with pytest.raises(FrozenInstanceError):
        result.attachments[0].dummy = 1


def test_reversed_core_map_changes_query_index_label():
    source, core = ch.from_smiles("CCOC"), ch.from_smiles("CC")
    a = ch.replace_core_checked(source, core, [0, 1])
    b = ch.replace_core_checked(source, core, (1, 0))
    assert a.attachments[0].core_atom == 1
    assert b.attachments[0].core_atom == 0
    assert a.molecule.smiles != b.molecule.smiles


@pytest.mark.parametrize("source,core,matching", [
    ("[13CH3:7]CO", "O", [2]),
    ("CC[NH3+]", "C", [0]),
    ("[2H]OC", "C", [2]),
    ("Cc1ccccc1", "c1ccccc1", list(range(1, 7))),
    ("CCO.CN", "O", [2]),
])
def test_retained_original_atom_attributes_and_h_counts(source, core, matching):
    mol = ch.from_smiles(source)
    before = json.loads(mol.to_moljson())["atoms"]
    result = ch.replace_core_checked(mol, ch.from_smiles(core), matching)
    after = json.loads(result.molecule.to_moljson())["atoms"]
    for old, new in enumerate(result.original_to_output):
        if new is not None:
            a, b = dict(before[old]), dict(after[new])
            a.pop("id"); b.pop("id")
            assert a == b


def test_whole_source_removal_is_empty_and_preserves_properties_independently():
    mol = ch.from_smiles("CCO")
    mol.SetProp("origin", "source")
    result = ch.replace_core_checked(mol, mol, range(3))
    assert result.molecule.atom_count == 0
    assert result.original_to_output == (None, None, None)
    assert result.attachments == ()
    assert result.molecule.GetProp("origin") == "source"
    result.molecule.SetProp("origin", "output")
    assert mol.GetProp("origin") == "source"


def test_whole_component_and_disconnected_core_are_supported():
    result = ch.replace_core_checked(ch.from_smiles("CCO.C"), ch.from_smiles("C"), [3])
    assert result.original_to_output == (0, 1, 2, None)
    assert result.attachments == ()
    assert result.molecule.smiles == ch.from_smiles("CCO").smiles
    result = ch.replace_core_checked(ch.from_smiles("CCOC"), ch.from_smiles("C.C"), [0, 3])
    assert result.original_to_output == (None, 0, 1, None)
    assert [a.core_atom for a in result.attachments] == [0, 1]


@pytest.mark.parametrize("source,core,matching", [
    ("CC", "C", []), ("CC", "C", [2]), ("CC", "CC", [0, 0]),
    ("CC", "O", [0]), ("C=C", "C", [0]), ("C1CC1", "C", [0]),
    ("*CC", "C", [1]), ("CC", "*", [0]),
    ("N[C@H](C)O", "N", [0]), ("C/C=C/CCC", "C", [5]),
    ("CC", "C=C", [0, 1]), ("COC", "CC", [0, 2]),
    ("[13CH3]C", "C", [0]), ("C[NH3+]", "N", [1]),
    ("c1ccccc1C", "C", [0]),
])
def test_invalid_graph_or_map_refusals_preserve_inputs(source, core, matching):
    mol, query = ch.from_smiles(source), ch.from_smiles(core)
    before = mol.smiles, query.smiles, list(matching)
    with pytest.raises(ValueError):
        ch.replace_core_checked(mol, query, matching)
    assert (mol.smiles, query.smiles, matching) == before


@pytest.mark.parametrize("value", [None, "0", b"0", bytearray(b"0"), 0, iter([0])])
def test_bad_match_container(value):
    with pytest.raises(TypeError):
        ch.replace_core_checked(ch.from_smiles("CC"), ch.from_smiles("C"), value)


@pytest.mark.parametrize("value", [True, False, 0.0, "0", None])
def test_match_index_type(value):
    with pytest.raises(TypeError):
        ch.replace_core_checked(ch.from_smiles("CC"), ch.from_smiles("C"), [value])


@pytest.mark.parametrize("value", [-1, 2**32])
def test_match_index_range(value):
    with pytest.raises(ValueError):
        ch.replace_core_checked(ch.from_smiles("CC"), ch.from_smiles("C"), [value])


@pytest.mark.parametrize("source,core", [(None, None), ("CC", None), (None, "C")])
def test_wrong_molecule_type(source, core):
    with pytest.raises(TypeError):
        ch.replace_core_checked(source, core, [0])


def test_empty_core_and_size_limits():
    with pytest.raises(ValueError):
        ch.replace_core_checked(ch.from_smiles("CC"), ch.RWMol().GetMol(), [])
    rw = ch.RWMol()
    for _ in range(1025):
        rw.AddAtom(6)
    with pytest.raises(ValueError, match="1024"):
        ch.replace_core_checked(rw.GetMol(), rw.GetMol(), range(1025))
    for _ in range(4097 - 1025):
        rw.AddAtom(6)
    with pytest.raises(ValueError, match="4096"):
        ch.replace_core_checked(rw.GetMol(), ch.from_smiles("C"), [0])


def test_native_aromatic_representation_must_agree():
    core = ch.from_smiles("c1ccccc1")
    with pytest.raises(ValueError, match="aromaticity"):
        ch.replace_core_checked(ch.from_smiles("CC1=CC=CC=C1"), core, range(1, 7))
    result = ch.replace_core_checked(ch.from_smiles("Cc1ccccc1"), core, range(1, 7))
    assert result.molecule.atom_count == 2
