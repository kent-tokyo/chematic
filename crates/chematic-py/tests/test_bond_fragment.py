"""Native fragment graphs and all-or-error caller counter updates."""

from dataclasses import FrozenInstanceError
import json
import pytest
import chematic as ch


def test_counters_replace_and_attachment_provenance():
    mol = ch.from_smiles("CCOC")
    before = mol.smiles
    counts = [99, -5, 99, 99]
    result = ch.fragment_on_bonds_checked(mol, [2, 0], cuts_per_atom=counts)
    assert counts == [1, 1, 1, 1]
    assert result.cuts_per_atom == tuple(counts)
    assert [
        (a.source_bond, a.anchor, a.dummy, a.isotope) for a in result.attachments
    ] == [
        (2, 3, 4, 2),
        (2, 2, 5, 3),
        (0, 1, 6, 0),
        (0, 0, 7, 1),
    ]
    assert result.molecule.atom_count == 8
    assert mol.smiles == before
    with pytest.raises(FrozenInstanceError):
        result.cuts_per_atom = ()
    with pytest.raises(FrozenInstanceError):
        result.attachments[0].isotope = 99


def test_multiple_cuts_at_one_atom_and_custom_labels():
    result = ch.fragment_on_bonds_checked(
        ch.from_smiles("CC(C)O"), [0, 1, 2], dummy_labels=[(7, 9), (0, 65535), (55, 55)]
    )
    assert result.cuts_per_atom == (1, 3, 1, 1)
    assert [a.isotope for a in result.attachments] == [7, 9, 0, 65535, 55, 55]
    assert "[65535*]" in result.molecule.smiles


@pytest.mark.parametrize(
    "smiles,cut",
    [
        ("CCO", 0),
        ("CC(=O)O", 0),
        ("c1ccccc1CO", 6),
        ("[13CH3:7]CO", 0),
        ("CC[NH3+]", 0),
        ("[2H]OC", 0),
    ],
)
def test_original_atoms_and_hydrogens_survive(smiles, cut):
    mol = ch.from_smiles(smiles)
    before = json.loads(mol.to_moljson())["atoms"]
    result = ch.fragment_on_bonds_checked(mol, [cut])
    after = json.loads(result.molecule.to_moljson())["atoms"]
    assert before == after[: len(before)]


def test_properties_are_copied_independently():
    mol = ch.from_smiles("CCO")
    mol.SetProp("source", "original")
    result = ch.fragment_on_bonds_checked(mol, [0])
    assert result.molecule.GetProp("source") == "original"
    result.molecule.SetProp("source", "result")
    assert mol.GetProp("source") == "original"


@pytest.mark.parametrize("smiles", ["CCO", "N[C@H](C)O", "*"])
def test_empty_cuts_copy_unchanged(smiles):
    mol = ch.from_smiles(smiles)
    counts = [17] * mol.atom_count
    result = ch.fragment_on_bonds_checked(mol, [], cuts_per_atom=counts)
    assert counts == [0] * mol.atom_count
    assert result.molecule.smiles == mol.smiles
    assert result.attachments == ()


@pytest.mark.parametrize(
    "smiles,indices",
    [
        ("CCO", [0, 2]),
        ("CCO", [0, 0]),
        ("C=C", [0]),
        ("C#N", [0]),
        ("C1CC1", [0]),
        ("c1ccccc1", [0]),
        ("N[C@H](C)O", [0]),
        ("C/C=C/CCC", [4]),
        ("*CC", [1]),
    ],
)
def test_native_refusals_leave_counters_and_source_unchanged(smiles, indices):
    mol = ch.from_smiles(smiles)
    before = mol.smiles
    counts = [17] * mol.atom_count
    with pytest.raises(ValueError):
        ch.fragment_on_bonds_checked(mol, indices, cuts_per_atom=counts)
    assert counts == [17] * mol.atom_count
    assert mol.smiles == before


@pytest.mark.parametrize("value", [None, "CC", object()])
def test_mol_type(value):
    with pytest.raises(TypeError):
        ch.fragment_on_bonds_checked(value, [0])


@pytest.mark.parametrize("value", [None, "0", b"0", 0, iter([0])])
def test_indices_container(value):
    with pytest.raises(TypeError):
        ch.fragment_on_bonds_checked(ch.from_smiles("CC"), value)


@pytest.mark.parametrize("value", [True, False, "0", 1.0, None])
def test_index_element_type(value):
    with pytest.raises(TypeError):
        ch.fragment_on_bonds_checked(ch.from_smiles("CC"), [value])


@pytest.mark.parametrize("value", [-1, 2**32])
def test_index_range(value):
    with pytest.raises(ValueError):
        ch.fragment_on_bonds_checked(ch.from_smiles("CC"), [value])


@pytest.mark.parametrize(
    "labels,error",
    [
        ([], ValueError),
        ([(1,)], ValueError),
        ([(1, 2, 3)], ValueError),
        ([(-1, 2)], ValueError),
        ([(65536, 2)], ValueError),
        ([(True, 2)], TypeError),
        ([(1.0, 2)], TypeError),
        ([("1", 2)], TypeError),
        ([None], TypeError),
        ("ab", TypeError),
    ],
)
def test_labels_validation(labels, error):
    counts = [17, 17]
    with pytest.raises(error):
        ch.fragment_on_bonds_checked(
            ch.from_smiles("CC"), [0], dummy_labels=labels, cuts_per_atom=counts
        )
    assert counts == [17, 17]


@pytest.mark.parametrize(
    "counts,error",
    [
        ((), TypeError),
        ([], ValueError),
        ([1], ValueError),
        ([True, 1], TypeError),
        ([1.0, 1], TypeError),
        ([None, 1], TypeError),
    ],
)
def test_counter_validation(counts, error):
    with pytest.raises(error):
        ch.fragment_on_bonds_checked(ch.from_smiles("CC"), [0], cuts_per_atom=counts)


def test_counter_list_subclass_is_refused():
    class Bad(list):
        def __setitem__(self, *args):
            pytest.fail("must not mutate a user-defined list subclass")

    with pytest.raises(TypeError):
        ch.fragment_on_bonds_checked(
            ch.from_smiles("CC"), [0], cuts_per_atom=Bad([0, 0])
        )


def test_cut_and_atom_caps():
    with pytest.raises(ValueError, match="1024"):
        ch.fragment_on_bonds_checked(ch.from_smiles("CC"), [0] * 1025)
    rw = ch.RWMol()
    for _ in range(4097):
        rw.AddAtom(6)
    with pytest.raises(ValueError, match="4096"):
        ch.fragment_on_bonds_checked(rw.GetMol(), [])


def test_index_and_label_sequences():
    result = ch.fragment_on_bonds_checked(
        ch.from_smiles("CCC"), range(2), dummy_labels=((0, 1), (2, 3))
    )
    assert result.cuts_per_atom == (1, 2, 1)
