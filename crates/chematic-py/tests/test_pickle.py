"""Versioned pickle round trips for Python Mol objects."""

import pickle

import pytest

import chematic


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)Oc1ccccc1C(=O)O",
        "N[C@@H](C)C(=O)O",
        "F/C=C\\Cl",
        "[13CH3:7][NH2+]C(=O)[O-]",
        "[Pt@SP2](Cl)(Cl)(N)N",
        "CC=[C@AL1]=CCC",
        "C1CCCCC1",
    ],
)
def test_mol_pickle_round_trip_preserves_chemistry(smiles):
    original = chematic.from_smiles(smiles)
    original.SetProp("name", "β-lactam")
    original.SetIntProp("row", 42)

    restored = pickle.loads(pickle.dumps(original, protocol=pickle.HIGHEST_PROTOCOL))

    assert restored.rdkit_smiles == original.rdkit_smiles
    assert restored.smiles == original.smiles
    assert restored.formula == original.formula
    assert restored.rdkit_ecfp4() == original.rdkit_ecfp4()
    assert restored.GetPropsAsDict() == original.GetPropsAsDict()


def test_mol_pickle_is_deterministic_for_property_insertion_order():
    left = chematic.from_smiles("CCO")
    left.SetProp("z", "last")
    left.SetProp("a", "first")
    right = chematic.from_smiles("CCO")
    right.SetProp("a", "first")
    right.SetProp("z", "last")

    assert pickle.dumps(left, protocol=5) == pickle.dumps(right, protocol=5)


def test_private_pickle_restore_rejects_malformed_state():
    with pytest.raises(ValueError, match="header|truncated"):
        chematic._restore_mol_pickle(b"not a molecule")

    valid = pickle.dumps(chematic.from_smiles("CCO"), protocol=5)
    assert pickle.loads(valid).smiles == chematic.from_smiles("CCO").smiles
