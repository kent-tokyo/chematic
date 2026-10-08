"""Parsers raise ChematicInputError with a category, code and format."""

import pytest

import chematic


@pytest.mark.parametrize(
    "call, fmt",
    [
        (lambda: chematic.from_smiles("C1CC("), "smiles"),
        (lambda: chematic.from_smiles("CC").has_substructure("[C;"), "smarts"),
        (lambda: chematic.from_mol_block("not a mol block"), "mol_block"),
        (lambda: chematic.from_mol_block_with_coords("not a mol block"), "mol_block"),
        (lambda: chematic.from_inchi("InChI=1S/garbage"), "inchi"),
        (lambda: chematic.parse_mmcif("data_x\nloop_\n_cell.length_a\n1\n"), "mmcif"),
    ],
)
def test_malformed_input(call, fmt):
    with pytest.raises(chematic.ChematicInputError) as info:
        call()
    error = info.value
    assert isinstance(error, ValueError)
    assert error.category == "malformed"
    assert error.code == f"{fmt}_parse"
    assert error.format == fmt


@pytest.mark.parametrize("smiles", ["O=O1C=CC=C1", "CO(C)C", "C=O=C", "[OH3]", "F(C)C"])
def test_over_valent_neutral_oxygen_or_fluorine_is_rejected(smiles):
    # #769: a neutral oxygen (or fluorine) with more bonds than it permits is
    # not a molecule; it used to be accepted and then aromatized.
    with pytest.raises(chematic.ChematicInputError) as info:
        chematic.from_smiles(smiles)
    assert info.value.category == "malformed"
    assert info.value.code == "smiles_valence"
    assert info.value.format == "smiles"
    assert "greater than permitted" in str(info.value)


@pytest.mark.parametrize(
    "smiles", ["C[O+](C)C", "c1cc[o+]cc1", "O=O", "c1ccoc1", "CS(=O)(=O)C", "CN(=O)=O", "[O]", "FC(F)(F)F"]
)
def test_charged_and_hypervalent_controls_still_parse(smiles):
    assert chematic.from_smiles(smiles).heavy_atoms > 0


def test_resource_limit():
    with pytest.raises(chematic.ChematicInputError) as info:
        chematic.parse_mmcif("data_x\n" * 100, max_input_bytes=10)
    assert info.value.category == "resource_limit"
    assert info.value.code == "input_too_large"
    assert info.value.format == "mmcif"


def test_still_a_value_error():
    with pytest.raises(ValueError, match="parenthesis"):
        chematic.from_smiles("C1CC(")
