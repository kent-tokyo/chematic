"""Explicit molecular bond geometry; no periodic bond perception."""
import math

import pytest
import chematic


def rows(text):
    return [line.split() for line in text.split("_ccdc_geom_bond_type\n")[1].splitlines()]


def test_explicit_water_graph_labels_and_existing_reader_contract():
    mol = chematic.from_smiles("[H]O[H]")
    coords = [[0.9, 0, 0], [0, 0, 0], [0, 1.1, 0]]
    text = chematic.write_cif_with_bonds_checked(mol, coords)
    assert rows(text) == [["H1", "O2", "0.9", ".", ".", "S"], ["O2", "H3", "1.1", ".", ".", "S"]]
    assert text.startswith(chematic.write_cif_checked(mol, coords))
    assert chematic.parse_cif(text)["coords"] == coords


@pytest.mark.parametrize("smiles,kind", [("CC", "S"), ("C=C", "D"), ("C#C", "T"), ("c1ccccc1", "A")])
def test_bond_types(smiles, kind):
    mol = chematic.from_smiles(smiles)
    coords = [[float(i), 0, 0] for i in range(6 if kind == "A" else 2)]
    assert {r[5] for r in rows(chematic.write_cif_with_bonds_checked(mol, coords))} == {kind}


def test_cell_does_not_wrap_bond_distance():
    mol = chematic.from_smiles("CC")
    coords = [[0.1, 0, 0], [9.9, 0, 0]]
    cell = dict(a=10., b=10., c=10., alpha=90., beta=90., gamma=90.)
    text = chematic.write_cif_with_bonds_checked(mol, coords, cell)
    assert "_atom_site_fract_x" in text
    assert float(rows(text)[0][2]) == 9.8
    assert rows(text) == rows(chematic.write_cif_with_bonds_checked(mol, coords))


@pytest.mark.parametrize("value", [1e308, 1e-300, math.ulp(0.)])
def test_extreme_representable_distances(value):
    text = chematic.write_cif_with_bonds_checked(chematic.from_smiles("CC"), [[0, 0, 0], [value, 0, 0]])
    assert float(rows(text)[0][2]) == value


@pytest.mark.parametrize("coords", [[], [[0,0,0]], [[0,0,0],[math.nan,0,0]], [[0,0,0],[math.inf,0,0]], [[0,0,0],[0,0,0]], [[1e308,0,0],[-1e308,0,0]]])
def test_invalid_coordinates_and_distances_refuse(coords):
    with pytest.raises(ValueError):
        chematic.write_cif_with_bonds_checked(chematic.from_smiles("CC"), coords)


def test_cell_fields_are_checked():
    mol = chematic.from_smiles("CC")
    with pytest.raises(ValueError, match="missing required key"):
        chematic.write_cif_with_bonds_checked(mol, [[0,0,0],[1,0,0]], {})
    with pytest.raises(ValueError, match="length b"):
        chematic.write_cif_with_bonds_checked(mol, [[0,0,0],[1,0,0]], dict(a=10., b=0., c=10., alpha=90., beta=90., gamma=90.))


def test_dative_bonds_are_not_downgraded():
    with pytest.raises(ValueError, match="unsupported order Dative"):
        chematic.write_cif_with_bonds_checked(chematic.from_smiles("N->[Cu]"), [[0,0,0],[1,0,0]])


def test_disconnected_ions_do_not_gain_bonds():
    mol = chematic.from_smiles("[Na+].[Cl-]")
    coords = [[0,0,0],[0,0,0]]
    assert chematic.write_cif_with_bonds_checked(mol, coords) == chematic.write_cif_checked(mol, coords)
