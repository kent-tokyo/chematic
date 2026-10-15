"""Checked numeric boundaries without changing the RDKit compatibility lane."""
import math

import chematic
import pytest


@pytest.mark.parametrize("weights", [False, True])
def test_normal_matrix_matches_legacy_exactly(weights):
    mol = chematic.from_smiles("CCO")
    coords = [[0, 0, 0], [3, 0, 0], [3, 4, 0]]
    assert mol.rdkit_distance_matrix_3d_checked(coords, weights) == mol.rdkit_distance_matrix_3d(coords, weights)


@pytest.mark.parametrize("value", [1e308, 1e-160, 1e-300, math.ulp(0.0)])
def test_representable_extremes(value):
    matrix = chematic.from_smiles("CC").rdkit_distance_matrix_3d_checked([[0, 0, 0], [value, 0, 0]])
    assert matrix == [[0, value], [value, 0]]


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_single_atom_coordinates_refuse(value):
    with pytest.raises(ValueError, match="non-finite coordinate"):
        chematic.from_smiles("[He]").rdkit_distance_matrix_3d_checked([[value, 0, 0]])


def test_unrepresentable_distance_and_undefined_diagonal_refuse():
    with pytest.raises(ValueError, match="non-finite 3D distance"):
        chematic.from_smiles("CC").rdkit_distance_matrix_3d_checked([[1e308, 0, 0], [-1e308, 0, 0]])
    dummy = chematic.from_smiles("*")
    assert dummy.rdkit_distance_matrix_3d_checked([[0, 0, 0]]) == [[0]]
    with pytest.raises(ValueError, match="atomic number zero"):
        dummy.rdkit_distance_matrix_3d_checked([[0, 0, 0]], True)


def test_prepared_hydrogen_order_and_coordinate_count():
    mol = chematic.from_smiles("[H]C([H])([H])[H]")
    assert mol.rdkit_distance_matrix_3d_checked([[0, 0, 0]]) == [[0]]
    with pytest.raises(ValueError, match="coordinates for 1 atoms"):
        mol.rdkit_distance_matrix_3d_checked([[0, 0, 0]] * 5)
    empty, _coords = chematic.from_cjson('{"atoms":{"elements":{"number":[]}}}')
    assert empty.rdkit_distance_matrix_3d_checked([]) == []


@pytest.mark.parametrize("n", [1025, 65536])
def test_oversized_inputs_refuse_before_coordinate_content(n):
    mol = chematic.from_smiles(".".join(["C"] * n))
    with pytest.raises(ValueError, match="distance matrix resource limit"):
        mol.rdkit_distance_matrix_3d_checked([[math.nan, 0, 0]])


def test_legacy_extreme_results_remain_available():
    mol = chematic.from_smiles("CC")
    assert math.isinf(mol.rdkit_distance_matrix_3d([[0, 0, 0], [1e308, 0, 0]])[0][1])
    assert mol.rdkit_distance_matrix_3d([[0, 0, 0], [1e-300, 0, 0]])[0][1] == 0


def test_added_hydrogens_keep_the_added_h_prepared_order():
    mol = chematic.from_smiles("C").add_hydrogens()
    coords = [[float(i), 0, 0] for i in range(5)]
    matrix = mol.rdkit_distance_matrix_3d_checked(coords)
    assert len(matrix) == 5
    assert matrix[0] == [0, 1, 2, 3, 4]
    assert matrix == mol.rdkit_distance_matrix_3d(coords)
