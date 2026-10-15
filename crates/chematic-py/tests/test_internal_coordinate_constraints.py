import math

import chematic
import pytest


def test_bounds_constraints_are_applied_before_smoothing():
    mol = chematic.from_smiles("CCCC")
    matrix = mol.rdkit_bounds_matrix_with_constraints(
        bond_lengths=[(0, 1, 1.8, 0.01)],
        bond_angles=[(0, 1, 2, 100.0, 1.0)],
        torsions=[(0, 1, 2, 3, 60.0, 2.0)],
    )
    assert 1.79 <= matrix[1][0] <= 1.81
    assert matrix[2][0] <= matrix[0][2]
    assert matrix[3][0] <= matrix[0][3]


def test_constrained_embedding_is_finite_and_honours_bond_length():
    mol = chematic.from_smiles("CCCC")
    coordinates = mol.rdkit_embed_with_constraints(
        bond_lengths=[(0, 1, 1.8, 0.01)],
        bond_angles=[(0, 1, 2, 100.0, 1.0)],
        torsions=[(0, 1, 2, 3, 60.0, 2.0)],
        random_seed=42,
    )
    assert len(coordinates) == mol.heavy_atoms
    assert all(math.isfinite(value) for row in coordinates for value in row)
    distance = math.dist(coordinates[0], coordinates[1])
    assert 1.75 <= distance <= 1.85


def test_invalid_or_inconsistent_constraints_fail_closed():
    mol = chematic.from_smiles("CCC")
    with pytest.raises(ValueError, match="outside molecule"):
        mol.rdkit_bounds_matrix_with_constraints(
            bond_lengths=[(0, 8, 1.5, 0.1)]
        )

    ring = chematic.from_smiles("C1CC1")
    with pytest.raises(ValueError, match="bounds inconsistent"):
        ring.rdkit_bounds_matrix_with_constraints(
            bond_lengths=[
                (0, 1, 1.0, 0.001),
                (1, 2, 1.0, 0.001),
                (2, 0, 3.0, 0.001),
            ]
        )
