"""RDKit #9647: Python receives errors instead of giant dense allocations."""
import chematic
import pytest


@pytest.mark.parametrize("n", [1025, 65536])
def test_oversized_graph_is_rejected_in_both_matrix_apis(n):
    mol = chematic.from_smiles(".".join(["C"] * n))
    for kwargs in [{}, {"use_bo": True}, {"use_atom_wts": True}]:
        with pytest.raises(ValueError, match=f"distance matrix resource limit: {n} atoms"):
            mol.rdkit_distance_matrix(**kwargs)
    with pytest.raises(ValueError, match="distance matrix resource limit"):
        mol.rdkit_distance_matrix_3d([])


def test_explicit_h_atoms_count_towards_the_preparation_limit():
    mol = chematic.from_smiles(".".join(["[H]C"] * 513))
    with pytest.raises(ValueError, match="1026 atoms"):
        mol.rdkit_distance_matrix()


def test_small_matrix_values_and_coordinate_errors_survive():
    mol = chematic.from_smiles("C=CO")
    assert mol.rdkit_distance_matrix() == [[0., 1., 2.], [1., 0., 1.], [2., 1., 0.]]
    assert mol.rdkit_distance_matrix(use_bo=True)[0] == [0., 0.5, 1.5]
    assert mol.rdkit_distance_matrix(use_atom_wts=True)[2][2] == 0.75
    assert mol.rdkit_distance_matrix_3d([[0, 0, 0], [3, 0, 0], [3, 4, 0]]) == [
        [0., 3., 5.], [3., 0., 4.], [5., 4., 0.]]
    with pytest.raises(ValueError, match="coordinates"):
        mol.rdkit_distance_matrix_3d([[0, 0, 0]])
