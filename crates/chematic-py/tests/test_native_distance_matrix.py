"""Native heavy-atom matrix boundaries; independent of the RDKit-prepared API."""
import chematic
import pytest


@pytest.mark.parametrize("n", [1025, 46341, 65535, 65536, 65537])
def test_issue_boundaries_reject_before_allocating(n):
    mol = chematic.from_smiles(".".join(["C"] * n))
    for kwargs in [{}, {"max_atoms": 2**63, "max_entries": 2**63}]:
        with pytest.raises(ValueError, match=f"{n} atoms exceeds 1024"):
            mol.topological_distance_matrix(**kwargs)


def test_exact_atom_cap_and_longest_path():
    matrix = chematic.from_smiles("C" * 1024).topological_distance_matrix()
    assert len(matrix) == 1024
    assert all(len(row) == 1024 for row in matrix)
    assert matrix[0][-1] == matrix[-1][0] == 1023
    assert matrix[512][512] == 0


@pytest.mark.parametrize("kwargs, message", [
    ({"max_atoms": 2}, "3 atoms exceeds 2"),
    ({"max_bonds": 1}, "2 bonds exceeds 1"),
    ({"max_entries": 8}, "9 entries exceeds 8"),
])
def test_tighter_budgets(kwargs, message):
    with pytest.raises(ValueError, match=message):
        chematic.from_smiles("CCC").topological_distance_matrix(**kwargs)


def test_exact_budgets_and_native_disconnected_sentinel():
    assert chematic.from_smiles("CCC").topological_distance_matrix(
        max_atoms=3, max_bonds=2, max_entries=9) == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
    assert chematic.from_smiles("C.O").topological_distance_matrix() == [[0, 4294967295], [4294967295, 0]]


def test_explicit_h_input_counts_before_filtering():
    mol = chematic.from_smiles(".".join(["[H]C"] * 513))
    with pytest.raises(ValueError, match="1026 atoms exceeds 1024"):
        mol.topological_distance_matrix()


def test_hydrogen_paths_and_wildcards_keep_native_semantics():
    assert chematic.from_smiles("C[H]O.*").topological_distance_matrix() == [
        [0, 4294967295, 4294967295], [4294967295, 0, 4294967295], [4294967295, 4294967295, 0]]
    assert chematic.from_smiles("[2H].[H][H]").topological_distance_matrix() == []


@pytest.mark.parametrize("argument", ["max_atoms", "max_bonds", "max_entries"])
def test_negative_limits_rejected(argument):
    with pytest.raises(OverflowError):
        chematic.from_smiles("CC").topological_distance_matrix(**{argument: -1})


def test_no_positional_budget_arguments():
    with pytest.raises(TypeError):
        chematic.from_smiles("CC").topological_distance_matrix(2)
