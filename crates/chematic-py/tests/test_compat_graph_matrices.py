"""Bounded RDKit-style NumPy entry points over existing native profiles."""

import numpy as np
import pytest

import chematic
from chematic import rdkit_compat as C


@pytest.mark.parametrize(
    "smiles,weighted",
    [
        ("C=C", [[0, 2], [2, 0]]),
        ("C#N", [[0, 3], [3, 0]]),
        ("N->[Cu]", [[0, 0], [1, 0]]),
        ("[Cu]<-N", [[0, 1], [0, 0]]),
        ("[Mo]$[Mo]", [[0, 4], [4, 0]]),
    ],
)
def test_adjacency_weights_and_binary_dtype(smiles, weighted):
    mol = C.MolFromSmiles(smiles)
    binary = C.GetAdjacencyMatrix(mol)
    actual = C.GetAdjacencyMatrix(mol, useBO=True)
    np.testing.assert_array_equal(actual, weighted)
    np.testing.assert_array_equal(binary, [[0, 1], [1, 0]])
    assert binary.dtype == np.int32
    assert actual.dtype == np.float64


def test_aromatic_edges_and_disconnected_sentinels():
    mol = C.MolFromSmiles("c1ccccc1.O")
    adjacency = C.GetAdjacencyMatrix(mol, useBO=True)
    assert adjacency.shape == (7, 7)
    assert adjacency[0, 1] == 1.5
    assert (adjacency[6] == 0).all()
    distance = C.GetDistanceMatrix(mol, useBO=True)
    assert distance[0, 3] == pytest.approx(2.0)
    assert distance[0, 6] == distance[6, 0] == 1e8


def test_distance_edge_lengths_and_atomic_weight_diagonal():
    mol = C.MolFromSmiles("C=CO")
    np.testing.assert_array_equal(
        C.GetDistanceMatrix(mol), [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
    )
    np.testing.assert_array_equal(
        C.GetDistanceMatrix(mol, True, True),
        [[1, 0.5, 1.5], [0.5, 1, 1], [1.5, 1, 0.75]],
    )
    assert C.GetDistanceMatrix(mol).dtype == np.float64


@pytest.mark.parametrize("weighted", [False, True])
def test_empty_arrays_own_zero_by_zero_storage(weighted):
    mol = C.RWMol().GetMol()
    for matrix in [
        C.GetAdjacencyMatrix(mol, weighted),
        C.GetDistanceMatrix(mol, weighted),
    ]:
        assert matrix.shape == (0, 0)
        assert matrix.flags.owndata
        assert matrix.flags.c_contiguous and matrix.flags.writeable


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
def test_repeated_force_calls_return_independent_mutable_results(function):
    mol = C.MolFromSmiles("CCO")
    mol.SetIntProp("sentinel", 9647)
    expected = function(mol)
    result = function(mol, force=True)
    assert result.flags.owndata and result.flags.c_contiguous and result.flags.writeable
    assert not np.shares_memory(result, expected)
    result[:] = 999
    np.testing.assert_array_equal(function(mol), expected)
    np.testing.assert_array_equal(function(mol, force=True), expected)
    assert mol.GetIntProp("sentinel") == 9647
    assert mol.GetPropNames() == ["sentinel"]


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
@pytest.mark.parametrize("count", [1025, 46341, 65535, 65536, 65537])
def test_issue_9647_boundaries_reject_before_dense_conversion(
    function, count, monkeypatch
):
    mol = C.Mol(chematic.from_smiles(".".join(["C"] * count)))

    # Prove no NumPy dense output allocation happens on these actual graphs.
    def forbidden(*args, **kwargs):
        raise AssertionError("dense array must not be allocated")

    monkeypatch.setattr(np, "array", forbidden)
    monkeypatch.setattr(np, "empty", forbidden)
    with pytest.raises(ValueError, match=f"{count} input atoms exceeds 1024"):
        function(mol)


def test_exact_adjacency_cap_and_nonshared_profiles():
    mol = C.Mol(chematic.from_smiles("C" * 1024))
    matrix = C.GetAdjacencyMatrix(mol)
    assert matrix.shape == (1024, 1024)
    assert matrix.sum() == 2046
    assert matrix[0, 1] == 1 and matrix[0, -1] == 0
    assert matrix[-1, -2] == 1


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
def test_raw_hydrogen_reindexing_is_refused(function):
    mol = C.Mol(chematic.from_smiles("[H]C([H])([H])[H]"))
    assert mol.GetNumAtoms() == 5
    with pytest.raises(ValueError, match="changed atom indices"):
        function(mol)
    # Native RDKit H addition explicitly retains the generated graph atoms.
    retained = C.Mol(chematic.from_smiles("C").rdkit_add_hydrogens())
    assert function(retained).shape == (5, 5)


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
def test_guard_applies_before_hydrogen_removal(function):
    mol = C.Mol(chematic.from_smiles(".".join(["[H]C"] * 513)))
    with pytest.raises(ValueError, match="1026 input atoms"):
        function(mol)


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
@pytest.mark.parametrize("bad", [None, "CC", False, chematic.from_smiles("CC")])
def test_molecule_type_boundary(function, bad):
    with pytest.raises(TypeError, match="compat.Mol"):
        function(bad)


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
@pytest.mark.parametrize(
    "option,value",
    [
        ("useBO", 0),
        ("useBO", np.bool_(True)),
        ("force", 1),
        ("prefix", b""),
        ("prefix", None),
    ],
)
def test_option_type_boundaries(function, option, value):
    with pytest.raises(TypeError):
        function(C.MolFromSmiles("CC"), **{option: value})


@pytest.mark.parametrize("value", [0, 1, np.bool_(False), None, "yes"])
def test_atom_weight_flag_requires_python_bool(value):
    with pytest.raises(TypeError):
        C.GetDistanceMatrix(C.MolFromSmiles("CC"), useAtomWts=value)


@pytest.mark.parametrize(
    "value,exception",
    [
        (1, ValueError),
        (-1, ValueError),
        (2**100, ValueError),
        (True, TypeError),
        (np.bool_(False), TypeError),
        (0.0, TypeError),
        ("0", TypeError),
    ],
)
def test_fill_value_boundary(value, exception):
    with pytest.raises(exception):
        C.GetAdjacencyMatrix(C.RWMol().GetMol(), emptyVal=value)
    assert C.GetAdjacencyMatrix(C.RWMol().GetMol(), emptyVal=np.int64(0)).shape == (
        0,
        0,
    )


@pytest.mark.parametrize("function", [C.GetAdjacencyMatrix, C.GetDistanceMatrix])
def test_nonempty_prefix_is_explicitly_unsupported(function):
    with pytest.raises(NotImplementedError, match="prefix"):
        function(C.MolFromSmiles("CC"), prefix="alternate")


def test_public_star_exports_and_native_api_are_both_available():
    namespace = {}
    exec("from chematic.rdkit_compat import *", namespace)
    assert namespace["GetAdjacencyMatrix"] is C.GetAdjacencyMatrix
    assert namespace["GetDistanceMatrix"] is C.GetDistanceMatrix
    native = chematic.from_smiles("CC")
    assert isinstance(native.rdkit_adjacency_matrix(), list)
    assert isinstance(native.rdkit_distance_matrix(), list)


def test_editor_snapshots_and_permutations_remain_independent():
    source = C.MolFromSmiles("CCO")
    matrix = C.GetAdjacencyMatrix(source)
    editor = C.RWMol(source)
    editor.RemoveAtom(2)
    changed = editor.GetMol()
    assert C.GetAdjacencyMatrix(changed).shape == (2, 2)
    assert C.GetAdjacencyMatrix(source).shape == (3, 3)
    np.testing.assert_array_equal(C.GetAdjacencyMatrix(source), matrix)
    reordered = C.RenumberAtoms(source, [2, 0, 1])
    order = [2, 0, 1]
    np.testing.assert_array_equal(
        C.GetDistanceMatrix(reordered),
        C.GetDistanceMatrix(source)[np.ix_(order, order)],
    )
