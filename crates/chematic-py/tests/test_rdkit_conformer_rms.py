import ast
from pathlib import Path

import pytest

import chematic


PROBE_CONFORMERS = [
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
    [[2.0, -1.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 2.0, 0.0]],
]


def test_all_conformer_best_rms_matches_repeated_calls_in_rdkit_order():
    mol = chematic.from_smiles("CCO")
    actual = mol.rdkit_all_conformer_best_rms(PROBE_CONFORMERS)
    expected = [
        mol.rdkit_best_rms(PROBE_CONFORMERS[1], PROBE_CONFORMERS[0]),
        mol.rdkit_best_rms(PROBE_CONFORMERS[2], PROBE_CONFORMERS[0]),
        mol.rdkit_best_rms(PROBE_CONFORMERS[2], PROBE_CONFORMERS[1]),
    ]
    assert actual == pytest.approx(expected, abs=1e-12)
    assert actual[0] < 1e-7


def test_all_conformer_best_rms_to_reference_is_reference_major():
    probe = chematic.from_smiles("CCO")
    reference = chematic.from_smiles("OCC")
    probe_conformers = [PROBE_CONFORMERS[0], PROBE_CONFORMERS[2]]
    reference_conformers = [
        [[1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[1.0, 3.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    ]
    actual = probe.rdkit_all_conformer_best_rms_to_reference(
        reference, probe_conformers, reference_conformers
    )
    assert len(actual) == 4
    assert actual[0] < 1e-7
    assert actual[2] > 0.1


def test_all_conformer_best_rms_rejects_bad_coordinate_count():
    mol = chematic.from_smiles("CCO")
    with pytest.raises(ValueError, match="1 positions for 3 atoms"):
        mol.rdkit_all_conformer_best_rms([[[0.0, 0.0, 0.0]]])


def test_symmetric_clustering_collapses_atom_permutations():
    mol = chematic.from_smiles("FC(F)(F)C")
    base = [
        [1.1, 0.1, -0.3],
        [0.0, 0.0, 0.0],
        [-0.5, 1.0, -0.3],
        [-0.5, -0.6, -0.9],
        [0.0, 0.0, 1.5],
    ]
    swapped = [point.copy() for point in base]
    swapped[0], swapped[2] = swapped[2], swapped[0]
    distinct = [point.copy() for point in base]
    distinct[0][0] += 4.0

    assert mol.rdkit_best_rms(base, swapped) < 1e-6
    assert mol.cluster_conformers_by_rms_symmetric(
        [base, swapped, distinct], 0.05
    ) == [0, 2]
    assert mol.cluster_conformers_by_rms_symmetric([base, swapped], 0.0) == [0, 1]


def test_symmetric_clustering_rejects_bad_coordinates_and_threshold():
    mol = chematic.from_smiles("CCO")
    with pytest.raises(ValueError, match="conformer 0 has 1 positions for 3 atoms"):
        mol.cluster_conformers_by_rms_symmetric([[[0.0, 0.0, 0.0]]])
    with pytest.raises(ValueError, match="non-finite coordinates"):
        mol.cluster_conformers_by_rms_symmetric(
            [[[0.0, 0.0, 0.0], [float("nan"), 0.0, 0.0], [1.0, 0.0, 0.0]]]
        )
    with pytest.raises(ValueError, match="rmsd_threshold must be finite"):
        mol.cluster_conformers_by_rms_symmetric([], float("nan"))


def test_symmetric_clustering_is_declared_on_mol_stub():
    stub = Path(__file__).resolve().parents[1] / "python/chematic/__init__.pyi"
    module = ast.parse(stub.read_text(encoding="utf-8"))
    mol = next(
        node for node in module.body if isinstance(node, ast.ClassDef) and node.name == "Mol"
    )
    assert any(
        isinstance(node, ast.FunctionDef) and node.name == "cluster_conformers_by_rms_symmetric"
        for node in mol.body
    )
