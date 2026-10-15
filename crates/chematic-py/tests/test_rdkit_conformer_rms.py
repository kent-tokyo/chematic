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
