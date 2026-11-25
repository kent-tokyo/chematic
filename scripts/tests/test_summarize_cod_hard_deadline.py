import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from summarize_cod_hard_deadline import recovered


def sample():
    return {
        'geometry_reference': {'bonds_supported': 1, 'bonds_total': 1,
                               'angles_supported': 1, 'angles_total': 1},
        'engines': {'chematic': {'rmsd': [0.4], 'statuses': [1],
            'stereo': [{'atoms': [], 'bonds': []}],
            'geometry': [{'finite': True, 'bond_violations': 0, 'angle_violations': 0,
                          'heavy_atom_clashes': 0}]}}
    }


def test_recovery_requires_complete_scoring_and_optional_convergence():
    row = sample()
    assert recovered(row, 'chematic')
    assert not recovered(row, 'chematic', converged=True)
    row['engines']['chematic']['statuses'][0] = 0
    assert recovered(row, 'chematic', converged=True)
    row['engines']['chematic']['stereo'][0]['atoms'] = [0]
    assert not recovered(row, 'chematic')
    row['engines']['chematic']['stereo'].clear()
    assert not recovered(row, 'chematic')
    assert not recovered(row, 'rdkit')
