"""Endpoint-aware 1-4 stereo inspection; never a whole-bond CIP label."""
import json
from pathlib import Path

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
CASES = json.loads((ROOT / 'validation/rdkit_issues/stereo-bounds/cases.json').read_text())['cases']

@pytest.mark.parametrize('case', CASES, ids=lambda c: c['smiles'])
def test_every_outer_pair_and_reverse_path_match_hand_derived_relations(case):
    mol = chematic.from_smiles(case['smiles'])
    if case['add_hs']:
        mol = mol.add_hydrogens()
    before = mol.smiles
    bounds = mol.rdkit_bounds_matrix(do_triangle_smoothing=False)
    for row in case['paths']:
        assert mol.rdkit_alkene_path_stereo(row['path']) == row['relation']
        assert mol.rdkit_alkene_path_stereo(tuple(reversed(row['path']))) == row['relation']
    assert mol.smiles == before
    assert mol.rdkit_bounds_matrix(do_triangle_smoothing=False) == bounds

@pytest.mark.parametrize('path', [[-1,1,3,5],[0,1,3,99],[0,1,3,0],[0,2,3,5],[0,1,2,3]])
def test_invalid_integer_paths_raise_value_error(path):
    with pytest.raises(ValueError, match='alkene path'):
        chematic.from_smiles('F/C(Cl)=C(Br)/I').rdkit_alkene_path_stereo(path)

@pytest.mark.parametrize('path', [[],[0,1,3],[0,1,3,5,2]])
def test_wrong_length_is_rejected(path):
    with pytest.raises(ValueError):
        chematic.from_smiles('F/C(Cl)=C(Br)/I').rdkit_alkene_path_stereo(path)

@pytest.mark.parametrize('path', [None,1,'0135',[0,1,3,'5'],[0,1,3,5.5]])
def test_noninteger_inputs_are_rejected(path):
    with pytest.raises(TypeError):
        chematic.from_smiles('F/C(Cl)=C(Br)/I').rdkit_alkene_path_stereo(path)

def test_single_centered_path_is_rejected():
    with pytest.raises(ValueError, match='double bond'):
        chematic.from_smiles('CCCC').rdkit_alkene_path_stereo([0,1,2,3])

def test_a_mentioned_hydrogen_must_be_explicit():
    mol = chematic.from_smiles('F/C=C/Cl')
    with pytest.raises(ValueError, match='out of range'):
        mol.rdkit_alkene_path_stereo([4,1,2,3])
    assert mol.add_hydrogens().rdkit_alkene_path_stereo([4,1,2,3]) == 'cis'

@pytest.mark.parametrize('order', [[5,4,3,2,1,0],[2,0,4,1,5,3],[1,2,3,4,5,0]])
def test_atom_renumbering_preserves_physical_path_relation(order):
    mol = chematic.from_smiles('F/C(Cl)=C(Br)/I').renumber_atoms(order)
    path = [order.index(i) for i in [0,1,3,5]]
    assert mol.rdkit_alkene_path_stereo(path) == 'trans'
    assert mol.rdkit_alkene_path_stereo(path[::-1]) == 'trans'
