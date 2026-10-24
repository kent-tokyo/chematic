"""Aromatic correction must not hide a carbon's fifth explicit valence."""
import json
from pathlib import Path

import chematic
import pytest

CASES=json.loads((Path(__file__).parents[3]/'validation/rdkit-9735-2026-10-24/cases.json').read_text(encoding="utf-8"))['cases']


@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_strict_aromatic_valence(case):
    mol=chematic.from_smiles(case['smiles'])
    problems=chematic.rdkit_detect_chemistry_problems(case['smiles'])
    if case['expected']=='accept':
        first=mol.rdkit_smiles
        assert chematic.from_smiles(first).rdkit_smiles==first
        assert chematic.rdkit_smiles_checked(mol)==first
        assert problems==[]
        return
    with pytest.raises(ValueError,match='greater than permitted'):
        _=mol.rdkit_smiles
    with pytest.raises(ValueError,match='greater than permitted'):
        chematic.rdkit_smiles_checked(mol)
    assert ('AtomValenceException',[case['atom']]) in problems
