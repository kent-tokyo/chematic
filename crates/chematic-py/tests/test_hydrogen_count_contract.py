"""Open Babel hydrogen-accounting cases; the formula implementation is unchanged."""
from pathlib import Path

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
CASES = [(s,int(n)) for row in (ROOT / 'validation/openbabel_issues/hydrogen-counts/cases.tsv').read_text(encoding='utf-8').splitlines() if not row.startswith('#') for s,n in [row.split('\t')]]


@pytest.mark.parametrize('smiles,expected', CASES)
def test_total_hydrogens_are_independent_of_their_representation(smiles,expected):
    mol = chematic.from_smiles(smiles)
    before = (mol.smiles,mol.atom_count,mol.formula)
    assert mol.num_hydrogens == expected
    assert mol.num_hydrogens == expected
    assert (mol.smiles,mol.atom_count,mol.formula) == before
    assert chematic.from_smiles(mol.smiles).num_hydrogens == expected


@pytest.mark.parametrize('smiles', ['C','[CH4]','[NH4+]','[OH2].[OH2]','CCO','c1ccccc1','[nH]1cccc1','N[C@@H](C)O','[2H]O[3H]','[H-][Al+3]([H-])[H-]'])
def test_expand_fold_keeps_count_and_formula(smiles):
    mol = chematic.from_smiles(smiles)
    expected = (mol.num_hydrogens,mol.formula)
    expanded = mol.add_hydrogens()
    folded = expanded.remove_hydrogens()
    assert (expanded.num_hydrogens,expanded.formula) == expected
    assert (folded.num_hydrogens,folded.formula) == expected
