"""Isotope labels, implicit/explicit H and total charge stay independent."""
import json
from pathlib import Path

import chematic
import pytest

CASES=json.loads((Path(__file__).parents[3]/'validation/openbabel-1464-2026-10-24/cases.json').read_text())['cases']


@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_authored_isotope_formula(case):
    mol=chematic.from_smiles(case['smiles'])
    before=mol.to_moljson();legacy=mol.formula
    assert mol.formula_with_isotopes==case['expected_formula']
    assert mol.formula==legacy
    assert mol.to_moljson()==before


def test_formula_does_not_encode_label_position_or_component_order():
    for a,b in [('C[13CH3]','[13CH3]C'),('[2H]O[3H]','[3H]O[2H]'),('[15NH4+].[37Cl-]','[37Cl-].[15NH4+]')]:
        assert chematic.from_smiles(a).formula_with_isotopes==chematic.from_smiles(b).formula_with_isotopes
