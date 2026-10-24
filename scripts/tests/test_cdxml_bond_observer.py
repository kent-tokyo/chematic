"""Counterfactual graph and refusal controls, with no RDKit dependency."""
from copy import deepcopy
import json
from pathlib import Path
import unittest
from scripts.validation.cdxml_bond_observer import cdxml_fields, compare, observe

CASES = json.loads((Path(__file__).parents[2] / 'validation/cdxml_bond_meaning/cases.json').read_text())['cases']


class CdxmlBondObserverTests(unittest.TestCase):
    def setUp(self):
        self.case = next(c for c in CASES if c['id'] == 'dative_forward')
        self.atoms, self.bond = cdxml_fields(self.case['input'])

    def test_authored_xml_fields_match_the_graph_oracle(self):
        for c in CASES:
            if c['group'] != 'refusal':
                atoms, bond = cdxml_fields(c['input'])
                self.assertTrue(all(v for v in compare(c['expected'], atoms, bond).values() if v is not None), c['id'])

    def test_dative_and_single_are_distinct_despite_equal_numeric_order(self):
        fields = compare(self.case['expected'], self.atoms, dict(self.bond, kind='single'))
        self.assertTrue(fields['endpoints'])
        self.assertFalse(fields['bond_class'])
        self.assertFalse(fields['dative_direction'])

    def test_reversed_dative_direction_is_loss(self):
        fields = compare(self.case['expected'], self.atoms, dict(self.bond, begin=1, end=0))
        self.assertTrue(fields['bond_class'])
        self.assertTrue(fields['endpoints'])
        self.assertFalse(fields['dative_direction'])

    def test_undirected_quadruple_and_atom_reordering_are_equivalent(self):
        c = next(c for c in CASES if c['id'] == 'quadruple')
        atoms, bond = cdxml_fields(c['input'])
        self.assertTrue(all(v for v in compare(c['expected'], atoms, dict(bond, begin=1, end=0)).values() if v is not None))
        self.assertTrue(all(v for v in compare(self.case['expected'], list(reversed(self.atoms)), dict(self.bond, begin=1, end=0)).values() if v is not None))

    def test_atom_labels_and_count_are_required(self):
        for key in ('number', 'isotope', 'charge', 'hydrogens'):
            atoms = deepcopy(self.atoms)
            atoms[0][key] += 1
            self.assertFalse(compare(self.case['expected'], atoms, self.bond)['atoms'])
        self.assertFalse(compare(self.case['expected'], self.atoms[:1], self.bond)['atoms'])

    def test_no_output_invalid_xml_and_missing_bonds_are_not_wins(self):
        for output in ('', '<CDXML>', self.case['input'].replace('<b B="1" E="2" Order="dative"/>', '')):
            self.assertFalse(observe(self.case, dict(output=output), 'cdxml')['preserved'])

    def test_output_with_a_warning_is_separate_from_typed_profile_refusal(self):
        c = next(c for c in CASES if c['group'] == 'refusal')
        rejected = observe(c, dict(error='unsupported CDXML bond order: 2.5'), 'smiles')
        self.assertTrue(rejected['profile_refused'])
        accepted = observe(c, dict(output='CO', stderr='warning: changed bond'), 'smiles')
        self.assertEqual(accepted['outcome'], 'output_outside_profile')
        self.assertFalse(accepted['profile_refused'])
        self.assertFalse(observe(c, dict(exit_code=0, output=''), 'smiles')['profile_refused'])


if __name__ == '__main__':
    unittest.main()
