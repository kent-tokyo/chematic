"""Negative controls for numeric and missing-field preservation observations."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.validation.mmcif_numeric_observer import bits, observe

CASES = json.loads((Path(__file__).parents[2] /
                    'validation/mmcif_numeric_preservation/cases.json').read_text())['cases']


def edit_row(source, index, replacement):
    lines = source.splitlines()
    at = next(i for i, line in enumerate(lines) if line.startswith('ATOM '))
    row = lines[at].split()
    row[index] = replacement
    lines[at] = ' '.join(row)
    return '\n'.join(lines) + '\n'


class MmcifNumericObserverTests(unittest.TestCase):
    def setUp(self):
        self.case = deepcopy(CASES[0])

    def test_authored_sources_and_ordinary_round_trip_profile(self):
        for case in CASES:
            self.assertTrue(observe(case['expected'], case['input'])['preserved'], case['id'])

    def test_coordinate_rounding_is_not_hidden_by_a_tolerance(self):
        result = observe(self.case['expected'], edit_row(self.case['input'], 10, '1.0000000000000002'))
        self.assertFalse(result['fields']['cartesian'])
        self.assertFalse(result['preserved'])

    def test_missing_occupancy_is_loss_even_when_the_default_would_match(self):
        source = self.case['input'].replace('_atom_site.occupancy\n', '')
        lines = source.splitlines()
        at = next(i for i, line in enumerate(lines) if line.startswith('ATOM '))
        row = lines[at].split()
        del row[13]
        lines[at] = ' '.join(row)
        result = observe(self.case['expected'], '\n'.join(lines))
        self.assertTrue(result['fields']['cartesian'])
        self.assertFalse(result['fields']['occupancy'])
        self.assertIn('_atom_site.occupancy', result['missing_fields'])

    def test_scalar_changes_and_atom_identity_are_required(self):
        for index, value, field in [(13, '0.5', 'occupancy'), (14, '9.999', 'b_iso'), (2, 'C', 'atom_identity')]:
            result = observe(self.case['expected'], edit_row(self.case['input'], index, value))
            self.assertFalse(result['fields'][field])
            self.assertFalse(result['preserved'])

    def test_cell_rounding_or_missing_parameters_are_losses(self):
        case = CASES[2]
        for replacement in ('_cell.length_a 61.25000000000001\n', ''):
            source = case['input'].replace('_cell.length_a 61.25\n', replacement)
            self.assertFalse(observe(case['expected'], source)['fields']['cell'])
        self.assertIsNone(observe(self.case['expected'], self.case['input'])['fields']['cell'])

    def test_nonfinite_and_wrong_numeric_types_are_not_numbers(self):
        for value in ('NaN', 'Infinity', '1e400'):
            self.assertFalse(observe(self.case['expected'], edit_row(self.case['input'], 10, value))['preserved'])
        for value in (True, None, float('inf')):
            with self.assertRaises(ValueError):
                bits([value])

    def test_empty_malformed_and_multiple_documents_are_not_wins(self):
        for source in (None, '', self.case['input'] + self.case['input'],
                       self.case['input'].replace('_atom_site.id\n', '_atom_site.type_symbol\n'),
                       self.case['input'].replace('ATOM 1 ', 'ATOM '), ';outside\nprofile\n;\n'):
            self.assertFalse(observe(self.case['expected'], source)['preserved'])


if __name__ == '__main__':
    unittest.main()
