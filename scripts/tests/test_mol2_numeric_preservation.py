"""Controls ensure the independent numeric observer detects real damage."""
import copy
import json
from pathlib import Path
import unittest

from scripts.validation.check_mol2_numeric_preservation import expected_snapshot, observe, raw_snapshot


class ObserverTests(unittest.TestCase):
    def setUp(self):
        packet = Path(__file__).resolve().parents[2] / 'validation/mol2_numeric_preservation/cases.json'
        self.cases = json.loads(packet.read_text())['cases']

    def test_source_tokens_match_authored_oracle(self):
        for case in self.cases:
            self.assertEqual(raw_snapshot(case['input']), expected_snapshot(case))

    def test_coordinate_rounding_is_detected(self):
        case = self.cases[1]
        damaged = case['input'].replace('1.2345678901234567', '1.2346')
        result = observe(expected_snapshot(case), damaged, raw_snapshot)
        self.assertFalse(result['coordinates_exact'])
        self.assertTrue(result['partial_charges_exact'])

    def test_charge_and_signed_zero_are_checked_independently(self):
        for index, old, new in [(3, '-0.12345678901234566', '-0.1235'), (4, '-0.0', '0.0')]:
            case = self.cases[index]
            damaged = case['input'].replace(old, new)
            result = observe(expected_snapshot(case), damaged, raw_snapshot)
            self.assertFalse(result['partial_charges_exact'])
        self.assertFalse(result['coordinates_exact'])

    def test_wrong_element_or_bond_order_fails_graph_check(self):
        case = self.cases[7]
        for text in [case['input'].replace('C.3', 'N.3', 1),
                     case['input'].replace('1 1 2 1\n', '1 1 2 2\n')]:
            result = observe(expected_snapshot(case), text, raw_snapshot)
            self.assertFalse(result['indexed_graph_exact'])

    def test_missing_or_extra_atom_rows_fail(self):
        case = self.cases[7]
        expected = expected_snapshot(case)
        row = next(line for line in case['input'].splitlines() if line.startswith('2 C2'))
        for text in [case['input'].replace(row+'\n', ''), case['input'].replace(row, row+'\n'+row), '']:
            result = observe(expected, text, raw_snapshot)
            self.assertFalse(result['coordinates_exact'])
            self.assertFalse(result['partial_charges_exact'])

    def test_changed_observed_value_cannot_pass_by_float_tolerance(self):
        case = self.cases[1]
        snapshot = raw_snapshot(case['input'])
        damaged = copy.deepcopy(snapshot)
        damaged['coordinates'][0][0] = '3ff3c0ca428c59fa'
        result = observe(expected_snapshot(case), case['input'], lambda _: damaged)
        self.assertFalse(result['coordinates_exact'])


if __name__ == '__main__':
    unittest.main()
