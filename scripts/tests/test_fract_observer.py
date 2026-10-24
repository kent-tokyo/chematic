"""Evidence controls: changed fields and crashes cannot become preservation."""
import importlib.util
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location('fract_observer', Path(__file__).resolve().parents[1] / 'validation/check_fract_preservation.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
SOURCE = 'title\n4 5 6 90 90 90\nC 0.125(02) -0.0 0.75\n'

class ObserverControls(unittest.TestCase):
    def test_equal_state(self):
        self.assertTrue(M.comparisons(M.observe(SOURCE), M.observe(SOURCE))['declared_state_exact'])

    def test_rounded_coordinate_fails(self):
        self.assertFalse(M.comparisons(M.observe(SOURCE.replace('0.125', '0.12')), M.observe(SOURCE))['nominal_numbers_exact'])

    def test_lost_uncertainty_fails(self):
        actual = M.observe(SOURCE.replace('(02)', ''))
        self.assertTrue(M.comparisons(actual, M.observe(SOURCE))['nominal_numbers_exact'])
        self.assertFalse(M.comparisons(actual, M.observe(SOURCE))['declared_state_exact'])

    def test_lost_signed_zero_fails(self):
        self.assertFalse(M.comparisons(M.observe(SOURCE.replace('-0.0', '0.0')), M.observe(SOURCE))['declared_state_exact'])

    def test_empty_output_is_unobservable(self):
        with self.assertRaises(ValueError):
            M.observe('')

    def test_crash_partial_output_and_wrong_prefix_are_not_typed_refusal(self):
        receipt = {'returncode': 2, 'stdout': '', 'stderr': 'chematic: fract line 1: invalid'}
        self.assertTrue(M.cli_refusal(receipt))
        for key, value in [('returncode', 101), ('returncode', -11), ('stdout', 'partial'), ('stderr', 'panic: bad')]:
            self.assertFalse(M.cli_refusal(dict(receipt, **{key: value})))

if __name__ == '__main__':
    unittest.main()
