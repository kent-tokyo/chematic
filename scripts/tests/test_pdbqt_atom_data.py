"""The external observer must reject truncated or malformed success output."""
import unittest
from copy import deepcopy

from scripts.check_interchange_preservation import ContractError
from scripts.validation.check_pdbqt_atom_data import read_probe, numeric_data_equal


class PdbqtAtomDataObserverTests(unittest.TestCase):
    def test_preserves_order_and_full_numeric_tokens(self):
        self.assertEqual(read_probe('2\n30 -0 1 2 0.123456789012345\n17 3 4 5 -1\n'),
                         dict(atomic_numbers=[30, 17], coords=[[-0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
                              partial_charges=[0.123456789012345, -1.0]))

    def test_observer_detects_zero_sign_precision_and_atom_order_loss(self):
        expected = dict(atomic_numbers=[6, 17], coords=[[-0.0, 0, 0], [1, 2, 3]],
                        partial_charges=[0.123456789012345, -1])
        self.assertTrue(numeric_data_equal(deepcopy(expected), expected))
        for field, value in [('coords', [[0.0, 0, 0], [1, 2, 3]]),
                             ('partial_charges', [0.123, -1]),
                             ('atomic_numbers', [17, 6])]:
            observed = deepcopy(expected)
            observed[field] = value
            self.assertFalse(numeric_data_equal(observed, expected))

    def test_rejects_partial_or_extra_atoms(self):
        for text in ['', '2\n30 0 0 0 2\n', '0\n30 0 0 0 2\n', '1\n30 0 0 0\n']:
            with self.subTest(text=text), self.assertRaises(ContractError):
                read_probe(text)


if __name__ == '__main__':
    unittest.main()
