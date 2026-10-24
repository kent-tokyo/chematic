"""Negative controls for the packet's unique-element correspondence rule."""
import copy
import unittest

from scripts.validation.check_v2000_property_precedence import normalize, observe


class PropertyPrecedenceObserverTests(unittest.TestCase):
    expected = dict(numbers=[6, 7, 8], charges=[0, 1, 0], isotopes=[0, 0, 18],
                    atom_maps=[7, 9, 11], bonds=[[0, 1, 1], [1, 2, 2]])

    def test_atom_reordering_and_reversed_bond_endpoints_are_equivalent(self):
        actual = {k: [self.expected[k][i] for i in [2, 0, 1]]
                  for k in ['numbers', 'charges', 'isotopes', 'atom_maps']}
        actual['bonds'] = [[2, 1, 1], [0, 2, 2]]
        self.assertEqual(normalize(self.expected), normalize(actual))
        self.assertTrue(observe(self.expected, '', lambda _: actual)['fields_exact'])

    def test_every_atom_label_is_observed_at_its_corresponding_element(self):
        for key, value in [('numbers', 9), ('charges', 1), ('isotopes', 13), ('atom_maps', 8)]:
            actual = copy.deepcopy(self.expected)
            actual[key][0] = value
            self.assertFalse(observe(self.expected, '', lambda _: actual)['fields_exact'])
        actual = dict(self.expected, atom_maps=[9, 7, 11])
        self.assertFalse(observe(self.expected, '', lambda _: actual)['fields_exact'])

    def test_connectivity_order_and_duplicate_bonds_cannot_pass(self):
        for bonds in [[[0, 2, 1], [1, 2, 2]], [[0, 1, 2], [1, 2, 2]],
                      [[0, 1, 1]], [[0, 1, 1], [0, 1, 1], [1, 2, 2]]]:
            self.assertFalse(observe(self.expected, '', lambda _: dict(self.expected, bonds=bonds))['fields_exact'])

    def test_packet_correspondence_limits_are_enforced(self):
        for changes in [dict(numbers=[6, 6, 8]), dict(charges=[0, 1]),
                        dict(isotopes=[False, 0, 18]), dict(bonds=[[0, 3, 1]]),
                        dict(bonds=[[0, 0, 1]]), dict(bonds=[[0, 1, 1.5]])]:
            result = observe(self.expected, '', lambda _: dict(self.expected, **changes))
            self.assertFalse(result['fields_exact'])
            self.assertFalse(result['readable'])

    def test_unreadable_output_never_counts_as_preserved(self):
        def broken(_):
            raise ValueError('unreadable')
        self.assertFalse(observe(self.expected, '', broken)['fields_exact'])

    def test_stale_unlisted_isotope_is_not_a_valid_oracle(self):
        stale = dict(self.expected, isotopes=[13, 15, 18])
        result = observe(self.expected, '', lambda _: stale)
        self.assertTrue(result['readable'])
        self.assertFalse(result['fields_exact'])


if __name__ == '__main__':
    unittest.main()
