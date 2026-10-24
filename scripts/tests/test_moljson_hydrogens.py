"""Damage controls for the H-count and external-graph observation gate."""
import json
import unittest

from scripts.validation.check_moljson_hydrogens import observe


class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.case = dict(expected_hydrogens=[1], expected_smiles='expected')
        self.row = dict(output=json.dumps(dict(atoms=[dict(hydrogens=1)])),
                        smiles='expected', common_smiles='expected')

    @staticmethod
    def reader(text):
        if text == 'invalid':
            raise ValueError('cannot read')
        return dict(graph=text)

    def observation(self):
        return observe(self.case, self.row, self.reader)

    def test_unchanged_output_passes(self):
        self.assertTrue(self.observation()['preserved'])

    def test_changed_h_count_fails_even_with_same_smiles(self):
        for count in [0, 2, True, 1.0, '1', None]:
            self.row['output'] = json.dumps(dict(atoms=[dict(hydrogens=count)]))
            self.assertFalse(self.observation()['non_node_h_exact'])
            self.assertFalse(self.observation()['preserved'])

    def test_missing_or_extra_atoms_fail(self):
        for atoms in [[], [dict(hydrogens=1), dict(hydrogens=0)], [{}]]:
            self.row['output'] = json.dumps(dict(atoms=atoms))
            self.assertFalse(self.observation()['preserved'])

    def test_changed_smiles_graph_fails_independently(self):
        self.row['smiles'] = 'different'
        result = self.observation()
        self.assertFalse(result['smiles_graph_exact'])
        self.assertTrue(result['non_node_h_exact'])
        self.assertFalse(result['preserved'])

    def test_changed_common_route_fails_independently(self):
        self.row['common_smiles'] = 'different'
        self.assertFalse(self.observation()['common_smiles_graph_exact'])
        self.assertFalse(self.observation()['preserved'])

    def test_unreadable_smiles_fails(self):
        for field in ['smiles', 'common_smiles']:
            self.row[field] = 'invalid'
            self.assertFalse(self.observation()['preserved'])


if __name__ == '__main__':
    unittest.main()
