"""Reject charge totals, reordered atoms and lazy-charge substitutions as proof."""
from copy import deepcopy
import unittest

from scripts.validation.check_cjson_mol2_formal_charges import observe, read_probe


class IndexedChargeObserverTests(unittest.TestCase):
    def setUp(self):
        self.expected=dict(atomic_numbers=[11,17],formal_charges=[1,-1],
                           partial_charges=[-0.25,0.25],coords=[[0,-0.0,0.125],[1,2,3]],bonds=[])

    def test_equal_totals_are_not_per_atom_preservation(self):
        actual=deepcopy(self.expected);actual['formal_charges']=[0,0]
        self.assertEqual(sum(actual['formal_charges']),sum(self.expected['formal_charges']))
        self.assertFalse(observe(self.expected,actual)['formal_charges'])
        actual['formal_charges']=[-1,1]
        self.assertFalse(observe(self.expected,actual)['formal_charges'])

    def test_partial_formal_and_geometry_are_separate_axes(self):
        actual=deepcopy(self.expected);actual['partial_charges']=[1,-1];actual['coords'][0][1]=0.0
        result=observe(self.expected,actual)
        self.assertTrue(result['formal_charges']);self.assertTrue(result['indexed_graph'])
        self.assertFalse(result['partial_charges_exact']);self.assertFalse(result['coordinates_exact'])
        actual['atomic_numbers'].reverse()
        self.assertFalse(observe(self.expected,actual)['formal_charges'])

    def test_probe_count_shapes_nonfinite_and_bond_indices_are_checked(self):
        good='2 0\n11 1 -0.25 0 -0 0.125\n17 -1 0.25 1 2 3\n'
        self.assertTrue(all(observe(self.expected,read_probe(good)).values()))
        for bad in ['', '2 0\n11 1 0 0 0 0\n',good.replace('-0.25','NaN'),
                    good.replace('2 0','2 1')+'0 2 1\n']:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                read_probe(bad)


if __name__=='__main__':
    unittest.main()
