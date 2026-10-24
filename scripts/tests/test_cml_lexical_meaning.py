"""Negative controls for the bounded comparison, without an RDKit dependency."""
from copy import deepcopy
import unittest

from scripts.validation.check_cml_lexical_meaning import compare_fields, observe, source_fields


class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.expected = dict(numbers=[6, 8], charges=[0, 0], isotopes=[13, 18],
                             hydrogens=[2, 0], atom_maps=[0, 0], bonds=[[0, 1, 2]])

    def test_identity_and_permuted_atoms_preserve(self):
        self.assertTrue(all(compare_fields(self.expected, self.expected).values()))
        actual = {k: list(reversed(v)) for k, v in self.expected.items() if k != 'bonds'}
        actual['bonds'] = [[1, 0, 2]]
        self.assertTrue(all(compare_fields(self.expected, actual).values()))

    def test_each_scalar_channel_is_required(self):
        for key in ('charges', 'isotopes', 'hydrogens', 'atom_maps'):
            actual = deepcopy(self.expected)
            actual[key][0] += 1
            self.assertFalse(compare_fields(self.expected, actual)[key], key)

    def test_atom_loss_or_injection_does_not_preserve(self):
        for numbers in ([6], [6, 8, 6], [6, 6], [7, 8]):
            actual = dict(self.expected, numbers=numbers)
            self.assertFalse(any(compare_fields(self.expected, actual).values()))

    def test_bond_order_loss_and_bad_endpoints_do_not_preserve(self):
        for bonds in ([], [[0, 1, 1]], [[0, 2, 2]]):
            self.assertFalse(compare_fields(self.expected, dict(self.expected, bonds=bonds))['bonds'])

    def test_no_output_and_multiple_records_are_not_wins(self):
        for text in (None, '', '   ', '[13CH2]=[18O]\nC\n'):
            self.assertFalse(observe(self.expected, text)['preserved'])

    def test_xml_source_observer_excludes_opaque_fake_atoms(self):
        source = '''<molecule><atomArray><atom id="a1" elementType="C" isotopeNumber="13" hydrogenCount="2"/>
          <atom id="a2" elementType="O" isotopeNumber="18" hydrogenCount="0"/></atomArray>
          <!-- <atom id="fake" elementType="N"/> -->
          <![CDATA[<atom id="fake" elementType="N"/>]]>
          <?memo <atom id="fake" elementType="N"/> ?>
          <bondArray><bond atomRefs2="a1 a2" order="2"/></bondArray></molecule>'''
        self.assertEqual(source_fields(source), self.expected)


if __name__ == '__main__':
    unittest.main()
