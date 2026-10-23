import copy
import json
from pathlib import Path
import unittest

from scripts.validation.isotope_formula_observer import observe,outcome

CASES=json.loads((Path(__file__).parents[2]/'validation/openbabel-1464-2026-10-24/cases.json').read_text())['cases']
CARBON=next(c for c in CASES if c['id']=='carbon-13')
WATER=next(c for c in CASES if c['id']=='heavy-water')


class IsotopeFormulaObservation(unittest.TestCase):
    def test_mass_labels_are_distinct_even_when_element_totals_match(self):
        self.assertEqual(outcome(CARBON,dict(ok=True,formula='[13C]H4'))['status'],'composition_match')
        for text in ['CH4','[12C]H4']:
            got=outcome(CARBON,dict(ok=True,formula=text))
            self.assertEqual(got['status'],'changed');self.assertTrue(got['element_totals_match'])

    def test_counts_and_net_charge_are_independent_checks(self):
        for text in ['[13C]2H4','[13C]H3']:
            got=outcome(CARBON,dict(ok=True,formula=text))
            self.assertFalse(got['element_totals_match']);self.assertTrue(got['charge_match'])
        case=copy.deepcopy(CARBON);case['expected_charge']=1
        got=outcome(case,dict(ok=True,formula='[13C]H4'))
        self.assertTrue(got['isotope_counts_match']);self.assertFalse(got['charge_match'])
        self.assertEqual(observe('Fe++'),observe('Fe+2'))
        self.assertEqual(observe('O4S--'),observe('O4S-2'))
        with self.assertRaises(ValueError):observe('Fe+-')

    def test_hydrogen_aliases_are_equivalent_and_dysprosium_is_not_deuterium(self):
        for text in ['D2O','[2H]2O']:
            self.assertEqual(outcome(WATER,dict(ok=True,formula=text))['status'],'composition_match')
        counts,_=observe('Dy');self.assertEqual(counts,[dict(element='Dy',isotope=None,count=1)])
        self.assertEqual(observe('T'),observe('[3H]'))

    def test_truncation_unknown_text_and_partial_refusals_never_match(self):
        for text in ['[13C','[13C]H4 garbage','', '[13C]H0','[2D]O','[65536C]H4']:
            self.assertEqual(outcome(CARBON,dict(ok=True,formula=text))['status'],'unreadable')
        with self.assertRaises(ValueError):outcome(CARBON,dict(ok=False,error='refused',formula='partial'))
        with self.assertRaises(ValueError):outcome(CARBON,dict(ok=1,formula='[13C]H4'))


if __name__=='__main__':unittest.main()
