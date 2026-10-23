"""Observer controls prevent rounded, truncated or relabelled output matching."""
import copy
import json
from pathlib import Path
import unittest

from scripts.validation.orca_checked_observer import observe, outcome

ROOT=Path(__file__).parents[2]
CASE=json.loads((ROOT/'validation/orca-checked-2026-10-24/cases.json').read_text())['cases'][0]
SOURCE='# authored ORCA geometry\n! Angs\n* xyz 0 1\nC 1.2345678901234567 -1.2345678901234567 0.41152263004115225\nO 2.5 0.0 3.0\n*\n'


class OrcaWriterObserver(unittest.TestCase):
    def result(self, text=SOURCE):
        return outcome(CASE,dict(ok=True,output=text))

    def test_numbers_and_coordinate_count_are_independent_of_metadata(self):
        self.assertEqual(self.result()['status'],'fields_match')
        got=self.result(SOURCE.replace('1.2345678901234567','1.234568'))
        self.assertFalse(got['numeric_match']);self.assertTrue(got['metadata_match'])
        self.assertEqual(self.result(SOURCE.replace('O 2.5 0.0 3.0\n',''))['status'],'changed')
        self.assertEqual(self.result(SOURCE.replace('*\n','O 2.5 0.0 3.0\n*\n'))['status'],'changed')

    def test_units_state_species_and_freeze_flags_are_not_numeric_matches(self):
        for before,after in [('Angs','Bohrs'),('xyz 0 1','xyz -1 2'),('C 1.','N 1.'),('3.0\n','3.0$\n')]:
            got=self.result(SOURCE.replace(before,after))
            self.assertEqual(got['status'],'changed');self.assertFalse(got['metadata_match'])

    def test_comments_keywords_opaque_blocks_and_extras_are_checked(self):
        for text in [SOURCE.replace('# authored ORCA geometry','# different'),SOURCE.replace('! Angs','! Angs TightSCF'),SOURCE.replace('! Angs\n','! Angs\n%scf\nMaxIter 3\nend\n'),SOURCE.replace('3.0\n','3.0 opaque\n')]:
            self.assertFalse(self.result(text)['metadata_match'])

    def test_signed_zero_nonfinite_and_partial_output(self):
        case=copy.deepcopy(CASE);case['coordinate_bits'][1][1]='8000000000000000'
        self.assertFalse(outcome(case,dict(ok=True,output=SOURCE))['numeric_match'])
        self.assertEqual(self.result(SOURCE.replace('2.5','NaN'))['status'],'unreadable')
        self.assertEqual(self.result(SOURCE[:-2])['status'],'unreadable')
        with self.assertRaises(ValueError):outcome(CASE,dict(ok=False,error='refused',output=SOURCE))
        with self.assertRaises(ValueError):observe(SOURCE+'unexpected\n')


if __name__=='__main__':unittest.main()
