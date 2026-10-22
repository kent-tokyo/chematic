"""The observer must not mistake altered identity/charge/radius for preservation."""
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).parents[2]
spec=importlib.util.spec_from_file_location('pqr_checked_validation',ROOT/'scripts/validation/validate_pqr_checked.py')
validation=importlib.util.module_from_spec(spec);spec.loader.exec_module(validation)
PACKET=json.loads((ROOT/'validation/pqr-checked-2026-10-22/cases.json').read_text())


class PqrObservation(unittest.TestCase):
    def case(self,name):return next(c for c in PACKET['cases'] if c['id']==name)

    def test_all_five_numeric_fields_are_checked(self):
        c=self.case('all-fields');source=c['pqr_input']
        self.assertEqual(validation.outcome(c,dict(ok=True,output=source))['status'],'numeric_match')
        for index in range(6,11):
            tokens=source.split();tokens[index]='0.0'
            self.assertEqual(validation.outcome(c,dict(ok=True,output=' '.join(tokens)))['status'],'numeric_mismatch')

    def test_signed_zero_is_not_numerically_equal_for_this_contract(self):
        c=self.case('signed-zero')
        result=validation.outcome(c,dict(ok=True,output=c['pqr_input'].replace('-0.0','0.0')))
        self.assertEqual(result['status'],'numeric_mismatch')
        self.assertFalse(result['scalar_bits']['charge']);self.assertFalse(result['scalar_bits']['radius'])

    def test_usage_truncation_nonfinite_and_extra_atoms_are_not_matches(self):
        c=self.case('all-fields')
        for text in ['Usage: obabel','', 'ATOM 1 CA ALA A 1 0 0 0 0', c['pqr_input'].replace('1.2345678901234567','NaN')]:
            self.assertEqual(validation.outcome(c,dict(ok=True,output=text))['status'],'unreadable')
        self.assertEqual(validation.outcome(c,dict(ok=True,output=c['pqr_input']*2))['status'],'mismatched_count')

    def test_metadata_changes_are_reported_separately_from_numeric_bits(self):
        c=self.case('all-fields');tokens=c['pqr_input'].split();tokens[1]='99'
        result=validation.outcome(c,dict(ok=True,output=' '.join(tokens)))
        self.assertEqual(result['status'],'numeric_match');self.assertFalse(result['metadata_match'])

    def test_element_extension_requires_explicit_dialect(self):
        c=PACKET['openbabel_element_extension_controls'][0]
        with self.assertRaises(ValueError):validation.observe(c['source'])
        r=validation.outcome(c,dict(ok=True,output=c['source']),True)
        self.assertEqual(r['status'],'numeric_match');self.assertTrue(r['element_column_match'])

    def test_refusal_has_no_partial_output_and_flags_are_boolean(self):
        c=self.case('all-fields')
        self.assertEqual(validation.outcome(c,dict(ok=False,error='refused'))['status'],'refused')
        for r in [dict(ok=False,error='refused',output='partial'),dict(ok=1,output=c['pqr_input']),dict(ok=True)]:
            with self.assertRaises(ValueError):validation.outcome(c,r)


if __name__=='__main__':unittest.main()
