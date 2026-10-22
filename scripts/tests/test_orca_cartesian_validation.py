import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).parents[2]
spec=importlib.util.spec_from_file_location('orca_validation',ROOT/'scripts/validation/validate_orca_cartesian.py');validation=importlib.util.module_from_spec(spec);spec.loader.exec_module(validation)
CASE=json.loads((ROOT/'validation/orca-cartesian-2026-10-22/cases.json').read_text())['cases'][0]


class OrcaGeometryObservation(unittest.TestCase):
    def result(self):
        fields=copy.deepcopy(CASE['fields']);fields.pop('coordinate_bits');fields['schema']='chematic.orca-cartesian.v1'
        return dict(ok=True,value=fields)

    def test_units_are_required_and_must_match_without_rescaling(self):
        result=self.result();self.assertEqual(validation.geometry_result(CASE,result),'verified')
        result['value']['coordinate_units']='bohr';self.assertEqual(validation.geometry_result(CASE,result),'mismatched')
        del result['value']['schema'];self.assertEqual(validation.geometry_result(CASE,result),'untagged')

    def test_zero_sign_and_source_state_are_checked(self):
        for field,value in [('charge',1),('multiplicity',0),('units_explicit',0),('frozen',[[1,False,True],[False,False,False]])]:
            result=self.result();result['value'][field]=value
            self.assertEqual(validation.geometry_result(CASE,result),'mismatched')
        result=self.result();result['value']['coords'][0][0]=0.
        self.assertEqual(validation.geometry_result(CASE,result),'mismatched')

    def test_no_partial_geometry_or_extra_rows_can_match(self):
        result=self.result();result['value']['coords']=result['value']['coords'][:1]
        self.assertEqual(validation.geometry_result(CASE,result),'mismatched')
        result=self.result();result['value']['coords']+=result['value']['coords']
        self.assertEqual(validation.geometry_result(CASE,result),'mismatched')
        with self.assertRaises(ValueError):validation.geometry_result(CASE,dict(ok=False,error='refused',value={'coords':[]}))

    def test_raw_json_transport_preserves_negative_zero(self):
        result=self.result();result['value']=json.dumps(result['value'])
        self.assertEqual(validation.geometry_result(CASE,result),'verified')
        result['value']=result['value'].replace('-0.0','0.0')
        self.assertEqual(validation.geometry_result(CASE,result),'mismatched')


if __name__=='__main__':unittest.main()
