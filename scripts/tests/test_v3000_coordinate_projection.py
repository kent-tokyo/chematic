"""Coordinate retention and indexed graph checks are separate observations."""
import json
import unittest
from scripts.validation.check_v3000_coordinate_projection import observe

class V3000ProjectionObserverTests(unittest.TestCase):
    def setUp(self):
        self.expected=dict(numbers=[6,8],formal_charges=[0,-1],isotopes=[13,0],bonds=[[0,1,1]],coordinates=[['0.123456789012345','-0.0','0.000012345678901'],['1.0','2.0','3.0']])
        self.value=dict(atoms=dict(elements=dict(number=[6,8]),formalCharges=[0,-1],isotopes=[13,0],coords={'3d':[float(v) for row in self.expected['coordinates'] for v in row]}),bonds=dict(connections=dict(index=[0,1]),order=[1]))
    def test_full_projection_and_reversed_bond_endpoints(self):
        self.assertTrue(all(observe(self.expected,json.dumps(self.value)).values()))
        self.value['bonds']['connections']['index']=[1,0]
        self.assertTrue(observe(self.expected,json.dumps(self.value))['indexed_graph'])
    def test_precision_signed_zero_and_z_loss(self):
        for axis,replacement in [(0,0.1235),(1,0.0),(2,0.0)]:
            value=json.loads(json.dumps(self.value));value['atoms']['coords']['3d'][axis]=replacement
            result=observe(self.expected,json.dumps(value))
            self.assertTrue(result['indexed_graph']);self.assertFalse(result['coordinate_bits'])
    def test_bond_charge_element_and_isotope_changes(self):
        for key,value in [('formalCharges',[0,0]),('elements',dict(number=[8,6]))]:
            out=json.loads(json.dumps(self.value));out['atoms'][key]=value
            self.assertFalse(observe(self.expected,json.dumps(out))['indexed_graph'])
        self.value['bonds']['order']=[2]
        self.assertFalse(observe(self.expected,json.dumps(self.value))['indexed_graph'])
        self.value['atoms']['isotopes']=[0,0]
        self.assertFalse(observe(self.expected,json.dumps(self.value))['isotopes'])
    def test_missing_coordinates_are_separate_from_valid_graph(self):
        del self.value['atoms']['coords']
        result=observe(self.expected,json.dumps(self.value))
        self.assertTrue(result['indexed_graph']);self.assertFalse(result['coordinates_present'])
        self.assertFalse(observe(self.expected,'')['json_object'])
    def test_absent_empty_bonds_are_not_scored_as_loss(self):
        self.expected['bonds']=[];del self.value['bonds']
        self.assertTrue(observe(self.expected,json.dumps(self.value))['indexed_graph'])
    def test_nonfinite_json_and_duplicate_keys_are_refused(self):
        self.assertFalse(observe(self.expected,'{"atoms":NaN}')['json_object'])
        self.assertFalse(observe(self.expected,'{"atoms":{},"atoms":{}}')['json_object'])
    def test_coordinate_strings_booleans_and_overflow_do_not_pass(self):
        for value in ['0.123456789012345', True, float('inf')]:
            output=json.loads(json.dumps(self.value));output['atoms']['coords']['3d'][0]=value
            self.assertFalse(observe(self.expected,json.dumps(output))['coordinate_bits'])
        output=json.dumps(self.value).replace('0.123456789012345', '1e999')
        self.assertFalse(observe(self.expected,output)['coordinate_bits'])

if __name__=='__main__':unittest.main()
