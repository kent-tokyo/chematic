import json
import unittest
from scripts.validation.check_mol2_cjson_common import observe
class Mol2CjsonObserverTests(unittest.TestCase):
    def setUp(self):
        self.expected=dict(numbers=[6,7],formal_charges=[0,1],bonds=[[0,1,1]],coordinates=[['0.123456789012345','-0.0','0.000012345678901'],['1','2','3']],name='ligand',method='USER_CHARGES',partial_charges=['0.7123456789012345','-0.0'])
        self.value=dict(atoms=dict(elements=dict(number=[6,7]),formalCharges=[0,1],coords={'3d':[float(v) for row in self.expected['coordinates'] for v in row]}),bonds=dict(connections=dict(index=[0,1]),order=[1]),properties=dict(name='ligand'),partialCharges=dict(USER_CHARGES=[0.7123456789012345,-0.0]))
    def test_all_fields_and_reversed_bond_endpoints(self):
        self.assertTrue(all(observe(self.expected,json.dumps(self.value)).values()))
        self.value['bonds']['connections']['index']=[1,0]
        self.assertTrue(observe(self.expected,json.dumps(self.value))['indexed_graph'])
    def test_method_and_numeric_values_are_independent(self):
        self.value['partialCharges']={'renamed':[0.7123456789012345,-0.0]}
        result=observe(self.expected,json.dumps(self.value))
        self.assertTrue(result['partial_values']);self.assertFalse(result['charge_method'])
        self.value['partialCharges']={'USER_CHARGES':[0.7123,0.0]}
        result=observe(self.expected,json.dumps(self.value))
        self.assertFalse(result['partial_values']);self.assertTrue(result['charge_method'])
    def test_partial_charges_are_not_formal_charges(self):
        self.value['atoms']['formalCharges']=[1,0]
        result=observe(self.expected,json.dumps(self.value))
        self.assertFalse(result['indexed_graph']);self.assertTrue(result['partial_values'])
    def test_coordinate_precision_zero_sign_z_and_name(self):
        for index,value in [(0,0.1235),(1,0.0),(2,0.0)]:
            output=json.loads(json.dumps(self.value));output['atoms']['coords']['3d'][index]=value
            self.assertFalse(observe(self.expected,json.dumps(output))['coordinates'])
        self.value['properties']['name']='changed'
        self.assertFalse(observe(self.expected,json.dumps(self.value))['name'])
    def test_absence_is_distinct_from_inventing_a_dataset(self):
        self.expected['method']=None;self.expected['partial_charges']=None
        self.assertFalse(observe(self.expected,json.dumps(self.value))['partial_values'])
        del self.value['partialCharges']
        self.assertTrue(observe(self.expected,json.dumps(self.value))['partial_values'])
    def test_malformed_json_and_nonnumeric_values_cannot_pass(self):
        for text in ['', '{"atoms":NaN}', '{"atoms":{},"atoms":{}}']:
            self.assertFalse(observe(self.expected,text)['json_object'])
        for bad in ['0.7123456789012345',True,float('inf')]:
            self.value['partialCharges']['USER_CHARGES'][0]=bad
            self.assertFalse(observe(self.expected,json.dumps(self.value))['partial_values'])
    def test_bonds_and_atom_correspondence_are_checked(self):
        self.value['bonds']['order']=[2]
        self.assertFalse(observe(self.expected,json.dumps(self.value))['indexed_graph'])
        self.value['atoms']['elements']['number']=[7,6]
        result=observe(self.expected,json.dumps(self.value))
        self.assertFalse(result['coordinates']);self.assertFalse(result['partial_values'])
if __name__=='__main__':unittest.main()
