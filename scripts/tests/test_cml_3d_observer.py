import copy
import json
import unittest
from scripts.validation.check_cml_3d import coordinate_bits, graph_matches, observe

CASE={'group':'valid','coords':[[1.25,2.5,-0.0]],'numbers':[8],'bond_indices':[],'bond_orders':[]}
DATA={'atoms':{'elements':{'number':[8]},'coords':{'3d':[1.25,2.5,-0.0]}},'bonds':{'connections':{'index':[]},'order':[]}}
class Cml3dObserverTests(unittest.TestCase):
    def test_source_order_and_signed_zero_are_observed(self):
        self.assertTrue(observe(CASE,{'output':json.dumps(DATA)})['coordinate_bits'])
        changed=copy.deepcopy(DATA);changed['atoms']['coords']['3d'][2]=0.0
        self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['coordinate_bits'])
        self.assertNotEqual(coordinate_bits([-0.0]),coordinate_bits([0.0]))
    def test_missing_and_reordered_coordinates_are_losses(self):
        changed=copy.deepcopy(DATA);changed['atoms']['coords']['3d']=[2.5,1.25,-0.0]
        self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['coordinate_bits'])
        del changed['atoms']['coords']
        self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['coordinate_bits'])
    def test_absence_is_not_zero_geometry(self):
        case={**CASE,'group':'absent','coords':None}
        self.assertFalse(observe(case,{'output':json.dumps(DATA)})['absence_retained'])
        changed=copy.deepcopy(DATA);del changed['atoms']['coords']
        self.assertTrue(observe(case,{'output':json.dumps(changed)})['absence_retained'])
    def test_nonempty_invalid_output_is_not_safe_refusal(self):
        result=observe({**CASE,'group':'invalid'},{'output':'{"atoms":'})
        self.assertEqual(result['outcome'],'invalid_output')
        self.assertNotIn('error_status',result)
    def test_error_empty_and_json_acceptance_are_separate(self):
        case={**CASE,'group':'invalid'}
        self.assertTrue(observe(case,{'error':'refused'})['error_status'])
        self.assertFalse(observe(case,{'output':'','exit_code':0})['error_status'])
        self.assertEqual(observe(case,{'output':json.dumps(DATA)})['outcome'],'accepted_outside_profile')
    def test_selected_graph_checks_elements_and_indexed_bonds(self):
        case={**CASE,'numbers':[8,6],'bond_indices':[1,0],'bond_orders':[2]}
        data={'atoms':{'elements':{'number':[8,6]}},'bonds':{'connections':{'index':[0,1]},'order':[2]}}
        self.assertTrue(graph_matches(data,case))
        data['bonds']['order']=[1]
        self.assertFalse(graph_matches(data,case))
        data['bonds']['order']=[]
        self.assertFalse(graph_matches(data,case))
    def test_invalid_coordinate_types_are_not_observed_as_numbers(self):
        for values in [[True], [None], [float('inf')]]:
            with self.assertRaises(ValueError): coordinate_bits(values)
if __name__=='__main__':unittest.main()
