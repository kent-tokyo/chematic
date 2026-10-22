import copy
import json
import unittest
from scripts.validation.check_cjson_checked import compact_json,observe,read_graph,read_cml_graph

EXPECTED={'numbers':[6,8],'charges':[0,-1],'isotopes':[0,0],'indices':[0,1],'orders':[1]}
DATA={'atoms':{'elements':{'number':[6,8]},'formalCharges':[0,-1]},'bonds':{'connections':{'index':[0,1]},'order':[1]}}
CASE={'group':'valid','expected':EXPECTED}
class CjsonCheckedObserverTests(unittest.TestCase):
    def test_compaction_keeps_strings_escapes_and_duplicate_keys(self):
        text=' { "value" : "a b \\\" quote ", "x" : 1, "x" : 2 } '
        compact=compact_json(text)
        self.assertIn('a b \\\" quote ',compact)
        self.assertEqual(compact.count('"x":'),2)
        self.assertEqual(json.loads(compact),json.loads(text))
    def test_selected_graph_and_undirected_endpoints(self):
        self.assertTrue(observe(CASE,{'output':json.dumps(DATA)})['selected_graph'])
        changed=copy.deepcopy(DATA);changed['bonds']['connections']['index']=[1,0]
        self.assertTrue(observe(CASE,{'output':json.dumps(changed)})['selected_graph'])
    def test_charge_element_isotope_and_order_changes_are_detected(self):
        for field,values in [('formalCharges',[0,0]),('isotopes',[13,0])]:
            changed=copy.deepcopy(DATA);changed['atoms'][field]=values
            self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['selected_graph'])
        changed=copy.deepcopy(DATA);changed['atoms']['elements']['number']=[8,6]
        self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['selected_graph'])
        changed=copy.deepcopy(DATA);changed['bonds']['order']=[2]
        self.assertFalse(observe(CASE,{'output':json.dumps(changed)})['selected_graph'])
    def test_invalid_output_and_refusal_are_not_preservation(self):
        invalid={**CASE,'group':'invalid'}
        self.assertEqual(observe(invalid,{'error':'refused'})['outcome'],'refused_no_output')
        self.assertEqual(observe(invalid,{'output':'','exit_code':0})['outcome'],'empty_output')
        self.assertEqual(observe(invalid,{'output':'{'})['outcome'],'invalid_output')
        self.assertEqual(observe(invalid,{'output':json.dumps(DATA)})['outcome'],'accepted_outside_profile')
    def test_process_and_resource_failure_have_separate_outcomes(self):
        self.assertEqual(observe(CASE,{'resource_status':'timeout','output':''})['outcome'],'resource_failure')
        self.assertEqual(observe(CASE,{'exit_code':1,'output':'partial'})['outcome'],'error_with_output')
    def test_wrong_types_and_lengths_are_not_numeric_matches(self):
        for value in [[True,0],[0],[0.0,-1]]:
            changed=copy.deepcopy(DATA);changed['atoms']['formalCharges']=value
            with self.assertRaises(ValueError):read_graph(changed)
    def test_namespaced_cml_references_and_source_order(self):
        text='<molecule xmlns="http://www.xml-cml.org/schema"><atomArray><atom id="C" elementType="C"/><atom id="O" elementType="O" formalCharge="-1"/></atomArray><bondArray><bond atomRefs2="O C" order="1"/></bondArray></molecule>'
        self.assertTrue(observe(CASE,{'output':text,'format':'cml'})['selected_graph'])
        self.assertFalse(observe(CASE,{'output':text.replace('formalCharge="-1"','formalCharge="0"'),'format':'cml'})['selected_graph'])
        with self.assertRaises(ValueError):read_cml_graph(text.replace('id="O"','id="C"'))
        self.assertEqual(observe(CASE,{'output':text[:-1],'format':'cml'})['outcome'],'invalid_output')
        self.assertEqual(observe(CASE,{'output':text.replace('order="1"','order="0"'),'format':'cml'})['outcome'],'graph_mismatch')
if __name__=='__main__':unittest.main()
