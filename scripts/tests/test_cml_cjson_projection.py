"""Independent observer controls; a successfully emitted document can still lose meaning."""
import importlib.util
from pathlib import Path
import json
import unittest
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('cml_projection',ROOT/'scripts/validation/check_cml_cjson_projection.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ObserverTests(unittest.TestCase):
    def setUp(self):
        s=json.loads((ROOT/'validation/cml_cjson_projection/cases.json').read_text())['valid'][1]['input']
        self.e=m.source_fields(s)
        self.d=dict(atoms=dict(elements=dict(number=self.e['numbers']),coords={'3d':[v for p in self.e['coordinates'] for v in p]}),bonds=dict(connections=dict(index=[0,1]),order=[1]),name=self.e['name'])
    def check(self,key):self.assertFalse(m.observe(self.e,json.dumps(self.d))['fields'][key])
    def test_control(self):self.assertTrue(m.observe(self.e,json.dumps(self.d))['full'])
    def test_wrong_bond(self):self.d['bonds']['order']=[2];self.check('graph')
    def test_negative_zero(self):self.d['atoms']['coords']['3d'][0]=0.;self.check('coordinates')
    def test_precision(self):self.d['atoms']['coords']['3d'][1]=1.23456789;self.check('coordinates')
    def test_missing_name(self):del self.d['name'];self.check('name')
    def test_isotope(self):self.d['atoms']['isotopes']=[13,0];self.check('isotopes')
    def test_boolean_is_not_an_atomic_number(self):self.d['atoms']['elements']['number']=[True,8];self.check('graph')
    def test_invalid_or_empty(self):
        for text in ['', '{', '[]']:self.assertFalse(m.observe(self.e,text)['full'])
if __name__=='__main__':unittest.main()
