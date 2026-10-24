"""The graph observer must detect changed chemistry, connectivity and keys."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('observer', ROOT/'scripts/validation/check_moljson_graph_input.py')
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.value = dict(atoms=[dict(id='a',element='Cr',charge=0,isotope=53,hydrogens=0,aromatic=False),
                                dict(id='b',element='Cr',charge=0,isotope=None,hydrogens=0,aromatic=False)],
                          bonds=[dict(source_id='a',target_id='b',order=4,aromatic=False)])
        self.text = json.dumps(self.value)
        self.expected = OBSERVER.snapshot(self.text)

    def test_unmodified_control(self):
        self.assertTrue(OBSERVER.observe(self.text,self.expected)['preserved'])

    def test_quadruple_to_aromatic_is_not_preserved(self):
        v=copy.deepcopy(self.value);v['bonds'][0].update(order=1.5,aromatic=True)
        self.assertFalse(OBSERVER.observe(json.dumps(v),self.expected)['preserved'])

    def test_charge_isotope_and_h_changes_are_not_preserved(self):
        for field,value in [('charge',-1),('isotope',1),('hydrogens',1),('aromatic',True),('charge',False)]:
            v=copy.deepcopy(self.value);v['atoms'][0][field]=value
            self.assertFalse(OBSERVER.observe(json.dumps(v),self.expected)['preserved'])

    def test_missing_bond_or_changed_endpoint_is_not_preserved(self):
        for bonds in [[],[dict(source_id='a',target_id='a',order=4,aromatic=False)]]:
            v=copy.deepcopy(self.value);v['bonds']=bonds
            self.assertFalse(OBSERVER.observe(json.dumps(v),self.expected)['preserved'])

    def test_panics_and_timeouts_are_not_safe_cli_refusals(self):
        valid=dict(resource_status='completed',exit_code=2,stdout='',stderr='chematic: MolJSON: invalid graph')
        self.assertEqual(OBSERVER.cli_status(valid),'refused')
        for change in [dict(exit_code=101,stderr='thread panicked'),dict(exit_code=-6),
                       dict(resource_status='timeout'),dict(stdout='partial molecule'),dict(stderr='other error')]:
            self.assertEqual(OBSERVER.cli_status(dict(valid,**change)),'error')

    def test_duplicate_emitted_json_keys_are_not_preserved(self):
        damaged=self.text.replace('"charge": 0','"charge": 1, "charge": 0',1)
        self.assertFalse(OBSERVER.observe(damaged,self.expected)['preserved'])


if __name__ == '__main__':
    unittest.main()
