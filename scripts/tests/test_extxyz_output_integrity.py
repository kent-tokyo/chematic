"""Exercise scoring on damaged output; ASE is independently run by the packet."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
P=Path(__file__).parents[2]
spec=importlib.util.spec_from_file_location('observer',P/'scripts/validation/check_extxyz_output_integrity.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.case={'elements':['C','O'],'coords':[[-0.,1.25,2.5],[3.,4.,5.]]}
        self.text='2\nProperties=species:S:1:pos:R:3\nC -0.0 1.25 2.5\nO 3 4 5\n'
        self.expected={'info':{'path':'alpha\\beta'},'arrays':{'weights':[.1,.2]}}
    def score(self,text,snapshot):
        with patch.object(m,'ase_snapshot',return_value=snapshot):return m.observe(self.case,text,self.expected)
    def test_unmodified_control(self):
        r=self.score(self.text,self.expected);self.assertTrue(r['geometry']);self.assertTrue(r['ase_full'])
    def test_signed_zero_is_not_normalized(self):
        self.assertFalse(self.score(self.text.replace('-0.0','0.0'),self.expected)['geometry'])
    def test_rounded_positions_or_swapped_elements_are_not_accepted(self):
        for text in [self.text.replace('1.25','1.25001'),self.text.replace('C -0.0','O -0.0')]:
            self.assertFalse(self.score(text,self.expected)['geometry'])
    def test_lost_metadata_or_atom_data_does_not_get_full_credit(self):
        for value in [{'info':{'path':'alphabeta'},'arrays':{'weights':[.1,.2]}},{'info':{'path':'alpha\\beta'},'arrays':{}}]:
            r=self.score(self.text,value);self.assertTrue(r['geometry']);self.assertFalse(r['ase_full'])
    def test_geometry_requires_exact_extent(self):
        with self.assertRaises(ValueError):m.raw_geometry(self.text+'extra\n')
    def test_unreadable_output_is_not_preserved(self):
        with patch.object(m,'ase_snapshot',side_effect=ValueError('missing columns')):
            r=m.observe(self.case,'not a frame',self.expected)
        self.assertFalse(r['geometry']);self.assertFalse(r['ase_readable']);self.assertFalse(r['ase_full'])
if __name__=='__main__':unittest.main()
