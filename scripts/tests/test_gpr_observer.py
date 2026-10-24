"""Evidence controls prevent data loss and crashes counting as preservation."""
import copy
import importlib.util
from pathlib import Path
import unittest
SPEC=importlib.util.spec_from_file_location('gpr_observer',Path(__file__).resolve().parents[1]/'validation/check_gpr_preservation.py')
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)
TEXT='!Header gpr 100\n!Info 1\n!Atoms 2\n0 6\n1 8\n!Bonds 1\n0 1 D\n!Coord\n0 .125 -0.0 0\n1 .25 .5 .75\n!Charges\n0 -.125\n1 .125\n!End\n'
class Controls(unittest.TestCase):
 def test_equal(self):self.assertTrue(M.comparisons(M.observe(TEXT),M.observe(TEXT))['declared_state_exact'])
 def test_rounded_position(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('0 .125','0 .12')),M.observe(TEXT))['coordinate_bits_exact'])
 def test_signed_zero(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('-0.0','0.0')),M.observe(TEXT))['declared_state_exact'])
 def test_rounded_charge(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('-.125','-.12')),M.observe(TEXT))['partial_charge_presence_and_bits_exact'])
 def test_invented_charges(self):
  r=M.observe(TEXT);missing=copy.deepcopy(r);missing['partial_charges']=None
  self.assertFalse(M.comparisons(r,missing)['declared_state_exact'])
 def test_changed_bond_is_loss_but_orientation_is_not(self):
  self.assertFalse(M.comparisons(M.observe(TEXT.replace('0 1 D','0 1 S')),M.observe(TEXT))['graph_exact'])
  self.assertTrue(M.comparisons(M.observe(TEXT.replace('0 1 D','1 0 D')),M.observe(TEXT))['declared_state_exact'])
 def test_empty_output_unobservable(self):
  with self.assertRaises(ValueError):M.observe('')
 def test_crashes_partial_output_wrong_prefix(self):
  r={'returncode':2,'stdout':'','stderr':'chematic: GPR line 1: invalid'}
  self.assertTrue(M.cli_refusal(r))
  for k,v in [('returncode',101),('returncode',-11),('stdout','partial'),('stderr','panic')]:self.assertFalse(M.cli_refusal(dict(r,**{k:v})))
if __name__=='__main__':unittest.main()
