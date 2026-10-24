"""Preservation observer controls distinguish semantic changes from formatting."""
import importlib.util
from pathlib import Path
import unittest
SP=importlib.util.spec_from_file_location('unixyz_observer',Path(__file__).resolve().parents[1]/'validation/check_unixyz_preservation.py')
M=importlib.util.module_from_spec(SP);SP.loader.exec_module(M)
TEXT='title\n1\n6 .125 -0.0 .75\n'
class Controls(unittest.TestCase):
 def test_equal(self):self.assertTrue(M.comparisons(M.observe(TEXT),M.observe(TEXT))['declared_state_exact'])
 def test_rounded_position(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('.125','.12')),M.observe(TEXT))['coordinate_bits_exact'])
 def test_signed_zero(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('-0.0','0.0')),M.observe(TEXT))['declared_state_exact'])
 def test_changed_title(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('title','other')),M.observe(TEXT))['title_exact'])
 def test_wrong_atomic_number(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('6 .125','8 .125')),M.observe(TEXT))['atomic_numbers_order_exact'])
 def test_spelling_is_not_numeric_loss(self):self.assertTrue(M.comparisons(M.observe(TEXT.replace('.125','0.12500')),M.observe(TEXT))['declared_state_exact'])
 def test_empty_truncated_extra_output_unobservable(self):
  for text in ['', 'title\n2\n6 0 0 0\n', TEXT+'8 0 0 0\n', TEXT.replace('.125 -0.0 .75','.125-0.0.75')]:
   with self.assertRaises(ValueError):M.observe(text)
 def test_crashes_partial_output_wrong_prefix(self):
  r={'returncode':2,'stdout':'','stderr':'chematic: UniChem XYZ line 1: invalid'}
  self.assertTrue(M.cli_refusal(r))
  for k,v in [('returncode',101),('returncode',-11),('stdout','partial'),('stderr','panic')]:self.assertFalse(M.cli_refusal(dict(r,**{k:v})))
if __name__=='__main__':unittest.main()
