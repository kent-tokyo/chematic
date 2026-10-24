"""Controls for unit normalization, metadata meaning and safe error accounting."""
import importlib.util
from pathlib import Path
import unittest
SP=importlib.util.spec_from_file_location('viewmol_observer',Path(__file__).resolve().parents[1]/'validation/check_viewmol_preservation.py')
M=importlib.util.module_from_spec(SP);SP.loader.exec_module(M)
TEXT='$title\ntitle\n$coord 2\n.125 -0.0 .5 C xz\n$end\n'
class Controls(unittest.TestCase):
 def test_equal(self):
  self.assertTrue(M.comparisons(M.observe(TEXT),M.observe(TEXT))['declared_state_exact'])
  self.assertEqual(M.observe(TEXT.replace('title\n$coord','a\u2028b\n$coord'))['title'],'a\u2028b')
 def test_normalization_is_not_physical_loss(self):
  actual=M.observe(TEXT.replace('$coord 2','$coord 1').replace('.125 -0.0 .5','.25 -0.0 1'))
  r=M.comparisons(actual,M.observe(TEXT));self.assertTrue(r['declared_state_exact']);self.assertFalse(r['raw_scale_coordinates_exact'])
 def test_lost_scale_changes_geometry(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('$coord 2','$coord 1')),M.observe(TEXT))['angstrom_positions_bits_exact'])
 def test_lost_axes_is_loss(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace(' C xz',' C')),M.observe(TEXT))['declared_state_exact'])
 def test_axis_order_is_not_physical_loss(self):self.assertTrue(M.comparisons(M.observe(TEXT.replace(' C xz',' C zx')),M.observe(TEXT))['declared_state_exact'])
 def test_lost_signed_zero(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('-0.0','0.0')),M.observe(TEXT))['declared_state_exact'])
 def test_changed_title(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace('title\n$coord','other\n$coord')),M.observe(TEXT))['declared_state_exact'])
 def test_changed_element(self):self.assertFalse(M.comparisons(M.observe(TEXT.replace(' C xz',' N xz')),M.observe(TEXT))['declared_state_exact'])
 def test_empty_output_unobservable(self):
  with self.assertRaises(ValueError):M.observe('')
 def test_crashes_partial_output_and_wrong_prefix(self):
  r={'returncode':2,'stdout':'','stderr':'chematic: ViewMol line 1: invalid'};self.assertTrue(M.cli_refusal(r))
  for k,v in [('returncode',101),('returncode',-11),('stdout','partial'),('stderr','panic')]:self.assertFalse(M.cli_refusal(dict(r,**{k:v})))
if __name__=='__main__':unittest.main()
