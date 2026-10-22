import unittest
from scripts.validation.check_v3000_source_roundtrip import atom_coordinates,observe

class V3000CoordinateObserverTests(unittest.TestCase):
    def setUp(self):
        self.text='M  V30 BEGIN ATOM\nM  V30 7 C 0.123456789012345 -0.0 0.000012345 0\nM  V30 END ATOM\n'
    def test_geometry_is_separate_from_header_and_source_ids(self):
        out='new header\n'+self.text.replace('7 C','1 C')
        self.assertEqual(observe(self.text,out),dict(atom_table_observed=True,exact_source=False,coordinates_exact=True,source_ids=False))
    def test_z_rounding_sign_and_precision_loss_are_detected(self):
        for before,after in [('0.000012345','0.0000'),('-0.0','0.0'),('0.123456789012345','0.1235')]:
            self.assertFalse(observe(self.text,self.text.replace(before,after))['coordinates_exact'])
    def test_continuation_and_atom_count_are_checked(self):
        continued=self.text.replace('-0.0 0.000012345','-0.0 -\nM  V30 0.000012345')
        self.assertEqual(atom_coordinates(self.text),atom_coordinates(continued))
        self.assertFalse(observe(self.text,self.text.replace('7 C','7 O'))['coordinates_exact'])
        with self.assertRaises(ValueError):atom_coordinates(self.text.replace('0.000012345','NaN'))

if __name__=='__main__':unittest.main()
