import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('nwchem_observer', ROOT/'scripts/validation/check_nwchem_geometry.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class NwchemObserverTests(unittest.TestCase):
    def test_same_indexed_xyz_and_export(self):
        expected = M.xyz_geometry('2\nname\nC 0.123456789 0 0\nO 1 0 0\n')
        self.assertEqual(M.observe('start job\ntitle text\ngeometry units angstroms print xyz autosym\nC 0.123456789 0 0\nO 1 0 0\nend\n'), expected)

    def test_reordering_and_rounding_differ(self):
        reference = M.observe('geometry\nC 0.123456789 0 0\nO 1 0 0\nend\n')
        for output in ('geometry\nO 1 0 0\nC 0.123456789 0 0\nend\n', 'geometry\nC 0.12346 0 0\nO 1 0 0\nend\n'):
            self.assertNotEqual(M.observe(output), reference)

    def test_wrong_units_and_invalid_output_refuse(self):
        for output in ('', 'geometry\nend\n', 'geometry\nC 0 0 0\n', 'geometry units bohr\nC 1 0 0\nend\n', 'geometry\nC NaN 0 0\nend\n'):
            with self.assertRaises(ValueError):
                M.observe(output)

    def test_multiple_geometry_blocks_refuse(self):
        with self.assertRaises(ValueError):
            M.observe('geometry\nC 0 0 0\nend\n'*2)


if __name__ == '__main__':
    unittest.main()
