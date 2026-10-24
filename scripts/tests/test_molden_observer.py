"""Observer controls distinguish physical data from source encoding changes."""
import copy
import importlib.util
from pathlib import Path
import unittest
SPEC=importlib.util.spec_from_file_location("molden_observer",Path(__file__).resolve().parents[1]/"validation/check_molden_preservation.py")
OBSERVER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(OBSERVER)


class ObserverTests(unittest.TestCase):
    def setUp(self):self.r=OBSERVER.observe("[Molden Format]\n[FREQ]\n-100\n[FR-COORD]\nH 0 0 0\n[FR-NORM-COORD]\nvibration 1\n.1 .2 .3\n[INT]\n2 0\n")

    def test_negative_frequency_and_zero_raman_are_not_absent(self):
        self.assertEqual(OBSERVER.axes(self.r)["frequencies"],[-100])
        self.assertEqual(OBSERVER.axes(self.r)["raman_values"],[[0,0]])
        changed=copy.deepcopy(self.r);changed["intensities"][0]["raman"]=None
        self.assertIsNone(OBSERVER.axes(changed)["raman_values"])

    def test_derived_atoms_do_not_mutate_original_fr_coordinates(self):
        changed=copy.deepcopy(self.r);changed["atoms"]={"unit":"angstrom","atoms":[{"symbol":"H","index":1,"atomic_number":1,"position":[0,0,0]}]}
        self.assertEqual(OBSERVER.axes(changed)["fr_coords_bohr"],OBSERVER.axes(self.r)["fr_coords_bohr"])
        self.assertIsNone(OBSERVER.axes(self.r)["atom_rows"])

    def test_each_vibration_axis_detects_changes(self):
        before=OBSERVER.axes(self.r)
        for field,axis,value in [("frequencies","frequencies",[100]),("normal_modes","normal_modes_bohr",None),("intensities","ir_values",None),("fr_coords","fr_coords_bohr",None)]:
            changed=copy.deepcopy(self.r);changed[field]=value
            self.assertNotEqual(OBSERVER.axes(changed)[axis],before[axis])

    def test_truncated_nonfinite_and_duplicate_sections_are_not_observable(self):
        for s in ["", "[Molden Format]\n[FR-COORD]\nH nan 0 0\n", "[Molden Format]\n[FR-COORD]\nH 0 0 0\n[FR-NORM-COORD]\nvibration 1\n", "[Molden Format]\n[FREQ]\n1\n[FREQ]\n2\n"]:
            with self.assertRaises(ValueError):OBSERVER.observe(s)


if __name__=="__main__":unittest.main()
