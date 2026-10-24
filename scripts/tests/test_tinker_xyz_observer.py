"""Negative controls for the comparison's field accounting."""
import copy
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("tinker_observer", Path(__file__).resolve().parents[1] / "validation/check_tinker_xyz_preservation.py")
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.record = OBSERVER.observe("3 title\n1 O 0 0 0 41 3 2\n2 H .9 0 0 42 1\n3 H -.2 .9 0 42 1\n")

    def test_neighbor_order_is_not_a_chemical_difference(self):
        changed = copy.deepcopy(self.record)
        changed["atoms"][0]["neighbors"].reverse()
        self.assertEqual(OBSERVER.axes(self.record), OBSERVER.axes(changed))

    def test_axes_detect_each_declared_field_change(self):
        mutations = {"indices":("index", 7), "names":("name", "OW"),
                     "atom_types":("atom_type",90001), "coordinates":("position_angstrom",[1e-12,0,0]),
                     "connectivity":("neighbors",[2])}
        expected = OBSERVER.axes(self.record)
        for axis, (field,value) in mutations.items():
            changed = copy.deepcopy(self.record);changed["atoms"][0][field] = value
            differences = [a for a in OBSERVER.AXES if OBSERVER.axes(changed)[a] != expected[a]]
            self.assertEqual(differences, [axis])
        changed = copy.deepcopy(self.record);changed["cell"] = {"lengths_angstrom":[10,10,10],"angles_degrees":[90,90,90]}
        self.assertNotEqual(OBSERVER.axes(changed)["cell"], expected["cell"])
        changed["title"] = "new title"
        self.assertNotEqual(OBSERVER.axes(changed)["title"], expected["title"])

    def test_missing_and_extra_records_are_not_observed(self):
        for text in ["", "1\n", "2\n1 O 0 0 0 41\n", "1\n1 O 0 0 0 41\n1\n1 O 0 0 0 41\n", "1\n1 O nan 0 0 41\n"]:
            with self.assertRaises((ValueError, IndexError)):
                OBSERVER.observe(text)


if __name__ == "__main__":
    unittest.main()
