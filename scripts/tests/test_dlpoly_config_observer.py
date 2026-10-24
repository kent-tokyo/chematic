"""Negative controls for raw CONFIG field comparisons."""
import copy
import importlib.util
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location("config_observer",Path(__file__).resolve().parents[1]/"validation/check_dlpoly_config_preservation.py")
OBSERVER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(OBSERVER)


class ObserverTests(unittest.TestCase):
    def setUp(self):
        self.record=OBSERVER.observe("title\n2 0 1\nOW 1 opaque\n1 2 3\n4 5 6\n7 8 9\n")

    def test_every_vector_kind_and_metadata_axis_is_independent(self):
        mapping={"names":("name","HW"),"indices":("index",2),"positions":("position",[.1,2,3]),
                 "velocities":("velocity",None),"forces":("force",None),"header_extra":("header_extra",[])}
        baseline=OBSERVER.axes(self.record)
        for axis,(field,value) in mapping.items():
            changed=copy.deepcopy(self.record);changed["atoms"][0][field]=value
            self.assertEqual([k for k in OBSERVER.AXES if OBSERVER.axes(changed)[k]!=baseline[k]],[axis])

    def test_missing_vectors_are_not_observed_as_zero_vectors(self):
        zero=OBSERVER.observe("title\n0 0\nOW 1\n1 2 3\n")
        self.assertIsNone(OBSERVER.axes(zero)["velocities"][0])
        self.assertNotEqual(OBSERVER.axes(zero)["velocities"],[[0,0,0]])
        self.assertNotEqual(OBSERVER.axes(zero)["forces"],OBSERVER.axes(self.record)["forces"])

    def test_truncation_and_nonfinite_output_are_unobservable(self):
        for text in ["", "title\n2 0 1\nOW 1\n1 2 3\n", "title\n0 0 2\nX\n0 0 0\n", "title\n0 0 1\nX\nnan 0 0\n"]:
            with self.assertRaises(ValueError):OBSERVER.observe(text)

    def test_cell_and_boundary_changes_remain_separate(self):
        baseline=OBSERVER.axes(self.record)
        changed=copy.deepcopy(self.record);changed["cell_vectors"]=[[10,1,0],[0,10,0],[0,0,10]]
        self.assertNotEqual(OBSERVER.axes(changed)["cell_vectors"],baseline["cell_vectors"])
        self.assertEqual(OBSERVER.axes(changed)["imcon"],baseline["imcon"])


if __name__=="__main__":unittest.main()
