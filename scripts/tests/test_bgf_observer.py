"""Observer controls: distinguish source encoding from bond meaning."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("bgf_observer", ROOT / "scripts/validation/check_bgf_preservation.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
PACKET = json.loads((ROOT / "validation/bgf/cases.json").read_text(encoding="utf-8"))


class BgfObserverTests(unittest.TestCase):
    def test_original_cases_independently_observed(self):
        for c in PACKET["valid"]:
            self.assertEqual(M.observe(c["input"]), c["expected"], c["id"])

    def test_zero_one_and_reciprocal_normalization(self):
        r = copy.deepcopy(PACKET["valid"][0]["expected"])
        for c in r["connectivity"]:
            c["orders"] = [0]
        before = M.axes(r)
        r["connectivity"] = r["connectivity"][:1]
        r["connectivity"][0]["orders"] = None
        after = M.axes(r)
        self.assertEqual(before["bond_graph"], after["bond_graph"])
        self.assertNotEqual(before["order_encoding"], after["order_encoding"])

    def test_charge_and_bond_changes_are_detected(self):
        r = copy.deepcopy(PACKET["valid"][0]["expected"])
        before = M.axes(r)
        r["atoms"][0]["partial_charge"] = 0.0
        for c in r["connectivity"]:
            c["orders"] = [1]
        after = M.axes(r)
        self.assertNotEqual(before["partial_charge"], after["partial_charge"])
        self.assertNotEqual(before["bond_graph"], after["bond_graph"])
        self.assertEqual(before["position_angstrom"], after["position_angstrom"])

    def test_invalid_outputs_not_counted(self):
        for text in ("", "BIOGRF 200\nEND\n", PACKET["valid"][0]["input"].replace("ORDER      1     2", "ORDER      1     4")):
            with self.assertRaises((ValueError, KeyError, IndexError)):
                M.observe(text)


if __name__ == "__main__":
    unittest.main()
