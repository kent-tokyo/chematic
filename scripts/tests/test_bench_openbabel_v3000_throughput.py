import copy
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[1]))
SPEC = importlib.util.spec_from_file_location(
    "v3000_throughput", Path(__file__).parents[1] / "bench_openbabel_v3000_throughput.py")
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


class V3000ThroughputGateTests(unittest.TestCase):
    def row(self):
        return {"graph": "CC", "atoms": [6, 6], "bonds": [(0, 1, 1)],
                "coordinates": [[0., 0., 0.], [1., 0., 0.]], "name": "molecule",
                "properties": {}, "explicit_unknown_stereo": []}

    def test_rejected_rows_keep_their_positions_and_fail_the_gate(self):
        reference = [self.row(), self.row(), self.row()]
        with patch.object(bench, "observe_rows", return_value=[self.row(), None, self.row()]):
            result = bench.rewrite_check(reference, b"output")
        self.assertFalse(result["passed"])
        self.assertEqual(result["observed_rows"], 3)
        self.assertEqual(result["accepted_rows"], 2)
        self.assertEqual(result["rejected_rows"], [1])

    def test_wrong_graph_after_a_rejected_row_keeps_the_original_index(self):
        wrong = copy.deepcopy(self.row())
        wrong["graph"] = "C=C"
        with patch.object(bench, "observe_rows", return_value=[None, wrong]):
            result = bench.rewrite_check([self.row(), self.row()], b"output")
        self.assertEqual(result["mismatch_examples"]["graph"], [1])
        self.assertFalse(result["passed"])

    def test_all_rejections_keep_zero_counters_and_fail(self):
        with patch.object(bench, "observe_rows", return_value=[None]):
            result = bench.rewrite_check([self.row()], b"output")
        self.assertEqual(set(result["mismatch_counts"]), set(bench.FIELDS))
        self.assertEqual(result["accepted_rows"], 0)
        self.assertEqual(result["rejected_rows"], [0])
        self.assertFalse(result["passed"])

    def test_missing_rows_are_not_treated_as_accepted_rejections(self):
        with patch.object(bench, "observe_rows", return_value=[self.row()]):
            result = bench.rewrite_check([self.row(), self.row()], b"output")
        self.assertTrue(result["record_count_mismatch"])
        self.assertFalse(result["passed"])

    def test_failed_rewrite_does_not_disable_a_separately_correct_parser(self):
        semantics = {name: {"passed": True} for name in
                     ("candidate", "candidate_parse", "baseline", "baseline_parse", "openbabel_parse")}
        semantics["openbabel"] = {"passed": False}
        self.assertTrue(bench.eligible(semantics, True, "parse", "openbabel"))
        self.assertFalse(bench.eligible(semantics, True, "write", "openbabel"))
        self.assertFalse(bench.eligible(semantics, True, "roundtrip", "openbabel"))


if __name__ == "__main__":
    unittest.main()
