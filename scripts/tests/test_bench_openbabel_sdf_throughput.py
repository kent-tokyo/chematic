import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "sdf_throughput", Path(__file__).parents[1] / "bench_openbabel_sdf_throughput.py")
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


class SdfThroughputGateTests(unittest.TestCase):
    def row(self):
        return {"graph": "CN=NC", "atoms": [6, 7, 7, 6], "bonds": [(1, 2, 2)],
                "coordinates": [[0.0, 0.0, 0.0]], "name": "molecule",
                "properties": {"source_index": "0"}, "explicit_unknown_stereo": [(1, 2)]}

    def test_wrong_stereo_missing_properties_and_nonfinite_coordinates_fail(self):
        reference = [self.row()]
        for field, value in (("graph", "C/N=N/C"), ("properties", {}),
                             ("coordinates", [[float("nan"), 0.0, 0.0]]),
                             ("coordinates", [[0.0, 0.0]]),
                             ("explicit_unknown_stereo", [])):
            output = copy.deepcopy(reference)
            output[0][field] = value
            self.assertFalse(bench.semantic_check(reference, output)["passed"])
        self.assertFalse(bench.semantic_check(reference, [])["passed"])

    def test_extra_either_marker_is_normalization_but_source_markers_must_survive(self):
        reference = [self.row()]
        output = copy.deepcopy(reference)
        output[0]["explicit_unknown_stereo"].append((2, 3))
        self.assertTrue(bench.semantic_check(reference, output)["passed"])

    def test_baseline_writer_defect_cannot_count_as_an_equivalent_write_win(self):
        semantics = {name: {"passed": True} for name in
                     ("candidate", "openbabel", "candidate_parse", "baseline_parse")}
        semantics["baseline"] = {"passed": False}
        self.assertTrue(bench.eligible(semantics, True, "parse", "baseline"))
        self.assertFalse(bench.eligible(semantics, True, "write", "baseline"))
        self.assertFalse(bench.eligible(semantics, True, "roundtrip", "baseline"))
        self.assertTrue(bench.eligible(semantics, True, "write", "openbabel"))

    def test_wrong_workload_and_nonpositive_native_time_are_rejected(self):
        good = {"elapsed_ns": 123, "records": 10, "repeats": 2, "mode": "parse"}
        for changes in ({"records": 9}, {"elapsed_ns": 0}, {"elapsed_ns": True}, {"elapsed_ns": float("nan")}):
            output = SimpleNamespace(stdout=json.dumps(good | changes).encode(),
                                     stderr=b"100 maximum resident set size\n")
            with patch.object(bench.subprocess, "run", return_value=output):
                with self.assertRaises(ValueError):
                    bench.timing_row(Path("binary"), Path("input"), "parse", 2, 5)

    def test_roundtrip_agreement_does_not_hide_a_parser_stereo_difference(self):
        semantics = {name: {"passed": True} for name in
                     ("candidate", "openbabel", "baseline", "candidate_parse", "baseline_parse", "openbabel_parse")}
        semantics["candidate_parse"]["passed"] = False
        self.assertFalse(bench.eligible(semantics, True, "parse", "openbabel"))
        self.assertFalse(bench.eligible(semantics, True, "parse", "baseline"))
        self.assertTrue(bench.eligible(semantics, True, "roundtrip", "openbabel"))


if __name__ == "__main__":
    unittest.main()
