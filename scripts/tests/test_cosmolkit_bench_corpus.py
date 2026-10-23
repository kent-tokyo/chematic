import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from bench_corpus import (
    _cosmolkit_fingerprint_method,
    operation_succeeded,
    paired_median_ratio_ci,
    pipeline_seconds,
    python_interpreters,
)


class CosmolkitBenchmarkReportTests(unittest.TestCase):
    def test_final_and_rc_fingerprint_names_resolve_before_measurement(self):
        class Final:
            def fingerprint_pattern(self):
                return "final"

            def pattern_fingerprint(self):
                raise AssertionError("final API must take precedence")

        class ReleaseCandidate:
            def pattern_fingerprint(self):
                return "rc"

        for cls, expected in ((Final, "final"), (ReleaseCandidate, "rc")):
            method = _cosmolkit_fingerprint_method(cls, "fingerprint_pattern", "pattern_fingerprint")
            self.assertEqual(method(cls()), expected)
        with self.assertRaises(AttributeError):
            _cosmolkit_fingerprint_method(object, "fingerprint_pattern", "pattern_fingerprint")

    def test_pipeline_time_includes_parse_or_3d_preparation(self):
        self.assertEqual(pipeline_seconds({"parse": 2.0}, "parse"), 2.0)
        self.assertEqual(
            pipeline_seconds({"parse": 2.0, "tpsa": 3.0}, "tpsa"), 5.0
        )
        self.assertEqual(
            pipeline_seconds(
                {
                    "parse": 2.0,
                    "mmff_energy_gradient__prepare": 4.0,
                    "mmff_energy_gradient": 5.0,
                },
                "mmff_energy_gradient",
            ),
            9.0,
        )

    def test_paired_bootstrap_report_is_deterministic(self):
        report = paired_median_ratio_ci([2.0] * 21, [1.0] * 21, seed=7)
        self.assertEqual(
            report,
            {"median": 2.0, "ci95": [2.0, 2.0], "winning_blocks": 21, "blocks": 21},
        )

    def test_all_python_override_does_not_create_a_fourth_engine(self):
        self.assertEqual(
            python_interpreters(["all=/tmp/python"], "/usr/bin/python"),
            {
                "rdkit": "/tmp/python",
                "chematic": "/tmp/python",
                "cosmolkit": "/tmp/python",
            },
        )

    def test_operation_with_any_row_error_is_not_a_valid_timing(self):
        samples = {
            "cosmolkit": [
                {"morgan2_chiral": 0.001, "morgan2_chiral__errors": 0},
                {"morgan2_chiral": 0.0001, "morgan2_chiral__errors": 1},
            ]
        }
        self.assertFalse(operation_succeeded(samples, "cosmolkit", "morgan2_chiral"))


if __name__ == "__main__":
    unittest.main()
