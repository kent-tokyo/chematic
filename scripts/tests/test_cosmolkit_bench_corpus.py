import unittest

from validation.cosmolkit_comparison.bench_corpus import (
    operation_succeeded,
    paired_median_ratio_ci,
    pipeline_seconds,
    python_interpreters,
)


class CosmolkitBenchmarkReportTests(unittest.TestCase):
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
