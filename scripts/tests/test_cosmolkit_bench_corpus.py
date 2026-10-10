import unittest

from validation.cosmolkit_comparison.bench_corpus import (
    paired_median_ratio_ci,
    pipeline_seconds,
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


if __name__ == "__main__":
    unittest.main()
