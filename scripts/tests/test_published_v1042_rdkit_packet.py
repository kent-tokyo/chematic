import unittest

from scripts.check_published_v1042_rdkit_packet import validate
from scripts.run_published_v1042_rdkit_python_packet import paired_ratio_ci


def record(equivalent: bool = True) -> dict:
    gate = {
        "compared": 10_000,
        "matches": 10_000 if equivalent else 9_996,
        "mismatch_count": 0 if equivalent else 4,
        "error_count": 0,
        "equivalent": equivalent,
    }
    interval = {
        "median": 2.0,
        "ci95": [2.0, 2.0],
        "winning_blocks": 21,
        "blocks": 21,
    }
    return {
        "schema": "published-v1.0.42-rdkit-packet/v1",
        "binding": "python",
        "versions": {"chematic": "1.0.42", "rdkit": "2026.9.1"},
        "artifacts": [
            {"file": "chematic.whl", "bytes": 1, "sha256": "a" * 64},
            {"file": "rdkit.whl", "bytes": 1, "sha256": "b" * 64},
        ],
        "corpus": {"rows": 10_000, "timing_rows": 1_000, "sha256": "c" * 64},
        "method": {"blocks": 21},
        "accuracy": {op: dict(gate) for op in ("tpsa", "labute_asa", "num_rings", "morgan2_chiral")},
        "speedup_rdkit_over_chematic": {
            op: {"eligible": equivalent, "prepared": interval if equivalent else None,
                 "pipeline": interval if equivalent else None}
            for op in ("tpsa", "labute_asa", "num_rings", "morgan2_chiral")
        },
        "raw_samples_seconds": {"chematic": [{}] * 21, "rdkit": [{}] * 21},
    }


class PublishedPacketTests(unittest.TestCase):
    def test_paired_bootstrap_is_deterministic(self):
        self.assertEqual(
            paired_ratio_ci([2.0] * 21, [1.0] * 21, 42),
            {"median": 2.0, "ci95": [2.0, 2.0], "winning_blocks": 21, "blocks": 21},
        )

    def test_valid_equivalent_record(self):
        self.assertEqual(validate(record()), [])

    def test_non_equivalent_lane_must_withhold_speed(self):
        report = record()
        report["accuracy"]["num_rings"] = record(False)["accuracy"]["num_rings"]
        report["speedup_rdkit_over_chematic"]["num_rings"] = {
            "eligible": False,
            "prepared": None,
            "pipeline": None,
        }
        self.assertEqual(validate(report), [])
        report["speedup_rdkit_over_chematic"]["num_rings"]["prepared"] = {
            "median": 2.0,
            "blocks": 21,
        }
        self.assertIn(
            "non-equivalent num_rings reports a speed interval", validate(report)
        )


if __name__ == "__main__":
    unittest.main()
