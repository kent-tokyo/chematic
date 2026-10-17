import copy
import gzip
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "morgan_numeric", ROOT / "scripts/check_morgan_numeric_packet.py"
)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class NumericCorruptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads(
            gzip.decompress(
                (
                    ROOT / "benchmarks/2026-10-17-morgan-inline-adjacency.json.gz"
                ).read_bytes()
            )
        )

    def test_changed_elapsed_sum_is_rejected(self):
        packet = copy.deepcopy(self.packet)
        row = packet["profiles"][0]["blocks"][0]["arms"]["baseline"]["morgan2_2048"][
            "first_use"
        ]
        row["batch_ns"][0] += 1
        with self.assertRaises(AssertionError):
            audit.check(packet)

    def test_missing_arm_is_rejected(self):
        packet = copy.deepcopy(self.packet)
        del packet["profiles"][0]["blocks"][0]["arms"]["candidate"]
        with self.assertRaises(AssertionError):
            audit.check(packet)

    def test_changed_output_digest_is_rejected(self):
        packet = copy.deepcopy(self.packet)
        row = packet["profiles"][0]["blocks"][0]["arms"]["baseline"]["morgan2_chiral"][
            "first_use"
        ]
        row["output_sha256_per_batch"] = ["0" * 64] * row["batches"]
        with self.assertRaises(AssertionError):
            audit.check(packet)


if __name__ == "__main__":
    unittest.main()
