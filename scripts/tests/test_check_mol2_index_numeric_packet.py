import copy
import contextlib
import io
import os
import tempfile
from unittest.mock import patch
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "mol2_numeric", ROOT / "scripts/check_mol2_index_numeric_packet.py"
)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class NumericPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packets = [
            audit.read_packet(ROOT / "benchmarks" / name)
            for name in (
                "2026-10-18-mol2-index-only-numeric.json.gz",
                "2026-10-18-mol2-index-fixed4-numeric.json.gz",
            )
        ]

    def test_both_epoch_decisions_recompute(self):
        for packet, expected in zip(self.packets, (False, True)):
            self.assertEqual(audit.check(packet)["source_gate_passed"], expected)

    def test_reproduction_prepares_all_frozen_input_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.dict(os.environ, {"CHEMATIC_MOL2_EVIDENCE": directory}):
                spec = importlib.util.spec_from_file_location(
                    "mol2_reproduction", ROOT / "scripts/benchmarks/bench_openbabel_mol2_index.py"
                )
                runner = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(runner)
                with contextlib.redirect_stdout(io.StringIO()):
                    runner.prepare_inputs(
                        ROOT / "benchmarks/2026-10-18-mol2-index-fixed4-numeric.json.gz"
                    )
                self.assertEqual(len(list(Path(directory).glob("*.mol2s"))), 9)

    def test_corrupt_accounting_order_gate_and_intervals_are_rejected(self):
        for original in self.packets:
            for mutation in (
                "records",
                "budget",
                "order",
                "gate",
                "interval",
                "metadata",
            ):
                with self.subTest(epoch=original["label"], mutation=mutation):
                    packet = copy.deepcopy(original)
                    block = packet["profiles"]["ethanol"]["blocks"][0]
                    if mutation == "records":
                        block["candidate"]["parse"]["records"] += 1
                    elif mutation == "budget":
                        block["candidate"]["parse"]["minimum_timed_ns"] = 1
                    elif mutation == "order":
                        block["order"] = list(reversed(block["order"]))
                    elif mutation == "gate":
                        packet["source_gate"]["all_lane_median_floor"] = 0.97
                    elif mutation == "interval":
                        packet["profiles"]["ethanol"]["summary"]["parse"]["source"][
                            "ci95"
                        ][0] += 0.001
                    else:
                        packet["preflight"]["profiles"]["ethanol"][
                            "external_eligible"
                        ] = True
                    with self.assertRaises(ValueError):
                        audit.check(packet)


if __name__ == "__main__":
    unittest.main()
