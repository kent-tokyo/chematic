"""Runner contract tests; no benchmarks or external libraries required."""
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "bench_hotpath_pair", Path(__file__).parents[1] / "bench_hotpath_pair.py"
)
BENCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCH)


class PairedBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.args = ["bench", "--pairs", "2", "--repeats", "1"]
        for name in ("baseline", "candidate", "smiles", "sdf"):
            path = root / name
            path.write_text(name)
            self.args.extend(["--" + name, str(path)])
        self.output = root / "report.json"
        self.args.extend(["--output", str(self.output)])
        self.order = []

    def fake_run(self, command, **kwargs):
        name = Path(command[0]).name
        if command[-1] == "--dump":
            return SimpleNamespace(stdout=b"same exact output")
        self.order.append(name)
        result = dict(output_fnv1a="same", smiles=10, records=3, repeats=1)
        result.update({lane: 2 if name == "baseline" else 1 for lane in
                       ("parse_us", "canonical_us", "sdf_read_us", "sdf_write_us")})
        return SimpleNamespace(stdout=json.dumps(result).encode())

    def run_main(self, runner=None):
        with patch.object(sys, "argv", self.args), patch.object(
            BENCH.subprocess, "run", side_effect=runner or self.fake_run
        ), patch.object(BENCH.platform, "platform", return_value="test-host"), patch("builtins.print"):
            BENCH.main()

    def test_alternating_order_summary_and_completed_resume(self):
        self.run_main()
        self.assertEqual(self.order, ["baseline", "candidate", "candidate", "baseline"])
        report = json.loads(self.output.read_text())
        self.assertTrue(report["complete"])
        self.assertEqual(report["summary"]["parse_us"]["paired_speedup_median"], 2)
        self.run_main()
        self.assertEqual(len(self.order), 4)
        Path(self.args[self.args.index("--candidate") + 1]).write_text("changed binary")
        with self.assertRaisesRegex(SystemExit, "refusing to resume"):
            self.run_main()

    def test_exact_output_mismatch_fails_before_timing(self):
        def different(command, **kwargs):
            return SimpleNamespace(stdout=Path(command[0]).name.encode())
        with self.assertRaisesRegex(SystemExit, "output bytes differ"):
            self.run_main(different)
        self.assertFalse(self.output.exists())

    def test_interruption_preserves_completed_pairs(self):
        def interrupt(command, **kwargs):
            if len(self.order) == 2 and command[-1] != "--dump":
                raise RuntimeError("interrupted")
            return self.fake_run(command, **kwargs)
        with self.assertRaisesRegex(RuntimeError, "interrupted"):
            self.run_main(interrupt)
        self.assertEqual(len(json.loads(self.output.read_text())["pairs"]), 1)
        self.run_main()
        self.assertEqual(self.order, ["baseline", "candidate", "candidate", "baseline"])

    def test_non_finite_timing_rejected(self):
        def invalid(command, **kwargs):
            output = self.fake_run(command, **kwargs)
            if command[-1] != "--dump":
                data = json.loads(output.stdout)
                data["parse_us"] = float("nan")
                output.stdout = json.dumps(data).encode()
            return output
        with self.assertRaisesRegex(SystemExit, "invalid timing"):
            self.run_main(invalid)
        self.assertFalse(json.loads(self.output.read_text())["complete"])


if __name__ == "__main__":
    unittest.main()
