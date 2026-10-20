"""Failure reports must retain provenance and remain visibly failed."""

import contextlib
import hashlib
import importlib.util
import io
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "same_process_report", SCRIPTS / "lib" / "same_process_report.py"
)
with patch.object(sys, "path", [str(SCRIPTS), *sys.path]):
    REPORT = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(REPORT)


class SameProcessReportTests(unittest.TestCase):
    def test_failed_report_keeps_fixture_identity_and_exit_code(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "sample.sdf"
            payload = b"fixture bytes\n"
            fixture.write_bytes(payload)
            rdkit = types.ModuleType("rdkit")
            rdkit.rdBase = types.SimpleNamespace(rdkitVersion="test-rdkit")
            with patch.dict(os.environ, {"SCHEMATIC_BENCHMARK_VERSION": "1.2.3"}), patch.dict(
                sys.modules, {"rdkit": rdkit}
            ):
                report = REPORT.make_report(
                    root, fixture, payload, 2, "sample_contract",
                    rows={"chematic": {"records": 0}},
                    comparison={"expected_records": 2, "speed_is_non_ranking_context": True},
                    errors=["missing records"],
                    extra={"input_blocks_per_repeat": 1},
                )
            output = Path("result.json")
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = REPORT.write_report(root, output, report)

            self.assertEqual(code, 1)
            self.assertEqual(stdout.getvalue(), (root / output).read_text())
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["target_version"], "1.2.3")
            self.assertEqual(report["fixture"], {
                "path": "sample.sdf", "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            })
            self.assertEqual(report["tool_versions"], {"rdkit": "test-rdkit"})
            self.assertEqual(report["input_blocks_per_repeat"], 1)


if __name__ == "__main__":
    unittest.main()
