"""Check the Windows pytest runner keeps every file and fails closed."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from scripts.checks.run_pytest_batches import main


class PytestBatchesTests(unittest.TestCase):
    def test_runs_all_files_in_deterministic_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tests = root / "tests"
            tests.mkdir()
            for name in ("test_c.py", "test_a.py", "test_b.py"):
                (tests / name).write_text("def test_ok(): pass\n")
            argv = [
                "run_pytest_batches.py", "--tests-dir", str(tests),
                "--report-dir", str(root / "reports"), "--batch-size", "2",
            ]
            with patch.object(sys, "argv", argv), patch("subprocess.run") as run:
                run.return_value = subprocess.CompletedProcess([], 0)
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(main(), 0)
            self.assertEqual(run.call_count, 2)
            first = run.call_args_list[0].args[0]
            second = run.call_args_list[1].args[0]
            self.assertIn(str(tests / "test_a.py"), first)
            self.assertIn(str(tests / "test_b.py"), first)
            self.assertIn(str(tests / "test_c.py"), second)
            self.assertNotIn(str(tests / "test_c.py"), first)

    def test_stops_at_first_failed_batch(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tests = root / "tests"
            tests.mkdir()
            for index in range(3):
                (tests / f"test_{index}.py").write_text("def test_ok(): pass\n")
            argv = [
                "run_pytest_batches.py", "--tests-dir", str(tests),
                "--report-dir", str(root / "reports"), "--batch-size", "1",
            ]
            with patch.object(sys, "argv", argv), patch("subprocess.run") as run:
                run.return_value = subprocess.CompletedProcess([], 124)
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(main(), 124)
            run.assert_called_once()


if __name__ == "__main__":
    unittest.main()
