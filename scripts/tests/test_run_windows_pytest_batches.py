"""The Windows batch supervisor must retain output and enforce its deadline."""

from pathlib import Path
import sys
import tempfile
import time
import unittest

from scripts.checks.run_windows_pytest_batches import run_with_deadline


class BatchSupervisorTests(unittest.TestCase):
    def test_child_output_is_captured(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "batch.log"
            result = run_with_deadline(
                [sys.executable, "-c", "print('captured')"], log, 5, 1
            )
            self.assertEqual(result, 0)
            self.assertEqual(log.read_text(), "captured\n")

    def test_stalled_child_exits_at_parent_deadline(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "batch.log"
            start = time.monotonic()
            result = run_with_deadline(
                [
                    sys.executable,
                    "-c",
                    "import time; print('started', flush=True); time.sleep(10)",
                ],
                log,
                0.3,
                0.1,
            )
            self.assertEqual(result, 124)
            self.assertLess(time.monotonic() - start, 3)
            self.assertIn("started", log.read_text())


if __name__ == "__main__":
    unittest.main()
