"""Check the outer CI deadline stops a command that ignores pytest limits."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "checks" / "run_bounded_command.py"


class BoundedCommandTests(unittest.TestCase):
    def test_timeout_exits_nonzero_and_reports_deadline(self) -> None:
        result = subprocess.run(
            [
                sys.executable, str(SCRIPT),
                "--timeout-seconds", "1", "--heartbeat-seconds", "1", "--",
                sys.executable, "-c", "import time; time.sleep(5)",
            ],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
        )
        self.assertEqual(result.returncode, 124)
        self.assertIn("CI command exceeded 1s", result.stderr)


if __name__ == "__main__":
    unittest.main()
