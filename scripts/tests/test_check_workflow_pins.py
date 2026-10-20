"""Local workflow calls are fixed to the caller commit; remote actions need SHAs."""

import tempfile
import unittest
from pathlib import Path

from scripts.check_workflow_pins import reference_error


class WorkflowPinTests(unittest.TestCase):
    def test_local_workflow_and_remote_sha(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workflow = root / ".github" / "workflows" / "ci.yml"
            workflow.parent.mkdir(parents=True)
            workflow.write_text("on: workflow_call\n")
            self.assertIsNone(reference_error("./.github/workflows/ci.yml", root))
            self.assertIsNone(reference_error("actions/checkout@" + "a" * 40, root))
            self.assertIsNotNone(reference_error("./.github/workflows/missing.yml", root))
            self.assertIsNotNone(reference_error("./../outside.yml", root))
            self.assertIsNotNone(reference_error("actions/checkout@main", root))


if __name__ == "__main__":
    unittest.main()
