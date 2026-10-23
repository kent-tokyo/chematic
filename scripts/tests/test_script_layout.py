from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from tools.check_script_layout import existing_files, validate_layout


class ScriptLayoutTests(unittest.TestCase):
    def test_ignores_tracked_paths_deleted_from_working_tree(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            present = root / "scripts" / "present.py"
            present.parent.mkdir()
            present.write_text("", encoding="utf-8")

            tracked = ["scripts/present.py", "scripts/deleted.py"]

            self.assertEqual(existing_files(root, tracked), ["scripts/present.py"])

    def test_accepts_matching_root_baseline(self) -> None:
        tracked = ["scripts/kept.py", "scripts/checks/current.py"]

        self.assertEqual(validate_layout(tracked, ["scripts/kept.py"]), [])


if __name__ == "__main__":
    unittest.main()
