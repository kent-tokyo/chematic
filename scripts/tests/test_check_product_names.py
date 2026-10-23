from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.checks.check_product_names import scan_paths


class ProductNameCheckTests(unittest.TestCase):
    def test_accepts_canonical_names(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "README.md"
            path.write_text("chematic and chematic-draw\n", encoding="utf-8")
            errors = scan_paths(root, [path])
        self.assertEqual(errors, [])

    def test_rejects_extra_leading_letter_case_insensitively(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "README.md"
            wrong = "s" + "chematic"
            path.write_text(
                f"{wrong}\n{wrong.title()}\n{wrong}-draw\n",
                encoding="utf-8",
            )
            errors = scan_paths(root, [path])
        self.assertEqual(len(errors), 3)


if __name__ == "__main__":
    unittest.main()
