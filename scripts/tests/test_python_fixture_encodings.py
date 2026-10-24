from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from tools.check_python_fixture_encodings import unencoded_read_text_calls


class PythonFixtureEncodingTests(unittest.TestCase):
    def test_accepts_explicit_utf8(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            tests = root / "crates" / "chematic-py" / "tests"
            tests.mkdir(parents=True)
            (tests / "test_fixture.py").write_text(
                'Path("fixture.json").read_text(encoding="utf-8")\n',
                encoding="utf-8",
            )

            self.assertEqual(unencoded_read_text_calls(root), [])

    def test_reports_host_locale_dependent_read(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            tests = root / "crates" / "chematic-py" / "tests"
            tests.mkdir(parents=True)
            (tests / "test_fixture.py").write_text(
                'Path("fixture.json").read_text()\n',
                encoding="utf-8",
            )

            self.assertEqual(
                unencoded_read_text_calls(root),
                ["crates/chematic-py/tests/test_fixture.py:1"],
            )


if __name__ == "__main__":
    unittest.main()
