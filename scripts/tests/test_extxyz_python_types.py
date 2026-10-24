"""The packet observer must catch type, precision, shape and comment losses."""
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "typed_xyz_observer", ROOT / "scripts/validation/check_extxyz_python_types.py"
)
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)
CASES = json.loads((ROOT / "validation/extxyz_python_types/cases.json").read_text(encoding="utf-8"))["valid"]


class ObserverTests(unittest.TestCase):
    def test_all_authored_controls_are_observable(self):
        for case in CASES:
            expected = OBSERVER.snapshot(case["text"])
            self.assertTrue(OBSERVER.observe(case["text"], expected)["preserved"])

    def test_rounded_wide_integer_is_not_preserved(self):
        text = CASES[1]["text"]
        expected = OBSERVER.snapshot(text)
        self.assertEqual(expected["rows"][0][2], [9007199254740993])
        damaged = text.replace("9007199254740993", "9007199254740992")
        self.assertFalse(OBSERVER.observe(damaged, expected)["preserved"])

    def test_integer_to_real_type_change_is_not_preserved(self):
        text = CASES[1]["text"]
        self.assertFalse(OBSERVER.observe(text.replace(":id:I:", ":id:R:"),
                                         OBSERVER.snapshot(text))["preserved"])

    def test_signed_zero_and_tiny_real_loss_are_not_preserved(self):
        text = CASES[9]["text"]
        for damaged in (text.replace("-0.0", "0.0"), text.replace("5e-324", "0.0")):
            self.assertFalse(OBSERVER.observe(damaged, OBSERVER.snapshot(text))["preserved"])

    def test_empty_frame_width_and_plain_comment_are_measured(self):
        for case, old, new in ((CASES[5], ":vec:R:4", ":vec:R:1"),
                               (CASES[6], "sample 甲", "")):
            text = case["text"]
            self.assertFalse(OBSERVER.observe(text.replace(old, new),
                                             OBSERVER.snapshot(text))["preserved"])

    def test_missing_row_or_data_is_not_a_success(self):
        text = CASES[4]["text"]
        for damaged in ("\n".join(text.splitlines()[:-1]) + "\n", text.replace(" -5 donor", " donor")):
            self.assertFalse(OBSERVER.observe(damaged, OBSERVER.snapshot(text))["preserved"])

    def test_equivalent_logical_spelling_is_allowed(self):
        text = CASES[3]["text"]
        self.assertTrue(OBSERVER.observe(text.replace("T false 1", "true F True"),
                                        OBSERVER.snapshot(text))["preserved"])


if __name__ == "__main__":
    unittest.main()
