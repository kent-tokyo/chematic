"""Keep finite preservation, safe refusals and missing observations distinct."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("xyz_observer", ROOT / "scripts/validation/check_xyz_numeric_read.py")
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


class ObserverTests(unittest.TestCase):
    def test_signed_zero_and_subnormal_bits_are_measured(self):
        for token, bits in [("-0", "8000000000000000"), ("5e-324", "0000000000000001")]:
            self.assertEqual(OBSERVER.output_scalar(f"1\nnumber\nC {token} 0 0\n")["bits"], bits)

    def test_nonfinite_output_is_an_unsafe_acceptance(self):
        for token in ["nan", "inf", "-inf"]:
            result = OBSERVER.output_scalar(f"1\nnumber\nC {token} 0 0\n")
            self.assertEqual(result["status"], "accepted")
            self.assertFalse(result["finite"])

    def test_missing_or_changed_atoms_are_not_numeric_observations(self):
        for text in ["", "0\nnumber\n", "1\nnumber\nO 0 0 0\n", "1\nnumber\nC 0 1 0\n",
                     "1\nnumber\nC 0 0 0\nC 0 0 0\n"]:
            with self.assertRaises((ValueError, IndexError)):
                OBSERVER.output_scalar(text)

    def test_refusals_do_not_count_as_preserved_values(self):
        records = [dict(expected="accept", lane=dict(status="refused", exact=False)),
                   dict(expected="refuse", lane=dict(status="refused", exact=False)),
                   dict(expected="refuse", lane=dict(status="no_record_or_unobservable", exact=False))]
        summary = OBSERVER.summarize(records, "lane")
        self.assertEqual(summary["exact_finite"], 0)
        self.assertEqual(summary["explicit_refusals"], 1)
        self.assertEqual(summary["no_record_or_unobservable"], 1)

    def test_authored_oracles_distinguish_real_zero_from_underflow(self):
        OBSERVER.validate_oracle(dict(token="-0e-999999", bits="8000000000000000"), True)
        OBSERVER.validate_oracle(dict(token="1e-9999", reason="underflow"), False)
        with self.assertRaises(ValueError):
            OBSERVER.validate_oracle(dict(token="1e-9999", bits="0000000000000000"), True)


if __name__ == "__main__":
    unittest.main()
