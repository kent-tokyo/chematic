import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).parents[1]))
from bench_python_scalar_reuse import validate_values


class ScalarReuseTests(unittest.TestCase):
    def test_every_retained_repeat_is_checked(self):
        validate_values([1.0, 2.0], [1.0, 2.0, 1.0, 2.0], 2)
        for values in ([1.0, 2.0], [1.0, 3.0, 1.0, 2.0], [1.0, 2.0, 1.0, 3.0],
                       [1.0, 2.0, 1.0, float("nan")], [1.0, 2.0, 1.0, 2]):
            with self.assertRaises(ValueError):
                validate_values([1.0, 2.0], values, 2)

    def test_nonfinite_or_empty_references_do_not_count_as_success(self):
        for reference in ([], [float("inf")], [float("nan")], [True]):
            with self.assertRaises(ValueError):
                validate_values(reference, reference, 1)
