#!/usr/bin/env python3
"""Semantic claim gates for the two-descriptor benchmark."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bench_python_scalar_bundle import agreement, comparable

class BundleGateTests(unittest.TestCase):
    def test_both_descriptors_must_agree(self):
        self.assertTrue(agreement([[1.0, 2.0]], [[1.0, 2.0]], exact=True)["passed"])
        for values in ([[1.0, 3.0]], [[3.0, 2.0]], [[1.0]], [], [[True, 2.0]]):
            self.assertFalse(agreement([[1.0, 2.0]], values, exact=True)["passed"])

    def test_nonfinite_and_signed_zero_are_not_exact_passes(self):
        for value in (math.nan, math.inf, -math.inf):
            self.assertFalse(agreement([[value, 2.0]], [[value, 2.0]], exact=True)["passed"])
        self.assertFalse(agreement([[0.0, 2.0]], [[-0.0, 2.0]], exact=True)["passed"])
        self.assertTrue(agreement([[0.0, 2.0]], [[-0.0, 2.0]], exact=False)["passed"])

    def test_tolerance_and_comparator_require_all_agreements(self):
        self.assertTrue(agreement([[1.0, 2.0]], [[1.0 + 1e-10, 2.0]], exact=False)["passed"])
        self.assertFalse(agreement([[1.0, 2.0]], [[1.0 + 1e-8, 2.0]], exact=False)["passed"])
        checks={"candidate_vs_baseline":{"passed":True}, "candidate_vs_rdkit":{"passed":True},
                "cosmolkit_vs_rdkit":{"passed":False}}
        self.assertTrue(comparable(checks,"baseline"))
        self.assertTrue(comparable(checks,"rdkit"))
        self.assertFalse(comparable(checks,"cosmolkit"))
        checks["candidate_vs_rdkit"]["passed"]=False
        self.assertFalse(comparable(checks,"rdkit"))

if __name__ == "__main__":
    unittest.main()
