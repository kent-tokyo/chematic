import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "speed_frontier", Path(__file__).parents[1] / "bench_python_speed_frontier.py"
)
frontier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(frontier)


class SpeedFrontierTests(unittest.TestCase):
    def test_balanced_permutations_cover_positions_and_pair_precedence(self):
        names = ("baseline", "candidate", "rdkit", "cosmolkit")
        orders = frontier.block_orders(names, 24, True)
        self.assertEqual(len(set(orders)), 24)
        self.assertEqual(orders, frontier.block_orders(names, 24, True))
        for name in names:
            for position in range(4):
                self.assertEqual(sum(order[position] == name for order in orders), 6)
            for other in names:
                if other != name:
                    self.assertEqual(sum(order.index(name) < order.index(other) for order in orders), 12)
        with self.assertRaises(ValueError):
            frontier.block_orders(names, 20, True)

    def test_fixed_slice_preserves_rows_and_does_not_filter_failures(self):
        data = b"CC first\n\ninvalid second\nO third\n"
        self.assertEqual(frontier.corpus_rows(data, 1, 2), ["invalid", "O"])
        for offset, limit in ((-1, 1), (0, 0), (2, 2)):
            with self.assertRaises(ValueError):
                frontier.corpus_rows(data, offset, limit)

    def test_rdkit_comparison_uses_reference_directly_and_cosmolkit_needs_both_gates(self):
        checks = {"candidate_vs_baseline": {"passed": True},
                  "candidate_vs_rdkit": {"passed": True},
                  "cosmolkit_vs_rdkit": {"passed": False}}
        self.assertTrue(frontier.comparable(checks, "rdkit"))
        self.assertTrue(frontier.comparable(checks, "baseline"))
        self.assertFalse(frontier.comparable(checks, "cosmolkit"))
        checks["candidate_vs_rdkit"]["passed"] = False
        self.assertFalse(frontier.comparable(checks, "rdkit"))

    def test_candidate_requires_exact_equality_but_oracle_allows_declared_tolerance(self):
        self.assertFalse(frontier.agreement([1.0], [1.0 + 1e-10], exact=True)["passed"])
        self.assertTrue(frontier.agreement([1.0], [1.0 + 1e-10], exact=False)["passed"])
        self.assertFalse(frontier.agreement([1.0], [1.0 + 1e-8], exact=False)["passed"])

    def test_missing_rows_nonfinite_values_and_wrong_types_never_count_as_agreement(self):
        for values in ([], [float("nan")], [float("inf")], ["1"]):
            self.assertFalse(frontier.agreement([1], values, exact=False)["passed"])
        self.assertFalse(frontier.agreement([1], [True], exact=True)["passed"])
        self.assertFalse(frontier.agreement([1.0], [True], exact=False)["passed"])

    def test_packed_fingerprints_keep_rdkit_bit_indices_and_maccs_offset(self):
        self.assertEqual(frontier.normal(bytes([129, 1]), "chematic", "morgan2_chiral"), [0, 7, 8])
        self.assertEqual(frontier.normal(bytes([129, 1]), "chematic", "maccs"), [1, 8, 9])

    def test_extended_fingerprints_keep_the_same_bit_numbering(self):
        class Bits:
            def GetOnBits(self):
                return [0, 7, 8]

            def on_bits(self):
                return [8, 0, 7]

        for op in ("fp_atom_pair", "fp_torsion", "fp_pattern", "fp_layered", "fp_rdkit"):
            self.assertEqual(frontier.normal(bytes([129, 1]), "chematic", op), [0, 7, 8])
            for engine in ("rdkit", "cosmolkit"):
                self.assertEqual(frontier.normal(Bits(), engine, op), [0, 7, 8])

    def test_mqn_requires_the_complete_ordered_vector(self):
        reference = [list(range(42))]
        self.assertTrue(frontier.agreement(reference, reference, exact=True)["passed"])
        for value in ([list(range(41))], [list(reversed(range(42)))]):
            self.assertFalse(frontier.agreement(reference, value, exact=False)["passed"])


if __name__ == "__main__":
    unittest.main()
