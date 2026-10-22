import copy
import unittest

from scripts.validation.check_cjson_mol2_formal_charges import (
    evaluate, mol2_snapshot, packet, same_bits, sdf_graph,
)


class ChargePacketOracleTests(unittest.TestCase):
    def test_explicit_charge_wins_over_type_and_partial_column(self):
        source = """@<TRIPOS>MOLECULE
original oracle example
2 0 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
7 N1 0 0 0 N.4 1 LIG 0.7
42 O2 1 0 0 O.2 1 LIG 0.7
@<TRIPOS>UNITY_ATOM_ATTR
42 1
charge -2
7 1
charge -1
@<TRIPOS>BOND
"""
        for value in [source, source.replace("@<TRIPOS>BOND\n", "").replace(
                "@<TRIPOS>UNITY_ATOM_ATTR", "@<TRIPOS>BOND\n@<TRIPOS>UNITY_ATOM_ATTR")]:
            snap = mol2_snapshot(value)
            self.assertEqual(snap["formal_charges"], [-1, -2])
            self.assertEqual(snap["partial_charges"], [0.7, 0.7])
        with self.assertRaises(ValueError):
            mol2_snapshot(source.replace("42 1", "999 1"))
        with self.assertRaises(ValueError):
            mol2_snapshot(source.replace("7 1\ncharge -1", "42 1\ncharge -1"))

    def test_exact_precision_does_not_ignore_signed_zero_or_ulp(self):
        self.assertFalse(same_bits([-0.0], [0.0]))
        self.assertFalse(same_bits([1.0000000000000002], [1.0]))
        self.assertFalse(same_bits([0], [None]))
        self.assertFalse(same_bits([], [0]))
        self.assertTrue(same_bits([5e-324, -0.0], [5e-324, -0.0]))

    def test_missing_charge_and_changed_bond_cannot_pass(self):
        _, source = packet()[-1]  # Authored connected nitrate graph.
        snap = {"elements": [7, 8, 8, 8], "formal_charges": [1, 0, -1, -1],
                "coordinates": source["atoms"]["coords"]["3d"],
                "bonds": [[0, 1, 2], [0, 2, 1], [0, 3, 1]],
                "partial_charges": [0.7, -0.5, -0.6, -0.6], "charge_method": "USER_CHARGES"}
        self.assertTrue(all(evaluate(source, snap).values()))
        bad = copy.deepcopy(snap)
        bad["formal_charges"][3] = 0
        self.assertFalse(evaluate(source, bad)["formal_charges"])
        bad = copy.deepcopy(snap)
        bad["bonds"][0][2] = 1
        self.assertFalse(evaluate(source, bad)["elements_and_indexed_bonds"])

    def test_external_sdf_charge_codes_and_m_chg_override(self):
        source = ("oracle\nprogram\n\n  2  0  0  0  0  0            999 V2000\n"
                  "    0.0000    0.0000    0.0000 Na  0  3  0  0  0  0\n"
                  "    1.0000    0.0000    0.0000 Cl  0  5  0  0  0  0\n"
                  "M  CHG  1   1   2\nM  END\n$$$$\n")
        self.assertEqual(sdf_graph(source), {
            "elements": [11, 17], "formal_charges": [2, -1], "bonds": []})
        with self.assertRaises(ValueError):
            sdf_graph(source.replace("M  CHG  1", "M  CHG  2"))
        with self.assertRaises(ValueError):
            sdf_graph(source.replace("M  END\n", ""))


if __name__ == "__main__":
    unittest.main()
