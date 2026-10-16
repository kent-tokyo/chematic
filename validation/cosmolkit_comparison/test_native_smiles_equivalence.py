"""False-positive controls for the native SMILES comparison oracle."""

import unittest

from rdkit import Chem

from check_native_smiles_roundtrip import classify


class NativeSmilesEquivalenceTest(unittest.TestCase):
    def check_pair(self, source: str, read_back: str, expected: str) -> None:
        source_mol = Chem.MolFromSmiles(source)
        read_back_mol = Chem.MolFromSmiles(read_back)
        self.assertIsNotNone(source_mol)
        self.assertIsNotNone(read_back_mol)
        self.assertEqual(classify(source_mol, read_back_mol), expected)

    def test_equivalent_permutations_and_resonance(self) -> None:
        self.check_pair("C[Pt@SP1](Cl)(N)N", "N[Pt@SP1](N)(Cl)C", "stereo_isomorphic")
        self.check_pair("F[P@TB1](F)(F)(F)F", "[P@TB10](F)(F)(F)(F)F", "stereo_isomorphic")
        self.check_pair("Cc1cccc[c+]1", "CC1=[C+]C=CC=C1", "resonance")

    def test_distinct_stereoisomers_and_connectivity(self) -> None:
        self.check_pair("F[Pt@SP1](Cl)(Br)I", "F[Pt@SP2](Cl)(Br)I", "unresolved")
        self.check_pair("F[P@TB1](Cl)(Br)(I)C", "F[P@TB2](Cl)(Br)(I)C", "unresolved")
        self.check_pair("F[S@OH1](Cl)(Br)(I)(N)O", "F[S@OH2](Cl)(Br)(I)(N)O", "unresolved")
        self.check_pair("F/C=C/F", "F/C=C\\F", "unresolved")
        self.check_pair("F[C@H](Cl)Br", "F[C@@H](Cl)Br", "unresolved")
        self.check_pair("Cc1ccccc1", "Cc1ccccn1", "unresolved")


if __name__ == "__main__":
    unittest.main()
