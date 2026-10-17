"""False-positive controls for the native SMILES comparison oracle."""

import unittest
from typing import Union

from rdkit import Chem

from check_native_smiles_roundtrip import classify


class NativeSmilesEquivalenceTest(unittest.TestCase):
    def check_pair(
        self, source: str, read_back: str, expected: Union[str, tuple[str, ...]]
    ) -> None:
        source_mol = Chem.MolFromSmiles(source)
        read_back_mol = Chem.MolFromSmiles(read_back)
        self.assertIsNotNone(source_mol)
        self.assertIsNotNone(read_back_mol)
        self.assertIn(
            classify(source_mol, read_back_mol, source, read_back),
            (expected,) if isinstance(expected, str) else expected,
        )

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

    def test_raw_aromatic_graph_survives_different_kekulization(self) -> None:
        self.check_pair(
            "Cc1c[c]([Pt])c(Oc2ccc(O)cc2CC(=O)O)s(Cl)c1",
            "c1(C)cs(c(Oc2ccc(O)cc2CC(=O)O)c([Pt])c1)Cl",
            "raw_aromatic_graph",
        )
        self.check_pair(
            "P(OC)(OC)(Nc1cc(OC)c(cc1)Nc1c2c(cccc2)nc2c1cc([as]c2)C)=O",
            "n2c1c(c(Nc4c(OC)cc(cc4)NP(=O)(OC)OC)c3cc(C)[as]cc23)cccc1",
            ("strict", "raw_aromatic_graph"),
        )

    def test_raw_aromatic_graph_preserves_remote_stereo(self) -> None:
        source = "c12ccccc2[as]ccc1C(N1CCN(C[C@H](CS)N)[C@@H](CCOC)C1)=O"
        reordered = "c23c(cc[as]c2cccc3)C(N1CCN(C[C@H](CS)N)[C@H](C1)CCOC)=O"
        self.check_pair(source, reordered, "raw_aromatic_graph")
        self.check_pair(
            source,
            reordered.replace("C[C@H](CS)N", "C[C@@H](CS)N"),
            "unresolved",
        )

    def test_raw_aromatic_graph_does_not_waive_other_differences(self) -> None:
        self.check_pair("F/C=C/c1ccccc1", "F/C=C\\c1ccccc1", "unresolved")
        self.check_pair("F[C@H](Cl)c1ccccc1", "F[C@@H](Cl)c1ccccc1", "unresolved")
        self.check_pair("Cc1ccccc1O", "Cc1cccc(O)c1", "unresolved")
        self.check_pair("c1cc([13CH3])ccc1", "Cc1ccccc1", "unresolved")
        self.check_pair("[c:7]1ccccc1", "c1ccccc1", "unresolved")
        self.check_pair("[H]/N=C(/N)N", "N=C(N)N", "unresolved")


if __name__ == "__main__":
    unittest.main()
