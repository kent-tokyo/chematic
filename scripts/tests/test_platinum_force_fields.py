import sys
import unittest
from pathlib import Path

import chematic
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_force_fields import compare


class ForceFieldComparisonTest(unittest.TestCase):
    def test_common_coordinates_and_finite_differences(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        coords = [
            list(mol.GetConformer().GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]
        other = chematic.from_smiles("CCO").add_hydrogens()
        for variant in ("MMFF94", "MMFF94s", "UFF"):
            with self.subTest(variant=variant):
                result = compare(mol, other, coords, variant)
                self.assertTrue(result["rdkit_supported"])
                self.assertTrue(result["chematic_supported"])
                self.assertLess(result["energy_absolute_error"], 1e-7)
                self.assertLess(max(result["gradient_absolute_errors"]), 1e-6)
                self.assertLess(max(result["finite_difference_absolute_errors"]), 1e-4)


if __name__ == "__main__":
    unittest.main()
