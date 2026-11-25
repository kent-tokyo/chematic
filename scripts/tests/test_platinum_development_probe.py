import sys
import unittest
from unittest.mock import patch
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_development_probe import rms_values, stereo_violations, score
import platinum_development_probe as probe
from summarize_platinum_probe import recovered


class PlatinumScoringTest(unittest.TestCase):
    def test_independent_scoring_runs_both_engines_when_precheck_fails(self):
        reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(reference, randomSeed=42), 0)
        with patch.object(AllChem, "MMFFHasAllMoleculeParams", return_value=False):
            skipped = score(reference, 1, 42, "compatible", False, 0)
            independent = score(reference, 1, 42, "compatible", False, 0, True)
        self.assertEqual(skipped["status"], "rdkit_mmff_unsupported")
        self.assertEqual(independent["status"], "ok")
        self.assertEqual(independent["chematic_count"], 1)
        self.assertEqual(independent["rdkit_count"], 1)

    def test_independent_scoring_continues_after_chematic_parse_error(self):
        reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(reference, randomSeed=42), 0)
        with patch.object(probe.chematic, "from_smiles", side_effect=ValueError("parse failed")):
            result = score(reference, 1, 42, "compatible", False, 0, True)
        self.assertIn("parse failed", result["chematic_error"])
        self.assertEqual(result["rdkit_count"], 1)

    def test_recovery_requires_valid_stereo_and_counts_failures(self):
        good = {"atoms": [], "bonds": []}
        bad = {"atoms": [1], "bonds": []}
        row = {
            "status": "ok",
            "chematic_rmsd": [0.5, 1.1],
            "chematic_stereo": [bad, good],
            "chematic_energies": [1.0, 2.0],
        }
        self.assertFalse(recovered(row, "chematic", "minimized"))
        row["chematic_stereo"] = [good, good]
        self.assertTrue(recovered(row, "chematic", "minimized", True))
        self.assertFalse(recovered(row, "chematic", "minimized", require_geometry=True))
        row["geometry_reference"] = {
            "bonds_supported": 1,
            "bonds_total": 1,
            "angles_supported": 0,
            "angles_total": 0,
        }
        screen = {
            "finite": True,
            "bond_violations": 0,
            "angle_violations": 0,
            "heavy_atom_clashes": 0,
        }
        row["chematic_geometry"] = [screen, screen]
        self.assertTrue(recovered(row, "chematic", "minimized", require_geometry=True))
        row["chematic_error"] = "failure"
        self.assertFalse(recovered(row, "chematic", "minimized"))

    def test_rmsd_aligns_translated_structure(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        reference = Chem.RemoveHs(mol)
        coordinates = [
            list(mol.GetConformer().GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]
        translated = [[x + 4, y - 2, z + 1] for x, y, z in coordinates]
        self.assertLess(rms_values(mol, reference, [translated])[0], 1e-6)

    def test_stereo_check_detects_inverted_tetrahedron(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("F[C@](Cl)(Br)I"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        coordinates = [
            list(mol.GetConformer().GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]
        reflected = [[-x, y, z] for x, y, z in coordinates]
        original, inverted = stereo_violations(mol, [coordinates, reflected])
        self.assertEqual(original["atoms"], [])
        self.assertEqual(inverted["atoms"], [1])

    def test_stereo_check_detects_double_bond_inversion(self):
        trans = Chem.AddHs(Chem.MolFromSmiles("F/C=C/F"))
        cis = Chem.AddHs(Chem.MolFromSmiles("F/C=C\\F"))
        self.assertEqual(AllChem.EmbedMolecule(trans, randomSeed=42), 0)
        self.assertEqual(AllChem.EmbedMolecule(cis, randomSeed=42), 0)
        coords = lambda mol: [
            list(mol.GetConformer().GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]
        original, inverted = stereo_violations(trans, [coords(trans), coords(cis)])
        self.assertEqual(original["bonds"], [])
        self.assertEqual(inverted["bonds"], [1])


if __name__ == "__main__":
    unittest.main()
