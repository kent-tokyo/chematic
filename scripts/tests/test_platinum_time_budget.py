import sys
import unittest
from pathlib import Path

import chematic
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_time_budget import attempts
from summarize_platinum_time_budget import recovered


class TimeBudgetTest(unittest.TestCase):
    def test_recovery_requires_screened_geometry_when_requested(self):
        row = {
            "chematic": {"rmsd": [0.5], "stereo": [{"atoms": [], "bonds": []}]},
            "geometry_reference": {
                "bonds_supported": 1,
                "bonds_total": 1,
                "angles_supported": 0,
                "angles_total": 0,
            },
        }
        self.assertTrue(recovered(row, "chematic"))
        self.assertFalse(recovered(row, "chematic", True))
        row["chematic"]["geometry"] = [
            {
                "finite": True,
                "bond_violations": 0,
                "angle_violations": 0,
                "heavy_atom_clashes": 0,
            }
        ]
        self.assertTrue(recovered(row, "chematic", True))
        self.assertFalse(recovered(row, "chematic", True, True))
        row["chematic"]["statuses"] = [1]
        self.assertFalse(recovered(row, "chematic", True, True))
        row["chematic"]["statuses"] = [0]
        self.assertTrue(recovered(row, "chematic", True, True))

    def test_one_attempt_per_engine(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        other = chematic.from_smiles("CCO").add_hydrogens()
        for engine in ("chematic", "rdkit"):
            coords, record = attempts(engine, mol, other, 5.0, 1)
            self.assertEqual(record["attempted"], 1)
            self.assertEqual(record["accepted"], 1)
            self.assertEqual(len(record["statuses"]), 1)
            self.assertEqual(len(coords[0]), mol.GetNumAtoms())


if __name__ == "__main__":
    unittest.main()
