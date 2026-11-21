import sys
import unittest
from pathlib import Path

import chematic
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_time_budget import attempts


class TimeBudgetTest(unittest.TestCase):
    def test_one_attempt_per_engine(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        other = chematic.from_smiles("CCO").add_hydrogens()
        for engine in ("chematic", "rdkit"):
            coords, record = attempts(engine, mol, other, 5.0, 1)
            self.assertEqual(record["attempted"], 1)
            self.assertEqual(record["accepted"], 1)
            self.assertEqual(len(coords[0]), mol.GetNumAtoms())


if __name__ == "__main__":
    unittest.main()
