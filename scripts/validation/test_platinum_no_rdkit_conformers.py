"""The fixed-count probe records an empty RDKit ensemble without optimizing it."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import platinum_development_probe as probe  # noqa: E402


class EmptyRdkitEnsembleTest(unittest.TestCase):
    def test_no_conformer_is_not_optimizer_error(self):
        reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(reference, randomSeed=42), 0)
        with patch.object(probe.AllChem, "EmbedMultipleConfs", return_value=[]), patch.object(
            probe.AllChem, "MMFFOptimizeMoleculeConfs", side_effect=AssertionError("must not optimize")
        ) as optimize:
            row = probe.score(reference, 1, 42, "corrected", False, 0)
        optimize.assert_not_called()
        self.assertNotIn("rdkit_error", row)
        self.assertEqual(row["rdkit_count"], 0)
        self.assertEqual(row["rdkit_rmsd"], [])
        self.assertEqual(row["rdkit_statuses"], [])


if __name__ == "__main__":
    unittest.main()
