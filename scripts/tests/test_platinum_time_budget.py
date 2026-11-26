import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

import chematic
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_time_budget import attempts, hard_attempts
from summarize_platinum_time_budget import recovered


class TimeBudgetTest(unittest.TestCase):
    def test_attempts_report_started_and_empty_results(self):
        other = Mock()
        other.etkdg_embed_multiple_confs_mmff_corrected.return_value = []
        events = []
        coords, record = attempts("chematic", None, other, 5.0, 2, on_event=events.append)
        self.assertEqual(coords, [])
        self.assertEqual(record["attempted"], 2)
        self.assertEqual(record["no_conformer"], 2)
        self.assertEqual(events, [("started", 1), ("no_conformer", 1),
                                  ("started", 2), ("no_conformer", 2)])

    def test_hard_deadline_terminates_an_attempt(self):
        import time
        start = time.perf_counter()
        coords, record = hard_attempts("chematic", "C" * 25, 0.001, 20)
        self.assertTrue(record["timed_out"])
        self.assertTrue(record["attempted_is_lower_bound"])
        self.assertEqual(record["accepted"], len(coords))
        self.assertLess(time.perf_counter() - start, 15)

    def test_hard_deadline_retains_completed_attempt(self):
        for engine in ("chematic", "rdkit"):
            coords, record = hard_attempts(engine, "CCO", 5.0, 1)
            self.assertFalse(record["timed_out"])
            self.assertEqual(record["attempted"], 1)
            self.assertEqual(record["accepted"], 1)
            self.assertEqual(len(coords), 1)

    def test_hard_deadline_keeps_force_field_failure_reason(self):
        coords, record = hard_attempts("chematic", "CB(C)C", 5.0, 1, profile="mixed-fallback")
        self.assertEqual(len(coords), 1)
        self.assertEqual(record["statuses"], [-1])
        self.assertIn("type assignment", record["force_field_errors"][0])

    def test_hard_deadline_reports_unevaluable_stereo(self):
        smiles = "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F"
        coords, record = hard_attempts("chematic", smiles, 5.0, 2, profile="mixed-fallback")
        self.assertEqual(coords, [])
        self.assertEqual(record["no_conformer"], 2)
        self.assertEqual(record["stereo_rejected"]["unevaluable"], 2)
        self.assertEqual(record["stereo_rejected"]["unclassified_attempts"], 0)

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

    def test_summary_includes_convergence_when_a_row_errors(self):
        complete = {"accepted": 1, "rmsd": [0.5], "stereo": [{"atoms": [], "bonds": []}], "statuses": [0]}
        report = {"seconds": 0.2, "max_candidates": 20, "results": [
            {"chematic": complete, "rdkit": complete},
            {"error": "failed after first engine", "chematic": {"accepted": 0, "statuses": []}},
        ]}
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.json"
            output_path = Path(directory) / "summary.json"
            input_path.write_text(json.dumps(report), encoding="utf-8")
            script = Path(__file__).resolve().parents[1] / "benchmarks" / "summarize_platinum_time_budget.py"
            subprocess.run([sys.executable, str(script), str(input_path), "--output", str(output_path)], check=True, capture_output=True)
            summary = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(summary["errors"], 1)
        self.assertEqual(summary["converged_recovery"]["chematic"], 1)
        self.assertEqual(summary["converged_recovery"]["rdkit"], 1)

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
