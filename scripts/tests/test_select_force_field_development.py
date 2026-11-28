import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem


class SelectForceFieldDevelopmentTest(unittest.TestCase):
    def test_selection_is_deterministic_and_never_uses_holdout(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sdf = directory / "input.sdf"
            writer = Chem.SDWriter(str(sdf))
            records = [("boron", "B(O)O", "development"),
                       ("selenium", "C[Se]C", "development"),
                       ("sulfur", "CS", "development"),
                       ("other", "CCO", "development"),
                       ("excluded", "CCN", "holdout")]
            rows = []
            for index, (name, smiles, partition) in enumerate(records):
                molecule = Chem.MolFromSmiles(smiles)
                self.assertIsNotNone(molecule)
                molecule.SetProp("_Name", name)
                writer.write(molecule)
                rows.append({"index": index, "name": name, "split": partition})
            writer.close()
            split = directory / "split.json"
            split.write_text(json.dumps({"source_sha256": hashlib.sha256(sdf.read_bytes()).hexdigest(),
                                        "rows": rows}), encoding="utf-8")
            script = Path(__file__).resolve().parents[1] / "validation" / "select_force_field_development.py"
            output = directory / "selected.json"
            command = [sys.executable, str(script), str(sdf), str(split), "--output", str(output)]
            subprocess.run(command, check=True, capture_output=True)
            first = json.loads(output.read_text(encoding="utf-8"))
            subprocess.run(command, check=True, capture_output=True)
            self.assertEqual(first, json.loads(output.read_text(encoding="utf-8")))
            self.assertEqual(first["counts"], {"boron": 1, "selenium": 1, "phosphorus_sulfur": 1, "other": 1})
            self.assertEqual({row["name"] for row in first["results"]}, {"boron", "selenium", "sulfur", "other"})
            split.write_text(json.dumps({"source_sha256": "wrong", "rows": rows}), encoding="utf-8")
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
