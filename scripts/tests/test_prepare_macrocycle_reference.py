"""Verify macrocycle reference input checks without scoring conformers."""
import importlib.util
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem


MODULE = Path(__file__).resolve().parents[1] / "validation" / "prepare_macrocycle_reference.py"
spec = importlib.util.spec_from_file_location("prepare_macrocycle_reference", MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def sdf(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=17) != 0:
        raise AssertionError("test molecule did not embed")
    return (Chem.MolToMolBlock(mol) + "$$$$\n").encode()


class ReferencePreparationTest(unittest.TestCase):
    def test_graph_identity_and_isotope_stereo_are_checked(self):
        raw = sdf(r"[2H]/C([H])=C(/F)Cl")
        details = mod.validate_reference(raw, r"[2H]/C([H])=C(/F)Cl")
        self.assertGreater(details["atoms"], details["heavy_atoms"])
        with self.assertRaisesRegex(ValueError, "graphs differ"):
            mod.validate_reference(raw, r"[2H]/C([H])=C(\F)Cl")
        with self.assertRaisesRegex(ValueError, "graphs differ"):
            mod.validate_reference(raw, r"[3H]/C([H])=C(/F)Cl")

    def test_title_replacement_preserves_all_other_sdf_bytes(self):
        raw = sdf("C/C=C/C")
        renamed = mod.title_replaced(raw, "case-01")
        self.assertEqual(renamed.splitlines(keepends=True)[1:], raw.splitlines(keepends=True)[1:])
        self.assertEqual(renamed.splitlines()[0], b"case-01")
        with self.assertRaises(ValueError):
            mod.title_replaced(raw + raw, "case-01")

    def test_archive_requires_exactly_one_reference_per_case(self):
        raw = sdf("C/C=C/C")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                info = tarfile.TarInfo("MacroCycle-Benchmark-abc/cases/case-01/ligand_unprotonated.sdf")
                info.size = len(raw)
                archive.addfile(info, io.BytesIO(raw))
            self.assertEqual(mod.read_cases(path, "abc", {"case-01"}), {"case-01": raw})
            with self.assertRaisesRegex(ValueError, "missing"):
                mod.read_cases(path, "abc", {"case-01", "case-02"})
            with self.assertRaisesRegex(ValueError, "missing"):
                mod.read_cases(path, "wrong-commit", {"case-01"})

    def test_preparation_records_graph_mismatches_and_reindexes_valid_cases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            index = [
                {"case_id": "dev", "unprotonated_isomeric_smiles": "C/C=C/C"},
                {"case_id": "holdout-mismatch", "unprotonated_isomeric_smiles": "C/C=C/C"},
                {"case_id": "excluded", "unprotonated_isomeric_smiles": "CC"},
            ]
            index_path = root / "index.csv"
            with index_path.open("w", newline="", encoding="utf-8") as out:
                writer = csv.DictWriter(out, fieldnames=index[0].keys())
                writer.writeheader()
                writer.writerows(index)
            index_bytes = index_path.read_bytes()
            partitions = ["development", "holdout", "excluded_known_scaffold"]
            rows = [{"index": i, "name": r["case_id"], "split": partitions[i],
                     "input_smiles_sha256": hashlib.sha256(r["unprotonated_isomeric_smiles"].encode()).hexdigest()}
                    for i, r in enumerate(index)]
            split = {"source_commit": "abc", "index_sha256": mod.sha(index_bytes),
                     "rdkit_version": mod.rdBase.rdkitVersion, "invalid_exclusions": [], "rows": rows}
            split_bytes = (json.dumps(split) + "\n").encode()
            (root / "split.json").write_bytes(split_bytes)
            audit = {"index_sha256": mod.sha(index_bytes), "split_sha256": mod.sha(split_bytes),
                     "invalid_exclusion_count": 0, "overlaps": []}
            (root / "audit.json").write_text(json.dumps(audit))
            raw = [sdf("C/C=C/C"), sdf(r"C/C=C\C"), sdf("CC")]
            with tarfile.open(root / "archive.tar.gz", "w:gz") as archive:
                for record, content in zip(index, raw):
                    info = tarfile.TarInfo(f"MacroCycle-Benchmark-abc/cases/{record['case_id']}/ligand_unprotonated.sdf")
                    info.size = len(content)
                    archive.addfile(info, io.BytesIO(content))
            command = [sys.executable, str(MODULE), str(root / "archive.tar.gz"), str(index_path),
                       str(root / "split.json"), str(root / "audit.json"), "--sdf", str(root / "prepared.sdf"),
                       "--prepared-split", str(root / "prepared.json"), "--report", str(root / "report.json")]
            subprocess.run(command, check=True, capture_output=True, text=True)
            report = json.loads((root / "report.json").read_text())
            result = json.loads((root / "prepared.json").read_text())
            self.assertEqual(report["prepared_records"], 2)
            self.assertEqual(report["failure_counts"], {"holdout": 1})
            self.assertIn("graphs differ", report["input_failures"][0]["error"])
            self.assertEqual([(r["index"], r["source_index"]) for r in result["rows"]], [(0, 0), (1, 2)])
            supplier = Chem.SDMolSupplier(str(root / "prepared.sdf"), removeHs=False)
            self.assertEqual([m.GetProp("_Name") for m in supplier], ["dev", "excluded"])


if __name__ == "__main__":
    unittest.main()
