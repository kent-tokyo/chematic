"""Declared-stereo preflight preserves its input boundary and reports CIP gaps."""
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

AUDIT = Path(__file__).with_name('audit_declared_stereo.py')


class DeclaredStereoAuditTest(unittest.TestCase):
    def test_cip_gap_and_source_hash(self):
        smiles = 'C[C@]1(C[NH2+]Cc2ccc[nH]2)CN(Cc2ccc[nH]2)[C@@H](c2ccc[nH]2)N(Cc2ccc[nH]2)C1'
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        mol.SetProp('_Name', 'synthetic-cip-gap')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sdf = root / 'sample.sdf'
            writer = Chem.SDWriter(str(sdf))
            writer.write(mol)
            writer.close()
            source_hash = hashlib.sha256(sdf.read_bytes()).hexdigest()
            split = {'source_sha256': source_hash, 'rows': [
                {'index': 0, 'name': 'synthetic-cip-gap', 'split': 'holdout'}]}
            split_path = root / 'split.json'
            split_path.write_text(json.dumps(split))
            output = root / 'audit.json'
            command = [sys.executable, str(AUDIT), str(sdf), str(split_path),
                       '--partition', 'holdout', '--output', str(output)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            result = json.loads(output.read_text())
            self.assertEqual(1, result['records'])
            self.assertEqual(0, result['errors'])
            self.assertGreater(result['tetrahedral_without_cip'], 0)
            split['source_sha256'] = '0' * 64
            split_path.write_text(json.dumps(split))
            self.assertNotEqual(0, subprocess.run(command, capture_output=True).returncode)


if __name__ == '__main__':
    unittest.main()
