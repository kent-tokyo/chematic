#!/usr/bin/env python3
"""Generate RDKit geometries for the seven exposed cases still without candidates."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_development_probe import graph_matches  # noqa: E402
from platinum_stereo_strict import stereo_violations_strict  # noqa: E402

RDLogger.DisableLog('rdApp.*')
source = gzip.decompress(Path('validation/results/cod-exposed-empty-linux-corrected-strict-2026-11-22.json.gz').read_bytes())
assert hashlib.sha256(source).hexdigest() == '83c31eae906112bc82840e713dbf3ad90bb9255e18eb4ec28dcbde8e5372e146'
rows = [r for r in json.loads(source)['results'] if r.get('chematic_count') == 0]
assert len(rows) == 7
output = []
for row in rows:
    mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    other = chematic.from_smiles(row['smiles']).add_hydrogens()
    assert graph_matches(other, mol)
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads = 42, 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))[:5]
    conformers = []
    for cid in ids:
        status = AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
        coords = [list(mol.GetConformer(cid).GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        external = stereo_violations_strict(mol, [coords])[0]
        conformers.append({'coords': coords, 'status': status, 'external': external})
    output.append({'name': row['name'], 'smiles': row['smiles'], 'conformers': conformers})
Path(sys.argv[1]).write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
print('cases:', len(output), 'conformers:', sum(len(r['conformers']) for r in output))
