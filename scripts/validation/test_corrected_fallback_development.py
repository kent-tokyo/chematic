#!/usr/bin/env python3
"""Verify a single already-exposed macrocycle through the Python API."""
import gzip
import hashlib
import json
import platform
import sys
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_stereo_strict import stereo_violations_strict  # noqa: E402

source = gzip.decompress(Path('validation/results/macrocycle-pdb-linux-development-v11-2026-11-22.json.gz').read_bytes())
assert hashlib.sha256(source).hexdigest() == '821e4cbb342bbb94db5f4a44a62560d17001ead9e1a7f1190f0cb58ab3d14708'
rows = [r for r in json.loads(source)['results'] if r['name'] == '5vk0__Ligand_phospholipids_ALA_L_4']
assert len(rows) == 1
smiles = rows[0]['smiles']
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
other = chematic.from_smiles(smiles).add_hydrogens()
kwargs = dict(random_seed=42, max_mmff_iters=200, repair_mmff_iters=400)
corrected = other.etkdg_embed_multiple_confs_mmff_corrected(20, **kwargs)
fallback = other.etkdg_embed_multiple_confs_mmff_corrected_fallback(20, **kwargs)
assert fallback
if corrected:
    assert corrected == fallback, 'fallback altered a successful corrected ensemble'
checks = stereo_violations_strict(mol, [r[2] for r in fallback])
assert all(not (x['violated_atoms'] or x['violated_bonds'] or x['unevaluable_atoms'] or x['unevaluable_bonds']) for x in checks)
result = {'case': rows[0]['name'], 'platform': platform.platform(), 'rdkit': rdBase.rdkitVersion,
          'corrected_count': len(corrected), 'fallback_count': len(fallback),
          'fallback_statuses': [r[0] for r in fallback], 'stereo_checked': len(checks),
          'stereo_contradicted': 0, 'stereo_unevaluable': 0}
Path('fallback-development-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
