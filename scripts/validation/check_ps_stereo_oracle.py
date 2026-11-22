#!/usr/bin/env python3
"""Check the exposed phosphorus/sulfur fixtures and their reflections with RDKit."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

from rdkit import Chem, rdBase

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_stereo_strict import stereo_violations_strict  # noqa: E402

assert rdBase.rdkitVersion == '2026.09.2'
raw = gzip.decompress(Path('validation/results/cod-exposed-stereo-rdkit-input-2026-11-22.json.gz').read_bytes())
assert hashlib.sha256(raw).hexdigest() == 'd6029b3350c3ea6f33f14b4713bf6b6be1830a10299ac730c4bb5a31045a167d'
source = json.loads(raw)
centers = {'OBNMTDLZGXAGLV-BEQSNLQOSA-N': 5, 'WZRCNRWIYMSYFJ-GHHQMCIASA-O': 9,
           'PMVWUWCKINFFIE-UHFFFAOYSA-N': 5, 'CKOAIHXEXMUXMP-MHZLTWQESA-N': 1}
expected = []
for row in source:
    if row['name'] in centers:
        for c in row['conformers']:
            if not c['external']['atoms'] and not c['external']['bonds']:
                expected.append({'name': row['name'], 'smiles': row['smiles'],
                                 'center': centers[row['name']], 'coords': c['coords']})
fixture = json.loads(Path('crates/chematic-3d/src/stereo_constraints_ps_oracle.json').read_text())
assert fixture['rows'] == expected and len(expected) == 17
controls = json.loads(Path('crates/chematic-3d/src/stereo_constraints_ps_controls.json').read_text())
assert controls['rdkit'] == '2026.09.2' and len(controls['rows']) == 20
assert {r['smiles'] for r in controls['rows']} == {
    'C[S@+](CC)CCC', 'C[S@](=NC)CC', 'C[S@](=O)(=NC)CC', 'N=[P@-](=S)(N)O'}
for row in expected + controls['rows']:
    mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    original = stereo_violations_strict(mol, [row['coords']])[0]
    reflected = stereo_violations_strict(mol, [[[-x, y, z] for x, y, z in row['coords']]])[0]
    assert not (original['atoms'] or original['bonds'])
    assert row['center'] in reflected['violated_atoms']
print('RDKit 2026.09.2 independently confirmed all 37 originals and reflected centers')
