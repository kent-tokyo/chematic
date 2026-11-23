#!/usr/bin/env python3
"""Compare exposed cage conformer energies and gradients with RDKit."""
import json
from pathlib import Path
import sys
import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
sys.path.insert(0, str(Path.cwd() / 'scripts/benchmarks'))
from platinum_development_probe import graph_matches
from platinum_force_fields import compare

names = {'XTULIPKREBMLNM-XHGUTPRISA-N', 'OXBMFCGRQJNAOY-TWXUZQDTSA-N'}
cases = json.loads(Path('validation/results/cod-second-mmff-outliers-v8-development-2026-11-22.json').read_text())['results']
rows = []
for case in cases:
    if case['name'] not in names:
        continue
    mol = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    other = chematic.from_smiles(case['smiles']).add_hydrogens()
    assert graph_matches(other, mol), case['name']
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    assert AllChem.EmbedMolecule(mol, params) == 0, case['name']
    coords = [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    result = compare(mol, other, coords, 'MMFF94')
    rows.append({'name': case['name'], 'energy_error': result['energy_absolute_error'],
                 'gradient_error_max': max(result['gradient_absolute_errors'])})
assert len(rows) == 2
Path('mmff-cage-angle-diagnostic.json').write_text(json.dumps({'rdkit': rdBase.rdkitVersion, 'results': rows},indent=2)+'\n')
print(json.dumps(rows))
assert all(x['energy_error'] < 1e-5 and x['gradient_error_max'] < 1e-4 for x in rows)
