#!/usr/bin/env python3
"""Compare force-field gradients at an exposed near-linear COD geometry."""
import json
import math
from pathlib import Path
import platform

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

case = next(x for x in json.loads(Path('validation/results/cod-second-field-outliers-development-2026-11-22.json').read_text())['results'] if x['name'] == 'YEBDEXSGHHFXMC-FNCQTZNRSA-N')
mol = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 1
if AllChem.EmbedMolecule(mol, params) != 0:
    raise RuntimeError('RDKit embedding failed')
coords = [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
other = chematic.from_smiles(case['smiles']).add_hydrogens()
u = [coords[10][i]-coords[9][i] for i in range(3)]
v = [coords[0][i]-coords[9][i] for i in range(3)]
dot = sum(a*b for a,b in zip(u,v)) / math.sqrt(sum(a*a for a in u)*sum(b*b for b in v))
angle = math.degrees(math.acos(max(-1.0,min(1.0,dot))))
out = {'name': case['name'], 'platform': platform.platform(), 'rdkit': rdBase.rdkitVersion,
       'angle_atoms': [10,9,0], 'angle_degrees': angle, 'atom': 9, 'axis': 1, 'fields': {}}
for variant in ('MMFF94','UFF'):
    ours = other.mmff_force_field(coords) if variant == 'MMFF94' else other.uff_force_field(coords)
    reference = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol)) if variant == 'MMFF94' else AllChem.UFFGetMoleculeForceField(mol)
    gradients = ours.gradient()
    numeric = {}
    for h in (1e-5,1e-6,1e-7,1e-8,1e-9):
        positions = [p[:] for p in coords]
        positions[9][1] += h
        plus = ours.set_positions_(positions).energy()
        positions[9][1] -= 2*h
        minus = ours.set_positions_(positions).energy()
        numeric[str(h)] = (plus-minus)/(2*h)
    out['fields'][variant] = {'chematic_gradient': gradients[9][1],
                              'rdkit_gradient': reference.CalcGrad()[28],
                              'finite_difference': numeric}
    if variant == 'MMFF94':
        ours.set_positions_(coords)
        out['fields'][variant]['consistent_plane_gradient'] = ours.gradient_with_consistent_planes()[9][1]
Path('near-linear-second-cod.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out))
