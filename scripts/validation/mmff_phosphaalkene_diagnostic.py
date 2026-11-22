#!/usr/bin/env python3
"""Check an exposed P=C field outlier against RDKit on common coordinates."""
import json
from pathlib import Path
import sys
import chematic
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
sys.path.insert(0,str(Path.cwd()/'scripts/benchmarks'))
from platinum_development_probe import graph_matches
from platinum_force_fields import compare

name='OZQDJSGRCNJHSF-UHFFFAOYSA-N'
rows=json.loads(Path('validation/results/cod-first-field-outliers-development-2026-11-22.json').read_text())['results']
row=next(r for r in rows if r['name']==name)
mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
other=chematic.from_smiles(row['smiles']).add_hydrogens()
assert graph_matches(other,mol)
props=AllChem.MMFFGetMoleculeProperties(mol)
assert props is not None
params=AllChem.ETKDGv3();params.randomSeed=42;params.numThreads=1
assert AllChem.EmbedMolecule(mol,params)==0
coords=[list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
field=compare(mol,other,coords,'MMFF94')
out={'rdkit':rdBase.rdkitVersion,'name':name,'energy_abs_error':field['energy_absolute_error'],
     'gradient_abs_error_max':max(field['gradient_absolute_errors'],default=0)}
Path('mmff-phosphaalkene-diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
assert out['energy_abs_error']<1e-5 and out['gradient_abs_error_max']<1e-4
