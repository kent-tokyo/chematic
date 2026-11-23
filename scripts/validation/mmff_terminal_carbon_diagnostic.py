#!/usr/bin/env python3
"""Check exposed terminal-carbon MMFF outliers against RDKit on common coordinates."""
import json
from pathlib import Path
import sys
import chematic
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
sys.path.insert(0,str(Path.cwd()/'scripts/benchmarks'))
from platinum_development_probe import graph_matches
from platinum_force_fields import compare

names={'OCRDFSOMVZHDRK-UHFFFAOYSA-N','BIHFVVRTTQNXMU-UHFFFAOYSA-N','CSSSZUOEJPLMJD-UHFFFAOYSA-M'}
rows=json.loads(Path('validation/results/cod-first-field-outliers-development-2026-11-22.json').read_text())['results']
out=[]
for row in rows:
    if row['name'] not in names:continue
    mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    other=chematic.from_smiles(row['smiles']).add_hydrogens()
    assert graph_matches(other,mol),row['name']
    props=AllChem.MMFFGetMoleculeProperties(mol)
    assert props is not None
    types=other.mmff94_numeric_atom_types(); charges=other.mmff94_charges_numeric()
    type_mismatches=[i for i,t in enumerate(types) if t!=props.GetMMFFAtomType(i)]
    max_charge_error=max(abs(q-props.GetMMFFPartialCharge(i)) for i,q in enumerate(charges))
    params=AllChem.ETKDGv3();params.randomSeed=42;params.numThreads=1
    assert AllChem.EmbedMolecule(mol,params)==0,row['name']
    coords=[list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    field=compare(mol,other,coords,'MMFF94')
    out.append({'name':row['name'],'type_mismatches':type_mismatches,'max_charge_error':max_charge_error,
                'energy_abs_error':field['energy_absolute_error'],
                'gradient_abs_error_max':max(field['gradient_absolute_errors'],default=0)})
assert len(out)==3
result={'rdkit':rdBase.rdkitVersion,'cases':out}
Path('mmff-terminal-carbon-diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
assert all(not x['type_mismatches'] and x['max_charge_error']<1e-8 and
           x['energy_abs_error']<1e-5 and x['gradient_abs_error_max']<1e-4 for x in out)
