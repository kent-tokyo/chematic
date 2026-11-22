#!/usr/bin/env python3
"""Diagnose typing and charges for exposed v8 force-field outliers."""
import json
from pathlib import Path
import sys
import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
sys.path.insert(0, str(Path.cwd() / 'scripts/benchmarks'))
from platinum_development_probe import graph_matches

cases = json.loads(Path('validation/results/cod-second-mmff-outliers-v8-development-2026-11-22.json').read_text())['results']
rows = []
for case in cases:
    row = {'name': case['name']}
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        native = chematic.from_smiles(case['smiles']).add_hydrogens()
        assert graph_matches(native, mol), 'atom or bond graph mismatch'
        props = AllChem.MMFFGetMoleculeProperties(mol)
        types = native.mmff94_numeric_atom_types()
        charges = native.mmff94_charges_numeric()
        assert len(types) == len(charges) == mol.GetNumAtoms()
        row['type_mismatches'] = [{'atom': i, 'element': mol.GetAtomWithIdx(i).GetSymbol(),
                                   'chematic': t, 'rdkit': props.GetMMFFAtomType(i)}
                                  for i,t in enumerate(types) if t != props.GetMMFFAtomType(i)]
        row['charge_mismatches'] = [{'atom': i, 'chematic': q, 'rdkit': props.GetMMFFPartialCharge(i)}
                                    for i,q in enumerate(charges) if abs(q-props.GetMMFFPartialCharge(i)) > 1e-8]
        row['formal_charge'] = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    except Exception as e:
        row['error'] = repr(e)
    rows.append(row)
output = {'rdkit': rdBase.rdkitVersion, 'classification': 'exposed development', 'results': rows}
Path('cod-second-mmff-types.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
print(json.dumps(output))
assert not any('error' in row for row in rows)
