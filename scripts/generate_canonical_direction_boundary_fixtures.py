"""Freeze canonical directional SMILES respellings with RDKit 2026.03.1."""
import json, random
from pathlib import Path
from rdkit import Chem, rdBase
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
sources = ['C/C=C/C=C/C', 'C/C=C\\C=C/C', 'F/C=C(/Cl)C(/Br)=C/I',
    'F/C=C(\\Cl)C(/Br)=C/I', 'F/C=C(/Cl)C(\\Br)=C/I',
    'C/C=C(/F)C=C/C', 'C/C=C(\\F)C=C/C', 'C/C=C/C=C/C=C/C',
    'C/C=C/C=C\\C=C/C', 'C/C=C\\C=C/C=C/C', 'C/C=C1/CCCCCCC1',
    'F/C=C1/CCCCCCC1', 'C1CCC/C=C/CCCC1', 'C1CCC/C=C\\CCCC1',
    'F/C=C(/Cl)C=C(/Br)I', 'Cl/C(F)=C(Br)/I', 'Cl/C(F)=C(Br)\\I',
    'N/C=C/C=N/O', 'N/C=C\\C=N/O', '*1=CC=CC=C1', '*1=CC=CC=C1*',
    'C1=*C=CC=C1', '*1=CC=C*1', 'O=C1C=CC(=O)C=C1',
    'C/C=N/c1ncc[nH]1', 'C/C=N\\c1ncc[nH]1',
    'F[C@H]1CC[C@@H](Cl)CC1', 'C[C@H]1CC[C@H](C)CC1']
randomizer = random.Random(940)
rows=[]
for source in sources:
    mol=Chem.MolFromSmiles(source)
    assert mol is not None, source
    canonical=Chem.MolToSmiles(mol)
    spellings={source}
    for _ in range(64):
        order=list(range(mol.GetNumAtoms()));randomizer.shuffle(order)
        spellings.add(Chem.MolToSmiles(Chem.RenumberAtoms(mol,order),canonical=False))
    rows.extend(dict(smiles=text,canonical=canonical) for text in sorted(spellings))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-canonical-direction-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print('wrote',len(rows),'cases')
