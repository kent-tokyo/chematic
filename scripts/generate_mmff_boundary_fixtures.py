"""Generate RDKit 2026.03.1 MMFF94 atom type/charge references."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
root=Path(__file__).resolve().parents[1]
parents=json.loads((root/'validation/rdkit-2026.03.1-descriptor-boundary.json').read_text())['rows'][::2]
smiles=[row['smiles'] for row in parents if row['smiles']]+[
    'C=S=O','S=P(S)(S)S','C[NH+]=C','C=N[O-]','c1cn[n-]n1','c1nn[nH]n1',
    'CC(=O)NC(=O)C','CO[O-]','N[O-]','CSS[O-]','[O-][S+](C)C','[OH3+]','[SH-]',
]
rows=[]
for text in smiles:
    mol=Chem.MolFromSmiles(text);assert mol is not None,text
    # The public Rust typer documents an organic atom-typing scope.
    # Free metal and halide ions need separate typing support.
    if any(a.GetAtomicNum() not in {1,6,7,8,9,14,15,16,17,35,53} or
           (a.GetDegree()==0 and a.GetAtomicNum() in {17,35} and a.GetFormalCharge()<0)
           for a in mol.GetAtoms()):
        continue
    mol=Chem.AddHs(mol)
    props=AllChem.MMFFGetMoleculeProperties(mol)
    if props is not None:
        rows.append(dict(smiles=text,types=[props.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())],charges=[props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms())]))
path=root/'validation/rdkit-2026.03.1-mmff-boundary.json'
path.write_text('{\n  "rdkit_version": "'+rdBase.rdkitVersion+'",\n  "rows": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print(f'wrote {len(rows)} supported cases')
