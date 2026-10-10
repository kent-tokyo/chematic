"""Generate RDKit 2026.03.1 distance-bound references for rings and amides."""
import itertools
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import rdDistGeom
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
smiles=['CC#CC','CC(=O)NC','CC(=O)N(C)C','CC(=O)OC','NCCN','OCCO',
        'C1CC1','C1CCC1','C1=CCC1','C1=CCCC1','C1CCCC1','C1CCCCCCC1',
        'C1CCCCCCCCCCC1','O=C1NCCCCCCCCCCC1','O=C1OCCCCCCCCCCC1',
        'O=C1N(C)CCCCCCCCCCC1','O=C1NCCNC(=O)NCCNCC1','C1CC2CCC1C2',
        'F/C=C/F','F/C=C\\F','CC=C=C(C)C','c1ccc2ccccc2c1',
        'P(F)(F)(F)(F)F','S(F)(F)(F)(F)(F)F']
rows=[]
for text in smiles:
    mol=Chem.MolFromSmiles(text);assert mol is not None,text
    for h in [False,True] if text in ['CC(=O)NC','NCCN','OCCO'] else [False]:
        source=Chem.AddHs(mol) if h else mol
        for set15,smooth,macro in itertools.product([False,True],repeat=3):
            bounds=rdDistGeom.GetMoleculeBoundsMatrix(source,set15bounds=set15,doTriangleSmoothing=smooth,useMacrocycle14config=macro)
            rows.append(dict(smiles=text,add_hs=h,set15=set15,smooth=smooth,macrocycle14=macro,bounds=bounds.tolist()))
root=Path(__file__).resolve().parents[1]
path=root/'validation/rdkit-2026.03.1-bounds-boundary.json'
path.write_text('{\n  "rdkit_version": "'+rdBase.rdkitVersion+'",\n  "rows": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print(f'wrote {len(rows)} cases')
