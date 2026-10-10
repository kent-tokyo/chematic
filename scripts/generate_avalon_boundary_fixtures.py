"""Generate fixed-version Avalon fingerprint references for boundary chemistry."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Avalon import pyAvalonTools
assert rdBase.rdkitVersion=="2026.03.1",rdBase.rdkitVersion
root=Path(__file__).resolve().parents[1]
parents=json.loads((root/'validation/rdkit-2026.03.1-descriptor-boundary.json').read_text())['rows'][::2]
smiles=[row['smiles'] for row in parents if row['smiles']]+[
    'NCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCNCCN',
    'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
    '*CC*','*C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)C(*)*',
    'N(CCO)(CCO)CCO','P(OCCN)(OCCN)(OCCN)(OCCN)OCCN','N(N(NN)NN)N(NN)NN',
    'NCC(O)(CN)CC(O)(CN)CC(O)(CN)CC(O)(CN)CC(O)(CN)CCN',
    'C1CC2(C1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCC1',
]
rows=[]
for text in smiles:
    mol=Chem.MolFromSmiles(text);assert mol is not None,text
    for flags in [15761407, (1 << 24) - 1] + [1 << i for i in range(24)]:
        rows.append(dict(smiles=text,bit_flags=flags,on_bits=list(pyAvalonTools.GetAvalonFP(mol,nBits=512,bitFlags=flags).GetOnBits())))
path=root/'validation/rdkit-2026.03.1-avalon-boundary.json'
path.write_text('{\n  "rdkit_version": "'+rdBase.rdkitVersion+'",\n  "rows": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print(f'wrote {len(rows)} cases')
