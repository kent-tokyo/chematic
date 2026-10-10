"""Freeze charge/bond edits of native sanitization cleanup, before perception."""
import json
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1'
rows=[]
sources=['CN(=O)=O','CN(=O)(C)C','CP(=O)=C','OP(=O)=C','O=Cl(=O)=O','O=Br(=O)=O','O=I(=O)=O','O=Cl(=O)(=O)O','CN(C)(C)[Cu]','CN(C)(C)[Fe]','CO(C)[Zn]','C[N+](=O)[O-]','CCO']
for source in sources:
    raw=Chem.MolFromSmiles(source,sanitize=False);assert raw is not None,source
    raw.UpdatePropertyCache(strict=False)
    before=Chem.Mol(raw)
    Chem.Cleanup(raw);Chem.CleanupOrganometallics(raw)
    charges=[[i,a.GetFormalCharge()] for i,a in enumerate(raw.GetAtoms()) if a.GetFormalCharge()!=before.GetAtomWithIdx(i).GetFormalCharge()]
    bonds=[]
    for b in raw.GetBonds():
        if b.GetBondType()!=before.GetBondWithIdx(b.GetIdx()).GetBondType():
            bonds.append([b.GetIdx(),str(b.GetBondType()),b.GetBeginAtomIdx() if b.GetBondType()==Chem.BondType.DATIVE else None])
    rows.append(dict(smiles=source,charges=charges,bonds=bonds))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-cleanup-edits-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'cleanup inputs',sum(bool(r['charges'] or r['bonds']) for r in rows),'changed')
