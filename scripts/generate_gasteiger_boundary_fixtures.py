"""Freeze default 12-step PEOE charges, including explicit parameter fallbacks."""
import json,math
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import rdPartialCharges
assert rdBase.rdkitVersion=='2026.03.1'
sources=['CCO','CC(=O)[O-]','C[NH+]=C(N)N','C[N+](C)(C)C','[NH4+]',
         '[Cl-]','[Na+]','[Fe+2].CC[O-]','*CC','C[Si](C)(C)C','C[SiH]=C',
         'C[Si]#C','[SiH4]','C[Be]C','C=[Be]','[BeH2]',
         'C[Mg]C','C=[Mg]','C#[Mg]','C[Al](C)C','C[Al]=C','[AlH3]',
         'CP(C)C','CP=C','C=S','CS','CS(=O)C','CS(=O)(=O)C',
         'C[N-][N+]#N','B(F)(F)F','[BH4-]','C#N','CN=O',
         'F[C@H](Cl)Br','CI','c1ccncc1','c1ncc[nH]1','[OH-]','[CH3]']
rows=[]
def encode(value):
    if math.isnan(value):return 'nan'
    if math.isinf(value):return '+inf' if value>0 else '-inf'
    return value
for source in sources:
    base=Chem.MolFromSmiles(source);assert base is not None,source
    for explicit in [False,True]:
        mol=Chem.AddHs(base) if explicit else Chem.Mol(base)
        rdPartialCharges.ComputeGasteigerCharges(mol,nIter=12)
        charges=[encode(float(a.GetProp('_GasteigerCharge'))) for a in mol.GetAtoms()]
        rows.append(dict(smiles=source,add_hs=explicit,charges=charges,
                         formal_charge=Chem.GetFormalCharge(mol)))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-gasteiger-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'charge references')
