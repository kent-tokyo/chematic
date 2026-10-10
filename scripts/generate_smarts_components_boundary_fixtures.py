"""Freeze rooted SMARTS/CXSMARTS references for disconnected molecules."""
import json,itertools
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1'
sources=['CC.CCC.O','[Na+].[Cl-]','[Mg+2].[O-]S(=O)(=O)[O-]',
'[Ca+2].[O-]C(=O)C.[O-]C(=O)C','[13CH3:2][NH2:4].[Cl-]',
'F[C@H](Cl)Br.[K+]','F/C=C/F.C/C=C\\C','[CH3].[OH]','[O].[N]',
'[Fe+3].[Cl-].[Cl-].[Cl-]','[NH3]->[Cu+2](<-[NH3])<-[NH3].[Cl-].[Cl-]',
'[H][H].[He]','[2H]C([3H])(O)C.[Ne]','C1CCCCC1.c1ccccc1','',
'[NH4+].[O-]C(=O)C','[O-2].[Zn+2]']
rows=[];cx=[]
for source in sources:
 m=Chem.MolFromSmiles(source);assert m is not None,source
 for iso,root in itertools.product([False,True],[-1]+list(range(m.GetNumAtoms()))):
  rows.append(dict(smiles=source,isomeric=iso,root=None if root<0 else root,smarts=Chem.MolToSmarts(m,isomericSmiles=iso,rootedAtAtom=root)))
 cx.append(dict(smiles=source,cx_smarts=Chem.MolToCXSmarts(m)))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-smarts-components-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows,cx=cx),indent=2)+'\n')
print(len(rows),'rooted SMARTS;',len(cx),'CX SMARTS')
