"""Freeze molecular hashes for charged, radical and unusual bond-order inputs."""
import json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import rdMolHash
assert rdBase.rdkitVersion=='2026.03.1'
sources=['','[Na+]','[Ca+2]','[Fe+3]','[O-2]','[P-2]','[NH4+]','[NH2+]','[NH+]','[N+2]',
         '[S]','[SH]','[SH2]','[SH3]','[S+2]','[P]','[PH]','[BH3]','[BH2]','[BH]','B(F)(F)F',
         'C#N','C#C','C$C','[Cr]$[Cr]','[H][H]','[2H]O[3H]','[13CH3][NH3+]',
         '[Na+].[Cl-]','[Mg+2].[O-]C(=O)C.[O-]C(=O)C','N=[N+]=[N-]',
         'O=NC=NO','NC(=N)C(=N)N','N=C1NC=NC(=O)N1','OC=C(O)C=O',
         'CC(=O)NC(=O)C','NNC(=O)NN','O=CC=CC=O','O=C(O)C=C(O)C=O','N=CC=CN',
         'S=C(N)N','O=C(N)N','[N-]=C=[N-]','C=[N+]=C','c1ncc[nH]1','c1cc[n+]([O-])cc1']
rows=[]
for source in sources:
 m=Chem.MolFromSmiles(source);assert m is not None,source
 for name in rdMolHash.HashFunction.names:
  # Default non-CX hashes avoid a native 2026.03.1 internal output-order
  # KeyError when a scaffold hash produces an empty graph.
  for cx in [False]:
   rows.append(dict(smiles=source,function=name,use_cx=cx,expected=rdMolHash.MolHash(Chem.Mol(m),rdMolHash.HashFunction.names[name],cx)))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-hash-ions-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'molecular hash references')
