"""Freeze Tripos phosphate/quaternary-N cleanup and explicit invalid type combinations."""
import itertools,json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import rdDepictor
assert rdBase.rdkitVersion=='2026.03.1'
cases=[('CP(=O)(O)O','phosphate'),('OP(=O)(O)O','phosphate'),
       ('[O-]P(=O)([O-])[O-]','phosphate'),('COP(=O)([O-])OC','phosphate'),
       ('[NH4+]','quat'),('C[N+](C)(C)C','quat'),
       ('CCOCC','bad-oxygen-degree'),('CO','bad-oxygen-neighbor'),('CC','bad-cation-neighbors')]
rows=[]
for source,mode in cases:
 mol=Chem.AddHs(Chem.MolFromSmiles(source));assert mol is not None
 rdDepictor.Compute2DCoords(mol)
 types=[]
 for a in mol.GetAtoms():
  symbol=a.GetSymbol();hybrid=str(a.GetHybridization())
  kind=symbol+'.'+{'SP':'1','SP2':'2','SP3':'3'}.get(hybrid,'3') if symbol in ['C','N','O','P','S'] else symbol
  if mode=='phosphate' and symbol=='P':kind='P.3'
  if mode=='phosphate' and symbol=='O' and a.GetDegree()==1:kind='O.co2'
  if mode=='quat' and symbol=='N':kind='N.4'
  if mode.startswith('bad-oxygen') and symbol=='O':kind='O.co2'
  if mode=='bad-cation-neighbors' and a.GetIdx()==0:kind='C.cat'
  types.append(kind)
 for reverse,cleanup in itertools.product([False,True],repeat=2):
  atoms=[];conf=mol.GetConformer()
  for i,a in enumerate(mol.GetAtoms()):
   p=conf.GetAtomPosition(i);atoms.append(f'{i+1} {a.GetSymbol()}{i+1} {p.x:.8f} {p.y:.8f} {p.z:.8f} {types[i]} 1 LIG 0')
  bonds=[]
  for i,b in enumerate(mol.GetBonds()):
   a,z=b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1
   if reverse:a,z=z,a
   bonds.append(f'{i+1} {a} {z} {int(b.GetBondTypeAsDouble())}')
  block=f'@<TRIPOS>MOLECULE\nfixture\n{len(atoms)} {len(bonds)} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n'+'\n'.join(atoms)+'\n@<TRIPOS>BOND\n'+'\n'.join(bonds)+'\n'
  rebuilt=Chem.MolFromMol2Block(block,cleanupSubstructures=cleanup)
  rows.append(dict(smiles=source,mode=mode,block=block,cleanup=cleanup,expected=Chem.MolToSmiles(rebuilt) if rebuilt is not None else None))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-mol2-ionic-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'Tripos inputs',sum(r['expected'] is None for r in rows),'native rejections')
