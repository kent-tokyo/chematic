"""Freeze canonical reading/writing of native-generated branched polyenes."""
import itertools,json,random
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1'
rows=[];randomizer=random.Random(9614)
for source in ['FC=CC(C=CCl)=C(C=CBr)C=CI',
               'FC=CC(C=CCl)=C1C=CC=CC=CC=C1',
               'FC=C1C=CC=CC(=CC=CCl)C=CC=C1',
               'CC=C(C=CF)C=C(C=CCl)C=CBr']:
 base=Chem.MolFromSmiles(source);assert base is not None,source
 indices=[]
 for b in base.GetBonds():
  if b.GetBondType()!=Chem.BondType.DOUBLE:continue
  a,z=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
  if base.GetAtomWithIdx(a).GetDegree()>1 and base.GetAtomWithIdx(z).GetDegree()>1:indices.append(b.GetIdx())
 for pattern in itertools.product([False,True],repeat=len(indices)):
  m=Chem.Mol(base)
  for idx,e in zip(indices,pattern):
   b=m.GetBondWithIdx(idx);a,z=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
   left=next(n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx()!=z)
   right=next(n.GetIdx() for n in m.GetAtomWithIdx(z).GetNeighbors() if n.GetIdx()!=a)
   b.SetStereoAtoms(left,right);b.SetStereo(Chem.BondStereo.STEREOE if e else Chem.BondStereo.STEREOZ)
  canonical=Chem.MolToSmiles(m)
  # The text is the actual input domain: reparse and freeze that molecule.
  m=Chem.MolFromSmiles(canonical);assert m is not None
  canonical=Chem.MolToSmiles(m)
  spellings={canonical}
  for _ in range(32):
   order=list(range(m.GetNumAtoms()));randomizer.shuffle(order)
   spellings.add(Chem.MolToSmiles(Chem.RenumberAtoms(m,order),canonical=False))
  for text in sorted(spellings):
   parsed=Chem.MolFromSmiles(text);assert parsed is not None,text
   # Freeze what the actual serialized input means to the native reader,
   # including any ambiguity cleanup. Do not assume assigned labels survive.
   expected=Chem.MolToSmiles(parsed)
   rows.append(dict(smiles=text,canonical=expected))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-branched-directions-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'branched direction spellings')
