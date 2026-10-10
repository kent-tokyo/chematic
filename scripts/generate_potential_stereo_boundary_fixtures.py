"""Freeze tetrahedral potential-center atom sets with RDKit 2026.03.1."""
import itertools,json
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
sources=[]
for left,right in itertools.product(['@','@@'],repeat=2):
    sources += [f'C[C{left}H](F)C(O)[C{right}H](C)F',
                f'O[C@H](C[C{left}H](F)Cl)C[C{right}H](F)Cl',
                f'O[C@H](C(=O)[C{left}H](F)Cl)C(=O)[C{right}H](F)Cl',
                f'O[C@H](C1C[C{left}H](F)CCC1)C1C[C{right}H](F)CCC1']
sources += ['CCC(F)Cl','CC(F)(Cl)Br','CC(C)C','C1CCCCC1', 'C[C@H]1CC[C@@H](C)CC1',
            'C[C@H]1CC[C@H](C)CC1', 'C[P@](F)Cl', 'C[As@](F)Cl', 'C[S@](=O)CC',
            'C[Se@](=O)CC','C[N@]1CC1','N1CC2CCC1C2','CC(=O)N(C)CC']
rows=[]
for source in sources:
    mol=Chem.MolFromSmiles(source);assert mol is not None,source
    centers=[s.centeredOn for s in Chem.FindPotentialStereo(mol) if str(s.type)=='Atom_Tetrahedral']
    rows.append(dict(smiles=source,centers=sorted(centers)))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-potential-stereo-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
