"""Freeze canonical directional SMILES respellings with RDKit 2026.03.1."""
import json, random
from pathlib import Path
from rdkit import Chem, rdBase
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
sources = ['C/C=C/C=C/C', 'C/C=C\\C=C/C', 'F/C=C(/Cl)C(/Br)=C/I',
    'F/C=C(\\Cl)C(/Br)=C/I', 'F/C=C(/Cl)C(\\Br)=C/I',
    'C/C=C(/F)C=C/C', 'C/C=C(\\F)C=C/C', 'C/C=C/C=C/C=C/C',
    'C/C=C/C=C\\C=C/C', 'C/C=C\\C=C/C=C/C', 'C/C=C1/CCCCCCC1',
    'F/C=C1/CCCCCCC1', 'C1CCC/C=C/CCCC1', 'C1CCC/C=C\\CCCC1',
    'F/C=C(/Cl)C=C(/Br)I', 'Cl/C(F)=C(Br)/I', 'Cl/C(F)=C(Br)\\I',
    'N/C=C/C=N/O', 'N/C=C\\C=N/O', '*1=CC=CC=C1', '*1=CC=CC=C1*',
    'C1=*C=CC=C1', '*1=CC=C*1', 'O=C1C=CC(=O)C=C1',
    'C/C=N/c1ncc[nH]1', 'C/C=N\\c1ncc[nH]1',
    'F[C@H]1CC[C@@H](Cl)CC1', 'C[C@H]1CC[C@H](C)CC1']
randomizer = random.Random(940)
# Multiple directional double bonds in a closed conjugated macrocycle share
# carrier bonds on both sides. These exercise direction reconciliation when
# the traversal reaches the last double bond from already marked neighbors.
for size in [12, 16]:
    base = Chem.MolFromSmiles('C1=' + 'CC=' * (size // 2 - 1) + 'C1')
    assert base is not None and base.GetNumAtoms() == size
    doubles = [b.GetIdx() for b in base.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
    pattern_randomizer = random.Random(960 + size)
    patterns = [0, (1 << len(doubles)) - 1] + [pattern_randomizer.randrange(1 << len(doubles)) for _ in range(14)]
    for pattern in patterns:
        mol = Chem.Mol(base)
        for i, index in enumerate(doubles):
            bond = mol.GetBondWithIdx(index)
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            left = next(n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != b)
            right = next(n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx() != a)
            bond.SetStereoAtoms(left, right)
            bond.SetStereo(Chem.BondStereo.STEREOE if pattern & (1 << i) else Chem.BondStereo.STEREOZ)
        sources.append(Chem.MolToSmiles(mol))
rows=[]
for source in sources:
    mol=Chem.MolFromSmiles(source)
    assert mol is not None, source
    canonical=Chem.MolToSmiles(mol)
    spellings={source}
    for _ in range(64):
        order=list(range(mol.GetNumAtoms()));randomizer.shuffle(order)
        spellings.add(Chem.MolToSmiles(Chem.RenumberAtoms(mol,order),canonical=False))
    rows.extend(dict(smiles=text,canonical=canonical) for text in sorted(spellings))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-canonical-direction-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print('wrote',len(rows),'cases')
