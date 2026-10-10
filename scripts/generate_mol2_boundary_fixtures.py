"""Generate RDKit 2026.03.1 references from explicit-H Tripos records."""
import json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import rdDepictor
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
smiles=['N=C(N)N','CN=C(N)N','NC(N)=NC','CN(C)C(N)=N','NC(=N)NC(=N)N',
        'N=C1NCCN1','CC(=O)[O-]','O=C(O)C','C[N+](=O)[O-]','[O-][N+](=O)c1ccccc1',
        'F/C=C/F','F/C=C\\F','CC/C=C(C)/C','Cl/C(F)=C(Br)/I',
        'C/C=C/C=C/C','C/C=C\\C=C/C','C/C=C(C)/C=C/C','CC/C(Cl)=C(Br)/CC',
        'C1=CCCCC1','CC(=O)NC','c1cc[nH]c1', 'NC=N', 'N(C)=CN(C)', 'N=C1NCC1', 'N=C1NCCC1', '[O-][N+](C)C=N', 'N=C(N)O', 'N=C(N)S', 'N=C1NC=CC=C1']
def sybyl(a):
    z=a.GetSymbol()
    if a.GetIsAromatic():
        return 'N.pl3' if z=='N' and any(n.GetAtomicNum()==1 for n in a.GetNeighbors()) else z+'.ar'
    if z in ['C','N','O','S','P']:
        h=str(a.GetHybridization());return z+'.'+{'SP':'1','SP2':'2','SP3':'3'}.get(h,'3')
    return z
rows=[]
for text in smiles:
    mol=Chem.MolFromSmiles(text);assert mol is not None,text
    mol=Chem.AddHs(mol);rdDepictor.Compute2DCoords(mol)
    types=[sybyl(a) for a in mol.GetAtoms()]
    for a in mol.GetAtoms():
        if a.GetAtomicNum()==6 and sum(n.GetAtomicNum()==7 for n in a.GetNeighbors())>=2:
            types[a.GetIdx()]='C.cat'
        if a.GetAtomicNum()==8 and any(n.GetAtomicNum()==6 and sum(nn.GetAtomicNum()==8 for nn in n.GetNeighbors())==2 and any(nn.GetFormalCharge()<0 for nn in n.GetNeighbors()) for n in a.GetNeighbors()):
            types[a.GetIdx()]='O.co2'
    for reverse in [False,True]:
        conf=mol.GetConformer();atoms=[]
        for i,a in enumerate(mol.GetAtoms()):
            p=conf.GetAtomPosition(i);atoms.append(f'{i+1} {a.GetSymbol()}{i+1} {p.x:.8f} {p.y:.8f} {p.z:.8f} {types[i]} 1 LIG 0')
        bonds=[]
        for i,b in enumerate(mol.GetBonds()):
            a,c=b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1
            if reverse:a,c=c,a
            bt='ar' if b.GetIsAromatic() else str(int(b.GetBondTypeAsDouble()))
            bonds.append(f'{i+1} {a} {c} {bt}')
        block=f'@<TRIPOS>MOLECULE\nfixture\n{len(atoms)} {len(bonds)} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n'+'\n'.join(atoms)+'\n@<TRIPOS>BOND\n'+'\n'.join(bonds)+'\n'
        rebuilt=Chem.MolFromMol2Block(block)
        assert rebuilt is not None,text
        variants=[]
        for sanitize in [True,False]:
            for remove_hs in [True,False]:
                for cleanup in [True,False]:
                    variant=Chem.MolFromMol2Block(block,sanitize=sanitize,removeHs=remove_hs,cleanupSubstructures=cleanup)
                    expected=dict(sanitize=sanitize,remove_hs=remove_hs,cleanup=cleanup,accepted=variant is not None)
                    if variant is not None:
                        expected.update(smiles=Chem.MolToSmiles(variant) if sanitize else '',
                            atom_count=variant.GetNumAtoms(),bond_count=variant.GetNumBonds(),
                            atomic_numbers=[a.GetAtomicNum() for a in variant.GetAtoms()],
                            charges=[a.GetFormalCharge() for a in variant.GetAtoms()])
                    variants.append(expected)
        rows.append(dict(source_smiles=text,block=block,smiles=Chem.MolToSmiles(rebuilt),variants=variants))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-mol2-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),separators=(',',':'))+'\n')
print(f'wrote {len(rows)} cases')
