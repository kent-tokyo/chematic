"""Freeze RDKit 2026.03.1 PDB stereo and malformed-record outcomes."""
import itertools, json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
blocks=[]
for center, ligands, coords in [
    (78,[9,17,35,53],[(2,0,0),(-2,0,0),(0,2,0),(0,-2,0)]),
    (15,[9,17,35,53,7],[(0,0,2),(0,0,-2),(2,0,0),(-1,1.732,0),(-1,-1.732,0)]),
    (27,[9,17,35,53,7,8],[(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2)]),
]:
    for order in itertools.permutations(range(len(ligands))):
        rw=Chem.RWMol(); rw.AddAtom(Chem.Atom(center))
        for number in ligands:
            j=rw.AddAtom(Chem.Atom(number));rw.AddBond(0,j,Chem.BondType.SINGLE)
        mol=rw.GetMol(); conf=Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
        for j,k in enumerate(order,1):conf.SetAtomPosition(j,coords[k])
        mol.AddConformer(conf)
        blocks.append((f'geometry-{center}-{order}',Chem.MolToPDBBlock(mol)))
for source in ['F[C@H](Cl)Br','F/C=C/Cl','N[C@@H](C)C(=O)O','[2H]O[3H]','c1ccncc1','C1CCCCC1','C#N','O=S(=O)(O)O']:
    mol=Chem.AddHs(Chem.MolFromSmiles(source));assert AllChem.EmbedMolecule(mol,randomSeed=12)==0
    blocks.append((source,Chem.MolToPDBBlock(mol)))
base=blocks[-1][1]
# Standard biological residue names restore bond orders when CONECT records
# carry connectivity only. Keep the full amino acid/base topology, including
# side chains and optional terminal phosphate caps, as the external oracle.
for sequence, flavor in [('ARNDCEQGHILKMFPSTWYV',0), ('ARNDCEQGHILKMFPSTWYV',1),
                         ('ACGU',2), ('ACGU',5), ('ACGT',6), ('ACGT',9)]:
    mol=Chem.MolFromFASTA(sequence,flavor=flavor)
    assert mol is not None,(sequence,flavor)
    connectivity=Chem.RWMol(mol)
    for bond in connectivity.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in connectivity.GetAtoms():
        atom.SetIsAromatic(False)
    blocks.append((f'biopolymer-{sequence}-{flavor}',Chem.MolToPDBBlock(connectivity)))
for off,n,value in [(6,5,'abcde'),(30,8,' 1.2e+03'),(38,8,'badcoord'),(22,4,'abcd'),(76,2,'Xx'),(78,2,'x+')]:
    lines=base.splitlines();line=lines[0].ljust(80);lines[0]=line[:off]+value.ljust(n)+line[off+n:]
    blocks.append((f'invalid-{off}','\n'.join(lines)+'\n'))
rows=[]
for label,block in blocks:
    for sanitize,remove_hs,flavor,proximity in ([(True,True,0,False)] if label.startswith('geometry-') else [(True,True,0,False),(True,False,0,False),(False,False,1,False),(True,True,8,False)]):
        mol=Chem.MolFromPDBBlock(block,sanitize=sanitize,removeHs=remove_hs,flavor=flavor,proximityBonding=proximity)
        row=dict(label=label,block=block,sanitize=sanitize,remove_hs=remove_hs,flavor=flavor,proximity=proximity,valid=mol is not None)
        if mol is not None:
            row.update(atoms=mol.GetNumAtoms(),bonds=mol.GetNumBonds(),smiles=Chem.MolToSmiles(mol) if sanitize else '')
        rows.append(row)
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-pdb-boundary.json').write_text('{\n  \"rdkit_version\": \"'+rdBase.rdkitVersion+'\",\n  \"rows\": [\n'+',\n'.join('    '+json.dumps(row) for row in rows)+'\n  ]\n}\n')
print('wrote',len(rows),'cases')
