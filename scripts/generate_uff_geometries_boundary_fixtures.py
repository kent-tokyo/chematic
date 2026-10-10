"""Freeze UFF small-ring and degenerate-coordinate reference calculations."""
import json,math
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion=='2026.03.1'
def encode(x):
    if math.isnan(x):return 'nan'
    if math.isinf(x):return '+inf' if x>0 else '-inf'
    return x
sources=['[H][H]','Cl','C.C','C1=CC1','C1=CCC1','C=C1CC1','C=C1CCC1',
         'C1=CC2CC12','C1=CC2CCC12','O=C1CC1','O=C1CCC1','N=C1CC1','C#N','CC=O','C=C']
rows=[]
for source in sources:
    mol=Chem.AddHs(Chem.MolFromSmiles(source));assert mol is not None
    assert AllChem.UFFHasAllMoleculeParams(mol),source
    n=mol.GetNumAtoms()
    for mode in ['distorted','collinear','coincident']:
        coords=[]
        for i in range(n):
            coords.append([math.cos(i*1.23)*(1+i*.17),math.sin(i*1.37)*(1+i*.11),i*.29] if mode=='distorted' else [i*1.1,0.,0.] if mode=='collinear' else [0.,0.,0.])
        c=Chem.Conformer(n)
        for i,xyz in enumerate(coords):c.SetAtomPosition(i,xyz)
        mol.RemoveAllConformers();mol.AddConformer(c)
        field=AllChem.UFFGetMoleculeForceField(mol,vdwThresh=100.,ignoreInterfragInteractions=False)
        assert field is not None,source
        rows.append(dict(smiles=source,mode=mode,coords=coords,energy=encode(field.CalcEnergy()),gradient=[encode(x) for x in field.CalcGrad()]))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-uff-geometries-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'UFF geometries',sum(r['energy'] in ['nan','+inf','-inf'] for r in rows),'nonfinite energies')
