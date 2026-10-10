"""Freeze RDKit UFF energies and analytic gradients for boundary types."""
import json,math
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion=='2026.03.1',rdBase.rdkitVersion
sources=['C=C','C#N','N=N','O=O','C1CC1','C1CCC1','C1=CC=CC=C1','C1CC1C1CC1','C1CCC1C1CCC1','CSSC','COOC','P(F)(F)(F)(F)F','S(F)(F)(F)(F)(F)F','O=[Sb](O)O','O=[Bi](O)O','[Si](Cl)(Cl)(Cl)Cl','[Ge](Cl)(Cl)(Cl)Cl','[Al](Cl)(Cl)Cl','[Cu]Cl','[Ag]Cl','[Be](F)F','[Hg](Cl)Cl','[Re](F)(F)(F)(F)(F)F','N[Pt](N)(Cl)Cl','F[Co](F)(F)(F)(F)F','F[Ni](F)(F)(F)(F)F','[Na+].[Cl-]','CC.O','[H][H]','[2H]O[3H]']
rows=[]
for source in sources:
    parent=Chem.MolFromSmiles(source);assert parent is not None,source
    mol=Chem.AddHs(parent)
    coords=[[i*1.25,0.8*(i%2),math.sin(i*0.7)*0.35] for i in range(mol.GetNumAtoms())]
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,p in enumerate(coords):conf.SetAtomPosition(i,p)
    mol.AddConformer(conf)
    for threshold,ignore in [(100.0,True),(100.0,False),(0.5,True)]:
        field=AllChem.UFFGetMoleculeForceField(mol,vdwThresh=threshold,ignoreInterfragInteractions=ignore)
        if field is None:continue
        rows.append(dict(smiles=source,coords=coords,threshold=threshold,ignore=ignore,energy=field.CalcEnergy(),gradient=list(field.CalcGrad()),all_params=AllChem.UFFHasAllMoleculeParams(mol)))
for number in range(1,104):
    mol=Chem.MolFromSmiles('['+Chem.GetPeriodicTable().GetElementSymbol(number)+']')
    if mol is not None:
        rows.append(dict(smiles=Chem.MolToSmiles(mol),coords=[[0.,0.,0.]],threshold=100.,ignore=True,energy=0.,gradient=[0.,0.,0.],all_params=AllChem.UFFHasAllMoleculeParams(mol)))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-uff-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print('wrote',len(rows),'cases')
