"""Freeze MMFF94 nonbonded energies and collision-guard gradients."""
import json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
assert rdBase.rdkitVersion=='2026.03.1'
rows=[]
for source in ['[Cl-].[Cl-]','[Cl-].[Na+]','[Br-].[Br-]'] + [x+'.[Cl-]' for x in ['[F-]','[Li+]','[K+]','[Mg+2]','[Ca+2]','[Zn+2]','[Fe+2]','[Fe+3]','[Cu+]','[Cu+2]']]:
    mol=Chem.MolFromSmiles(source);assert mol is not None
    props=AllChem.MMFFGetMoleculeProperties(mol);assert props is not None
    for distance in [0.,1e-8,1.,4.]:
        conformer=Chem.Conformer(2)
        conformer.SetAtomPosition(0,(0.,0.,0.))
        conformer.SetAtomPosition(1,(distance,0.,0.))
        mol.RemoveAllConformers();mol.AddConformer(conformer)
        for ignore in [False,True]:
            field=AllChem.MMFFGetMoleculeForceField(mol,props,nonBondedThresh=100.,ignoreInterfragInteractions=ignore)
            assert field is not None
            rows.append(dict(smiles=source,distance=distance,ignore_interfrag=ignore,
                             energy=field.CalcEnergy(),gradient=list(field.CalcGrad()),
                             atom_types=[props.GetMMFFAtomType(i) for i in range(2)],
                             charges=[props.GetMMFFPartialCharge(i) for i in range(2)]))
root=Path(__file__).resolve().parents[1]
(root/'validation/rdkit-2026.03.1-mmff-nonbonded-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'nonbonded references')
