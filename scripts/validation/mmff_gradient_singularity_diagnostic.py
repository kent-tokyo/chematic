#!/usr/bin/env python3
"""Separate native and RDKit MMFF gradients from finite differences on an exposed input."""
import json
from pathlib import Path
from rdkit import Chem,rdBase,RDLogger
from rdkit.Chem import AllChem,rdMolTransforms
import chematic
RDLogger.DisableLog('rdApp.*')
name='PCCNTHRNUKLHFF-XDHOZWIPSA-N'
rows=json.loads(Path('validation/results/cod-first-field-outliers-development-2026-11-22.json').read_text())['results']
row=next(r for r in rows if r['name']==name)
mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
params=AllChem.ETKDGv3();params.randomSeed=42;params.numThreads=1
assert AllChem.EmbedMolecule(mol,params)==0
coords=[list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
flat=[v for p in coords for v in p]
other=chematic.from_smiles(row['smiles']).add_hydrogens()
props=AllChem.MMFFGetMoleculeProperties(mol)
reference=AllChem.MMFFGetMoleculeForceField(mol,props,nonBondedThresh=100.0)
native=other.mmff_force_field(coords)
native_grad=[v for p in native.gradient() for v in p];rdkit_grad=list(reference.CalcGrad())
axes=sorted({round(i*(len(flat)-1)/5) for i in range(6)})
def fd(field_energy,axis,step=1e-5):
    plus=flat[:];minus=flat[:];plus[axis]+=step;minus[axis]-=step
    return (field_energy(plus)-field_energy(minus))/(2*step)
native_fd=lambda x:other.mmff_force_field([x[i:i+3] for i in range(0,len(x),3)]).energy()
checks=[{'axis':i,'native_gradient':native_grad[i],'rdkit_gradient':rdkit_grad[i],
         'native_fd':fd(native_fd,i),'rdkit_fd':fd(reference.CalcEnergy,i)} for i in axes]
worst=max(checks,key=lambda x:abs(x['native_gradient']-x['native_fd']))['axis']
consistent_grad=[v for point in native.gradient_with_consistent_planes() for v in point]
step_sweep=[{'step':h,'native_fd':fd(native_fd,worst,h),
             'rdkit_fd':fd(reference.CalcEnergy,worst,h)}
            for h in (1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)]
methods={'bond':'SetMMFFBondTerm','angle':'SetMMFFAngleTerm','stretch_bend':'SetMMFFStretchBendTerm',
         'oop':'SetMMFFOopTerm','torsion':'SetMMFFTorsionTerm','vdw':'SetMMFFVdWTerm','ele':'SetMMFFEleTerm'}
terms=[]
for term in methods:
    p=AllChem.MMFFGetMoleculeProperties(mol)
    for name2,method in methods.items():getattr(p,method)(name2==term)
    rf=AllChem.MMFFGetMoleculeForceField(mol,p,nonBondedThresh=100.0)
    energy,grad=other._rdkit_mmff_terms(coords,terms=[term])
    def isolated(x):return other._rdkit_mmff_terms([x[i:i+3] for i in range(0,len(x),3)],terms=[term])[0]
    terms.append({'term':term,'native_gradient':grad[worst],'rdkit_gradient':rf.CalcGrad()[worst],
                  'native_fd':fd(isolated,worst),'rdkit_fd':fd(rf.CalcEnergy,worst)})
angles=[]
for a in mol.GetAtoms():
    n=[b.GetIdx() for b in a.GetNeighbors()]
    for i in range(len(n)):
        for j in range(i+1,len(n)):
            degree=rdMolTransforms.GetAngleDeg(mol.GetConformer(),n[i],a.GetIdx(),n[j])
            if degree>175:angles.append({'atoms':[n[i],a.GetIdx(),n[j]],'degrees':degree})
out={'name':name,'rdkit':rdBase.rdkitVersion,'worst_sampled_axis':worst,
     'consistent_plane_gradient':consistent_grad[worst],'step_sweep':step_sweep,
     'sampled_axes':checks,'terms':terms,'angles_over_175':angles}
Path('mmff-gradient-singularity-diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
