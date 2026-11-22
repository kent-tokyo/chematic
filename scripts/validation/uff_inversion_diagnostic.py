#!/usr/bin/env python3
"""Compare exposed P-containing UFF development cases against RDKit and energy differences."""
import json
from pathlib import Path
import sys
import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path.cwd() / 'scripts/benchmarks'))
from platinum_development_probe import graph_matches

# These cases were exposed in the first COD holdout; they are development
# diagnostics and must never be described as independent test cases again.
cases = {
    'LQYQIJSGWHWKOT-UUWRZZSWSA-N': 'C=C1C=C(C)N(c2c(C(C)C)cccc2C(C)C)[Si@@](P)(Cl)N1c1c(C(C)C)cccc1C(C)C',
    'YRMADJRSLOCNRA-UHFFFAOYSA-N': 'CC(C)(C)c1cc(P(CCP(c2cc(C(C)(C)C)cc(C(C)(C)C)c2)c2cc(C(C)(C)C)cc(C(C)(C)C)c2)c2cc(C(C)(C)C)cc(C(C)(C)C)c2)cc(C(C)(C)C)c1',
    'VXFCZRMFSLOXNB-UHFFFAOYSA-N': 'CC(C)n1nc(-c2ccccc2P(c2ccccc2)c2ccccc2)n[n+](C(C)C)c1=O',
}
observations=[]
for name, smiles in cases.items():
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    other=chematic.from_smiles(smiles).add_hydrogens()
    assert graph_matches(other,mol), name
    params=AllChem.ETKDGv3(); params.randomSeed=42; params.numThreads=1
    assert AllChem.EmbedMolecule(mol,params)==0, name
    coords=[list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    reference=AllChem.UFFGetMoleculeForceField(mol)
    native=other.uff_force_field(coords)
    rd_grad=list(reference.CalcGrad())
    native_grad=[v for atom in native.gradient() for v in atom]
    energy_error=abs(reference.CalcEnergy()-native.energy())
    gradient_error=max(abs(a-b) for a,b in zip(rd_grad,native_grad))
    phosphorus=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==15]
    fd_errors=[]
    for atom in phosphorus:
        for axis in range(3):
            coords[atom][axis]+=1e-5; plus=native.set_positions_(coords).energy()
            coords[atom][axis]-=2e-5; minus=native.set_positions_(coords).energy()
            coords[atom][axis]+=1e-5
            fd_errors.append(abs((plus-minus)/2e-5-native_grad[3*atom+axis]))
    item={'name':name,'energy_abs_error':energy_error,'gradient_abs_error_max':gradient_error,
          'p_finite_difference_error_max':max(fd_errors,default=0.0)}
    observations.append(item)
out={'rdkit':rdBase.rdkitVersion,'cases':observations}
Path('uff-inversion-diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
assert all(v['energy_abs_error']<1e-6 and v['gradient_abs_error_max']<1e-5
           and v['p_finite_difference_error_max']<0.002 for v in observations)
