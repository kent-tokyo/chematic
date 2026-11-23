#!/usr/bin/env python3
"""Inspect the exposed COD development case with a charged E/Z bond."""
import gzip
import json
from pathlib import Path
import sys

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolTransforms
import chematic

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_development_probe import stereo_violations

with gzip.open('validation/results/cod-linux-development-v2-2026-11-22.json.gz', 'rt') as source:
    record = next(row for row in json.load(source)['results'] if row['name'] == 'MMHMSBIEZRBUKF-AIGMYPEUSA-N')
reference = Chem.AddHs(Chem.MolFromSmiles(record['smiles']))
bond = reference.GetBondWithIdx(13)
assert str(bond.GetStereo()) == 'STEREOZ'
assert [bond.GetBeginAtom().GetFormalCharge(), bond.GetEndAtom().GetFormalCharge()] == [-1, -1]
indices = [bond.GetStereoAtoms()[0], bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetStereoAtoms()[1]]

def inspect(coordinates, status):
    molecule = Chem.Mol(reference)
    conf = Chem.Conformer(molecule.GetNumAtoms())
    for i, xyz in enumerate(coordinates):
        conf.SetAtomPosition(i, xyz)
    molecule.AddConformer(conf)
    angle = rdMolTransforms.GetDihedralDeg(conf, *indices)
    Chem.AssignStereochemistryFrom3D(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return {'status': status, 'dihedral_degrees': angle,
            'observed': str(molecule.GetBondWithIdx(13).GetStereo()),
            'violations': stereo_violations(reference, [coordinates])[0]}

other = chematic.from_smiles(record['smiles']).add_hydrogens()
candidates = other.etkdg_embed_multiple_confs_mmff_corrected(20, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400)
rdmol = Chem.Mol(reference)
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 1
ids = list(AllChem.EmbedMultipleConfs(rdmol, numConfs=20, params=params))
def coordinates(conf):
    return [list(conf.GetAtomPosition(i)) for i in range(rdmol.GetNumAtoms())]
raw = [coordinates(rdmol.GetConformer(i)) for i in ids]
optimized = AllChem.MMFFOptimizeMoleculeConfs(rdmol, numThreads=1, maxIters=200)
result = {'name': record['name'], 'rdkit': rdBase.rdkitVersion, 'indices': indices,
          'chematic': [inspect(c, s) for s, _, c in candidates],
          'rdkit_raw': [inspect(c, None) for c in raw],
          'rdkit_optimized': [inspect(coordinates(rdmol.GetConformer(i)), status) for i, (status, _) in zip(ids, optimized)],
          'chematic_first_coordinates': candidates[0][2] if candidates else None,
          'rdkit_first_raw_coordinates': raw[0] if raw else None,
          'rdkit_first_optimized_coordinates': coordinates(rdmol.GetConformer(ids[0])) if ids else None}
Path('cod-charged-stereo-diagnostic.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if not k.endswith('_coordinates')}))
