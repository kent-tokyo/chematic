#!/usr/bin/env python3
"""Reproduce an exposed Prime development-set double-bond stereo mismatch."""
import gzip
import json
from pathlib import Path
import sys

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolTransforms
import chematic

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_development_probe import stereo_violations

with gzip.open('validation/results/prime-linux-development-100-2026-11-22.json.gz', 'rt') as source:
    rows = json.load(source)['results']
record = next(row for row in rows if row['name'] == '1YXQ_SWI')
molecule = Chem.AddHs(Chem.MolFromSmiles(record['smiles']))
other = chematic.from_smiles(record['smiles']).add_hydrogens()
candidates = other.etkdg_embed_multiple_confs_mmff_corrected(44, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400)
observations = []
for index, (status, energy, coordinates) in enumerate(candidates):
    violations = stereo_violations(molecule, [coordinates])[0]
    if not (violations['atoms'] or violations['bonds']):
        continue
    probe = Chem.Mol(molecule)
    conf = Chem.Conformer(probe.GetNumAtoms())
    for atom_index, xyz in enumerate(coordinates):
        conf.SetAtomPosition(atom_index, xyz)
    probe.AddConformer(conf)
    torsions = []
    for bond_index in violations['bonds']:
        bond = probe.GetBondWithIdx(bond_index)
        neighbours = list(bond.GetStereoAtoms())
        atom_indices = [neighbours[0], bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), neighbours[1]]
        torsions.append({'bond': bond_index, 'expected': str(bond.GetStereo()), 'atoms': atom_indices,
                         'dihedral_degrees': rdMolTransforms.GetDihedralDeg(probe.GetConformer(), *atom_indices)})
    observations.append({'index': index, 'status': status, 'energy': energy, 'violations': violations, 'torsions': torsions})
result = {'name': record['name'], 'rdkit': rdBase.rdkitVersion, 'count': len(candidates), 'observations': observations}
Path('prime-stereo-diagnostic.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
