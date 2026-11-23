#!/usr/bin/env python3
"""Inspect the exposed Prime embedded conformer under the frozen v3 wheel."""
import gzip
import json
from pathlib import Path
import sys
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import chematic

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_development_probe import stereo_violations
with gzip.open('validation/results/prime-linux-development-100-2026-11-22.json.gz', 'rt') as source:
    record = next(r for r in json.load(source)['results'] if r['name'] == '1YXQ_SWI')
reference = Chem.AddHs(Chem.MolFromSmiles(record['smiles']))
coords = chematic.from_smiles(record['smiles']).add_hydrogens().etkdg_embed_corrected(random_seed=1848)
probe = Chem.Mol(reference)
conf = Chem.Conformer(probe.GetNumAtoms())
for i, xyz in enumerate(coords):
    conf.SetAtomPosition(i, xyz)
probe.AddConformer(conf)
bonds = []
for bond in probe.GetBonds():
    if bond.GetStereo() not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
        continue
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    angles = []
    for end, other in [(a,b), (b,a)]:
        for atom in probe.GetAtomWithIdx(end).GetNeighbors():
            if atom.GetIdx() != other:
                angles.append({'end': end, 'other': other, 'sub': atom.GetIdx(),
                               'degrees': rdMolTransforms.GetAngleDeg(conf, atom.GetIdx(), end, other)})
    bonds.append({'bond': bond.GetIdx(), 'stereo': str(bond.GetStereo()), 'angles': angles})
out = {'name': record['name'], 'violations': stereo_violations(reference, [coords])[0], 'bonds': bonds, 'coordinates': coords}
Path('prime-linear-stereo-diagnostic.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k:v for k,v in out.items() if k != 'coordinates'}))
