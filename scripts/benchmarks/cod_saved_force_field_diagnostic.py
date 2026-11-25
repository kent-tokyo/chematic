#!/usr/bin/env python3
"""Full-axis exposed-coordinate force-field finite-difference diagnostic."""
import argparse
import gzip
import hashlib
import json
import math
import importlib
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from platinum_development_probe import graph_matches


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flatten(values):
    return [v for point in values for v in point]


def measure(analytic, numerical):
    errors = [abs(a - b) for a, b in zip(analytic, numerical)]
    i = max(range(len(errors)), key=errors.__getitem__)
    return {'axis': i, 'absolute_error': errors[i], 'analytic': analytic[i],
            'numerical': numerical[i], 'scaled_error': errors[i] / max(1., abs(analytic[i]), abs(numerical[i]))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', required=True)
    parser.add_argument('--coordinates', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    data = json.load(gzip.open(args.coordinates, 'rt'))
    supplier = Chem.SDMolSupplier(args.sdf, removeHs=False)
    out = {'source_sdf_sha256': sha(args.sdf), 'coordinate_file_sha256': sha(args.coordinates),
           'script_sha256': sha(__file__), 'rdkit': rdBase.rdkitVersion,
           'chematic': chematic.__version__, 'chematic_extension_sha256': sha(importlib.import_module('chematic.chematic').__file__), 'platform': platform.platform(),
           'steps': [1e-6, 1e-7, 1e-8], 'cases': []}
    for case in data['cases']:
        if 'coordinates' not in case:
            continue
        ref = supplier[case['index']]
        if ref is None or ref.GetProp('_Name') != case['name']:
            raise ValueError('reference identity mismatch')
        smiles = Chem.MolToSmiles(Chem.RemoveHs(ref), isomericSmiles=True)
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        ours = chematic.from_smiles(smiles).add_hydrogens()
        if not graph_matches(ours, mol):
            raise ValueError('graph mismatch')
        coords = case['coordinates']
        if len(coords) != mol.GetNumAtoms() * 3:
            raise ValueError('coordinate count mismatch')
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, coords[3*i:3*i+3])
        mol.AddConformer(conf, assignId=True)
        points = [coords[i:i+3] for i in range(0, len(coords), 3)]
        row = {'index': case['index'], 'name': case['name'], 'smiles': smiles, 'fields': {}}
        for variant in ('MMFF94', 'UFF'):
            if variant == 'UFF':
                field = AllChem.UFFGetMoleculeForceField(mol)
                own = ours.uff_force_field(points)
            else:
                props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
                field = AllChem.MMFFGetMoleculeForceField(mol, props) if props else None
                own = ours.mmff_force_field(points, variant=variant)
            if field is None:
                row['fields'][variant] = {'status': 'unsupported'}
                continue
            gradients = {'rdkit': list(field.CalcGrad()), 'compatible': flatten(own.gradient())}
            if variant != 'UFF':
                gradients['consistent_planes'] = flatten(own.gradient_with_consistent_planes())
                gradients['consistent_torsions'] = flatten(own.gradient_with_consistent_torsions())
            result = {'rdkit_energy': field.CalcEnergy(), 'chematic_energy': own.energy(), 'steps': {}}
            for step in out['steps']:
                rd, ch = [], []
                for axis in range(len(coords)):
                    hi = list(coords); lo = list(coords)
                    hi[axis] += step; lo[axis] -= step
                    rd.append((field.CalcEnergy(hi) - field.CalcEnergy(lo)) / (2 * step))
                    hp = [hi[i:i+3] for i in range(0, len(hi), 3)]
                    lp = [lo[i:i+3] for i in range(0, len(lo), 3)]
                    ch.append((own.set_positions_(hp).energy() - own.set_positions_(lp).energy()) / (2 * step))
                own.set_positions_(points)
                result['steps'][str(step)] = {'rdkit': measure(gradients['rdkit'], rd),
                    'compatible': measure(gradients['compatible'], ch),
                    'energy_derivative_cross_engine': measure(rd, ch)}
                for key in ('consistent_planes', 'consistent_torsions'):
                    if key in gradients:
                        result['steps'][str(step)][key] = measure(gradients[key], ch)
            if not all(math.isfinite(v) for v in (result['rdkit_energy'], result['chematic_energy'])):
                raise ValueError('nonfinite energy')
            row['fields'][variant] = result
        out['cases'].append(row)
        print(case['name'], flush=True)
    Path(args.output).write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
