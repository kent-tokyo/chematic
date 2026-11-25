#!/usr/bin/env python3
"""Evaluate opt-in unoptimized fallback on previously exposed B/Se development molecules."""
import argparse
import hashlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from cod_hard_deadline import annotate, run_engine
from platinum_geometry import prepare_geometry_checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', type=Path, required=True)
    parser.add_argument('--split', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != split['source_sha256']:
        raise ValueError('SDF hash mismatch')
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    output = {'source_sha256': split['source_sha256'], 'split_sha256': hashlib.sha256(split_bytes).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'api_sha256': hashlib.sha256(Path(chematic._timed_conformers.__file__).read_bytes()).hexdigest(),
              'rdkit': rdBase.rdkitVersion, 'count': 20, 'seconds': 30, 'rows': []}
    for row in split['rows']:
        if row.get('split') != 'development':
            continue
        ref = supplier[row['index']]
        if ref is None:
            continue
        elements = {atom.GetAtomicNum() for atom in ref.GetAtoms()}
        if not elements.intersection({5, 34}):
            continue
        entry = {'index': row['index'], 'name': row['name'], 'elements': sorted(elements.intersection({5, 34}))}
        try:
            if ref.GetProp('_Name') != row['name']:
                raise ValueError('reference identity mismatch')
            ref = Chem.RemoveHs(ref)
            smiles = Chem.MolToSmiles(ref, isomericSmiles=True)
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            other = chematic.from_smiles(smiles).add_hydrogens()
            from platinum_development_probe import graph_matches
            if not graph_matches(other, mol):
                raise ValueError('indexed graphs differ')
            checks = prepare_geometry_checks(mol)
            entry['geometry_reference'] = {'bonds_supported': len(checks[0]), 'bonds_total': mol.GetNumBonds(),
                'angles_supported': len(checks[1]), 'angles_total': sum(a.GetDegree() * (a.GetDegree()-1)//2 for a in mol.GetAtoms())}
            result = chematic.timed_conformer_batch([smiles], seconds=30, num_confs=20,
                                                     repair_mmff_iters=400, allow_unoptimized=True)[0]
            entry['fallback_reason'] = result['fallback_reason']
            candidate = {'status': result['status'], 'elapsed': result['elapsed'],
                         'error': result['error'],
                         'candidates': [(status, coords) for status, _, coords in result['conformers']]}
            entry['engines'] = {'chematic': annotate(candidate, mol, ref, checks),
                                'rdkit': annotate(run_engine('rdkit', smiles, 20, 30, 400), mol, ref, checks)}
        except Exception as exc:
            entry['input_error'] = repr(exc)
        output['rows'].append(entry)
        args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
        print(entry['index'], entry.get('input_error', 'ok'), flush=True)


if __name__ == '__main__':
    main()
