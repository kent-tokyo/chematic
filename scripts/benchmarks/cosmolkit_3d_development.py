#!/usr/bin/env python3
"""Score COSMolKit on the same exposed development rows as the fixed-count probe."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import cosmolkit as ck
from rdkit import Chem, rdBase
from platinum_development_probe import rms_values
from platinum_geometry import check_geometry, prepare_geometry_checks
from platinum_stereo_strict import stereo_violations_strict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--sdf', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    baseline = json.loads(args.baseline.read_text())
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_hash != baseline['source_sha256']:
        raise ValueError('reference SDF differs from baseline')
    result = {'cosmolkit': ck.__version__, 'rdkit_scorer': rdBase.rdkitVersion,
              'baseline_sha256': hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
              'source_sha256': source_hash, 'count': baseline['candidate_limit'],
              'seed': baseline['seed'], 'mmff_iterations': 200, 'rows': []}
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    for row in baseline['results']:
        item = {'index': row['index'], 'name': row['name']}
        try:
            reference = supplier[row['index']]
            if reference is None or reference.GetProp('_Name') != row['name']:
                raise ValueError('reference identity mismatch')
            smiles = row['smiles']
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            other = ck.Molecule.from_smiles(smiles).with_hydrogens()
            atoms = [atom.atomic_number for atom in other.atoms]
            bonds = sorted((min(b.begin, b.end), max(b.begin, b.end), b.order_name.upper()) for b in other.bonds)
            reference_bonds = sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()).upper()) for b in mol.GetBonds())
            if atoms != [atom.GetAtomicNum() for atom in mol.GetAtoms()] or bonds != reference_bonds:
                raise ValueError('indexed molecular graphs differ')
            start = time.perf_counter()
            other.embed_3d_conformers_(result['count'], random_seed=result['seed'], num_threads=1)
            optimized = other.with_mmff_optimized_conformers(num_threads=1, max_iterations=200)
            item['seconds'] = time.perf_counter() - start
            item['statuses'] = [outcome.status_code for outcome in optimized.conformer_results]
            coords = [conformer.coordinates for conformer in optimized.molecule.conformers_3d]
            checks = prepare_geometry_checks(mol)
            item['geometry_reference'] = {'bonds_supported': len(checks[0]), 'bonds_total': mol.GetNumBonds(),
                'angles_supported': len(checks[1]), 'angles_total': sum(a.GetDegree() * (a.GetDegree() - 1) // 2 for a in mol.GetAtoms())}
            item['rmsd'] = rms_values(mol, Chem.RemoveHs(reference), coords)
            item['stereo'] = stereo_violations_strict(mol, coords)
            item['geometry'] = [check_geometry(conf, checks) for conf in coords]
            item['status'] = 'ok'
        except Exception as exc:
            item['status'] = 'error'
            item['error'] = repr(exc)
        result['rows'].append(item)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        print(item['index'], item['status'], flush=True)


if __name__ == '__main__':
    main()
