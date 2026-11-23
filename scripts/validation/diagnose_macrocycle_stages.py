#!/usr/bin/env python3
"""Diagnose one exposed macrocycle through embedding and MMFF stages."""
import argparse
import gzip
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_stereo_strict import stereo_violations_strict  # noqa: E402


def stereo_counts(mol, conformers):
    checks = stereo_violations_strict(mol, conformers)
    return {'conformers': len(checks),
            'contradicted': sum(bool(x['violated_atoms'] or x['violated_bonds']) for x in checks),
            'unevaluable': sum(bool(x['unevaluable_atoms'] or x['unevaluable_bonds']) for x in checks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('development', type=Path)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    source = gzip.decompress(args.development.read_bytes())
    if hashlib.sha256(source).hexdigest() != '821e4cbb342bbb94db5f4a44a62560d17001ead9e1a7f1190f0cb58ab3d14708':
        raise ValueError('development result fingerprint mismatch')
    rows = [r for r in json.loads(source)['results']
            if r['name'] == '5vk0__Ligand_phospholipids_ALA_L_4']
    if len(rows) != 1:
        raise ValueError('expected exactly one exposed development case')
    smiles = rows[0]['smiles']
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    other = chematic.from_smiles(smiles).add_hydrogens()
    result = {'name': rows[0]['name'], 'smiles_sha256': hashlib.sha256(smiles.encode()).hexdigest(),
              'atoms': mol.GetNumAtoms(), 'rdkit_version': rdBase.rdkitVersion,
              'chematic_version': importlib.metadata.version('chematic'),
              'platform': platform.platform(), 'wheel_sha256': hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seed': 42, 'attempts': 20, 'mmff_iterations': 200, 'repair_iterations': 400,
              'stages': {}}
    methods = (
        ('chematic_raw', lambda: other.rdkit_embed_multiple_confs(20, random_seed=42)),
        ('chematic_compatible', lambda: other.rdkit_embed_multiple_confs_mmff(
            20, random_seed=42, max_mmff_iters=200)),
        ('chematic_corrected', lambda: other.etkdg_embed_multiple_confs_mmff_corrected(
            20, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400)),
    )
    for name, method in methods:
        start = time.perf_counter()
        try:
            candidates = method()
            coords = candidates if name == 'chematic_raw' else [candidate[2] for candidate in candidates]
            stage = {'count': len(candidates), 'seconds': time.perf_counter() - start,
                     'stereo': stereo_counts(mol, coords)}
            if name != 'chematic_raw':
                stage['statuses'] = [candidate[0] for candidate in candidates]
        except Exception as exc:
            stage = {'error': repr(exc), 'seconds': time.perf_counter() - start}
        result['stages'][name] = stage
    start = time.perf_counter()
    try:
        ref = Chem.Mol(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed, params.numThreads = 42, 1
        ids = list(AllChem.EmbedMultipleConfs(ref, numConfs=20, params=params))
        raw = [[list(ref.GetConformer(cid).GetAtomPosition(i)) for i in range(ref.GetNumAtoms())]
               for cid in ids]
        raw_seconds = time.perf_counter() - start
        result['stages']['rdkit_raw'] = {'count': len(ids), 'seconds': raw_seconds,
                                         'stereo': stereo_counts(mol, raw)}
        statuses = [AllChem.MMFFOptimizeMolecule(ref, confId=cid, maxIters=200) for cid in ids]
        coords = [[list(ref.GetConformer(cid).GetAtomPosition(i)) for i in range(ref.GetNumAtoms())]
                  for cid in ids]
        result['stages']['rdkit_optimized'] = {'count': len(ids), 'statuses': statuses,
                                               'seconds_including_embedding': time.perf_counter() - start,
                                               'stereo': stereo_counts(mol, coords)}
    except Exception as exc:
        result['stages']['rdkit_error'] = repr(exc)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result['stages'], indent=2))
    if any('error' in stage for stage in result['stages'].values() if isinstance(stage, dict)) or 'rdkit_error' in result['stages']:
        raise SystemExit('diagnostic encountered an error')


if __name__ == '__main__':
    main()
