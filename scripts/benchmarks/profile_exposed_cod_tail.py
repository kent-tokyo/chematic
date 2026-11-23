#!/usr/bin/env python3
"""Profile embedding and MMFF separately on one exposed COD development case."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import time
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from platinum_development_probe import graph_matches


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sdf', required=True, type=Path)
    p.add_argument('--baseline', required=True, type=Path)
    p.add_argument('--wheel', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    baseline = json.loads(args.baseline.read_text())
    source_sha = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_sha != baseline['source_sha256']:
        raise ValueError('reference SDF hash mismatch')
    selected = [r for r in baseline['results'] if r['index'] == 3425]
    if len(selected) != 1 or selected[0]['name'] != 'AIAOXFULFQJBOK-ONALTULYSA-N':
        raise ValueError('exposed diagnostic identity mismatch')
    row = selected[0]
    reference = Chem.SDMolSupplier(str(args.sdf), removeHs=False)[row['index']]
    if reference is None or reference.GetProp('_Name') != row['name']:
        raise ValueError('reference identity mismatch')
    mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    other = chematic.from_smiles(row['smiles']).add_hydrogens()
    if not graph_matches(other, mol):
        raise ValueError('atom or bond identity mismatch')
    result = {'index': row['index'], 'name': row['name'], 'source_sha256': source_sha,
              'baseline_sha256': hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
              'wheel_sha256': hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'rdkit': rdBase.rdkitVersion, 'chematic': importlib.metadata.version('chematic'),
              'platform': platform.platform(), 'seed': 42, 'max_mmff_iters': 200,
              'scope': 'exposed development case; compatibility embedding and MMFF only; no stereo checks or repair',
              'stages': {}}

    def checkpoint():
        tmp = args.output.with_suffix(args.output.suffix + '.tmp')
        tmp.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        tmp.replace(args.output)

    def stage(name, operation):
        started = time.perf_counter()
        try:
            value = operation()
            result['stages'][name] = {'seconds': time.perf_counter() - started, 'status': 'ok'}
            checkpoint()
            return value
        except Exception as exc:
            result['stages'][name] = {'seconds': time.perf_counter() - started, 'status': 'error', 'error': repr(exc)}
            checkpoint()
            return None

    # Warm both implementations on the same small molecule, outside measured phases.
    small = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    small_other = chematic.from_smiles('CCO').add_hydrogens()
    AllChem.EmbedMolecule(small, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(small, maxIters=200)
    small_coords = small_other.rdkit_embed_multiple_confs(1, random_seed=42)
    if not small_coords:
        raise RuntimeError('warmup failed')
    small_other.rdkit_mmff_optimize(small_coords[0], max_iters=200)
    checkpoint()
    coords = stage('chematic_embed', lambda: other.rdkit_embed_multiple_confs(1, random_seed=42))
    if coords:
        optimized = stage('chematic_mmff', lambda: other.rdkit_mmff_optimize(coords[0], max_iters=200))
        if optimized is not None:
            result['stages']['chematic_mmff']['optimizer_status'] = optimized[0]
            checkpoint()
    generated = Chem.Mol(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = 42; params.numThreads = 1
    embedded = stage('rdkit_embed', lambda: AllChem.EmbedMolecule(generated, params))
    if embedded == 0:
        status = stage('rdkit_mmff', lambda: AllChem.MMFFOptimizeMolecule(generated, maxIters=200))
        if status is not None:
            result['stages']['rdkit_mmff']['optimizer_status'] = status
            checkpoint()
    elif embedded is not None:
        result['stages']['rdkit_embed']['embed_status'] = embedded
        checkpoint()


if __name__ == '__main__':
    main()
