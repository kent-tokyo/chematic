#!/usr/bin/env python3
"""Evaluate the opt-in unoptimized fallback using the packaged Python API."""
import argparse
import hashlib
import json
import platform
import resource
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from platinum_geometry import prepare_geometry_checks


from cod_hard_deadline import annotate, run_engine

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', required=True, type=Path)
    parser.add_argument('--split', required=True, type=Path)
    parser.add_argument('--audit', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--mode', choices=['fixed', 'timed'], required=True)
    parser.add_argument('--partition', choices=['development', 'holdout'], default='development')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--fixed-timeout', type=float, default=30.0)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    source = args.sdf.read_bytes()
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    audit = json.loads(args.audit.read_text())
    if hashlib.sha256(source).hexdigest() != split['source_sha256']:
        raise ValueError('SDF hash mismatch')
    if (audit['split_sha256'] != hashlib.sha256(split_bytes).hexdigest()
            or audit['source_sha256'] != split['source_sha256'] or audit['overlaps']):
        raise ValueError('invalid exclusion audit')
    selected = [r for r in split['rows'] if r.get('split') == args.partition]
    selected.sort(key=lambda r: hashlib.sha256(r['name'].encode()).digest())
    if args.limit is not None:
        selected = selected[:args.limit]
    if not selected:
        raise ValueError('no records selected')
    seconds = 0.2 if args.mode == 'timed' else args.fixed_timeout
    hash_file = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    header = {'type': 'metadata', 'mode': args.mode, 'seconds': seconds, 'count': 20,
              'seed_schedule': '42 * attempt, one-based', 'repair_iters': 400,
              'rdkit': rdBase.rdkitVersion, 'platform': platform.platform(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'chematic_version': chematic.__version__,
              'chematic_extension_sha256': hash_file(chematic.chematic.__file__),
              'candidate_api_sha256': hash_file(chematic._timed_conformers.__file__),
              'rdkit_extension_sha256': hash_file(rdBase.__file__),
              'scorers_sha256': {name: hash_file(Path(__file__).with_name(name)) for name in
                                   ('cod_hard_deadline.py', 'platinum_development_probe.py', 'platinum_stereo_strict.py', 'platinum_geometry.py')},
              'split_sha256': hashlib.sha256(split_bytes).hexdigest(),
              'source_sha256': split['source_sha256'], 'total': len(selected),
              'selection': args.partition, 'timeout_policy': 'retain complete candidates delivered by deadline',
              'process_peak_scope': 'RDKit worker and parent high water RSS, bytes on macOS and KiB on Linux; candidate child RSS unavailable',
              'timing_scope': 'candidate budget includes parsing and hydrogen addition, RDKit budget begins after preparation',
              'candidate_policy': 'opt-in stereo-checked unoptimized fallback for MMFF type assignment errors'}
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    completed = set()
    if args.resume:
        with args.output.open() as previous:
            saved = [json.loads(line) for line in previous]
        if not saved or saved[0] != header:
            raise ValueError('resume metadata differs')
        if any(row.get('type') == 'complete' for row in saved[1:]):
            raise ValueError('run already complete')
        completed = {row['index'] for row in saved[1:] if row.get('type') == 'result'}
        expected = {row['index'] for row in selected}
        if len(completed) != len(saved) - 1 or not completed <= expected:
            raise ValueError('resume records contain duplicates or unexpected indices')
    with args.output.open('a' if args.resume else 'x') as f:
        if not args.resume:
            f.write(json.dumps(header, allow_nan=False) + '\n'); f.flush()
        for number, row in enumerate(selected):
            if row['index'] in completed:
                continue
            result = {'type': 'result', 'index': row['index'], 'name': row['name'], 'engines': {}}
            try:
                reference = supplier[row['index']]
                if reference is None or reference.GetProp('_Name') != row['name']:
                    raise ValueError('reference identity mismatch')
                reference = Chem.RemoveHs(reference)
                smiles = Chem.MolToSmiles(reference, isomericSmiles=True)
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                checks = prepare_geometry_checks(mol)
                result['geometry_reference'] = {'bonds_supported': len(checks[0]), 'bonds_total': mol.GetNumBonds(),
                    'angles_supported': len(checks[1]), 'angles_total': sum(a.GetDegree() * (a.GetDegree()-1) // 2 for a in mol.GetAtoms())}
                order = ['chematic', 'rdkit'] if number % 2 == 0 else ['rdkit', 'chematic']
                result['order'] = order
                for engine in order:
                    try:
                        if engine == 'chematic':
                            candidate = chematic.timed_conformer_batch(
                                [smiles], seconds=seconds, num_confs=20,
                                repair_mmff_iters=400, allow_unoptimized=True)[0]
                            outcome = {'status': candidate['status'], 'error': candidate['error'],
                                       'fallback_reason': candidate['fallback_reason'],
                                       'elapsed': candidate['elapsed'],
                                       'elapsed_with_cleanup': candidate['elapsed'],
                                       'candidates': [(status, coords) for status, _, coords
                                                      in candidate['conformers']]}
                        else:
                            outcome = run_engine(engine, smiles, 20, seconds, 400)
                    except Exception as exc:
                        result['engines'][engine] = {'status': 'error', 'error': repr(exc)}
                        continue
                    try:
                        result['engines'][engine] = annotate(outcome, mol, reference, checks)
                    except Exception as exc:
                        outcome['scoring_error'] = repr(exc)
                        result['engines'][engine] = outcome
            except Exception as exc:
                result['input_error'] = repr(exc)
            result['parent_peak_rss'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            f.write(json.dumps(result, allow_nan=False) + '\n'); f.flush()
            print(number + 1, '/', len(selected), row['index'], flush=True)
        f.write(json.dumps({'type': 'complete', 'rows': len(selected)}) + '\n')


if __name__ == '__main__':
    main()
