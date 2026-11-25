#!/usr/bin/env python3
"""Paired conformer evaluation with independently interruptible engines."""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import platform
import queue
import resource
import threading
import time
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from platinum_development_probe import graph_matches, rms_values
from platinum_geometry import check_geometry, prepare_geometry_checks
from platinum_stereo_strict import stereo_violations_strict


def engine_worker(conn, engine, smiles, count, repair_iters):
    try:
        import chematic
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        other = None
        if engine == 'chematic':
            other = chematic.from_smiles(smiles).add_hydrogens()
            if not graph_matches(other, mol):
                raise ValueError('indexed atom or bond graph mismatch')
        conn.send(('ready', None))
        if conn.recv() != 'start':
            return
        for attempt in range(count):
            seed = 42 * (attempt + 1)
            if engine == 'chematic':
                candidates = other.etkdg_embed_multiple_confs_mmff_corrected(
                    1, random_seed=seed, max_mmff_iters=200, repair_mmff_iters=repair_iters)
                candidate = None if not candidates else (candidates[0][0], candidates[0][2])
            else:
                generated = Chem.Mol(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed, params.numThreads = seed, 1
                candidate = None
                if AllChem.EmbedMolecule(generated, params) == 0:
                    status = AllChem.MMFFOptimizeMolecule(generated, maxIters=200)
                    if status == 0:
                        props = AllChem.MMFFGetMoleculeProperties(generated)
                        field = AllChem.MMFFGetMoleculeForceField(generated, props)
                        if any(not math.isfinite(g) or abs(g) > 0.01 for g in field.CalcGrad()):
                            status = 1
                    coords = [list(generated.GetConformer().GetAtomPosition(i))
                              for i in range(generated.GetNumAtoms())]
                    candidate = (status, coords)
            conn.send(('candidate', (candidate, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
        conn.send(('done', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    except BaseException as exc:
        try:
            conn.send(('error', f'{type(exc).__name__}: {exc}'))
        except (EOFError, OSError):
            pass
    finally:
        conn.close()


def run_engine(engine, smiles, count, seconds, repair_iters):
    ctx = mp.get_context('spawn')
    parent, child = ctx.Pipe()
    process = ctx.Process(target=engine_worker, args=(child, engine, smiles, count, repair_iters))
    process.start()
    child.close()
    events = queue.Queue()
    def reader():
        try:
            while True:
                events.put((time.monotonic(), *parent.recv()))
        except (EOFError, OSError):
            events.put((time.monotonic(), 'closed', None))
    thread = threading.Thread(target=reader, daemon=True)
    thread.start()
    output = {'status': 'error', 'attempted': 0, 'candidates': [], 'error': None}
    try:
        try:
            _, kind, value = events.get(timeout=30)
        except queue.Empty:
            output['error'] = 'worker startup timeout'
            return output
        if kind != 'ready':
            output['error'] = value or 'worker exited before ready'
            return output
        start = time.monotonic()
        deadline = start + seconds
        parent.send('start')
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                output['status'] = 'timeout'
                break
            try:
                received, kind, value = events.get(timeout=remaining)
            except queue.Empty:
                output['status'] = 'timeout'
                break
            if received > deadline:
                output['status'] = 'timeout'
                break
            if kind == 'candidate':
                candidate, peak_rss = value
                output['worker_peak_rss'] = peak_rss
                output['attempted'] += 1
                if candidate is not None:
                    output['candidates'].append(candidate)
            elif kind == 'done':
                output['worker_peak_rss'] = value
                output['status'] = 'ok'
                break
            else:
                output['error'] = value or 'worker exited before completion'
                break
        output['elapsed'] = time.monotonic() - start
    finally:
        if process.is_alive():
            process.terminate()
        process.join(timeout=1)
        if process.is_alive():
            process.kill()
            process.join()
        thread.join(timeout=1)
        parent.close()
    output['elapsed_with_cleanup'] = time.monotonic() - start
    return output


def annotate(out, mol, reference, checks):
    candidates = out.pop('candidates')
    coords = [coords for _, coords in candidates]
    out['statuses'] = [status for status, _ in candidates]
    out['accepted'] = len(coords)
    out['rmsd'] = rms_values(mol, reference, coords)
    out['stereo'] = stereo_violations_strict(mol, coords)
    out['geometry'] = [check_geometry(c, checks) for c in coords]
    return out


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
              'rdkit_extension_sha256': hash_file(rdBase.__file__),
              'scorers_sha256': {name: hash_file(Path(__file__).with_name(name)) for name in
                                   ('platinum_development_probe.py', 'platinum_stereo_strict.py', 'platinum_geometry.py')},
              'split_sha256': hashlib.sha256(split_bytes).hexdigest(),
              'source_sha256': split['source_sha256'], 'total': len(selected),
              'selection': args.partition, 'timeout_policy': 'retain complete candidates delivered by deadline',
              'process_peak_scope': 'per-worker and parent process high water RSS; bytes on macOS, KiB on Linux; timed-out workers are sampled only at completed attempts'}
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
