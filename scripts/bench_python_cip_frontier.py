#!/usr/bin/env python3
"""Fresh-process CIP comparison with two source wheels and RDKit.

Native results are retained before normalization. RDKit mutates molecule
properties; chematic returns diagnostic dictionaries. Both public API calls
are timed as offered, and label extraction is outside timing. Consequently
these are Python API timings, not a comparison of identical internal kernels.
All input rows and all blocks are retained. A mismatching lane is unranked.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import importlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
import zipfile

from bench_python_speed_frontier import agreement, corpus_rows, normal, table
from verify_python_cip_frontier import rdkit_labels, source_labels

OPS = ('cip', 'canonical_smiles', 'smarts_write', 'qed', 'morgan2_chiral', 'maccs', 'fp_layered')
MODES = ('first_use', 'repeated_use', 'parse_inclusive')
ARMS = ('baseline', 'candidate', 'rdkit')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def orders(blocks):
    if blocks < 24 or blocks % 6:
        raise ValueError('at least 24 blocks, a multiple of six')
    values = list(itertools.permutations(ARMS)) * (blocks // 6)
    random.Random(20261019).shuffle(values)
    return values


def interval(ratios):
    rng = random.Random(20261019)
    boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10_000))
    return {'median': statistics.median(ratios), 'ci95': [boot[249], boot[9749]],
            'winning_blocks': sum(x > 1 for x in ratios), 'blocks': len(ratios)}


def worker(a):
    engine = 'rdkit' if a.worker == 'rdkit' else 'chematic'
    parse, funcs = table(engine)
    data = a.corpus.read_bytes()
    rows = corpus_rows(data, a.offset, a.limit)
    def parsed():
        molecules = [parse(s) for s in rows]
        if any(m is None for m in molecules):
            raise ValueError('parse refusal; lane cannot be ranked')
        return molecules
    result = {'operations': {}, 'corpus_sha256': sha(data)}
    for op in a.ops.split(','):
        fn = funcs[op]
        fn(parse(rows[0]))  # disjoint warm-up object
        times = {}
        reference = None
        for mode in MODES:
            molecules = None if mode == 'parse_inclusive' else parsed()
            if mode == 'repeated_use':
                for m in molecules:
                    fn(m)
            start = time.perf_counter_ns()
            retained = ([(m := parse(s), fn(m)) for s in rows] if mode == 'parse_inclusive'
                        else [(m, fn(m)) for m in molecules])
            times[mode] = time.perf_counter_ns() - start
            if op == 'cip':
                values = [rdkit_labels(m) if engine == 'rdkit' else source_labels(v)
                          for m, v in retained]
            else:
                values = [normal(v, engine, op) for _, v in retained]
            if reference is not None and values != reference:
                raise ValueError(f'{op}: outputs differ across preparation modes')
            reference = values
            del retained, molecules
        result['operations'][op] = {**times, 'values': reference}
    module = __import__(engine)
    native = {p.name: sha(p.read_bytes()) for p in Path(module.__file__).parent.glob('*.so')}
    if engine == 'rdkit':
        native = {name: sha(Path(importlib.import_module(name).__file__).read_bytes())
                  for name in ('rdkit.rdBase', 'rdkit.Chem.rdCIPLabeler', 'rdkit.Chem.rdchem',
                               'rdkit.Chem.rdmolfiles', 'rdkit.Chem.rdMolDescriptors',
                               'rdkit.Chem.rdFingerprintGenerator')}
    result.update(version=module.__version__, extension_sha256=native,
        peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if sys.platform == 'darwin' else 1024))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker', choices=ARMS)
    p.add_argument('--corpus', type=Path, required=True)
    p.add_argument('--offset', type=int, default=0)
    p.add_argument('--limit', type=int, default=1000)
    p.add_argument('--ops', default=','.join(OPS))
    p.add_argument('--blocks', type=int, default=24)
    p.add_argument('--output', type=Path)
    p.add_argument('--baseline-identity', type=Path)
    p.add_argument('--candidate-identity', type=Path)
    for arm in ARMS[:2]:
        p.add_argument('--'+arm+'-wheel', type=Path)
        p.add_argument('--'+arm+'-path', type=Path)
    a = p.parse_args()
    ops = a.ops.split(',')
    if set(ops) - set(OPS):
        p.error('unsupported operation')
    if a.worker:
        print(json.dumps(worker(a), sort_keys=True, allow_nan=False))
        return
    if any(getattr(a, x) is None for x in ('output', 'baseline_identity', 'candidate_identity',
                                         'baseline_wheel', 'candidate_wheel', 'baseline_path', 'candidate_path')):
        p.error('both wheel identities, wheels, installed paths and output required')
    if a.output.exists():
        raise ValueError('refusing to overwrite existing observations')
    order = orders(a.blocks)
    plan = {'schema': 'python-cip-frontier/v1', 'ops': ops, 'modes': list(MODES),
            'orders': order, 'rows': a.limit, 'offset': a.offset,
            'corpus_sha256': sha(a.corpus.read_bytes()),
            'identities': {arm: json.loads(getattr(a, arm+'_identity').read_text()) for arm in ARMS[:2]},
            'wheel_sha256': {arm: sha(getattr(a, arm+'_wheel').read_bytes()) for arm in ARMS[:2]},
            'script_sha256': {name: sha((Path(__file__).parent.parent / name).read_bytes())
                              for name in ('scripts/bench_python_cip_frontier.py',
                                           'scripts/bench_python_speed_frontier.py',
                                           'scripts/verify_python_cip_frontier.py',
                                           'validation/cosmolkit_comparison/bench_corpus.py')},
            'adoption_policy': 'Exact old/new CIP dictionaries and unresolved results on all 14072 rows; source CIP first-use and parse-inclusive CI95 lower bounds above one. Retain and disclose every control lane and every comparator mismatch; no toolkit-wide speed claim.',
            'protocol': __doc__, 'created_utc': datetime.now(timezone.utc).isoformat(),
            'host': {'platform': platform.platform(), 'python': sys.version, 'cpus': os.cpu_count()}}
    plan_path = a.output.with_suffix('.plan.json')
    if plan_path.exists():
        raise ValueError('measurement plan already exists')
    plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    report = {'plan': plan, 'plan_sha256': sha(plan_path.read_bytes()), 'blocks': [],
              'agreement': {}, 'summary': {}, 'complete': False}
    expected = {}
    native = {}
    for arm in ARMS[:2]:
        with zipfile.ZipFile(getattr(a, arm+'_wheel')) as z:
            native[arm] = {Path(n).name: sha(z.read(n)) for n in z.namelist() if n.endswith('.so')}
    def save():
        a.output.write_text(json.dumps(report, indent=2)+'\n')
    save()
    for i, names in enumerate(order):
        block = {'order': names, 'load_average': list(os.getloadavg())}
        for arm in names:
            env = os.environ.copy()
            env['PYTHONDONTWRITEBYTECODE'] = '1'
            if arm == 'rdkit':
                env.pop('PYTHONPATH', None)
            else:
                env['PYTHONPATH'] = str(getattr(a, arm+'_path'))
            cmd = [sys.executable, __file__, '--worker', arm, '--corpus', str(a.corpus),
                   '--limit', str(a.limit), '--offset', str(a.offset), '--ops', a.ops]
            proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=300)
            if proc.returncode:
                raise RuntimeError(f'{arm}: worker failed; no speed claim\n{proc.stderr[-2000:]}')
            sample = json.loads(proc.stdout)
            if sample['corpus_sha256'] != plan['corpus_sha256']:
                raise ValueError('input changed')
            if arm in native and sample['extension_sha256'] != native[arm]:
                raise ValueError('imported native binary differs from supplied wheel')
            outputs = {op: v.pop('values') for op, v in sample['operations'].items()}
            if arm not in expected:
                expected[arm] = outputs
            if expected[arm] != outputs:
                raise ValueError(f'{arm}: nondeterministic output')
            sample['output_sha256'] = sha(json.dumps(outputs, sort_keys=True).encode())
            block[arm] = sample
        report['blocks'].append(block)
        save()
        print('block', i+1, '/', a.blocks, flush=True)
    for op in ops:
        checks = {arm: agreement(expected[arm][op], expected['candidate'][op], exact=arm=='baseline')
                  for arm in ('baseline', 'rdkit')}
        report['agreement'][op] = checks
        for mode in MODES:
            report['summary'][op+'/'+mode] = {
                arm: interval([b[arm]['operations'][op][mode] / b['candidate']['operations'][op][mode]
                               for b in report['blocks']]) if checks[arm]['passed'] else None
                for arm in ('baseline', 'rdkit')}
    report['complete'] = True
    save()


if __name__ == '__main__':
    main()
