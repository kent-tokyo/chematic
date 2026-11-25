#!/usr/bin/env python3
"""Summarize independently interrupted COD results with scaffold-cluster bootstrap."""
import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path


def recovered(row, engine, converged=False):
    result = row.get('engines', {}).get(engine, {})
    if result.get('scoring_error'):
        return False
    ref = row.get('geometry_reference', {})
    if (ref.get('bonds_supported') != ref.get('bonds_total')
            or ref.get('angles_supported') != ref.get('angles_total') or not ref):
        return False
    arrays = [result.get(name, []) for name in ('rmsd', 'stereo', 'geometry', 'statuses')]
    if len({len(x) for x in arrays}) != 1:
        return False
    for rms, stereo, geom, status in zip(*arrays):
        if (math.isfinite(rms) and rms <= 1 and not stereo['atoms'] and not stereo['bonds']
                and geom['finite'] and geom['bond_violations'] == 0
                and geom['angle_violations'] == 0 and geom['heavy_atom_clashes'] == 0
                and (not converged or status == 0)):
            return True
    return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input', type=Path)
    p.add_argument('--split', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--bootstrap', type=int, default=10000)
    args = p.parse_args()
    records = [json.loads(line) for line in args.input.read_text().splitlines()]
    header, marker = records[0], records[-1]
    rows = records[1:-1]
    if marker.get('type') != 'complete' or marker['rows'] != len(rows) or header['total'] != len(rows):
        raise ValueError('incomplete results')
    split = json.loads(args.split.read_text())
    import hashlib
    if hashlib.sha256(args.split.read_bytes()).hexdigest() != header['split_sha256']:
        raise ValueError('split fingerprint mismatch')
    by_index = {r['index']: r for r in split['rows'] if r.get('split') == header['selection']}
    if len({r['index'] for r in rows}) != len(rows) or not set(r['index'] for r in rows) <= set(by_index):
        raise ValueError('unexpected or duplicate record')
    result = {'rows': len(rows), 'mode': header['mode'], 'source_script_sha256': header['script_sha256'],
              'bootstrap_seed': 42, 'bootstrap_samples': args.bootstrap, 'engines': {}, 'recovery': {}}
    for engine in ('chematic', 'rdkit'):
        values = [r.get('engines', {}).get(engine, {}) for r in rows]
        stereo = [v for item in values for v in item.get('stereo', [])]
        elapsed = sorted(v['elapsed_with_cleanup'] for v in values if 'elapsed_with_cleanup' in v)
        memory = sorted(v['worker_peak_rss'] for v in values if 'worker_peak_rss' in v)
        result['engines'][engine] = {'statuses': dict(Counter(v.get('status', 'missing') for v in values)),
            'input_errors': sum('input_error' in row for row in rows),
            'scoring_errors': sum('scoring_error' in v for v in values),
            'conformers': sum(v.get('accepted', 0) for v in values),
            'converged': sum(s == 0 for v in values for s in v.get('statuses', [])),
            'stereo_violated': sum(bool(v.get('violated_atoms') or v.get('violated_bonds')) for v in stereo),
            'stereo_unevaluable': sum(bool(v.get('unevaluable_atoms') or v.get('unevaluable_bonds')) for v in stereo),
            'stereo_coverage_missing': sum(v.get('accepted', 0) - len(v.get('stereo', [])) for v in values),
            'elapsed_seconds': {'count': len(elapsed), 'median': elapsed[len(elapsed)//2] if elapsed else None,
                                'p95': elapsed[math.ceil(len(elapsed)*.95)-1] if elapsed else None,
                                'max': elapsed[-1] if elapsed else None},
            'worker_peak_rss': {'count': len(memory), 'median': memory[len(memory)//2] if memory else None,
                                'max': memory[-1] if memory else None}}
    for converged in (False, True):
        groups = {}
        for row in rows:
            pair = []
            for engine in ('chematic', 'rdkit'):
                complete = header['mode'] != 'fixed' or row.get('engines', {}).get(engine, {}).get('status') == 'ok'
                pair.append(int(complete and recovered(row, engine, converged)))
            key = by_index[row['index']]['scaffold_sha256']
            groups.setdefault(key, []).append(pair)
        aggregate = [(sum(x for x, _ in group), sum(y for _, y in group), len(group)) for group in groups.values()]
        rng = random.Random(42)
        samples = []
        for _ in range(args.bootstrap):
            draw = [aggregate[rng.randrange(len(aggregate))] for _ in aggregate]
            samples.append(sum(a-b for a,b,n in draw) / sum(n for a,b,n in draw))
        samples.sort()
        total_a = sum(a for a,b,n in aggregate)
        total_b = sum(b for a,b,n in aggregate)
        result['recovery']['converged_only' if converged else 'any'] = {
            'chematic': total_a, 'rdkit': total_b, 'difference': (total_a-total_b)/len(rows),
            'scaffolds': len(aggregate), 'interval_95': [samples[int(.025*(len(samples)-1))], samples[int(.975*(len(samples)-1))]]}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
