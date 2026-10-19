#!/usr/bin/env python3
"""Recompute the 24-order source experiment, including every control and non-win."""
import argparse
import gzip
import hashlib
import itertools
import json
import random
import statistics
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_python_speed_frontier import comparable


def exact_interval(numerator, denominator, seed):
    assert len(numerator) == len(denominator) == 24
    ratios = [n / d for n, d in zip(numerator, denominator, strict=True)]
    rng = random.Random(seed)
    estimates = sorted(statistics.median(ratios[rng.randrange(24)] for _ in ratios)
                       for _ in range(10000))
    return {'median': statistics.median(ratios), 'ci95': [estimates[249], estimates[9749]],
            'winning_blocks': sum(r > 1 for r in ratios), 'blocks': 24}


def precise_summary(report):
    output = {}
    for lane in report['summary']:
        op, mode = lane.split('/')
        comparisons = {}
        for seed, name in enumerate(('baseline','rdkit','cosmolkit')):
            interval = exact_interval(
                [b[name]['operations'][op][mode] for b in report['blocks']],
                [b['candidate']['operations'][op][mode] for b in report['blocks']], seed
            ) if comparable(report['agreement'][op], name) else None
            comparisons[name] = {'speedup': interval, 'demonstrated_win': bool(
                interval and interval['ci95'][0] > 1 and interval['winning_blocks'] == 24)}
        output[lane] = comparisons
    return output


def source_gate(primary, controls):
    return bool(primary and primary['ci95'][0] > 1
                and all(r and r['median'] >= 0.99 for r in controls))


def check(packet):
    assert packet['schema'] == 'python-speed-frontier-numeric/v1'
    report = packet['benchmark']
    plan = packet['measurement_plan']
    assert hashlib.sha256((json.dumps(plan,indent=2)+'\n').encode()).hexdigest() == packet['measurement_plan_sha256']
    assert plan['primary'] == ['full-10000/morgan2_chiral/first_use']
    assert plan['base'] == report['base_revision']
    assert report['complete'] is True
    assert report['corpus']['limit'] == 10000 and report['corpus']['offset'] == 0
    names = ('baseline', 'candidate', 'rdkit', 'cosmolkit')
    ops = ('tpsa', 'labute_asa', 'num_rings', 'morgan2_2048', 'morgan2_chiral', 'morgan3_2048', 'maccs')
    modes = ('first_use', 'repeated_use', 'parse_inclusive')
    blocks = report['blocks']
    assert len(blocks) == 24
    assert {tuple(b['order']) for b in blocks} == set(itertools.permutations(names))
    first = blocks[0]
    for arm in ('baseline','candidate'):
        assert {Path(k).name:v for k,v in packet['build'][arm]['native_sha256'].items()} == first[arm]['extension_sha256']
    for b in blocks:
        for name in names:
            arm = b[name]
            assert arm['extension_sha256'] == first[name]['extension_sha256']
            assert arm['output_sha256'] == first[name]['output_sha256']
            assert arm['version'] == first[name]['version']
            assert type(arm['parse_ns']) is int and arm['parse_ns'] > 0
            assert set(arm['operations']) == set(ops)
            for op in ops:
                assert set(arm['operations'][op]) == set(modes)
                assert all(type(arm['operations'][op][m]) is int and arm['operations'][op][m] > 0 for m in modes)
    assert set(report['summary']) == {op+'/'+mode for op in ops for mode in modes}
    assert set(report['agreement']) == set(ops)
    source = []
    exact = precise_summary(report)
    assert packet['full_precision_summary'] == exact
    for op in ops:
        gates = report['agreement'][op]
        assert set(gates) == {'candidate_vs_baseline','baseline_vs_rdkit','candidate_vs_rdkit','cosmolkit_vs_rdkit'}
        for v in gates.values():
            assert v['rows'] == 10000 and type(v['mismatches']) is int and 0 <= v['mismatches'] <= 10000
            assert type(v['passed']) is bool and v['passed'] == (v['mismatches'] == 0)
        for mode in modes:
            saved = report['summary'][op+'/'+mode]
            assert set(saved) == {'baseline','rdkit','cosmolkit'}
            for seed, name in enumerate(('baseline','rdkit','cosmolkit')):
                interval = exact[op+'/'+mode][name]['speedup']
                rounded = ({**interval, 'median': round(interval['median'], 3), 'ci95': [round(v,3) for v in interval['ci95']]} if interval else None)
                assert saved[name] == {'speedup':rounded,'demonstrated_win': bool(rounded and rounded['ci95'][0]>1 and rounded['winning_blocks']==24)}
                if name == 'baseline':source.append(interval)
    assert packet['validation']['rows'] == 14072
    assert packet['validation']['source_and_rdkit_exact_cells'] == 56288
    assert packet['validation']['refusals'] == 0
    primary = exact['morgan2_chiral/first_use']['baseline']['speedup']
    adopted = source_gate(primary, source)
    assert packet['source_adoption_passed'] is adopted
    return {'source_adoption_passed':adopted,'blocks':24,'source_lanes':21,'comparisons':63,'exact_cells':56288}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet',type=Path)
    args=parser.parse_args()
    raw=args.packet.read_bytes()
    if args.packet.suffix == '.gz':raw=gzip.decompress(raw)
    print(json.dumps(check(json.loads(raw)),indent=2))


if __name__ == '__main__':main()
