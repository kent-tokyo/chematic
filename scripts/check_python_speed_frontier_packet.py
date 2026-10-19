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


def check_rejected_pilot(packet):
    """Verify a stopped diagnostic without treating it as formal evidence."""
    assert packet['source_adoption_passed'] is False
    assert packet['proceeded_to_formal'] is False
    plan = packet['measurement_plan']
    assert hashlib.sha256((json.dumps(plan, indent=2)+'\n').encode()).hexdigest() == packet['measurement_plan_sha256']
    assert plan['python_hash_seed'] == 20261019
    assert plan['primary'] == ['full-10000/morgan2_chiral/first_use',
                               'full-10000/tpsa/first_use',
                               'full-10000/num_rings/repeated_use']
    pilot = packet['pilot']
    assert pilot['diagnostic_only'] is True and pilot['not_adoption_evidence'] is True
    assert len(pilot['blocks']) == 8
    assert sum(b['order'] == ['baseline','candidate'] for b in pilot['blocks']) == 4
    assert sum(b['order'] == ['candidate','baseline'] for b in pilot['blocks']) == 4
    ops = ('tpsa','labute_asa','num_rings','morgan2_2048','morgan2_chiral','morgan3_2048','maccs')
    modes = ('first_use','repeated_use','parse_inclusive')
    for block in pilot['blocks']:
        for arm in ('baseline','candidate'):
            data = block[arm]
            assert data['extension_sha256'] == {Path(k).name:v for k,v in packet['build'][arm]['native_sha256'].items()}
            assert data['output_sha256'] == pilot['blocks'][0]['baseline']['output_sha256']
            assert set(data['operations']) == set(ops)
            for op in ops:
                assert set(data['operations'][op]) == set(modes)
                assert all(type(data['operations'][op][m]) is int and data['operations'][op][m] > 0 for m in modes)
    expected = {op+'/'+mode: statistics.median(
        b['baseline']['operations'][op][mode]/b['candidate']['operations'][op][mode]
        for b in pilot['blocks']) for op in ops for mode in modes}
    assert pilot['source_medians'] == expected
    assert any(value < 0.99 for value in expected.values())
    assert packet['validation']['source_and_rdkit_exact_cells'] == 56288
    assert packet['scalar_validation']['source_exact_cells_including_refusals'] == 42216
    return {'source_adoption_passed':False,'proceeded_to_formal':False,'pilot_blocks':8,'source_lanes':21}


def check(packet):
    if packet['schema'] == 'python-speed-frontier-pilot/v1':
        return check_rejected_pilot(packet)
    assert packet['schema'] in ('python-speed-frontier-numeric/v1', 'python-speed-frontier-numeric/v2', 'python-speed-frontier-numeric/v3')
    report = packet['benchmark']
    plan = packet['measurement_plan']
    assert hashlib.sha256((json.dumps(plan,indent=2)+'\n').encode()).hexdigest() == packet['measurement_plan_sha256']
    primary_lanes = ['morgan2_chiral/first_use'] if packet['schema'].endswith('/v1') else ['morgan2_chiral/first_use','tpsa/first_use']
    if packet['schema'].endswith('/v3'):
        primary_lanes.append('num_rings/repeated_use')
    assert plan['primary'] == ['full-10000/'+lane for lane in primary_lanes]
    if not packet['schema'].endswith('/v1'):
        assert plan['python_hash_seed'] == 20261019
        assert report['source_manifest']['python_hash_seed'] == 20261019
    measured = report['source_manifest']
    assert measured['measurement_plan_sha256'] == packet['measurement_plan_sha256']
    for arm in ('baseline','candidate'):
        assert measured[arm+'_core_source_sha256'] == packet['build'][arm]['source_sha256']
    assert report['protocol']['blocks'] == 24
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
    normal_header = next(h for h in packet['validation']['headers']
                         if h['corpus'] == 'rdkit-js-browser-10k-v1.smi')
    assert normal_header['rows'] == 10000
    assert normal_header['sha256'] == report['corpus']['sha256']
    assert packet['validation']['rows'] == 14072
    assert packet['validation']['source_and_rdkit_exact_cells'] == 56288
    assert packet['validation']['refusals'] == 0
    adopted = all(source_gate(exact[lane]['baseline']['speedup'], source) for lane in primary_lanes)
    if not packet['schema'].endswith('/v1'):
        scalar = packet['scalar_validation']
        assert scalar['source_rows_equal'] is True
        assert scalar['source_exact_cells_including_refusals'] == 42216
        assert scalar['headers'] == packet['validation']['headers']
        corpora = {'rdkit-js-browser-10k-v1.smi': 10000, 'stress-v1.smi': 4072}
        assert {h['corpus']:h['rows'] for h in scalar['headers']} == corpora
        expected = {f'{arm}/{corpus}/{op}' for arm in ('candidate','cosmolkit')
                    for corpus in corpora for op in ('tpsa','labute_asa','num_rings')}
        assert set(scalar['comparison']) == expected
        states = {'equal','mismatch','both_unavailable','source_unavailable','rdkit_unavailable'}
        for key, counts in scalar['comparison'].items():
            assert counts and set(counts) <= states
            assert all(type(v) is int and v >= 0 for v in counts.values())
            assert sum(counts.values()) == corpora[key.split('/')[1]]
        for corpus, rows in corpora.items():
            assert scalar['comparison'][f'candidate/{corpus}/tpsa'] == {'equal': rows}
        pilot = packet['pilot']
        assert pilot['diagnostic_only'] is True and pilot['not_adoption_evidence'] is True
        assert len(pilot['blocks']) == 8
        assert sum(v['order'] == ['baseline','candidate'] for v in pilot['blocks']) == 4
        assert sum(v['order'] == ['candidate','baseline'] for v in pilot['blocks']) == 4
        pilot_ops = ops if packet['schema'].endswith('/v3') else ('tpsa','morgan2_chiral')
        for block in pilot['blocks']:
            for arm in ('baseline','candidate'):
                data = block[arm]
                assert data['extension_sha256'] == first[arm]['extension_sha256']
                assert data['output_sha256'] == pilot['blocks'][0]['baseline']['output_sha256']
                assert set(data['operations']) == set(pilot_ops)
                for op in pilot_ops:
                    assert set(data['operations'][op]) == set(modes)
                    assert all(type(data['operations'][op][m]) is int and data['operations'][op][m] > 0 for m in modes)
        expected_pilot = {op+'/'+mode: statistics.median(
            v['baseline']['operations'][op][mode] / v['candidate']['operations'][op][mode]
            for v in pilot['blocks']) for op in pilot_ops for mode in modes}
        assert pilot['source_medians'] == expected_pilot
    assert packet['source_adoption_passed'] is adopted
    return {'source_adoption_passed':adopted,'blocks':24,'source_lanes':21,'comparisons':63,'exact_cells':56288}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet',type=Path)
    parser.add_argument('--require-adoption',action='store_true',help='fail unless the measured source candidate passes its registered gates')
    args=parser.parse_args()
    raw=args.packet.read_bytes()
    if args.packet.suffix == '.gz':raw=gzip.decompress(raw)
    result = check(json.loads(raw))
    print(json.dumps(result,indent=2))
    if args.require_adoption and not result['source_adoption_passed']:
        parser.error('source adoption gate failed')


if __name__ == '__main__':main()
