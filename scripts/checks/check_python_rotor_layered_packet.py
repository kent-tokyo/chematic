#!/usr/bin/env python3
"""Recompute the bounded rotor/Layered timing packet and its adoption decision."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_python_speed_frontier import block_orders, comparable, MODES
from check_python_speed_frontier_packet import exact_interval

ARMS = ('baseline', 'candidate', 'rdkit', 'cosmolkit')
OPS = ('rotatable_bonds', 'mqn', 'fp_layered', 'qed', 'tpsa', 'labute_asa', 'morgan2_chiral', 'maccs')
PRIMARY = ('rotatable_bonds/first_use', 'mqn/first_use', 'fp_layered/parse_inclusive')


def check(packet):
    assert packet['schema'] == 'python-rotor-layered-numeric/v1'
    plan = packet['measurement_plan']
    assert json.loads(packet['measurement_plan_text']) == plan
    assert packet['measurement_plan_sha256'] == hashlib.sha256(
        packet['measurement_plan_text'].encode()).hexdigest()
    benchmark = packet['benchmark']
    manifest = benchmark['source_manifest']
    assert manifest['measurement_plan_sha256'] == packet['measurement_plan_sha256']
    assert manifest['build'] == packet['build']
    for arm in ('baseline', 'candidate'):
        source_files = packet['source_files'][arm]
        assert source_files and all(not name.startswith('/') and '..' not in name.split('/') for name in source_files)
        assert all(re.fullmatch('[0-9a-f]{64}', digest) for digest in source_files.values())
        assert hashlib.sha256(json.dumps(source_files, sort_keys=True).encode()).hexdigest() == packet['build'][arm]['source_sha256']
    assert plan['operations'] == list(OPS) and plan['modes'] == list(MODES)
    assert plan['primary_source_improvements'] == list(PRIMARY)
    assert plan['offset'] == 0 and plan['limit'] == 1000 and plan['blocks'] == 24
    assert benchmark['complete'] is True and len(benchmark['blocks']) == 24
    assert benchmark['base_revision'] == plan['source_baseline'] == packet['build']['baseline']['revision']
    assert packet['build']['candidate']['revision'].startswith(plan['candidate'])
    assert benchmark['corpus']['sha256'] == 'f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f'
    assert benchmark['corpus']['limit'] == 1000 and benchmark['corpus']['offset'] == 0
    expected_orders = block_orders(ARMS, 24, True)
    first = benchmark['blocks'][0]
    assert first['baseline']['output_sha256'] == first['candidate']['output_sha256']
    for index, block in enumerate(benchmark['blocks']):
        assert tuple(block['order']) == expected_orders[index]
        for arm in ARMS:
            sample = block[arm]
            for key in ('extension_sha256', 'version', 'output_sha256'):
                assert sample[key] == first[arm][key]
            assert re.fullmatch('[0-9a-f]{64}', sample['output_sha256'])
            assert sample['corpus_sha256'] == benchmark['corpus']['sha256']
            assert type(sample['parse_ns']) is int and sample['parse_ns'] > 0
            assert set(sample['operations']) == set(OPS)
            assert all(set(sample['operations'][op]) == set(MODES) for op in OPS)
            assert all(type(sample['operations'][op][mode]) is int and sample['operations'][op][mode] > 0
                       for op in OPS for mode in MODES)
            if arm in packet['build']:
                assert sample['extension_sha256'] == packet['build'][arm]['native_sha256']
    source = {}
    for op in OPS:
        checks = benchmark['agreement'][op]
        assert checks['candidate_vs_baseline'] == {'passed': True, 'rows': 1000, 'mismatches': 0}
        for mode in MODES:
            lane = op + '/' + mode
            candidate = [b['candidate']['operations'][op][mode] for b in benchmark['blocks']]
            for seed, arm in enumerate(ARMS[:1] + ARMS[2:]):
                numerator = [b[arm]['operations'][op][mode] for b in benchmark['blocks']]
                allowed = comparable(checks, arm)
                interval = exact_interval(numerator, candidate, seed) if allowed else None
                assert benchmark['summary'][lane][arm] == {
                    'speedup': interval,
                    'demonstrated_win': bool(interval and interval['ci95'][0] > 1 and interval['winning_blocks'] == 24),
                }
                if arm == 'baseline':
                    source[lane] = interval
    correctness = packet['correctness']
    assert correctness['dump_sha256']['baseline'] == correctness['dump_sha256']['candidate']
    assert correctness['source_exact_states'] == correctness['source_successful_values'] == 56288
    assert [h['rows'] for h in correctness['corpora']] == [10000, 4072]
    for arm in ('candidate', 'cosmolkit'):
        for op in ('rotatable_bonds', 'fp_layered', 'qed'):
            assert correctness['comparisons'][arm+'/rdkit-js-browser-10k-v1.smi/'+op]['counts'] == {'equal': 10000}
    for key, comparison in correctness['comparisons'].items():
        rows = 10000 if '/rdkit-js-browser-10k-v1.smi/' in key else 4072
        counts = comparison['counts']
        assert set(counts) <= {'equal', 'different', 'both_unavailable', 'source_unavailable', 'oracle_unavailable'}
        assert all(type(n) is int and n > 0 for n in counts.values())
        assert sum(counts.values()) == rows
        indices = comparison['non_equal_row_indices']
        assert indices == sorted(set(indices)) and all(0 <= n < rows for n in indices)
        assert len(indices) == rows - counts.get('equal', 0)
    for digest in correctness['dump_sha256'].values():
        assert re.fullmatch('[0-9a-f]{64}', digest)
    assert packet['source_adoption_passed'] == all(source[lane]['ci95'][0] > 1 for lane in PRIMARY)
    return {'source_adoption_passed': packet['source_adoption_passed'], 'comparisons': len(OPS)*len(MODES)*3,
            'source_intervals': source}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet', type=Path)
    parser.add_argument('--require-adoption', action='store_true')
    args = parser.parse_args()
    result = check(json.loads(gzip.decompress(args.packet.read_bytes())))
    if args.require_adoption:
        assert result['source_adoption_passed'], 'primary improvement not demonstrated'
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
