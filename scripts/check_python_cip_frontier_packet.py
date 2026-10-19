#!/usr/bin/env python3
"""Verify the retained CIP measurements, input identities and speed summaries."""
import argparse
from collections import Counter
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import re
import statistics

ROOT = Path(__file__).resolve().parents[1]
OPS = ('cip', 'canonical_smiles', 'smarts_write', 'qed', 'morgan2_chiral', 'maccs', 'fp_layered')
MODES = ('first_use', 'repeated_use', 'parse_inclusive')
ARMS = ('baseline', 'candidate', 'rdkit')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def check(packet, *, verify_current_source=False):
    report = packet['measurement']
    plan = report['plan']
    assert report['complete'] is True
    assert plan['schema'] == 'python-cip-frontier/v1'
    assert plan['ops'] == list(OPS) and plan['modes'] == list(MODES)
    assert plan['rows'] == 1000 and plan['offset'] == 0
    expected_orders = list(itertools.permutations(ARMS)) * 4
    random.Random(20261019).shuffle(expected_orders)
    assert plan['orders'] == [list(x) for x in expected_orders]
    assert len(report['blocks']) == 24
    assert Counter(tuple(b['order']) for b in report['blocks']) == Counter(expected_orders)
    assert report['plan_sha256'] == digest((json.dumps(plan, indent=2)+'\n').encode())
    corpus = ROOT / 'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'
    assert plan['corpus_sha256'] == digest(corpus.read_bytes())
    for name, value in plan['script_sha256'].items():
        assert not name.startswith('/') and '..' not in name.split('/')
        assert re.fullmatch('[0-9a-f]{64}', value)
        if verify_current_source:
            assert digest((ROOT / name).read_bytes()) == value, name
    for arm in ARMS[:2]:
        identity = plan['identities'][arm]
        assert digest(json.dumps(identity['source_files'], sort_keys=True).encode()) == identity['source_sha256']
        assert identity['wheel_sha256'] == plan['wheel_sha256'][arm]
        assert identity['flags'] == {'CARGO_BUILD_JOBS': '1', 'CARGO_INCREMENTAL': '0',
                                     'CARGO_PROFILE_RELEASE_DEBUG': '0', 'CARGO_PROFILE_RELEASE_STRIP': 'none'}
        for name, value in identity['source_files'].items():
            assert not name.startswith('/') and '..' not in name.split('/')
            assert re.fullmatch('[0-9a-f]{64}', value)
            if arm == 'candidate' and verify_current_source:
                assert digest((ROOT / name).read_bytes()) == value, name
    for i, block in enumerate(report['blocks']):
        assert block['order'] == list(expected_orders[i])
        for arm in ARMS:
            sample = block[arm]
            assert sample['corpus_sha256'] == plan['corpus_sha256']
            assert set(sample['operations']) == set(OPS)
            for op in OPS:
                assert set(sample['operations'][op]) == set(MODES)
                assert all(type(x) is int and x > 0 for x in sample['operations'][op].values())
            assert sample['output_sha256'] == report['blocks'][0][arm]['output_sha256']
            assert sample['extension_sha256'] == report['blocks'][0][arm]['extension_sha256']
            assert sample['version'] == report['blocks'][0][arm]['version']
            if arm != 'rdkit':
                assert sample['extension_sha256'] == plan['identities'][arm]['native_sha256']
            else:
                assert 'rdkit.Chem.rdCIPLabeler' in sample['extension_sha256']
    screening = packet['screening']
    assert screening['diagnostic_only'] is True and screening['limit'] == 1000
    assert len(screening['blocks']) == 3
    for block in screening['blocks']:
        for engine in ('chematic', 'rdkit'):
            assert block[engine]['parsed'] == 1000
            for op in screening['ops']:
                assert math.isfinite(block[engine][op]) and block[engine][op] > 0
                assert block[engine].get(op+'__errors', 0) == 0
    correctness = packet['correctness']
    assert correctness['source_exact_rows'] == 14072 and correctness['source_changed_rows'] == 0
    assert correctness['normal_rdkit_raw_matches'] == correctness['normal_rows'] == 10000
    assert correctness['stress_rows'] == 4072
    assert correctness['stress_rdkit_raw_matches'] == 4014
    assert correctness['stress_rdkit_model_matches'] == 4070
    assert correctness['stress_ranked'] is False
    assert correctness['raw_file_sha256']['baseline'] == correctness['raw_file_sha256']['candidate']
    for name, value in correctness['corpus_sha256'].items():
        assert digest((ROOT / name).read_bytes()) == value
    for op in OPS:
        for arm in ('baseline', 'rdkit'):
            gate = report['agreement'][op][arm]
            assert gate['rows'] == 1000
            assert gate['passed'] is (gate['mismatches'] == 0)
            for mode in MODES:
                saved = report['summary'][op+'/'+mode][arm]
                if not gate['passed']:
                    assert saved is None
                    continue
                ratios = [b[arm]['operations'][op][mode] / b['candidate']['operations'][op][mode]
                          for b in report['blocks']]
                rng = random.Random(20261019)
                estimates = sorted(statistics.median(rng.choices(ratios, k=24)) for _ in range(10000))
                expected = {'median': statistics.median(ratios), 'ci95': [estimates[249], estimates[9749]],
                            'blocks': 24, 'winning_blocks': sum(r > 1 for r in ratios)}
                assert saved == expected, (op, mode, arm)
                assert all(math.isfinite(x) for x in ratios)
    adoption = all(report['summary']['cip/'+mode]['baseline']['ci95'][0] > 1
                   for mode in ('first_use', 'parse_inclusive'))
    assert adoption is True and packet['source_adoption_passed'] is True
    return len(report['blocks'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('packet', type=Path)
    p.add_argument('--verify-current-source', action='store_true',
                   help='also compare live source files; omit for immutable historical packets')
    a = p.parse_args()
    packet = json.loads(gzip.decompress(a.packet.read_bytes()))
    print('verified', check(packet, verify_current_source=a.verify_current_source), 'balanced blocks; all seven operations and three timing modes')


if __name__ == '__main__':
    main()
