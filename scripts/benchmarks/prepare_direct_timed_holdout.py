#!/usr/bin/env python3
"""Prepare timed shards from a sealed audited population, without a fixed-count run."""
import argparse
import json
from pathlib import Path

import shard_platinum_holdout as fixed
import shard_platinum_time as timed


def prepare(split_path, audit_path, input_path, output_dir, count, wheel_sha256):
    split_raw, input_raw = split_path.read_bytes(), input_path.read_bytes()
    split, input_data = json.loads(split_raw), json.loads(input_raw)
    if (input_data.get('split_sha256') != fixed.digest(split_raw)
            or input_data.get('source_sha256') != split['source_sha256']
            or input_data.get('partition') != 'holdout'):
        raise ValueError('input population fingerprint mismatch')
    expected = {(r['index'], r['name']) for r in split['rows'] if r.get('split') == 'holdout'}
    records = input_data['results']
    actual = [(r['index'], r['name']) for r in records]
    if not expected or len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('input population incomplete, duplicated or unexpected')
    if any(not isinstance(r.get('smiles'), str) or not r['smiles'] or r.get('input_error') for r in records):
        raise ValueError('every selected record must have a parsed SMILES')
    selection_dir = output_dir / 'selection'
    selection = fixed.prepare(split_path, audit_path, selection_dir, count, wheel_sha256, 'mixed_checked')
    by_index = {r['index']: r for r in records}
    selection_raw = (selection_dir / 'manifest.json').read_bytes()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        'selection_manifest_sha256': fixed.digest(selection_raw),
        'baseline_sha256': fixed.digest(input_raw),
        'source_sha256': selection['source_sha256'],
        'wheel_sha256': wheel_sha256,
        'shard_count': count,
        'protocol': {'rdkit': '2026.09.2', 'seconds': 0.2,
                     'stereo_policy': 'rdkit-declared-stereo-strict-v1',
                     'max_candidates': 20, 'repair_iters': 400,
                     'geometry_checked': True,
                     'engine_error_policy': 'independent generation and scoring per engine',
                     'seed_schedule': '42 * (attempt + 1)',
                     'profile': 'mixed_checked',
                     'gradient_schedule': 'first floor(max_candidates / 2) compatible, remainder corrected'},
        'profile': 'mixed_checked',
        'selection': 'audited holdout; SHA256(name) round-robin shards, same order within each shard',
        'fingerprints': timed.fingerprints('mixed_checked'),
        'shards': [],
    }
    for shard in selection['shards']:
        subset = {'source_sha256': selection['source_sha256'],
                  'split_sha256': fixed.digest(split_raw),
                  'input_sha256': fixed.digest(input_raw),
                  'results': [by_index[r['index']] for r in shard['records']]}
        raw = fixed.encoded(subset)
        (output_dir / f"baseline-{shard['index']:03d}.json").write_bytes(raw)
        manifest['shards'].append({'index': shard['index'], 'baseline_sha256': fixed.digest(raw),
                                   'count': shard['count'], 'records': shard['records']})
    (output_dir / 'manifest.json').write_bytes(fixed.encoded(manifest))
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--split', required=True, type=Path)
    p.add_argument('--audit', required=True, type=Path)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--shards', required=True, type=int)
    p.add_argument('--wheel-sha256', required=True)
    a = p.parse_args()
    m = prepare(a.split, a.audit, a.input, a.output_dir, a.shards, a.wheel_sha256)
    print(f"Prepared {sum(s['count'] for s in m['shards'])} direct timed records")


if __name__ == '__main__':
    main()
