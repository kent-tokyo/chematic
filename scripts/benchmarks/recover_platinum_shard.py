#!/usr/bin/env python3
"""Recover a checkpointed fixed-count shard without changing its scorer.

Run with the frozen scorer and sharder from the original manifest. The remaining
records are evaluated as deterministic audited subshards, then joined to the
checkpoint after all subshards pass the original merge and stereo checks.
"""
import argparse
import hashlib
import json
from pathlib import Path

from shard_platinum_holdout import encoded, digest, merge, prepare, seal
from platinum_stereo_strict import annotate_result


def context(manifest_path, index, split_path, audit_path, checkpoint_path):
    manifest = json.loads(manifest_path.read_bytes())
    shard = manifest['shards'][index]
    split_raw, audit_raw = split_path.read_bytes(), audit_path.read_bytes()
    if digest(split_raw) != shard['split_sha256'] or digest(audit_raw) != shard['audit_sha256']:
        raise ValueError('original split or audit fingerprint mismatch')
    split, audit = json.loads(split_raw), json.loads(audit_raw)
    raw = checkpoint_path.read_bytes()
    checkpoint = json.loads(raw)
    expected = sorted((r['index'], r['name']) for r in shard['records'])
    completed = [(r['index'], r['name']) for r in checkpoint['results']]
    if not completed or len(completed) >= len(expected) or completed != expected[:len(completed)]:
        raise ValueError('checkpoint is not an incomplete prefix of this shard')
    checks = {'split_sha256': shard['split_sha256'], 'source_sha256': manifest['source_sha256'],
              'script_sha256': manifest['scorer_sha256'], 'rdkit': manifest['protocol']['rdkit'],
              'profile': 'corrected', 'geometry_checked': True, 'repair_iters': 400,
              'candidate_limit': 20, 'mmff_iterations': 200, 'seed': 42}
    if any(checkpoint.get(k) != v for k, v in checks.items()):
        raise ValueError('checkpoint metadata mismatch')
    return manifest, shard, split, audit, checkpoint, completed, digest(raw)


def make(args):
    manifest, _, split, audit, _, completed, checkpoint_hash = context(
        args.manifest, args.index, args.split, args.audit, args.checkpoint)
    done = {i for i, _ in completed}
    subset = dict(split)
    subset['rows'] = [dict(r, split='deferred_holdout') if r.get('split') == 'holdout' and r['index'] in done else r for r in split['rows']]
    args.directory.mkdir(parents=True, exist_ok=True)
    subset_path, audit_path = args.directory / 'remaining-split.json', args.directory / 'remaining-audit.json'
    subset_path.write_bytes(encoded(subset))
    audit = dict(audit, split_sha256=digest(subset_path.read_bytes()),
                 parent_split_sha256=digest(args.split.read_bytes()),
                 parent_audit_sha256=digest(args.audit.read_bytes()),
                 checkpoint_sha256=checkpoint_hash,
                 derivation='unscored suffix of the audited original shard')
    audit_path.write_bytes(encoded(audit))
    prepared = prepare(subset_path, audit_path, args.directory / 'plan', args.parts,
                       manifest['wheel_sha256'], 'mixed_checked')
    print(f"Preserved {len(completed)} checkpoint rows; planned {sum(s['count'] for s in prepared['shards'])} remaining rows")


def finish(args):
    manifest, shard, original_split, original_audit, checkpoint, completed, checkpoint_hash = context(
        args.manifest, args.index, args.split, args.audit, args.checkpoint)
    remaining_split = json.loads((args.directory / 'remaining-split.json').read_bytes())
    remaining_audit = json.loads((args.directory / 'remaining-audit.json').read_bytes())
    done = {index for index, _ in completed}
    expected_split = dict(original_split)
    expected_split['rows'] = [dict(r, split='deferred_holdout') if r.get('split') == 'holdout' and r['index'] in done else r for r in original_split['rows']]
    if remaining_split != expected_split:
        raise ValueError('remaining split is not the checkpoint suffix')
    expected_audit = dict(original_audit, split_sha256=digest(encoded(expected_split)),
                          parent_split_sha256=digest(args.split.read_bytes()),
                          parent_audit_sha256=digest(args.audit.read_bytes()),
                          checkpoint_sha256=checkpoint_hash,
                          derivation='unscored suffix of the audited original shard')
    if remaining_audit != expected_audit:
        raise ValueError('remaining audit provenance mismatch')
    remaining_path = args.directory / 'remaining-merged.json'
    remaining = merge(args.directory / 'plan/manifest.json', args.results, remaining_path)
    if remaining['split_sha256'] != digest((args.directory / 'remaining-split.json').read_bytes()):
        raise ValueError('remaining split mismatch')
    expected = {(r['index'], r['name']) for r in shard['records']}
    actual = completed + [(r['index'], r['name']) for r in remaining['results']]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('checkpoint and remaining rows do not partition original shard')
    for key in ('rdkit', 'chematic', 'script_sha256', 'source_sha256', 'geometry_checked',
                'repair_iters', 'candidate_limit', 'mmff_iterations', 'seed'):
        if checkpoint[key] != remaining[key]:
            raise ValueError(f'checkpoint and remaining metadata differ: {key}')
    checkpoint['results'] += remaining['results']
    checkpoint['results'].sort(key=lambda r: r['index'])
    checkpoint['profile'] = 'mixed_checked'
    checkpoint['recovery_provenance'] = {
        'checkpoint_sha256': checkpoint_hash,
        'checkpoint_platform': checkpoint['platform'],
        'remaining_platforms': remaining['platform'],
        'remaining_result_sha256': digest(remaining_path.read_bytes()),
        'remaining_shards': remaining['shard_provenance'],
    }
    args.output.write_bytes(encoded(checkpoint))
    annotate_result(args.output, Path(__file__).with_name('platinum_mixed_checked_development.py'))
    seal(args.manifest, args.index, args.wheel, args.audit, args.output)
    print(f"Recovered and sealed {len(actual)} rows")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'finish'))
    for name in ('manifest', 'split', 'audit', 'checkpoint', 'directory'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--parts', type=int, default=8)
    parser.add_argument('--wheel', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('results', type=Path, nargs='*')
    args = parser.parse_args()
    if args.action == 'prepare':
        make(args)
    else:
        if not args.wheel or not args.output or not args.results:
            parser.error('finish requires --wheel, --output, and subshard results')
        finish(args)


if __name__ == '__main__':
    main()
