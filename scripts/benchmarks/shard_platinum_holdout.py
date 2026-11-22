#!/usr/bin/env python3
"""Prepare and verify deterministic shards without inspecting holdout outcomes."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_platinum_stereo_strict import summarize as summarize_strict_stereo  # noqa: E402


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def prepare(split_path, audit_path, output_dir, count, wheel_sha256):
    split_bytes = split_path.read_bytes()
    audit_bytes = audit_path.read_bytes()
    split, audit = json.loads(split_bytes), json.loads(audit_bytes)
    if (audit['split_sha256'] != digest(split_bytes)
            or audit['source_sha256'] != split['source_sha256']
            or audit['invalid_exclusion_count'] != len(split['invalid_exclusions'])
            or audit['overlaps']):
        raise ValueError('original invalid-exclusion audit does not match the split')
    if len(wheel_sha256) != 64 or any(c not in '0123456789abcdef' for c in wheel_sha256):
        raise ValueError('wheel SHA-256 must be 64 lowercase hex characters')
    rows = [r for r in split['rows'] if r.get('split') == 'holdout']
    rows.sort(key=lambda r: hashlib.sha256(r['name'].encode()).digest())
    indices = [r['index'] for r in rows]
    names = [r['name'] for r in rows]
    if not rows or len(set(indices)) != len(rows) or len(set(names)) != len(rows):
        raise ValueError('holdout must be nonempty and have unique indices and names')
    if any(not r.get('scaffold_sha256') for r in rows):
        raise ValueError('every holdout record must have a scaffold identity')
    if count < 1 or count > len(rows):
        raise ValueError('invalid shard count')
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {'original_split_sha256': digest(split_bytes),
                'original_audit_sha256': digest(audit_bytes),
                'source_sha256': split['source_sha256'],
                'wheel_sha256': wheel_sha256, 'shard_count': count,
                'sharder_sha256': digest(Path(__file__).read_bytes()),
                'scorer_sha256': digest(Path(__file__).with_name('platinum_development_probe.py').read_bytes()),
                'geometry_sha256': digest(Path(__file__).with_name('platinum_geometry.py').read_bytes()),
                'summarizer_sha256': digest(Path(__file__).with_name('summarize_platinum_holdout.py').read_bytes()),
                'summary_helper_sha256': digest(Path(__file__).with_name('summarize_platinum_probe.py').read_bytes()),
                'stereo_verifier_sha256': digest(Path(__file__).with_name('platinum_stereo_strict.py').read_bytes()),
                'stereo_wrapper_sha256': digest(Path(__file__).with_name('platinum_development_strict.py').read_bytes()),
                'strict_summary_sha256': digest(Path(__file__).with_name('summarize_platinum_stereo_strict.py').read_bytes()),
                'selection': 'SHA256(name) order, round-robin shards',
                'protocol': {'rdkit': '2026.09.2', 'profile': 'corrected', 'geometry_checked': True,
                             'stereo_policy': 'rdkit-declared-stereo-strict-v1',
                             'repair_iters': 400, 'candidate_limit': 20,
                             'mmff_iterations': 200, 'seed': 42}, 'shards': []}
    for shard_index in range(count):
        selected = rows[shard_index::count]
        selected_indices = {r['index'] for r in selected}
        shard_split = dict(split)
        shard_split['rows'] = [dict(r, split='deferred_holdout')
                               if r.get('split') == 'holdout' and r['index'] not in selected_indices
                               else r for r in split['rows']]
        shard_bytes = encoded(shard_split)
        shard_audit = dict(audit, split_sha256=digest(shard_bytes),
                           parent_split_sha256=digest(split_bytes),
                           parent_audit_sha256=digest(audit_bytes),
                           derivation='subset of the audited holdout')
        audit_data = encoded(shard_audit)
        (output_dir / f'split-{shard_index:03d}.json').write_bytes(shard_bytes)
        (output_dir / f'audit-{shard_index:03d}.json').write_bytes(audit_data)
        manifest['shards'].append({'index': shard_index, 'count': len(selected),
                                   'split_sha256': digest(shard_bytes),
                                   'audit_sha256': digest(audit_data),
                                   'records': [{'index': r['index'], 'name': r['name']} for r in selected]})
    (output_dir / 'manifest.json').write_bytes(encoded(manifest))
    return manifest


def seal(manifest_path, shard_index, wheel_path, audit_path, result_path):
    manifest = json.loads(manifest_path.read_bytes())
    if shard_index < 0 or shard_index >= manifest['shard_count']:
        raise ValueError('invalid shard index')
    shard = manifest['shards'][shard_index]
    if shard['index'] != shard_index:
        raise ValueError('manifest shard index mismatch')
    for key, path in (('sharder_sha256', Path(__file__)),
                      ('scorer_sha256', Path(__file__).with_name('platinum_development_probe.py')),
                      ('geometry_sha256', Path(__file__).with_name('platinum_geometry.py')),
                      ('stereo_verifier_sha256', Path(__file__).with_name('platinum_stereo_strict.py')),
                      ('stereo_wrapper_sha256', Path(__file__).with_name('platinum_development_strict.py')),
                      ('strict_summary_sha256', Path(__file__).with_name('summarize_platinum_stereo_strict.py'))):
        if digest(path.read_bytes()) != manifest[key]:
            raise ValueError(f'{key} mismatch')
    wheel_hash = digest(wheel_path.read_bytes())
    audit_hash = digest(audit_path.read_bytes())
    if wheel_hash != manifest['wheel_sha256'] or audit_hash != shard['audit_sha256']:
        raise ValueError('actual wheel or audit fingerprint mismatch')
    data = json.loads(result_path.read_bytes())
    if data.get('split_sha256') != shard['split_sha256']:
        raise ValueError('result split mismatch')
    if (data.get('stereo_policy') != manifest['protocol']['stereo_policy']
            or data.get('stereo_verifier_sha256') != manifest['stereo_verifier_sha256']
            or data.get('stereo_wrapper_sha256') != manifest['stereo_wrapper_sha256']):
        raise ValueError('strict stereo policy or fingerprint mismatch')
    expected = {(r['index'], r['name']) for r in shard['records']}
    actual = [(r['index'], r['name']) for r in data['results']]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('cannot seal an incomplete or incorrect shard')
    summarize_strict_stereo(data)
    data['wheel_sha256'] = wheel_hash
    data['invalid_exclusion_audit_sha256'] = audit_hash
    temporary = result_path.with_suffix(result_path.suffix + '.seal.tmp')
    temporary.write_bytes(encoded(data))
    temporary.replace(result_path)


def merge(manifest_path, shard_paths, output_path):
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    if digest(Path(__file__).read_bytes()) != manifest['sharder_sha256']:
        raise ValueError('sharder fingerprint mismatch')
    if digest(Path(__file__).with_name('summarize_platinum_holdout.py').read_bytes()) != manifest['summarizer_sha256']:
        raise ValueError('summarizer fingerprint mismatch')
    if digest(Path(__file__).with_name('summarize_platinum_probe.py').read_bytes()) != manifest['summary_helper_sha256']:
        raise ValueError('summary helper fingerprint mismatch')
    if digest(Path(__file__).with_name('summarize_platinum_stereo_strict.py').read_bytes()) != manifest['strict_summary_sha256']:
        raise ValueError('strict stereo summary fingerprint mismatch')
    shards = manifest['shards']
    if len(shards) != manifest['shard_count'] or len(shard_paths) != len(shards):
        raise ValueError('missing or extra shard output')
    by_split = {s['split_sha256']: s for s in shards}
    if len(by_split) != len(shards):
        raise ValueError('duplicate shard split hash')
    common_keys = ('rdkit', 'chematic', 'script_sha256', 'profile', 'geometry_checked',
                   'stereo_policy', 'stereo_verifier_sha256', 'stereo_wrapper_sha256',
                   'repair_iters', 'candidate_limit', 'mmff_iterations', 'seed', 'source_sha256')
    expected_common = None
    seen, combined, platforms = set(), [], []
    shard_hashes = []
    for path in shard_paths:
        raw = path.read_bytes()
        data = json.loads(raw)
        shard = by_split.get(data.get('split_sha256'))
        if shard is None or shard['index'] in seen:
            raise ValueError('unexpected or duplicate shard')
        seen.add(shard['index'])
        if (data.get('wheel_sha256') != manifest['wheel_sha256']
                or data.get('invalid_exclusion_audit_sha256') != shard['audit_sha256']):
            raise ValueError('wheel or shard audit fingerprint mismatch')
        common = {key: data[key] for key in common_keys}
        if (common['script_sha256'] != manifest['scorer_sha256']
                or common['stereo_verifier_sha256'] != manifest['stereo_verifier_sha256']
                or common['stereo_wrapper_sha256'] != manifest['stereo_wrapper_sha256']):
            raise ValueError('scorer or strict-stereo fingerprint mismatch')
        if common['source_sha256'] != manifest['source_sha256']:
            raise ValueError('source mismatch')
        if any(common[key] != value for key, value in manifest['protocol'].items()):
            raise ValueError('evaluation protocol mismatch')
        if expected_common is None:
            expected_common = common
        elif common != expected_common:
            raise ValueError('scorer, version or evaluation parameters differ')
        expected = {(r['index'], r['name']) for r in shard['records']}
        actual = [(r['index'], r['name']) for r in data['results']]
        if len(actual) != len(expected) or set(actual) != expected:
            raise ValueError('shard records are incomplete, duplicated or unexpected')
        summarize_strict_stereo(data)
        combined.extend(data['results'])
        platforms.append(data['platform'])
        shard_hashes.append({'index': shard['index'], 'result_sha256': digest(raw)})
    if seen != set(range(len(shards))):
        raise ValueError('shard indices are incomplete')
    combined.sort(key=lambda r: r['index'])
    result = dict(expected_common, split_sha256=manifest['original_split_sha256'],
                  platform=sorted(set(platforms)),
                  selection=manifest['selection'], results=combined,
                  shard_provenance={'manifest_sha256': digest(manifest_bytes),
                                    'wheel_sha256': manifest['wheel_sha256'],
                                    'shards': sorted(shard_hashes, key=lambda x: x['index'])})
    output_path.write_bytes(encoded(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('prepare')
    p.add_argument('--split', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--shards', type=int, required=True)
    p.add_argument('--wheel-sha256', required=True)
    p = commands.add_parser('seal')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--wheel', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True)
    p = commands.add_parser('merge')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('shard_results', type=Path, nargs='+')
    args = parser.parse_args()
    if args.command == 'prepare':
        result = prepare(args.split, args.audit, args.output_dir, args.shards, args.wheel_sha256)
        print(f"Prepared {sum(s['count'] for s in result['shards'])} records in {result['shard_count']} shards")
    elif args.command == 'seal':
        seal(args.manifest, args.index, args.wheel, args.audit, args.result)
        print(f'Sealed shard {args.index}')
    else:
        result = merge(args.manifest, args.shard_results, args.output)
        print(f"Verified and merged {len(result['results'])} records")


if __name__ == '__main__':
    main()
