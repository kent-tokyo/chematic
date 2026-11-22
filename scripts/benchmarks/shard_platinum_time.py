#!/usr/bin/env python3
"""Prepare and verify timed shards from an already verified fixed-count run."""
import argparse
import json
from pathlib import Path

from shard_platinum_holdout import digest, encoded

HERE = Path(__file__).resolve().parent


def fingerprints():
    return {name: digest((HERE / filename).read_bytes()) for name, filename in (
        ('timed_sharder_sha256', 'shard_platinum_time.py'),
        ('fixed_sharder_sha256', 'shard_platinum_holdout.py'),
        ('time_scorer_sha256', 'platinum_time_budget.py'),
        ('fixed_scorer_sha256', 'platinum_development_probe.py'),
        ('geometry_sha256', 'platinum_geometry.py'),
        ('time_summary_sha256', 'summarize_platinum_time_holdout.py'),
        ('time_summary_helper_sha256', 'summarize_platinum_time_budget.py'),
    )}


def prepare(fixed_manifest_path, baseline_path, output_dir):
    manifest_bytes = fixed_manifest_path.read_bytes()
    fixed = json.loads(manifest_bytes)
    baseline_bytes = baseline_path.read_bytes()
    baseline = json.loads(baseline_bytes)
    fingerprints_now = fingerprints()
    if (fingerprints_now['fixed_sharder_sha256'] != fixed['sharder_sha256']
            or fingerprints_now['fixed_scorer_sha256'] != fixed['scorer_sha256']
            or fingerprints_now['geometry_sha256'] != fixed['geometry_sha256']):
        raise ValueError('fixed scorer fingerprint mismatch')
    provenance = baseline.get('shard_provenance', {})
    if (provenance.get('manifest_sha256') != digest(manifest_bytes)
            or provenance.get('wheel_sha256') != fixed['wheel_sha256']
            or baseline.get('split_sha256') != fixed['original_split_sha256']
            or baseline.get('source_sha256') != fixed['source_sha256']):
        raise ValueError('baseline is not the verified fixed-count result')
    if (baseline.get('script_sha256') != fixed['scorer_sha256']
            or any(baseline.get(k) != v for k, v in fixed['protocol'].items())):
        raise ValueError('fixed-count protocol mismatch')
    expected = {(r['index'], r['name']) for shard in fixed['shards'] for r in shard['records']}
    actual = [(r['index'], r['name']) for r in baseline['results']]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('fixed-count baseline is incomplete, duplicated or unexpected')
    by_index = {r['index']: r for r in baseline['results']}
    output_dir.mkdir(parents=True, exist_ok=True)
    timed = {'fixed_manifest_sha256': digest(manifest_bytes),
             'baseline_sha256': digest(baseline_bytes),
             'source_sha256': fixed['source_sha256'],
             'wheel_sha256': fixed['wheel_sha256'],
             'shard_count': fixed['shard_count'],
             'protocol': {'rdkit': '2026.09.2', 'seconds': 0.2,
                          'max_candidates': 20, 'repair_iters': 400,
                          'geometry_checked': True,
                          'seed_schedule': '42 * (attempt + 1)'},
             'selection': 'fixed-count shard membership; SHA256(name) order within each shard',
             'fingerprints': fingerprints_now, 'shards': []}
    for shard in fixed['shards']:
        subset = dict(baseline)
        subset['results'] = [by_index[r['index']] for r in shard['records']]
        subset_bytes = encoded(subset)
        (output_dir / f"baseline-{shard['index']:03d}.json").write_bytes(subset_bytes)
        timed['shards'].append({'index': shard['index'], 'baseline_sha256': digest(subset_bytes),
                                'count': shard['count'], 'records': shard['records']})
    (output_dir / 'manifest.json').write_bytes(encoded(timed))
    return timed


def seal(manifest_path, index, wheel_path, baseline_path, result_path):
    manifest = json.loads(manifest_path.read_bytes())
    if fingerprints() != manifest['fingerprints']:
        raise ValueError('scorer or helper fingerprint mismatch')
    if index < 0 or index >= manifest['shard_count']:
        raise ValueError('invalid shard index')
    shard = manifest['shards'][index]
    if shard['index'] != index:
        raise ValueError('manifest shard index mismatch')
    wheel_hash = digest(wheel_path.read_bytes())
    baseline_hash = digest(baseline_path.read_bytes())
    if wheel_hash != manifest['wheel_sha256'] or baseline_hash != shard['baseline_sha256']:
        raise ValueError('wheel or baseline fingerprint mismatch')
    data = json.loads(result_path.read_bytes())
    if (data.get('source_sha256') != manifest['source_sha256']
            or data.get('script_sha256') != manifest['fingerprints']['time_scorer_sha256']
            or any(data.get(k) != v for k, v in manifest['protocol'].items())):
        raise ValueError('timed scoring protocol mismatch')
    expected = {(r['index'], r['name']) for r in shard['records']}
    actual = [(r['index'], r['name']) for r in data['results']]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('cannot seal incomplete or incorrect timed shard')
    data['wheel_sha256'] = wheel_hash
    data['baseline_sha256'] = baseline_hash
    temporary = result_path.with_suffix(result_path.suffix + '.seal.tmp')
    temporary.write_bytes(encoded(data))
    temporary.replace(result_path)


def merge(manifest_path, shard_paths, output_path):
    raw_manifest = manifest_path.read_bytes()
    manifest = json.loads(raw_manifest)
    if fingerprints() != manifest['fingerprints']:
        raise ValueError('scorer or helper fingerprint mismatch')
    shards = manifest['shards']
    if len(shards) != manifest['shard_count'] or len(shard_paths) != len(shards):
        raise ValueError('missing or extra timed shard')
    by_baseline = {s['baseline_sha256']: s for s in shards}
    if len(by_baseline) != len(shards):
        raise ValueError('duplicate baseline hash')
    keys = ('rdkit', 'chematic', 'script_sha256', 'source_sha256', 'seconds',
            'max_candidates', 'repair_iters', 'geometry_checked', 'seed_schedule', 'timing')
    common, seen, results, platforms, hashes = None, set(), [], [], []
    for path in shard_paths:
        raw = path.read_bytes()
        data = json.loads(raw)
        shard = by_baseline.get(data.get('baseline_sha256'))
        if shard is None or shard['index'] in seen or data.get('wheel_sha256') != manifest['wheel_sha256']:
            raise ValueError('unexpected, duplicate or mixed-wheel timed shard')
        seen.add(shard['index'])
        current = {key: data[key] for key in keys}
        if (current['source_sha256'] != manifest['source_sha256']
                or current['script_sha256'] != manifest['fingerprints']['time_scorer_sha256']
                or any(current[k] != v for k, v in manifest['protocol'].items())):
            raise ValueError('timed scoring protocol mismatch')
        if common is None:
            common = current
        elif common != current:
            raise ValueError('timed shard versions or settings differ')
        expected = {(r['index'], r['name']) for r in shard['records']}
        actual = [(r['index'], r['name']) for r in data['results']]
        if len(actual) != len(expected) or set(actual) != expected:
            raise ValueError('timed shard records are incomplete, duplicated or unexpected')
        results.extend(data['results'])
        platforms.append(data['platform'])
        hashes.append({'index': shard['index'], 'result_sha256': digest(raw)})
    if seen != set(range(len(shards))):
        raise ValueError('timed shard indices are incomplete')
    results.sort(key=lambda r: r['index'])
    output = dict(common, platform=sorted(set(platforms)), results=results,
                  shard_provenance={'manifest_sha256': digest(raw_manifest),
                                    'baseline_sha256': manifest['baseline_sha256'],
                                    'wheel_sha256': manifest['wheel_sha256'],
                                    'shards': sorted(hashes, key=lambda x: x['index'])})
    output_path.write_bytes(encoded(output))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--fixed-manifest', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p = sub.add_parser('seal')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--wheel', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True)
    p = sub.add_parser('merge')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('shard_results', type=Path, nargs='+')
    args = parser.parse_args()
    if args.command == 'prepare':
        m = prepare(args.fixed_manifest, args.baseline, args.output_dir)
        print(f"Prepared {sum(s['count'] for s in m['shards'])} timed records")
    elif args.command == 'seal':
        seal(args.manifest, args.index, args.wheel, args.baseline, args.result)
        print(f'Sealed timed shard {args.index}')
    else:
        result = merge(args.manifest, args.shard_results, args.output)
        print(f"Verified and merged {len(result['results'])} timed records")


if __name__ == '__main__':
    main()
