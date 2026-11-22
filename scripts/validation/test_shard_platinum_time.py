"""Verify provenance, completeness and policy for timed holdout shards."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / 'benchmarks'
sys.path.insert(0, str(BENCH))
import shard_platinum_holdout as fixed  # noqa: E402
import shard_platinum_time as timed  # noqa: E402


class TimedShardingTest(unittest.TestCase):
    def test_verified_baseline_and_timed_shards(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split = {'source_sha256': 'a' * 64, 'invalid_exclusions': [],
                     'rows': [{'index': i, 'name': f'mol-{i}', 'split': 'holdout',
                               'scaffold_sha256': f'scaffold-{i}'} for i in range(3)]}
            split_path = root / 'split.json'
            split_path.write_bytes(fixed.encoded(split))
            audit = {'source_sha256': split['source_sha256'],
                     'split_sha256': fixed.digest(split_path.read_bytes()),
                     'invalid_exclusion_count': 0, 'overlaps': []}
            audit_path = root / 'audit.json'
            audit_path.write_bytes(fixed.encoded(audit))
            wheel = root / 'candidate.whl'
            wheel.write_bytes(b'synthetic wheel')
            wheel_hash = fixed.digest(wheel.read_bytes())
            fixed_dir = root / 'fixed'
            m = fixed.prepare(split_path, audit_path, fixed_dir, 2, wheel_hash)
            fixed_paths = []
            for shard in m['shards']:
                result = dict(m['protocol'], chematic='test', script_sha256=m['scorer_sha256'],
                              source_sha256=m['source_sha256'], split_sha256=shard['split_sha256'],
                              platform='test', results=[dict(r, status='ok') for r in shard['records']])
                path = root / f"fixed-result-{shard['index']}.json"
                path.write_bytes(fixed.encoded(result))
                fixed.seal(fixed_dir / 'manifest.json', shard['index'], wheel,
                           fixed_dir / f"audit-{shard['index']:03d}.json", path)
                fixed_paths.append(path)
            baseline = root / 'baseline.json'
            fixed.merge(fixed_dir / 'manifest.json', fixed_paths, baseline)
            timed_dir = root / 'timed'
            tm = timed.prepare(fixed_dir / 'manifest.json', baseline, timed_dir)
            results = []
            for shard in tm['shards']:
                data = dict(tm['protocol'], chematic='test', source_sha256=tm['source_sha256'],
                            script_sha256=tm['fingerprints']['time_scorer_sha256'],
                            platform='test', timing='same timed policy',
                            results=[dict(r) for r in shard['records']])
                path = root / f"timed-result-{shard['index']}.json"
                path.write_bytes(fixed.encoded(data))
                timed.seal(timed_dir / 'manifest.json', shard['index'], wheel,
                           timed_dir / f"baseline-{shard['index']:03d}.json", path)
                results.append(path)
            merged = timed.merge(timed_dir / 'manifest.json', results, root / 'merged.json')
            self.assertEqual([0, 1, 2], [r['index'] for r in merged['results']])
            with self.assertRaises(ValueError):
                timed.merge(timed_dir / 'manifest.json', results[:1], root / 'bad.json')
            with self.assertRaises(ValueError):
                timed.merge(timed_dir / 'manifest.json', [results[0], results[0]], root / 'bad.json')
            altered = json.loads(results[0].read_text())
            altered['seconds'] = 1.0
            results[0].write_bytes(fixed.encoded(altered))
            with self.assertRaises(ValueError):
                timed.merge(timed_dir / 'manifest.json', results, root / 'bad.json')
            invalid_wheel = root / 'wrong.whl'
            invalid_wheel.write_bytes(b'wrong')
            with self.assertRaises(ValueError):
                timed.seal(timed_dir / 'manifest.json', 0, invalid_wheel,
                           timed_dir / 'baseline-000.json', results[0])


if __name__ == '__main__':
    unittest.main()
