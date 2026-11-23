"""Verify complete, provenance-checked assembly of scaffold holdout shards."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / 'benchmarks' / 'shard_platinum_holdout.py'
spec = importlib.util.spec_from_file_location('shard_platinum_holdout', MODULE)
sharding = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sharding)


class HoldoutShardingTest(unittest.TestCase):
    def test_mixed_profile_fingerprints_and_protocol(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split = {'source_sha256': 'a' * 64, 'invalid_exclusions': [],
                     'rows': [{'index': 0, 'name': 'mol', 'split': 'holdout',
                               'scaffold_sha256': 'b' * 64}]}
            split_path = root / 'split.json'
            split_path.write_bytes(sharding.encoded(split))
            audit_path = root / 'audit.json'
            audit_path.write_bytes(sharding.encoded({
                'split_sha256': sharding.digest(split_path.read_bytes()),
                'source_sha256': split['source_sha256'],
                'invalid_exclusion_count': 0, 'overlaps': [],
            }))
            wheel = root / 'candidate.whl'
            wheel.write_bytes(b'synthetic wheel')
            manifest = sharding.prepare(split_path, audit_path, root / 'shards', 1,
                                        sharding.digest(wheel.read_bytes()), 'mixed_checked')
            wrapper = MODULE.with_name('platinum_mixed_checked_development.py')
            self.assertEqual(manifest['stereo_wrapper_sha256'], sharding.digest(wrapper.read_bytes()))
            self.assertEqual(manifest['protocol']['profile'], 'mixed_checked')
            shard = manifest['shards'][0]
            result = dict(manifest['protocol'], chematic='test', platform='test',
                          script_sha256=manifest['scorer_sha256'], source_sha256=split['source_sha256'],
                          stereo_verifier_sha256=manifest['stereo_verifier_sha256'],
                          stereo_wrapper_sha256=manifest['stereo_wrapper_sha256'],
                          split_sha256=shard['split_sha256'],
                          results=[dict(r, status='ok') for r in shard['records']])
            result_path = root / 'result.json'
            result_path.write_bytes(sharding.encoded(result))
            sharding.seal(root / 'shards' / 'manifest.json', 0, wheel,
                          root / 'shards' / 'audit-000.json', result_path)
            merged = sharding.merge(root / 'shards' / 'manifest.json', [result_path], root / 'merged.json')
            self.assertEqual(merged['profile'], 'mixed_checked')
            altered = json.loads(result_path.read_text())
            altered['profile'] = 'corrected'
            result_path.write_bytes(sharding.encoded(altered))
            with self.assertRaises(ValueError):
                sharding.merge(root / 'shards' / 'manifest.json', [result_path], root / 'bad.json')

    def test_complete_and_tampered_shards(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split = {'source_sha256': 'a' * 64, 'invalid_exclusions': ['x'],
                     'rows': [{'index': i, 'name': f'mol-{i}', 'split': 'holdout',
                               'scaffold_sha256': str(i // 2)} for i in range(7)]}
            split_path = root / 'split.json'
            split_path.write_bytes(sharding.encoded(split))
            audit = {'split_sha256': sharding.digest(split_path.read_bytes()),
                     'source_sha256': split['source_sha256'],
                     'invalid_exclusion_count': 1, 'overlaps': []}
            audit_path = root / 'audit.json'
            audit_path.write_bytes(sharding.encoded(audit))
            wheel_path = root / 'candidate.whl'
            wheel_path.write_bytes(b'synthetic wheel')
            wheel_hash = sharding.digest(wheel_path.read_bytes())
            manifest = sharding.prepare(split_path, audit_path, root / 'shards', 3, wheel_hash)
            self.assertEqual([3, 2, 2], [s['count'] for s in manifest['shards']])
            altered_audit = dict(audit, overlaps=[{'name': 'conflict'}])
            bad_audit = root / 'bad-audit.json'
            bad_audit.write_bytes(sharding.encoded(altered_audit))
            with self.assertRaises(ValueError):
                sharding.prepare(split_path, bad_audit, root / 'bad-shards', 3, wheel_hash)
            result_paths = []
            for shard in manifest['shards']:
                result = dict(manifest['protocol'], chematic='test',
                              script_sha256=manifest['scorer_sha256'], source_sha256=split['source_sha256'],
                              stereo_verifier_sha256=manifest['stereo_verifier_sha256'],
                              stereo_wrapper_sha256=manifest['stereo_wrapper_sha256'],
                              split_sha256=shard['split_sha256'],
                              platform='test',
                              results=[dict(r, status='ok') for r in shard['records']])
                path = root / f"result-{shard['index']}.json"
                path.write_bytes(sharding.encoded(result))
                sharding.seal(root / 'shards' / 'manifest.json', shard['index'], wheel_path,
                              root / 'shards' / f"audit-{shard['index']:03d}.json", path)
                result_paths.append(path)
            wrong_wheel = root / 'wrong.whl'
            wrong_wheel.write_bytes(b'wrong wheel')
            with self.assertRaises(ValueError):
                sharding.seal(root / 'shards' / 'manifest.json', 0, wrong_wheel,
                              root / 'shards' / 'audit-000.json', result_paths[0])
            with self.assertRaises(ValueError):
                sharding.seal(root / 'shards' / 'manifest.json', 0, wheel_path,
                              root / 'shards' / 'audit-001.json', result_paths[0])
            manifest_path = root / 'shards' / 'manifest.json'
            result = sharding.merge(manifest_path, result_paths, root / 'merged.json')
            self.assertEqual(list(range(7)), [r['index'] for r in result['results']])
            self.assertEqual(sharding.digest(split_path.read_bytes()), result['split_sha256'])
            with self.assertRaises(ValueError):
                sharding.merge(manifest_path, result_paths[:-1], root / 'bad.json')
            with self.assertRaises(ValueError):
                sharding.merge(manifest_path, [result_paths[0]] * 3, root / 'bad.json')
            altered = json.loads(result_paths[0].read_text())
            altered['wheel_sha256'] = 'd' * 64
            result_paths[0].write_bytes(sharding.encoded(altered))
            with self.assertRaises(ValueError):
                sharding.merge(manifest_path, result_paths, root / 'bad.json')
            altered['wheel_sha256'] = wheel_hash
            altered['results'][0]['name'] = 'different'
            result_paths[0].write_bytes(sharding.encoded(altered))
            with self.assertRaises(ValueError):
                sharding.merge(manifest_path, result_paths, root / 'bad.json')


if __name__ == '__main__':
    unittest.main()
