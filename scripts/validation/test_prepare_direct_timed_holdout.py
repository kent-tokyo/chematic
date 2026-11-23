"""Validate direct timed plans without fabricating fixed-count results."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / 'benchmarks'
sys.path.insert(0, str(BENCH))
import shard_platinum_holdout as fixed  # noqa: E402
import shard_platinum_time as timed  # noqa: E402
import prepare_direct_timed_holdout as direct  # noqa: E402


class DirectTimedTest(unittest.TestCase):
    def test_complete_population_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            split = {'source_sha256': 'a' * 64, 'invalid_exclusions': [],
                     'rows': [{'index': i, 'name': f'mol-{i}', 'split': 'holdout',
                               'scaffold_sha256': f'scaffold-{i}'} for i in range(5)]}
            split_path = root / 'split.json'; split_path.write_bytes(fixed.encoded(split))
            audit_path = root / 'audit.json'
            audit_path.write_bytes(fixed.encoded({'source_sha256': split['source_sha256'],
                'split_sha256': fixed.digest(split_path.read_bytes()), 'invalid_exclusion_count': 0, 'overlaps': []}))
            data = {'source_sha256': split['source_sha256'], 'split_sha256': fixed.digest(split_path.read_bytes()),
                    'partition': 'holdout', 'results': [{'index': i, 'name': f'mol-{i}', 'smiles': 'CC'} for i in range(5)]}
            inputs = root / 'input.json'; inputs.write_bytes(fixed.encoded(data))
            wheel = root / 'wheel.whl'; wheel.write_bytes(b'candidate')
            out = root / 'plan'; m = direct.prepare(split_path, audit_path, inputs, out, 2, fixed.digest(wheel.read_bytes()))
            self.assertEqual(5, sum(s['count'] for s in m['shards']))
            shards = []
            for s in m['shards']:
                i = s['index']; baseline = out / f'baseline-{i:03d}.json'
                result = dict(m['protocol'], chematic='test', source_sha256=m['source_sha256'],
                              script_sha256=m['fingerprints']['time_scorer_sha256'],
                              stereo_verifier_sha256=m['fingerprints']['stereo_verifier_sha256'],
                              stereo_wrapper_sha256=m['fingerprints']['time_stereo_wrapper_sha256'],
                              platform='test', timing='policy', results=s['records'])
                path = out / f'result-{i}.json'; path.write_bytes(fixed.encoded(result))
                timed.seal(out / 'manifest.json', i, wheel, baseline, path); shards.append(path)
            merged = timed.merge(out / 'manifest.json', shards, root / 'merged.json')
            self.assertEqual([0, 1, 2, 3, 4], [r['index'] for r in merged['results']])
            self.assertEqual(m['protocol']['engine_error_policy'], merged['engine_error_policy'])
            data['results'][0]['name'] = 'wrong'; inputs.write_bytes(fixed.encoded(data))
            with self.assertRaises(ValueError):
                direct.prepare(split_path, audit_path, inputs, root / 'bad', 2, fixed.digest(wheel.read_bytes()))


if __name__ == '__main__':
    unittest.main()
