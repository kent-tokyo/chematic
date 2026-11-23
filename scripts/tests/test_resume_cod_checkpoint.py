import json
from pathlib import Path
import tempfile
import unittest
from argparse import Namespace

from scripts.validation.resume_cod_checkpoint import digest, merge, prepare


class ResumeCheckpointTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.split = self.root / 'split.json'
        self.audit = self.root / 'audit.json'
        self.checkpoint = self.root / 'checkpoint.json'
        self.resume_split = self.root / 'resume-split.json'
        self.resume_audit = self.root / 'resume-audit.json'
        self.resumed = self.root / 'resumed.json'
        self.output = self.root / 'combined.json'
        self.write(self.split, {'source_sha256': 'source', 'invalid_exclusions': [1],
                                'rows': [{'index': i, 'name': str(i), 'split': 'holdout'} for i in range(4)]})
        self.write(self.audit, {'split_sha256': digest(self.split), 'source_sha256': 'source',
                                'invalid_exclusion_count': 1, 'overlaps': []})
        self.meta = {'rdkit': '2026.09.2', 'chematic': '1.2.2', 'source_sha256': 'source',
                     'script_sha256': 'script', 'profile': 'corrected', 'geometry_checked': True,
                     'repair_iters': 400, 'candidate_limit': 20, 'mmff_iterations': 200, 'seed': 42,
                     'selection': 'holdout only sorted by SHA256(name)', 'platform': 'macOS-test-arm64'}
        self.write(self.checkpoint, {**self.meta, 'split_sha256': digest(self.split),
                                     'results': [{'index': i, 'name': str(i)} for i in range(2)]})
        prepare(Namespace(split=self.split, audit=self.audit, checkpoint=self.checkpoint,
                          out_split=self.resume_split, out_audit=self.resume_audit))
        self.write(self.resumed, {**self.meta, 'split_sha256': digest(self.resume_split),
                                  'results': [{'index': i, 'name': str(i)} for i in range(2, 4)]})
        self.merge_args = Namespace(split=self.split, audit=self.audit, checkpoint=self.checkpoint,
                                    resume_split=self.resume_split, resumed=self.resumed, output=self.output)

    @staticmethod
    def write(path, obj):
        path.write_text(json.dumps(obj))

    def test_complete_resume_preserves_all_records(self):
        merge(self.merge_args)
        out = json.loads(self.output.read_text())
        self.assertEqual([r['index'] for r in out['results']], [0, 1, 2, 3])
        self.assertEqual(out['resume_provenance']['completed_before_resume'], 2)
        self.assertEqual(out['split_sha256'], digest(self.split))

    def test_missing_or_duplicate_records_rejected(self):
        obj = json.loads(self.resumed.read_text())
        for bad in ([obj['results'][0]], [obj['results'][0], obj['results'][0]]):
            obj['results'] = bad
            self.write(self.resumed, obj)
            with self.assertRaisesRegex(ValueError, 'missing, duplicated, or misordered'):
                merge(self.merge_args)

    def test_wrong_platform_and_settings_rejected(self):
        obj = json.loads(self.resumed.read_text())
        obj['platform'] = 'Linux-x86_64'
        self.write(self.resumed, obj)
        with self.assertRaisesRegex(ValueError, 'macOS arm64'):
            merge(self.merge_args)
        obj['platform'] = 'macOS-test-arm64'
        obj['seed'] = 17
        self.write(self.resumed, obj)
        with self.assertRaisesRegex(ValueError, 'metadata: seed'):
            merge(self.merge_args)


if __name__ == '__main__':
    unittest.main()
