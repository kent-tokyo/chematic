"""Corruption controls for release and measured wheel provenance."""
import copy
import json
import unittest

from scripts.checks.evidence import check_v122_artifact_consistency as gate


class ArtifactConsistencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = [json.loads(path.read_text()) for path in
                       (gate.REGISTRY, gate.PYPI, gate.CRATES, gate.NPM, gate.CHANNELS, gate.GATES, gate.MCS, gate.SPEED)]

    def test_current_records_agree(self):
        gate.check(*self.records)

    def test_current_cosmolkit_wheel_records_agree(self):
        records = copy.deepcopy(self.records)
        records[0] = json.loads(gate.CURRENT_REGISTRY.read_text())
        records[6] = json.loads(gate.CURRENT_MCS.read_text())
        records[7] = json.loads(gate.CURRENT_SPEED.read_text())
        versions = dict(gate.VERSIONS, cosmolkit=('cosmolkit', '0.6.3'))
        gate.check(*records, versions)
        records[0]['wheels'][-1]['version'] = '0.6.2'
        with self.assertRaisesRegex(ValueError, 'benchmark registry/MCS'):
            gate.check(*records, versions)

    def test_candidate_benchmark_wheel_hash_cannot_change(self):
        records = copy.deepcopy(self.records)
        next(row for row in records[0]['wheels'] if row['arm'] == 'v122')['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'benchmark registry/MCS'):
            gate.check(*records)

    def test_missing_published_distribution_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[1]['files'].pop(next(iter(records[1]['files'])))
        with self.assertRaises(ValueError):
            gate.check(*records)

    def test_changed_crate_tag_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[2]['files'][0]['vcs']['git']['sha1'] = '0' * 40
        with self.assertRaisesRegex(ValueError, 'crate provenance'):
            gate.check(*records)

    def test_unverified_npm_integrity_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[3]['computed_integrity_matches'] = False
        with self.assertRaisesRegex(ValueError, 'npm release'):
            gate.check(*records)

    def test_mismatched_speed_comparator_is_rejected(self):
        records = copy.deepcopy(self.records)
        name = next(row['filename'] for row in records[0]['wheels'] if row['arm'] == 'rdkit')
        records[7]['comparator_wheels'][name] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'speed comparator'):
            gate.check(*records)

    def test_missing_chemistry_job_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[5]['jobs'].pop()
        with self.assertRaisesRegex(ValueError, 'runtime or chemistry gate'):
            gate.check(*records)

    def test_missing_measurement_block_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[6]['blocks'].pop()
        with self.assertRaisesRegex(ValueError, 'block count'):
            gate.check(*records)

    def test_mcs_query_failure_cannot_be_ranked(self):
        records = copy.deepcopy(self.records)
        records[6]['blocks'][0]['v122']['validation'][0]['matches'][0] = False
        with self.assertRaisesRegex(ValueError, 'MCS output contract'):
            gate.check(*records)

    def test_empty_mcs_matches_cannot_be_ranked(self):
        records = copy.deepcopy(self.records)
        records[6]['blocks'][0]['v122']['validation'][0]['matches'] = []
        with self.assertRaisesRegex(ValueError, 'MCS output contract'):
            gate.check(*records)


if __name__ == '__main__':
    unittest.main()
