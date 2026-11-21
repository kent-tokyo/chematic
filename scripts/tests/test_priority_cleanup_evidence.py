"""Cross-record controls for the prioritized Cleanup audit."""
import copy
import json
import unittest

from scripts.checks.evidence import check_standardization_source_evidence as gate


class PriorityCleanupEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = [json.loads(path.read_text()) for path in (
            gate.PRIORITY_BASE, gate.RECORD, gate.PRIORITY_AUDIT, gate.PRIORITY_ROOTED,
            gate.PRIORITY_WHEEL, gate.PRIORITY_BASE_AUDIT,
        )]

    def test_current_records_agree(self):
        gate.check_priority_audit(*self.records)

    def test_other_wheel_cannot_be_substituted(self):
        records = copy.deepcopy(self.records)
        records[4]['wheel_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'wheel identities'):
            gate.check_priority_audit(*records)

    def test_different_candidate_output_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[4]['rows'][0]['output'] = 'C'
        with self.assertRaisesRegex(ValueError, 'wheel/source output'):
            gate.check_priority_audit(*records)

    def test_formula_change_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[4]['rows'][0]['reparsed_formula'] = 'H2O'
        with self.assertRaisesRegex(ValueError, 'wheel formula or charge'):
            gate.check_priority_audit(*records)

    def test_rooted_variant_change_is_rejected(self):
        records = copy.deepcopy(self.records)
        records[3]['cases'][0]['variants'][0]['chematic'] = 'C'
        with self.assertRaisesRegex(ValueError, 'rooted output'):
            gate.check_priority_audit(*records)

    def test_repeat_claim_must_match_observed_rdkit_strings(self):
        records = copy.deepcopy(self.records)
        row = records[2]['rows'][0]
        row['rdkit_second_cleanup_matches'] = not row['rdkit_second_cleanup_matches']
        with self.assertRaisesRegex(ValueError, 'repeat status'):
            gate.check_priority_audit(*records)

    def test_historical_graph_difference_cannot_be_reclassified_silently(self):
        records = copy.deepcopy(self.records)
        row = next(x for x in records[5]['rows'] if x['index'] == 65)
        row['charge_h_connectivity_isomorphic'] = True
        row['bond_isotope_graph_isomorphic'] = True
        with self.assertRaises(ValueError):
            gate.check_priority_audit(*records)


if __name__ == '__main__':
    unittest.main()
