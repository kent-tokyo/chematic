"""Independent observer negative controls for the original CML packet."""
import json
import unittest

from scripts.validation.check_cml_source_conversion import ROOT, check_report, compare, snapshot


class CmlSourceObserverTests(unittest.TestCase):
    def setUp(self):
        self.cases = json.loads((ROOT/'validation/cml_source_preservation/cases.json').read_text())['cases']

    def test_dual_dimensionality_cannot_pass_with_only_one_coordinate_set(self):
        source = self.cases[1]['source']
        changed = source.replace('x2="1.25" y2="-2.5"', '')
        result = compare(snapshot(source), snapshot(changed))
        self.assertFalse(result['all_selected_fields_retained'])
        self.assertFalse(result['fields']['C_coordinate_presence'])

    def test_missing_unknown_atom_attribute_is_still_a_loss(self):
        source = self.cases[3]['source']
        changed = source.replace(' v:label="site α"', '')
        self.assertFalse(compare(snapshot(source), snapshot(changed))['all_selected_fields_retained'])

    def test_repeated_property_order_is_not_silently_collapsed(self):
        source = self.cases[3]['source']
        first = '<property dictRef="v:note"><scalar>first α &amp; β</scalar></property>'
        second = '<property dictRef="v:note"><scalar>second</scalar></property>'
        changed = source.replace(first+'\n'+second, second+'\n'+first)
        self.assertFalse(compare(snapshot(source), snapshot(changed))['fields']['property_extension_payload'])

    def test_coordinate_equality_does_not_hide_changed_graph(self):
        with self.assertRaises(ValueError):
            snapshot(self.cases[0]['source'].replace('order="1"', 'order="2"'))

    def test_added_unscored_metadata_is_not_source_information_loss(self):
        source = self.cases[0]['source']
        changed = source.replace('</molecule>', '<metadataList><metadata name="generated"/></metadataList></molecule>')
        self.assertTrue(compare(snapshot(source), snapshot(changed))['all_selected_fields_retained'])

    def test_source_retention_does_not_claim_graph_interpretation(self):
        source = self.cases[0]['source']
        report = {'schema': 'chematic.conversion-report.v1', 'policy': 'strict',
                  'status': 'completed', 'verification_scope': 'exact_source_bytes',
                  'candidate_graph_verified': True, 'output': source,
                  'candidate_atom_mapping': [], 'candidate_bond_mapping': [],
                  'observations': []}
        with self.assertRaises(ValueError):
            check_report(report, source, 'strict')


if __name__ == '__main__':
    unittest.main()
