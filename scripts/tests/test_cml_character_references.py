import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

from scripts.validation.check_cml_character_references import SOURCE, packet, snapshot, failed_process


class CmlReferenceOracleTests(unittest.TestCase):
    def test_signals_panics_and_runtime_errors_are_not_rejections(self):
        for engine in ['chematic', 'openbabel']:
            for code in [-11, -9, 101, 127]:
                self.assertEqual(failed_process(engine, {'exit_code':code, 'stderr':'Error: failure'}, True), 'process_failure')
        self.assertEqual(failed_process('chematic', {'exit_code':1,'stderr':'Error: MalformedXml'}, True), 'rejected')
        self.assertEqual(failed_process('chematic', {'exit_code':1,'stderr':'Error: MalformedXml'}, False), 'read_write_failure')
        self.assertEqual(failed_process('openbabel', {'exit_code':0,'stderr':'0 molecules converted'}, True), 'rejected')

    def test_valid_lexical_variants_have_the_same_independent_snapshot(self):
        expected = snapshot(SOURCE)
        cases = packet()
        self.assertEqual(len(cases), 14)
        self.assertEqual(sum(c["expected_refusal"] for c in cases), 7)
        for case in cases:
            if not case["expected_refusal"]:
                self.assertEqual(snapshot(case["source"]), expected, case["id"])

    def test_changed_charge_isotope_or_coordinate_is_detected(self):
        expected = snapshot(SOURCE)
        for before, after in [('formalCharge="1"','formalCharge="0"'),
                              ('isotopeNumber="15"','isotopeNumber="14"'),
                              ('x2="1.25"','x2="0"')]:
            self.assertNotEqual(snapshot(SOURCE.replace(before, after)), expected)

    def test_duplicate_decoded_ids_and_invalid_xml_are_not_oracle_graphs(self):
        duplicate = next(c for c in packet() if c["id"] == 'duplicate-decoded-id')
        with self.assertRaises(ValueError): snapshot(duplicate['source'])
        with self.assertRaises(ET.ParseError): snapshot(SOURCE.replace('formalCharge="1"', 'formalCharge="&#0;"'))

    def test_checked_fixture_matches_the_recipe(self):
        data = json.loads((Path(__file__).resolve().parents[2] /
                           'validation/cml_character_references/cases.json').read_text())
        self.assertEqual(data['cases'], packet())
        self.assertEqual(data['expected'], snapshot(SOURCE))


if __name__ == '__main__': unittest.main()
