import unittest
from scripts.validation.check_cml_common_retention import observe

SOURCE = '<molecule xmlns="http://www.xml-cml.org/schema"><atomArray><atom id="c" elementType="C"/><atom id="o" elementType="O"/></atomArray><bondArray><bond atomRefs2="c o" order="1"/></bondArray></molecule>'


class CommonCmlObserverTests(unittest.TestCase):
    def test_removed_bond_is_a_graph_mismatch_even_when_xml_is_valid(self):
        row = dict(output=SOURCE.replace('<bond atomRefs2="c o" order="1"/>', ''))
        observe(SOURCE, row)
        self.assertEqual(row['outcome'], 'selected_graph_mismatch')
        self.assertFalse(row['exact_source_bytes'])

    def test_xml_failures_and_absent_outputs_are_distinct(self):
        for output, expected in [('', 'no_output'), ('<broken', 'unresolved_output')]:
            row = dict(output=output)
            observe(SOURCE, row)
            self.assertEqual(row['outcome'], expected)

    def test_identical_output_is_retained_without_a_chemistry_claim(self):
        row = dict(output=SOURCE)
        observe(SOURCE, row)
        self.assertTrue(row['exact_source_bytes'])
        self.assertEqual(row['outcome'], 'selected_fields_retained')


if __name__ == '__main__':
    unittest.main()
