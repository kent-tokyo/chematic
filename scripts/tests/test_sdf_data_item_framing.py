"""Independent SD block observer controls; no native reader in these assertions."""
import json
import unittest

from scripts.validation.check_sdf_data_item_framing import ROOT, graph, properties


class SdfFieldObserverTests(unittest.TestCase):
    def setUp(self):
        self.packet = json.loads((ROOT/'validation/sdf_data_item_framing/cases.json').read_text())

    def test_header_shaped_lines_are_data_until_empty_separator(self):
        source = self.packet['ctabs']['v2000']+'> <NOTE>\n> <LITERAL>\nlast\n\n$$$$\n'
        self.assertEqual(properties(source), {'NOTE':'> <LITERAL>\nlast'})

    def test_empty_values_do_not_hide_fake_extra_fields(self):
        source = self.packet['ctabs']['v2000']+'> <NOTE>\n\n> <LITERAL>\nlast\n\n$$$$\n'
        self.assertEqual(properties(source), {'NOTE':'','LITERAL':'last'})

    def test_whitespace_is_not_normalized_by_the_observer(self):
        source = self.packet['ctabs']['v2000']+'> 17 <NOTE>\n  text\t \n\n$$$$\n'
        self.assertEqual(properties(source)['NOTE'], '  text\t ')

    def test_changed_graph_is_not_a_field_preservation_success(self):
        for fmt in ['v2000','v3000']:
            ctab = self.packet['ctabs'][fmt]
            graph(ctab)
            if fmt == 'v2000':
                changed = ctab.replace('  1  2  1', '  1  2  2')
            else:
                changed = ctab.replace('M  V30 1 1 1 2','M  V30 1 2 1 2')
            with self.assertRaises(ValueError): graph(changed)

    def test_record_count_changes_are_not_silently_skipped(self):
        source = self.packet['ctabs']['v2000']+'> <NOTE>\ntext\n\n$$$$\n'
        with self.assertRaises(ValueError): properties(source+source)


if __name__ == '__main__': unittest.main()
