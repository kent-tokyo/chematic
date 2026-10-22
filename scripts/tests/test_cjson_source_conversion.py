"""Protect the independent source-field observer from false wins/losses."""
import unittest

from scripts.validation.check_cjson_source_conversion import observe, source_fields_equal, changed_source_paths


class CjsonSourceObserverTests(unittest.TestCase):
    def test_added_fields_do_not_penalize_retention(self):
        self.assertTrue(source_fields_equal({'atoms': {'n': [6, 8]}}, {'atoms': {'n': [6, 8], 'extra': 0}, 'extra': 1}))

    def test_missing_charge_precision_and_order_are_losses(self):
        for a, b in [({'q': [0.123456789012345]}, {}),
                     ({'q': [0.123456789012345]}, {'q': [0.123]}),
                     ([6, 8], [8, 6]), ([6, 8], [6]),
                     (9007199254740993, 9007199254740992.0),
                     (-0.0, 0.0), (True, 1), (None, 'null')]:
            self.assertFalse(source_fields_equal(a, b))

    def test_changed_paths_escape_opaque_keys(self):
        self.assertEqual(changed_source_paths({'vendor/a~b': {'id': 7}}, {'vendor/a~b': {'id': 8}}),
                         ['/vendor~1a~0b/id'])

    def test_json_fields_and_exact_source_are_separate(self):
        self.assertEqual(observe('{"n":1}', '{ "n": 1, "new": 2 }'),
                         dict(exact_source=False, source_fields_retained=True, changed_source_paths=[]))
        self.assertEqual(observe('{"n":1}', None), dict(exact_source=False, source_fields_retained=False))
        self.assertFalse(observe('{"n":1}', 'bad')['source_fields_retained'])


if __name__ == '__main__':
    unittest.main()
