"""Negative controls for the four-field observer, without an external reader."""
import json
import unittest

from scripts.validation.check_v3000_atom_scalars import (
    capture, moljson_fields_exact, observe, process_outcome,
)


class AtomScalarObserverTests(unittest.TestCase):
    expected = dict(number=6, charge=-1, isotope=13, atom_map=7)

    def test_each_field_changes_the_score(self):
        self.assertTrue(observe(self.expected, '', lambda _: self.expected)['fields_exact'])
        for key in self.expected:
            changed = dict(self.expected, **{key: self.expected[key] + 1})
            self.assertFalse(observe(self.expected, '', lambda _: changed)['fields_exact'])
        changed = dict(self.expected, charge=True)
        self.assertFalse(observe(dict(self.expected, charge=1), '', lambda _: changed)['fields_exact'])

    def test_unreadable_is_not_preserved(self):
        def broken(_):
            raise ValueError('unreadable')
        result = observe(self.expected, '', broken)
        self.assertFalse(result['readable'])
        self.assertFalse(result['fields_exact'])

    def test_moljson_types_and_singleton_shape(self):
        value = dict(atoms=[dict(element='C', charge=-1, isotope=13)], bonds=[])
        self.assertTrue(moljson_fields_exact(self.expected, json.dumps(value)))
        for key, replacement in [('element', 'N'), ('charge', 0), ('charge', '-1'),
                                 ('charge', True), ('isotope', 12), ('isotope', '13')]:
            changed = json.loads(json.dumps(value))
            changed['atoms'][0][key] = replacement
            self.assertFalse(moljson_fields_exact(self.expected, json.dumps(changed)))
        for changed in [dict(value, atoms=[]), dict(value, atoms=value['atoms'] * 2),
                        dict(value, bonds=[dict()]), dict(), dict(atoms=[dict()], bonds=[])]:
            self.assertFalse(moljson_fields_exact(self.expected, json.dumps(changed)))
        self.assertFalse(moljson_fields_exact(self.expected, '{"atoms":NaN}'))

    def test_unspecified_isotope_requires_null(self):
        value = dict(atoms=[dict(element='C', charge=-1, isotope=None)], bonds=[])
        expected = dict(self.expected, isotope=0)
        self.assertTrue(moljson_fields_exact(expected, json.dumps(value)))
        value['atoms'][0]['isotope'] = 0
        self.assertFalse(moljson_fields_exact(expected, json.dumps(value)))

    def test_capture_separates_errors_from_output(self):
        self.assertEqual(capture(lambda: 'C'), dict(accepted=True, output='C'))
        def broken():
            raise ValueError('MASS invalid')
        self.assertEqual(capture(broken), dict(accepted=False, error='MASS invalid'))

    def test_terminal_outcomes_are_distinct(self):
        run = dict(resource_status='completed', exit_code=0, stdout='C\n')
        self.assertEqual(process_outcome(run), 'converted')
        for changes, outcome in [({'stdout': ' \n'}, 'no_output'),
                                 ({'exit_code': 1}, 'error'),
                                 ({'exit_code': -11}, 'signal'),
                                 ({'resource_status': 'timeout'}, 'resource_limit')]:
            self.assertEqual(process_outcome(dict(run, **changes)), outcome)


if __name__ == '__main__':
    unittest.main()
