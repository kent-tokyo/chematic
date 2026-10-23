import json
import unittest

from scripts.validation.check_cjson_output import classify


class CjsonOutputObservationTests(unittest.TestCase):
    def test_supported_case_requires_graph_and_geometry(self):
        expected = dict(numbers=[6, 7], formal_charges=[0, 0], isotopes=[0, 0],
                        coordinates=[[0, 0, 1], [1, 2, 3]], bonds=[[0, 1, 2]])
        case = dict(accepted=True, expected=expected)
        for order, coords in [(1, [0, 0, 1, 1, 2, 3]), (2, [0, 0, 0, 0, 0, 0]), (2, None)]:
            atoms = dict(elements=dict(number=[6, 7]))
            if coords is not None:
                atoms['coords'] = {'3d': coords}
            row = dict(output=json.dumps(dict(atoms=atoms, bonds=dict(connections=dict(index=[0, 1]), order=[order]))))
            classify(case, row)
            self.assertEqual(row['outcome'], 'changed')

    def test_output_plus_error_is_not_a_safe_refusal(self):
        row = dict(error='cannot preserve bond order Zero',
                   output='{"atoms":{"elements":{"number":[6,7]}},"bonds":{"order":[1]}}')
        classify(dict(accepted=False, error_order='Zero'), row)
        self.assertEqual(row['outcome'], 'unsupported_output')
        self.assertFalse(row['typed_refusal'])

    def test_wrong_refusal_reason_is_not_accepted(self):
        row = dict(error='CJSON: invalid coordinate count')
        classify(dict(accepted=False, error_order='Zero'), row)
        self.assertFalse(row['typed_refusal'])

    def test_empty_exit_zero_is_not_a_refusal_receipt(self):
        row = dict(output='', exit_code=0, stderr='0 molecules converted')
        classify(dict(accepted=False, error_order='Zero'), row)
        self.assertEqual(row['outcome'], 'no_output')
        self.assertFalse(row['typed_refusal'])


if __name__ == '__main__':
    unittest.main()
