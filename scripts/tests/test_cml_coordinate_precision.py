"""Negative controls for the independent coordinate/graph observer."""
import unittest

from scripts.validation.check_cml_coordinate_precision import compare, coordinates, packet


class CmlCoordinateObserverTests(unittest.TestCase):
    def test_missing_coordinates_are_not_zero(self):
        source = packet()[0]['source'].replace('x2="1.25"', '')
        with self.assertRaises(KeyError):
            coordinates(source)

    def test_correct_coordinates_do_not_hide_changed_graph(self):
        for before, after in [('order="1"', 'order="2"'),
                              ('elementType="O"', 'elementType="N"'),
                              ('elementType="O"', 'elementType="O" formalCharge="1"')]:
            with self.assertRaises(ValueError):
                coordinates(packet()[0]['source'].replace(before, after))

    def test_atom_serialization_order_does_not_define_correspondence(self):
        source = packet()[0]['source']
        lines = source.splitlines()
        lines[2], lines[3] = lines[3], lines[2]
        self.assertEqual(coordinates('\n'.join(lines)), coordinates(source))

    def test_loss_has_no_relative_tolerance_or_translation_escape(self):
        expected = [[1000000.0, 0.0], [1000001.4, 0.0]]
        result = compare(expected, [[1000000.000001, 0.0], [1000001.400001, 0.0]])
        self.assertFalse(result['within_tolerance'])
        self.assertFalse(result['exact_binary64_values'])

    def test_exact_and_tolerance_are_separate(self):
        result = compare([[0.0, 0.0], [1.4, 0.1]], [[1e-11, 0.0], [1.4, 0.1]])
        self.assertTrue(result['within_tolerance'])
        self.assertFalse(result['exact_binary64_values'])


if __name__ == '__main__':
    unittest.main()
