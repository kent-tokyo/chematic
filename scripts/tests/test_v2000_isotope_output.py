"""The isotope observer must detect loss, field shifts and malformed M ISO."""
import copy
import unittest

from scripts.validation.check_v2000_isotope_output import isotope_layout, observe, passed


def block(mass=12):
    atom = '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0'
    return f'example\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{atom}\nM  ISO  1   1{mass:>4}\nM  END\n'


class IsotopeOutputObserverTests(unittest.TestCase):
    expected = dict(numbers=[6], charges=[0], isotopes=[12], atom_maps=[7], bonds=[])

    def test_exact_fields_and_layout_are_required(self):
        self.assertTrue(passed(observe(self.expected, block(), lambda _: self.expected)))
        for key in ['numbers', 'charges', 'isotopes', 'atom_maps']:
            actual = copy.deepcopy(self.expected)
            actual[key][0] += 1
            result = observe(self.expected, block(), lambda _: actual)
            self.assertTrue(result['readable'])
            self.assertFalse(passed(result))
        actual = dict(self.expected, bonds=[[0, 0, 1]])
        self.assertFalse(passed(observe(self.expected, block(), lambda _: actual)))

    def test_unreadable_never_passes(self):
        def broken(_):
            raise ValueError('unreadable')
        result = observe(self.expected, block(), broken)
        self.assertFalse(result['readable'])
        self.assertFalse(passed(result))

    def test_four_digit_mass_and_profile_endpoints(self):
        for mass in [1, 12, 999, 9999]:
            self.assertTrue(isotope_layout(block(mass)))
        for mass in [0, -1, 65535]:
            self.assertFalse(isotope_layout(block(mass)))

    def test_legacy_mass_difference_range_and_alignment(self):
        for value in [-3, 0, 4]:
            self.assertTrue(isotope_layout(block().replace('C   0', f'C  {value:>2}', 1)))
        for value in [-4, 5, 100]:
            self.assertFalse(isotope_layout(block().replace('C   0', f'C  {value:>2}', 1)))

    def test_iso_assignment_count_target_and_duplicates(self):
        text = block()
        for changed in [text.replace('M  ISO  1', 'M  ISO  2'),
                        text.replace('M  ISO  1   1', 'M  ISO  1   2'),
                        text.replace('M  END', 'M  ISO  1   1  13\nM  END'),
                        text.replace('M  ISO  1', 'M  ISO  9'),
                        text.replace('  12\n', 'oops\n')]:
            self.assertFalse(isotope_layout(changed))

    def test_layout_alone_does_not_establish_isotope_retention(self):
        text = block().replace('M  ISO  1   1  12\n', '')
        self.assertTrue(isotope_layout(text))
        actual = dict(self.expected, isotopes=[0])
        self.assertFalse(passed(observe(self.expected, text, lambda _: actual)))


if __name__ == '__main__':
    unittest.main()
