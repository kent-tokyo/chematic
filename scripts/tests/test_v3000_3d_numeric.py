"""Damage controls use the real observation path and independent numeric lexer."""
import json
from pathlib import Path
import unittest

from scripts.validation.check_v3000_3d_numeric import atom_rows, bits, observe, raw_coordinates


class ObserverTests(unittest.TestCase):
    def setUp(self):
        path=Path(__file__).resolve().parents[2]/'validation/v3000_3d_numeric/cases.json'
        self.cases=json.loads(path.read_text())['cases']
        self.case=self.cases[1]
        self.expected=[[bits(v) for v in row] for row in self.case['expected_coordinates']]
        self.reference=self.reader(self.case['input'])

    @staticmethod
    def reader(text):
        rows=atom_rows(text)
        return dict(coordinates=raw_coordinates(text),graph=[row[1] for row in rows],stereo=[0]*len(rows))

    def observation(self,text):
        return observe(self.expected,self.reference,text,self.reader)

    def test_unchanged_output_passes(self):
        self.assertTrue(self.observation(self.case['input'])['joint_preserved'])

    def test_rounding_is_detected(self):
        text=self.case['input'].replace('1.2345678901234567','1.2346')
        result=self.observation(text)
        self.assertTrue(result['indexed_graph_exact'])
        self.assertFalse(result['coordinates_exact'])
        self.assertFalse(result['joint_preserved'])

    def test_signed_zero_is_not_ordinary_zero(self):
        case=next(c for c in self.cases if c['id']=='signed_zero')
        expected=[[bits(v) for v in row] for row in case['expected_coordinates']]
        result=observe(expected,self.reader(case['input']),case['input'].replace('-0.0','0.0'),self.reader)
        self.assertFalse(result['coordinates_exact'])

    def test_changed_element_or_extra_atom_fails(self):
        for text in [self.case['input'].replace('1 C ','1 N '),
                     self.case['input'].replace('M  V30 END ATOM','M  V30 2 C 0 0 0 0\nM  V30 END ATOM')]:
            self.assertFalse(self.observation(text)['indexed_graph_exact'])
            self.assertFalse(self.observation(text)['joint_preserved'])

    def test_stereo_is_checked_independently(self):
        actual=self.reader(self.case['input']);actual['stereo']=[1]
        result=observe(self.expected,self.reference,self.case['input'],lambda _:actual)
        self.assertTrue(result['coordinates_exact'])
        self.assertFalse(result['source_stereo_exact'])
        self.assertFalse(result['joint_preserved'])

    def test_continuations_and_incomplete_output(self):
        text=self.case['input'].replace('3.456789012345678 0','-\nM  V30 3.456789012345678 0')
        self.assertEqual(raw_coordinates(text),self.expected)
        self.assertFalse(self.observation('')['readable'])
        self.assertFalse(self.observation('M  V30 BEGIN ATOM\nM  V30 1 C 0 -\n')['readable'])


if __name__=='__main__':unittest.main()
