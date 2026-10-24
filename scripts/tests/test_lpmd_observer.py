"""Negative controls for the independent LPMD preservation observer."""
import copy
import json
from pathlib import Path
import unittest
from scripts.validation.check_lpmd_preservation import read_layout, physical, geometry, bits

PACKET=json.loads((Path(__file__).resolve().parents[2]/'validation/lpmd_preservation/cases.json').read_text(encoding='utf-8'))

class LpmdObserverTests(unittest.TestCase):
    def test_layout_all_frames_and_fields(self):
        for c in PACKET['valid']: self.assertEqual(read_layout(c['text']),c['record'])
    def test_missing_symbol_header_is_invalid_even_when_process_exits_zero(self):
        with self.assertRaises(ValueError):read_layout(PACKET['valid'][0]['text'].replace('HDR SYM','HDR'))
    def test_lost_column_refuses_instead_of_becoming_observable(self):
        with self.assertRaises(ValueError):read_layout(PACKET['valid'][0]['text'].replace('SYM X Y Z','SYM X Y'))
    def test_truncated_frame_refuses(self):
        with self.assertRaises(ValueError):read_layout(PACKET['valid'][0]['text']+'2\n')
    def test_row_vector_orientation(self):
        c=next(c for c in PACKET['valid'] if c['id']=='triclinic')
        self.assertEqual(physical(c['record'])[0][0],[.75,2.125,8.5])
    def test_basis_and_fraction_equivalence_is_not_physical_loss(self):
        a=PACKET['valid'][0]['record'];b=copy.deepcopy(a)
        for frame in b['frames']:
            frame['cell']=[[str(float(v)*2) for v in row] for row in frame['cell']]
            for row in frame['atoms']:
                for i in (1,2,3): row[i]=str(float(row[i])/2)
        self.assertNotEqual(a,b);self.assertEqual(bits(physical(a)),bits(physical(b)))
    def test_empty_external_output_unobservable(self):
        self.assertFalse(geometry(PACKET['valid'][0]['record'],'')['observable'])
    def test_wrong_elements_and_coordinates_detected(self):
        result=geometry(PACKET['valid'][0]['record'],'2\n\n* 0 1 4\n* 0 6 8\n')
        self.assertTrue(result['observable']);self.assertFalse(result['elements']);self.assertFalse(result['cartesian_bits'])
    def test_missing_frames_detected(self):
        c=next(c for c in PACKET['valid'] if c['id']=='three-frames-variable-count')
        self.assertFalse(geometry(c['record'],'2\n\nC .5 2 8\nCl 3 4 0\n')['frame_and_atom_counts'])
    def test_unsafe_projection_not_silent_zero_or_infinity(self):
        for c in PACKET['projection_refused']:
            with self.assertRaises(ValueError):physical(c['record'])

if __name__=='__main__':unittest.main()
