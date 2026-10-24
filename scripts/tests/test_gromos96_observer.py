"""Negative controls for independent G96 framing, unit and geometry observations."""
import copy
import json
from pathlib import Path
import unittest
from scripts.validation.check_gromos96_preservation import project,bits,inspect_text,decimal_positions,writer_record

P=json.loads((Path(__file__).resolve().parents[2]/'validation/gromos96_preservation/cases.json').read_text(encoding='utf-8'))

class Gromos96ObserverTests(unittest.TestCase):
    def test_inspector_retains_all_blocks_comments_and_fields(self):
        for c in P['valid']:self.assertEqual(inspect_text(c['text'],c['record']['length_unit']),c['record'])
    def test_explicit_units_scale_cell_and_positions_together(self):
        g=project(P['valid'][0]['record'])[0]
        self.assertEqual(g['positions_angstrom'][0],[1.25,2.5,5.])
        self.assertEqual(g['box_vectors_angstrom'],[[20.,0.,0.],[0.,40.,0.],[0.,0.,80.]])
    def test_unknown_units_never_get_a_default_scale(self):
        r=copy.deepcopy(P['valid'][0]['record']);r['length_unit']='unspecified'
        with self.assertRaises(ValueError):project(r)
        with self.assertRaises(ValueError):decimal_positions(r)
    def test_wrong_unit_changes_geometry(self):
        a=P['valid'][0]['record'];b=copy.deepcopy(a);b['length_unit']='angstrom'
        self.assertNotEqual(bits(project(a)),bits(project(b)))
    def test_spaced_fields_change_layout_without_physical_loss(self):
        c=next(c for c in P['valid'] if c['id']=='ob-spaced-fields')
        self.assertNotEqual(c['record'],P['valid'][0]['record']);self.assertEqual(bits(project(c['record'])),bits(project(P['valid'][0]['record'])))
    def test_names_are_never_element_truth(self):
        c=next(c for c in P['valid'] if c['id']=='labels-are-not-elements')
        g=project(c['record'])[0];self.assertEqual([a['atom_name'] for a in g['identities']],['CA','QQ7']);self.assertNotIn('elements',g)
    def test_dropped_frame_detected(self):
        c=next(c for c in P['valid'] if c['id']=='two-frames-full');r=copy.deepcopy(c['record']);r['blocks']=r['blocks'][:3]
        self.assertEqual(len(project(c['record'])),2);self.assertEqual(len(project(r)),1)
    def test_missing_end_or_positions_are_unobservable(self):
        for text in ('','TITLE\nx\nEND\n','TITLE\nx\n'):
            with self.assertRaises(ValueError):inspect_text(text,'nm')
    def test_overflow_is_not_an_infinite_success(self):
        for c in P['projection_refused']:
            with self.assertRaises(ValueError):project(c['record'])
    def test_writer_probe_unit_changes_are_physically_equal(self):
        for probe in P['writer_probes']:self.assertEqual(decimal_positions(writer_record(probe,'nm')),decimal_positions(writer_record(probe,'angstrom')))
    def test_signed_zero_and_comment_loss_are_visible(self):
        self.assertNotEqual(bits(-0.),bits(0.))
        c=next(c for c in P['valid'] if c['id']=='comments-layout');r=copy.deepcopy(c['record']);r['trailing']=[]
        self.assertNotEqual(r,c['record']);self.assertEqual(bits(project(r)),bits(project(c['record'])))

if __name__=='__main__':unittest.main()
