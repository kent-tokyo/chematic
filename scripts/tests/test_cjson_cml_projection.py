import unittest
from scripts.validation.check_cjson_cml_projection import observe

EXPECTED=dict(numbers=[6,8],charges=[0,0],isotopes=[0,0],bonds=[[0,1,1]],coordinates=[[-0.0,0,1.123456789012345],[1,0,0]],name='source')
XML='<molecule xmlns="http://www.xml-cml.org/schema" title="source"><atomArray><atom id="a" elementType="C" x3="-0.0" y3="0" z3="1.123456789012345"/><atom id="b" elementType="O" x3="1" y3="0" z3="0"/></atomArray><bondArray><bond atomRefs2="a b" order="1"/></bondArray></molecule>'


class ProjectionObserverTests(unittest.TestCase):
    def test_matching_fields_accept_title_as_a_name(self):
        self.assertTrue(all(observe(EXPECTED,XML).values()))

    def test_zero_sign_and_rounding_loss_are_detected(self):
        for output in [XML.replace('x3="-0.0"','x3="0.0"'),XML.replace('1.123456789012345','1.12346')]:
            self.assertFalse(observe(EXPECTED,output)['coordinate_bits'])

    def test_absent_coordinates_and_graph_changes_do_not_pass(self):
        self.assertFalse(observe(EXPECTED,XML.replace('order="1"','order="2"'))['indexed_graph'])
        self.assertFalse(observe(EXPECTED,XML.replace(' z3="0"',''))['coordinate_presence'])

    def test_name_loss_is_separate_from_graph_preservation(self):
        result=observe(EXPECTED,XML.replace(' title="source"',''))
        self.assertTrue(result['indexed_graph']);self.assertFalse(result['name'])

    def test_bad_xml_and_no_output_never_pass(self):
        for output in ['',None,'<broken']:
            self.assertFalse(any(observe(EXPECTED,output).values()))


if __name__=='__main__':unittest.main()
