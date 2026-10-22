import unittest
from scripts.validation.check_sdf_source_document import framing,observe,framed_carrier
A='first\n\n\n  0  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n>  <DATA>\nprefix $$$$ suffix\nM  END\n\n$$$$\n'
B=A.replace('first\n','second\n',1)
class SdfFramingObserverTests(unittest.TestCase):
    def test_exact_count_order_and_crlf(self):
        for text in [A+B,(A+B).replace('\n','\r\n')]:
            self.assertTrue(all(observe(text,text,['first','second']).values()))
    def test_first_record_prefix_and_reordering_are_detected(self):
        self.assertFalse(observe(A+B,A,['first','second'])['record_count'])
        result=observe(A+B,B+A,['first','second'])
        self.assertTrue(result['record_count']);self.assertFalse(result['ordered_titles']);self.assertFalse(result['exact_source'])
    def test_title_and_values_are_not_false_terminators_or_delimiters(self):
        text=A.replace('first\n','M  END\n',1)
        self.assertEqual(framing(text)['record_count'],1)
        self.assertEqual(framing(text)['titles'],['M  END'])
    def test_blank_title_and_unterminated_last_record(self):
        self.assertEqual(framing(A.replace('first\n','\n',1))['titles'],[''])
        self.assertEqual(framing(A+B.removesuffix('$$$$\n'))['titles'],['first','second'])
    def test_trailing_nuls_are_separate_from_observed_count(self):
        result=observe(A,A+'\x00\x00',['first'])
        self.assertTrue(result['record_count']);self.assertFalse(result['clean_framing']);self.assertFalse(result['exact_source'])
    def test_empty_output_and_missing_terminator_do_not_count_as_clean(self):
        self.assertFalse(observe(A,'',['first'])['clean_framing'])
        text='first\n\n\n  0  0  0  0  0  0  0  0  0  0  0999 V2000\n$$$$\n'
        self.assertFalse(framing(text)['complete_terminators'])
    def test_copy_carrier_is_explicitly_framed(self):
        self.assertEqual(framed_carrier(A.rstrip('\n')),A)
        self.assertTrue(framed_carrier(A.removesuffix('$$$$\n')).endswith('$$$$\n'))
        self.assertTrue(framed_carrier(A.removesuffix('$$$$\n').rstrip('\n')).endswith('\n\n$$$$\n'))
        self.assertTrue(framed_carrier(A.replace('\n','\r\n')).endswith('$$$$\r\n'))
    def test_unicode_separator_is_not_a_physical_header_line(self):
        text=A.replace('first\n','first\u2028title\n',1)
        self.assertEqual(framing(text)['titles'],['first\u2028title'])
        self.assertTrue(observe(text,text,['first\u2028title'])['clean_framing'])
if __name__=='__main__':unittest.main()
