import copy
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('history_observer',ROOT/'scripts/validation/check_dlpoly_history.py')
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class HistoryObserverTests(unittest.TestCase):
    def frame(self):
        return {'atoms':{'elements':{'number':[6,8]},'coords':{'3d':[0.123456,0.,0.,1.,0.,0.]}}}

    def test_all_frames_counted(self):
        a=self.frame();b=copy.deepcopy(a);b['atoms']['coords']['3d'][0]=2.0
        self.assertEqual(len(M.observe_cjson([a,b])),2)
        self.assertNotEqual(M.observe_cjson([a,b]),M.observe_cjson([a]))

    def test_changes_and_reordering_detected(self):
        a=self.frame();b=copy.deepcopy(a);b['atoms']['coords']['3d'][0]=0.12346
        self.assertNotEqual(M.observe_cjson([a]),M.observe_cjson([b]))
        b=copy.deepcopy(a);b['atoms']['elements']['number'].reverse()
        self.assertNotEqual(M.observe_cjson([a]),M.observe_cjson([b]))

    def test_invalid_output_not_scored(self):
        for values in ([],[0.0]*5,[math.nan]*6):
            a=self.frame();a['atoms']['coords']['3d']=values
            with self.assertRaises(ValueError):M.observe_cjson([a])
        with self.assertRaises(ValueError):M.observe_cjson([])

    def test_matching_record_and_cjson_projection(self):
        r={'frames':[{'atoms':[{'name':'C','position':[0.123456,0.,0.]},{'name':'O','position':[1.,0.,0.]}]}]}
        self.assertEqual(M.observe_record(r),M.observe_cjson([self.frame()]))


if __name__=='__main__':unittest.main()
