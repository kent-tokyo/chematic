import copy
import unittest

from scripts.validation.validate_rdkit_issue9735 import status, verify_candidate

CASE=dict(id='invalid',expected='valence_error',atom=4)
ERROR=dict(ok=False,error='Explicit valence for atom # 4 C, 5, is greater than permitted')


class AromaticValenceObservation(unittest.TestCase):
    def row(self):
        return dict(canonical=copy.deepcopy(ERROR),checked=copy.deepcopy(ERROR),problems=[['AtomValenceException',[4]]])

    def test_acceptance_and_unrelated_refusal_do_not_count_as_valence_rejection(self):
        self.assertEqual(verify_candidate(CASE,self.row(),True),'valence_error')
        for result in [dict(ok=True,value='returned SMILES'),dict(ok=False,error="Can't kekulize mol")]:
            row=self.row();row['canonical']=result
            with self.assertRaises(ValueError):verify_candidate(CASE,row,True)

    def test_diagnosis_must_name_the_affected_atom(self):
        for problems in [[],[['AtomValenceException',[5]]],[['KekulizeException',[4]]]]:
            row=self.row();row['problems']=problems
            with self.assertRaises(ValueError):verify_candidate(CASE,row,True)

    def test_partial_outputs_and_invalid_envelopes_are_refused(self):
        for result in [dict(ERROR,value='partial'),dict(ok=1,value='accepted'),dict(ok=True,value='')]:
            with self.assertRaises(ValueError):status(result)

    def test_valid_controls_require_clean_diagnostics_and_matching_checked_output(self):
        case=dict(id='valid',expected='accept');value=dict(ok=True,value='c1ccccc1')
        row=dict(canonical=value,checked=copy.deepcopy(value),problems=[])
        self.assertEqual(verify_candidate(case,row,True),'accepted')
        row['checked']=copy.deepcopy(ERROR)
        with self.assertRaises(ValueError):verify_candidate(case,row,True)
        row['checked']=value;row['problems']=[['AtomValenceException',[0]]]
        with self.assertRaises(ValueError):verify_candidate(case,row,True)


if __name__=='__main__':unittest.main()
