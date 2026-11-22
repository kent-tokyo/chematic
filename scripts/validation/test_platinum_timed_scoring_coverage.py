"""Check scoring coverage when engines fail independently or jointly."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from summarize_platinum_timed_scoring_coverage import summarize  # noqa: E402


class CoverageTest(unittest.TestCase):
    def test_unscored_accepted_and_errors(self):
        data = {'results': [
            {'error': 'legacy failure', 'rdkit': {'accepted': 2}},
            {'engine_errors': {'chematic': {'stage': 'generation', 'error': 'failure'}},
             'rdkit': {'accepted': 1, 'stereo': [{}]}},
            {'engine_errors': {'rdkit': {'stage': 'scoring', 'error': 'failure'}},
             'rdkit': {'accepted': 2, 'stereo': [{}]},
             'chematic': {'accepted': 0, 'stereo': []}},
        ]}
        s = summarize(data)
        self.assertEqual(3, s['engines']['rdkit']['unscored_accepted_conformers'])
        self.assertEqual(2, s['engines']['rdkit']['records_with_unscored_accepted'])
        self.assertEqual(1, s['engines']['rdkit']['scoring_errors'])
        self.assertEqual(1, s['engines']['chematic']['generation_errors'])
        self.assertEqual(2, s['engines']['chematic']['missing_engine_records'])

    def test_reject_more_scores_than_accepted(self):
        with self.assertRaises(ValueError):
            summarize({'results': [{'rdkit': {'accepted': 0, 'stereo': [{}]}}]})


if __name__ == '__main__':
    unittest.main()
