#!/usr/bin/env python3
"""Regression for asymmetric timed-engine failures and recovery scoring."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
import platinum_time_budget as probe  # noqa: E402
from summarize_platinum_time_budget import recovered  # noqa: E402


class EngineIsolationTest(unittest.TestCase):
    def test_other_engine_is_scored_for_either_order(self):
        calls = []

        def fake_attempts(engine, *_args):
            calls.append(engine)
            if engine == 'chematic':
                raise ValueError('unsupported atom')
            return ['coords'], {'attempted': 1, 'accepted': 1, 'elapsed': 0.1, 'statuses': [0]}

        stereo = {'atoms': [], 'bonds': [], 'violated_atoms': [], 'violated_bonds': [],
                  'unevaluable_atoms': [], 'unevaluable_bonds': []}
        geometry = {'finite': True, 'bond_violations': 0, 'angle_violations': 0, 'heavy_atom_clashes': 0}
        for order in (('chematic', 'rdkit'), ('rdkit', 'chematic')):
            with self.subTest(order=order):
                calls.clear()
                result = {'geometry_reference': {'bonds_supported': 1, 'bonds_total': 1,
                                                  'angles_supported': 1, 'angles_total': 1}}
                with patch.object(probe, 'attempts', side_effect=fake_attempts), \
                     patch.object(probe, 'rms_values', return_value=[0.5]), \
                     patch.object(probe, 'stereo_violations', return_value=[stereo]), \
                     patch.object(probe, 'check_geometry', return_value=geometry):
                    probe.score_engines(result, order, None, None, None, object(), 0.2, 20, 400)
                self.assertEqual(calls, list(order))
                self.assertEqual(result['engine_errors']['chematic']['stage'], 'generation')
                self.assertTrue(recovered(result, 'rdkit', True))
                self.assertFalse(recovered(result, 'chematic', True))
                self.assertTrue(recovered(result, 'rdkit', True, True))

    def test_scoring_failure_is_local_to_engine(self):
        result = {}
        def fake_attempts(engine, *_args):
            return [engine], {'attempted': 1, 'accepted': 1, 'elapsed': 0.1, 'statuses': [0]}
        def fake_rms(_mol, _ref, conformers):
            if conformers == ['chematic']:
                raise ValueError('cannot score')
            return [0.5]
        with patch.object(probe, 'attempts', side_effect=fake_attempts), \
             patch.object(probe, 'rms_values', side_effect=fake_rms), \
             patch.object(probe, 'stereo_violations', return_value=[{'atoms': [], 'bonds': []}]):
            probe.score_engines(result, ('rdkit', 'chematic'), None, None, None, None, 0.2, 20, 0)
        self.assertEqual(result['engine_errors']['chematic']['stage'], 'scoring')
        self.assertTrue(recovered(result, 'rdkit'))
        self.assertFalse(recovered(result, 'chematic'))


if __name__ == '__main__':
    unittest.main()
