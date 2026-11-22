"""Verify that individual timed attempts follow the mixed gradient schedule."""
import importlib.util
import sys
import unittest
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / 'benchmarks'
sys.path.insert(0, str(BENCH))
spec = importlib.util.spec_from_file_location('platinum_time_mixed_strict', BENCH / 'platinum_time_mixed_strict.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TimedMixedScheduleTest(unittest.TestCase):
    def test_first_half_uses_compatible_and_second_half_corrected(self):
        class Molecule:
            def etkdg_embed_multiple_confs_mmff_compatible_checked(self, count, **kwargs):
                return ('compatible', count, kwargs['random_seed'])

            def etkdg_embed_multiple_confs_mmff_corrected(self, count, **kwargs):
                return ('corrected', count, kwargs['random_seed'])

        molecule = module.MixedAttempts(Molecule(), 10)
        results = [molecule.etkdg_embed_multiple_confs_mmff_corrected(1, random_seed=42 * (i + 1))
                   for i in range(20)]
        self.assertEqual([r[0] for r in results], ['compatible'] * 10 + ['corrected'] * 10)
        self.assertEqual([r[2] for r in results], [42 * (i + 1) for i in range(20)])


if __name__ == '__main__':
    unittest.main()
