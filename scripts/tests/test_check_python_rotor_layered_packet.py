import copy
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'benchmarks'))
SPEC = importlib.util.spec_from_file_location('rotor_layered_packet', Path(__file__).parents[1] / 'checks' / 'check_python_rotor_layered_packet.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class RotorLayeredPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-rotor-layered-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))

    def test_real_packet_recomputes_all_lanes(self):
        result = audit.check(self.original)
        self.assertEqual(result['comparisons'], 72)
        self.assertTrue(result['source_adoption_passed'])

    def reject(self, mutate):
        packet = copy.deepcopy(self.original)
        mutate(packet)
        with self.assertRaises(AssertionError):
            audit.check(packet)

    def test_duplicate_order_is_rejected(self):
        self.reject(lambda p: p['benchmark']['blocks'][1].update(order=p['benchmark']['blocks'][0]['order']))

    def test_native_identity_mismatch_is_rejected(self):
        self.reject(lambda p: p['benchmark']['blocks'][0]['candidate'].update(extension_sha256={}))

    def test_source_digest_mutation_is_rejected(self):
        self.reject(lambda p: p['source_files']['candidate'].update({'Cargo.toml': '0'*64}))

    def test_nonpositive_sample_is_rejected(self):
        self.reject(lambda p: p['benchmark']['blocks'][0]['candidate']['operations']['fp_layered'].update(first_use=0))

    def test_plan_mutation_is_rejected(self):
        self.reject(lambda p: p['measurement_plan'].update(limit=999))

    def test_source_state_loss_is_rejected(self):
        self.reject(lambda p: p['correctness'].update(source_exact_states=56287))

    def test_adoption_cannot_be_changed(self):
        self.reject(lambda p: p.update(source_adoption_passed=not p['source_adoption_passed']))

    def test_mqn_output_difference_cannot_be_turned_into_a_speed_claim(self):
        self.reject(lambda p: p['benchmark']['summary']['mqn/first_use'].update(cosmolkit={'speedup': {'median': 99}, 'demonstrated_win': True}))


if __name__ == '__main__':
    unittest.main()
