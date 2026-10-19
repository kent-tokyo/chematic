import copy
import gzip
import importlib.util
import json
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('frontier_packet', Path(__file__).parents[1] / 'check_python_speed_frontier_packet.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class ExactIntervalTests(unittest.TestCase):
    def test_adoption_uses_unrounded_improvement_and_control_floor(self):
        primary = audit.exact_interval([1000001] * 24, [1000000] * 24, 0)
        self.assertEqual(round(primary['ci95'][0], 3), 1.0)
        self.assertTrue(audit.source_gate(primary, [primary]))
        control = audit.exact_interval([989999] * 24, [1000000] * 24, 0)
        self.assertEqual(round(control['median'], 3), 0.99)
        self.assertFalse(audit.source_gate(primary, [control]))


class NumericPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-morgan-direct-connectivity-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))

    def test_real_packet_is_recomputed(self):
        self.assertEqual(audit.check(self.original)['comparisons'], 63)

    def test_duplicate_order_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['benchmark']['blocks'][0]['order'] = packet['benchmark']['blocks'][1]['order']
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_native_hash_tampering_is_rejected(self):
        packet = copy.deepcopy(self.original)
        hashes = packet['build']['candidate']['native_sha256']
        hashes[next(iter(hashes))] = '0' * 64
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_nonpositive_duration_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['benchmark']['blocks'][0]['candidate']['operations']['morgan2_chiral']['first_use'] = 0
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_changed_adoption_decision_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['source_adoption_passed'] = not packet['source_adoption_passed']
        with self.assertRaises(AssertionError): audit.check(packet)


if __name__ == '__main__': unittest.main()
