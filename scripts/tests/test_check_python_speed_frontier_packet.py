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

    def test_source_hash_tampering_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['build']['candidate']['source_sha256'] = '0' * 64
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_plan_hash_must_match_the_measured_manifest(self):
        import hashlib
        packet = copy.deepcopy(self.original)
        packet['measurement_plan']['hypothesis'] += ' modified after measurement'
        packet['measurement_plan_sha256'] = hashlib.sha256(
            (json.dumps(packet['measurement_plan'], indent=2)+'\n').encode()).hexdigest()
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_nonpositive_duration_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['benchmark']['blocks'][0]['candidate']['operations']['morgan2_chiral']['first_use'] = 0
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_changed_adoption_decision_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['source_adoption_passed'] = not packet['source_adoption_passed']
        with self.assertRaises(AssertionError): audit.check(packet)


class NumericPacketV2Tests(NumericPacketTests):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-morgan-tpsa-environments-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))

    def test_scalar_accounting_corruption_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['scalar_validation']['source_exact_cells_including_refusals'] += 1
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_second_primary_cannot_be_removed(self):
        import hashlib
        packet = copy.deepcopy(self.original)
        packet['measurement_plan']['primary'].pop()
        packet['measurement_plan_sha256'] = hashlib.sha256(
            (json.dumps(packet['measurement_plan'], indent=2)+'\n').encode()).hexdigest()
        packet['benchmark']['source_manifest']['measurement_plan_sha256'] = packet['measurement_plan_sha256']
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_pilot_observations_do_not_decide_formal_adoption(self):
        import statistics
        packet = copy.deepcopy(self.original)
        for block in packet['pilot']['blocks']:
            block['candidate']['operations']['tpsa']['first_use'] *= 100
        packet['pilot']['source_medians']['tpsa/first_use'] = statistics.median(
            v['baseline']['operations']['tpsa']['first_use'] / v['candidate']['operations']['tpsa']['first_use']
            for v in packet['pilot']['blocks'])
        self.assertEqual(audit.check(packet)['source_adoption_passed'], self.original['source_adoption_passed'])


class NumericPacketV3Tests(NumericPacketV2Tests):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-python-ring-cache-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))


class RejectedPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-morgan-tpsa-hot-scalars-pilot-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))

    def test_stopped_pilot_is_never_adopted(self):
        self.assertFalse(audit.check(self.original)['source_adoption_passed'])
        packet = copy.deepcopy(self.original)
        packet['source_adoption_passed'] = True
        with self.assertRaises(AssertionError): audit.check(packet)

    def test_changed_pilot_duration_is_rejected(self):
        packet = copy.deepcopy(self.original)
        packet['pilot']['blocks'][0]['candidate']['operations']['tpsa']['first_use'] *= 100
        with self.assertRaises(AssertionError): audit.check(packet)


class RejectedCompactPilotTests(RejectedPilotTests):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parents[2] / 'benchmarks/2026-10-19-morgan-tpsa-compact-cache-pilot-numeric.json.gz'
        cls.original = json.loads(gzip.decompress(path.read_bytes()))


class AdoptionCliTests(unittest.TestCase):
    def test_rejected_source_cannot_pass_required_adoption_cli(self):
        import subprocess
        import sys
        root = Path(__file__).parents[2]
        result = subprocess.run([
            sys.executable, str(root/'scripts/check_python_speed_frontier_packet.py'),
            str(root/'benchmarks/2026-10-19-morgan-direct-connectivity-numeric.json.gz'),
            '--require-adoption',
        ], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('source adoption gate failed', result.stderr)


if __name__ == '__main__': unittest.main()
