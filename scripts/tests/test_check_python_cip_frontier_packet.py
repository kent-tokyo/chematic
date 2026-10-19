import copy
import gzip
import importlib.util
import json
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('cip_packet', Path(__file__).parents[1] / 'check_python_cip_frontier_packet.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class CipPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads(gzip.decompress((Path(__file__).parents[2] / 'benchmarks/2026-10-19-cip-lazy-preparation-numeric.json.gz').read_bytes()))

    def reject(self, mutate):
        packet = copy.deepcopy(self.packet)
        mutate(packet)
        with self.assertRaises(AssertionError):
            audit.check(packet)

    def test_real_packet(self):
        self.assertEqual(audit.check(self.packet), 24)

    def test_partial_run(self):
        self.reject(lambda p: p['measurement'].update(complete=False))

    def test_reordered_block(self):
        self.reject(lambda p: p['measurement']['blocks'][0].update(order=['candidate', 'candidate', 'rdkit']))

    def test_changed_native_binary(self):
        self.reject(lambda p: p['measurement']['blocks'][0]['candidate'].update(extension_sha256={}))

    def test_stale_source(self):
        self.reject(lambda p: p['measurement']['plan']['identities']['candidate']['source_files'].update({'Cargo.toml': '0'*64}))

    def test_nonpositive_timing(self):
        self.reject(lambda p: p['measurement']['blocks'][0]['candidate']['operations']['cip'].update(first_use=0))

    def test_dropped_correctness_row(self):
        self.reject(lambda p: p['correctness'].update(source_exact_rows=14071))

    def test_stress_difference_cannot_be_ranked(self):
        self.reject(lambda p: p['correctness'].update(stress_ranked=True))

    def test_invented_speedup(self):
        self.reject(lambda p: p['measurement']['summary']['cip/first_use']['rdkit'].update(median=999))


if __name__ == '__main__':
    unittest.main()
