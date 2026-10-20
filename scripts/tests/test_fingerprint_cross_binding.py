"""Contract checks for the consolidated fingerprint parity runner."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from scripts.validation import check_fingerprint_cross_binding as gate


class FingerprintCrossBindingTests(unittest.TestCase):
    def test_current_runtime_versions_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = root / "crates/chematic-wasm/pkg-node/package.json"
            package.parent.mkdir(parents=True)
            package.write_text('{"version": "1.2.3"}', encoding="utf-8")
            (package.parent / "chematic_wasm.js").write_text("export {};", encoding="utf-8")
            (package.parent / "chematic_wasm_bg.wasm").write_bytes(b"wasm")
            with patch.object(gate, "ROOT", root), patch.object(
                gate.importlib, "import_module", return_value=SimpleNamespace(__version__="1.2.3", __file__=__file__)
            ):
                provenance = gate.runtime_provenance("1.2.3")
                self.assertEqual(provenance["node_package_version"], "1.2.3")
            with patch.object(gate, "ROOT", root), patch.object(
                gate.importlib, "import_module", return_value=SimpleNamespace(__version__="1.2.4", __file__=__file__)
            ):
                with self.assertRaisesRegex(ValueError, "Node/WASM package version"):
                    gate.runtime_provenance("1.2.4")
            with patch.object(gate, "ROOT", root), patch.object(
                gate.importlib, "import_module", return_value=SimpleNamespace(__version__="1.2.2", __file__=__file__)
            ):
                with self.assertRaisesRegex(ValueError, "Python chematic version"):
                    gate.runtime_provenance("1.2.3")

    def test_profile_payloads_keep_the_old_fields_and_order(self) -> None:
        class Molecule:
            def rdkit_ecfp4_detail(self):
                return None, {42: 3, 4: 1}, {42: [(3, 2), (1, 0)]}, {7: [(2, 1), (0, 0)]}

            def path_fp(self):
                return b"\x01\xff"

            ecfp4 = path_fp
            rdkit_ecfp4 = path_fp
            rdkit_torsion_fp = path_fp
            topo_path_fp = path_fp
            torsion_fp = path_fp

        molecule = Molecule()
        expected = {
            "rdkit-ecfp4-bitinfo": [[7, [[0, 0], [2, 1]]]],
            "rdkit-ecfp4-raw-bitinfo": [[42, [[1, 0], [3, 2]]]],
            "rdkit-ecfp4-sparse": [[4, 1], [42, 3]],
        }
        for name, value in expected.items():
            with self.subTest(profile=name):
                self.assertEqual(gate.payload(molecule, gate.PROFILES[name]), value)
        for name in ("ecfp4", "rdkit-ecfp4", "rdkit-path", "rdkit-torsion", "topo-path", "torsion"):
            with self.subTest(profile=name):
                self.assertEqual(gate.payload(molecule, gate.PROFILES[name]), "01ff")

    def test_row_identity_and_payload_are_required(self) -> None:
        profile = gate.PROFILES["rdkit-path"]
        valid = [{"index": 2, "smiles": "CC", "status": "ok", "path_hex": "01"}]
        gate.validate_rows(valid, [(2, "CC")], profile.field, "rust")
        for invalid in (
            [],
            [{**valid[0], "index": 1}],
            [{**valid[0], "smiles": "CO"}],
            [{**valid[0], "status": "unknown"}],
            [{key: value for key, value in valid[0].items() if key != "path_hex"}],
        ):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                gate.validate_rows(invalid, [(2, "CC")], profile.field, "rust")

    def test_error_comparison_retains_profile_boundary(self) -> None:
        left = [{"index": 0, "smiles": "bad", "status": "error", "error": "Rust parse"}]
        right = [{**left[0], "error": "Python parse"}]
        self.assertEqual(gate.compare(left, right, gate.PROFILES["rdkit-path"])["mismatches"], 0)
        self.assertEqual(gate.compare(left, right, gate.PROFILES["topo-path"])["mismatches"], 1)
        right[0]["smiles"] = "different"
        self.assertEqual(gate.compare(left, right, gate.PROFILES["rdkit-path"])["mismatches"], 1)

    def test_report_accounts_for_all_lanes_and_fails_on_a_mismatch(self) -> None:
        profile = gate.PROFILES["rdkit-ecfp4-sparse"]
        rust = [{"index": 0, "smiles": "CC", "status": "ok", "sparse_counts": [[1, 2]]}]
        python = [{**rust[0], "sparse_counts": [[1, 1]]}]
        with tempfile.TemporaryDirectory() as directory:
            corpus = Path(directory) / "input.smi"
            corpus.write_text("CC\n", encoding="utf-8")
            with patch.object(gate, "external_records", side_effect=[rust, rust]) as external, patch.object(gate, "python_records", return_value=python):
                report = gate.build_report(profile, corpus, "CC\n", "1.2.3")
        self.assertEqual(report["schema_version"], 1)
        self.assertEqual(report["target_version"], "1.2.3")
        self.assertEqual(report["corpus"]["rows"], 1)
        self.assertEqual(report["binding_status_counts"]["rust"]["ok"], 1)
        self.assertFalse(report["gate_passed"])
        self.assertEqual(report["pairwise"]["rust_vs_python"]["mismatches"], 1)
        self.assertEqual(report["pairwise"]["rust_vs_node_wasm"]["mismatches"], 0)
        self.assertEqual(external.call_count, 2)


if __name__ == "__main__":
    unittest.main()
