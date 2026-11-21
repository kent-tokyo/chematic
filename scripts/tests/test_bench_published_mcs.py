import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

SPEC = importlib.util.spec_from_file_location(
    "bench_published_mcs", Path(__file__).parents[1] / "benchmarks" / "bench_published_mcs.py"
)
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


class PublishedMcsTests(unittest.TestCase):
    def test_query_graph_predicates_must_match(self):
        from rdkit import Chem

        first = Chem.MolFromSmarts("[#6]-[#6]-[#7]")
        reordered = Chem.MolFromSmarts("[#7]-[#6]-[#6]")
        different_atom = Chem.MolFromSmarts("[#8]-[#6]-[#6]")
        different_bond = Chem.MolFromSmarts("[#7]=[#6]-[#6]")
        self.assertTrue(bench.same_query_predicates(first, reordered, Chem))
        self.assertFalse(bench.same_query_predicates(first, different_atom, Chem))
        self.assertFalse(bench.same_query_predicates(first, different_bond, Chem))

    def test_record_suppresses_ratios_for_invalid_or_different_results(self):
        path = Path(__file__).parents[2] / "benchmarks/2026-11-21-v122-published-mcs.json"
        report = json.loads(path.read_text())
        self.assertIsNone(report["paired_repeated_search"]["v120"][1]["competitor_over_v122"])
        for arm in ("v120", "v121", "rdkit", "cosmolkit"):
            self.assertIsNone(report["paired_repeated_search"][arm][4]["competitor_over_v122"])
        self.assertIsNotNone(report["paired_repeated_search"]["v120"][0]["competitor_over_v122"])

    def test_uniform_smarts_record_excludes_conversion_errors(self):
        path = Path(__file__).parents[2] / "benchmarks/2026-11-21-v122-published-mcs-uniform.json"
        report = json.loads(path.read_text())
        self.assertTrue(report["uniform_smarts"])
        self.assertTrue(all(row["passes_all_blocks"] for row in report["case_summary"]["v122"]))
        for arm in ("v120", "v121"):
            self.assertIn("output error:", report["blocks"][0][arm]["validation"][4]["error"])
            self.assertIsNone(report["paired_repeated_search"][arm][4]["competitor_over_v122"])

    def test_verify_wheel_rejects_an_installed_file_change(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            wheel_dir = root / "wheels"
            wheel_dir.mkdir()
            installed = root / "installed" / "chematic"
            installed.mkdir(parents=True)
            wheel = wheel_dir / "chematic-1.2.2-cp313-cp313-macosx_11_0_arm64.whl"
            with zipfile.ZipFile(wheel, "w") as archive:
                archive.writestr("chematic/__init__.py", b"original")
            target = installed / "__init__.py"
            target.write_bytes(b"original")
            result = bench.verify_wheel("v122", wheel_dir, installed.parent)
            self.assertEqual(result["installed_files_verified"], 1)
            target.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "installed file differs"):
                bench.verify_wheel("v122", wheel_dir, installed.parent)

    def test_verify_wheel_pins_the_requested_cosmolkit_version(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            installed = root / "cosmolkit"
            installed.mkdir()
            (installed / "__init__.py").write_bytes(b"new")
            for version, content in (("0.6.2", b"old"), ("0.6.3", b"new")):
                with zipfile.ZipFile(root / f"cosmolkit-{version}-test.whl", "w") as archive:
                    archive.writestr("cosmolkit/__init__.py", content)
            self.assertIn("0.6.3", bench.verify_wheel("cosmolkit", root, root, "0.6.3")["filename"])
            with self.assertRaisesRegex(ValueError, "installed file differs"):
                bench.verify_wheel("cosmolkit", root, root, "0.6.2")


if __name__ == "__main__":
    unittest.main()
