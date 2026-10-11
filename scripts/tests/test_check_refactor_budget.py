from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

import check_refactor_budget as budget  # noqa: E402

function_metrics = budget.function_metrics
mask_rust = budget.mask_rust


class RefactorBudgetScannerTests(unittest.TestCase):
    def test_masks_literals_and_comments_without_losing_braces(self) -> None:
        source = '''
fn sample(value: char) {
    let byte = b'a';
    let quote = '\\'';
    let raw = r###"{ not a body }"###;
    let text = "also { not a body }";
    // }
    /* { nested /* } */ comment } */
    if value == 'x' { println!("ok"); }
}
'''
        masked = mask_rust(source)
        self.assertEqual(masked.count("\n"), source.count("\n"))
        self.assertEqual(masked.count("{"), 2)
        self.assertEqual(masked.count("}"), 2)

    def test_scans_impl_methods_and_ignores_cfg_test_module(self) -> None:
        source = '''
impl Demo {
    pub(crate) async fn calculate(&self) -> usize {
        if true { 1 } else { 0 }
    }
}

#[cfg(test)]
mod tests {
    fn intentionally_long_test_helper() {
        if true { loop { break; } }
    }
}
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.rs"
            path.write_text(source, encoding="utf-8")
            metrics = function_metrics(path, "sample.rs")
        self.assertEqual([metric.name for metric in metrics], ["calculate"])
        self.assertEqual(metrics[0].decisions, 1)

    def test_duplicate_function_names_receive_stable_occurrences(self) -> None:
        source = "fn new() {}\nmod nested {\n    fn new() {}\n}\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.rs"
            path.write_text(source, encoding="utf-8")
            metrics = function_metrics(path, "sample.rs")
        self.assertEqual(
            [metric.key for metric in metrics],
            ["sample.rs::new#1", "sample.rs::new#2"],
        )

    def test_ratchet_rejects_growth_and_requires_resolved_allowance_cleanup(self) -> None:
        seed = {
            "schema_version": 1,
            "limits": {
                "max_file_lines": 30,
                "max_function_lines": 3,
                "max_function_decisions": 1,
            },
            "excluded_paths": [],
            "file_allowances": {},
            "function_allowances": {},
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "crates" / "demo" / "src" / "lib.rs"
            source.parent.mkdir(parents=True)
            source.write_text(
                "fn calculate() {\n    if true {\n        work();\n    }\n}\n",
                encoding="utf-8",
            )
            with patch.object(budget, "ROOT", root):
                baseline = budget.generated_baseline(seed)
                self.assertEqual(budget.check(baseline), [])

                source.write_text(
                    "fn calculate() {\n"
                    "    if true {\n"
                    "        work();\n"
                    "        work();\n"
                    "        work();\n"
                    "    }\n"
                    "}\n",
                    encoding="utf-8",
                )
                self.assertTrue(
                    any("function grew" in error for error in budget.check(baseline))
                )

                source.write_text("fn calculate() {}\n", encoding="utf-8")
                self.assertTrue(
                    any(
                        "resolved function allowance" in error
                        for error in budget.check(baseline)
                    )
                )

    def test_dedicated_test_modules_are_not_production_sources(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_dir = root / "crates" / "demo" / "src"
            source_dir.mkdir(parents=True)
            (source_dir / "lib.rs").write_text("fn production() {}\n", encoding="utf-8")
            (source_dir / "boundary_tests.rs").write_text(
                "fn test_helper() {}\n", encoding="utf-8"
            )
            nested = source_dir / "feature" / "tests.rs"
            nested.parent.mkdir()
            nested.write_text("fn test_helper() {}\n", encoding="utf-8")
            sources = budget.rust_sources(root, set())
        self.assertEqual(
            [path.relative_to(root).as_posix() for path in sources],
            ["crates/demo/src/lib.rs"],
        )

    def test_exclusions_must_be_generated_crate_sources(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = root / "crates" / "demo" / "src" / "table.rs"
            generated.parent.mkdir(parents=True)
            generated.write_text("// Generated from the reference table.\n", encoding="utf-8")
            self.assertEqual(
                budget.validate_excluded_paths(root, ["crates/demo/src/table.rs"]),
                [],
            )

            generated.write_text("// Maintained by hand.\n", encoding="utf-8")
            errors = budget.validate_excluded_paths(
                root,
                ["crates/demo/src/table.rs", "crates/demo/src/table.rs", "missing.rs"],
            )
        self.assertTrue(any("duplicate excluded path" in error for error in errors))
        self.assertTrue(any("lacks a generated header" in error for error in errors))
        self.assertTrue(any("outside crate production Rust" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
