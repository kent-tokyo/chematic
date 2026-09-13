#!/usr/bin/env python3
"""Focused tests for the fail-closed descriptor evidence validators."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROMOTION = load_script("check_descriptor_promotion_gate.py")
ACCURACY = load_script("check_rdkit_accuracy_gate.py")
BASELINE_CANDIDATE = load_script("check_descriptor_baseline_candidate.py")


class AccuracyGateTests(unittest.TestCase):
    def complete_fields(self, rows: int = 2) -> dict:
        return {
            name: {
                "parsed": rows,
                "matches": rows,
                "strict_matches": rows,
                "mismatches": 0,
            }
            for name in PROMOTION.REQUIRED_FIELDS
        }

    def complete_holdout(self, rows: int = 2) -> dict:
        return {
            "rows": rows,
            "fields": {
                name: {
                    "rows": rows,
                    "strict_passed": rows,
                    "failed": 0,
                    "checks": [
                        {
                            "id": f"case-{i}",
                            "chematic": 1.0,
                            "rdkit": 1.0,
                            "absolute_error": 0.0,
                        }
                        for i in range(rows)
                    ],
                }
                for name in PROMOTION.REQUIRED_FIELDS
            },
        }

    def test_descriptor_schema_rejects_missing_or_extra_fields(self):
        fields = self.complete_fields()
        self.assertTrue(PROMOTION.strict_field_checks(fields, 2))
        fields.pop("hba")
        fields["unexpected"] = {}
        self.assertFalse(PROMOTION.strict_field_checks(fields, 2))

    def test_holdout_requires_all_fields_and_unique_rows(self):
        holdout = self.complete_holdout()
        self.assertEqual(ACCURACY.validate_holdout(holdout, list(PROMOTION.REQUIRED_FIELDS)), [])
        holdout["fields"]["hba"]["checks"][1]["id"] = "case-0"
        errors = ACCURACY.validate_holdout(holdout, list(PROMOTION.REQUIRED_FIELDS))
        self.assertTrue(any("duplicate row id" in error for error in errors))

    def test_holdout_rejects_empty_and_legacy_shape(self):
        self.assertTrue(ACCURACY.validate_holdout({}, list(PROMOTION.REQUIRED_FIELDS)))
        legacy = {"rows": 2, "passed": 2, "failed": 0, "checks": []}
        self.assertTrue(ACCURACY.validate_holdout(legacy, list(PROMOTION.REQUIRED_FIELDS)))

    def test_holdout_rejects_boolean_row_counts_and_unhashable_ids(self):
        holdout = self.complete_holdout()
        holdout["rows"] = True
        self.assertTrue(ACCURACY.validate_holdout(holdout, list(PROMOTION.REQUIRED_FIELDS)))
        holdout = self.complete_holdout()
        holdout["fields"]["hba"]["checks"][1]["id"] = []
        errors = ACCURACY.validate_holdout(holdout, list(PROMOTION.REQUIRED_FIELDS))
        self.assertTrue(any("missing or duplicate row id" in error for error in errors))

    def test_holdout_rejects_non_finite_values(self):
        holdout = self.complete_holdout()
        holdout["fields"]["hba"]["checks"][0]["absolute_error"] = float("nan")
        errors = ACCURACY.validate_holdout(holdout, list(PROMOTION.REQUIRED_FIELDS))
        self.assertTrue(any("absolute_error is non-finite" in error for error in errors))

    def test_diagnostics_rejects_bool_counts_extra_fields_and_missing_provenance(self):
        manifest = {
            "descriptor_contract": {
                "required_fields": list(PROMOTION.REQUIRED_FIELDS),
                "strict_tolerances": {name: 0.0 for name in PROMOTION.REQUIRED_FIELDS},
                "corpus": "fixture.smi",
            }
        }
        diagnostics = {
            "schema_version": 1,
            "profile": "test",
            "chematic_version": "1.0.13",
            "corpus": "fixture.smi",
            "rows": 2,
            "parsed": 2,
            "parse_failures": 0,
            "rdkit_version": "2025.09.3",
            "unsupported_values": 0,
            "fields": self.complete_fields(),
        }
        diagnostics["fields"]["hba"]["parsed"] = True
        diagnostics["fields"]["extra"] = {}
        diagnostics["chematic_version"] = None
        errors = ACCURACY.validate(manifest, diagnostics, self.complete_holdout())
        self.assertTrue(any("chematic version is missing" in error for error in errors))
        self.assertTrue(any("fields must exactly match" in error for error in errors))
        self.assertTrue(any("parsed must be a non-negative integer" in error for error in errors))

    def test_diagnostics_rejects_non_finite_aggregate_and_missing_raw_row(self):
        manifest = {
            "descriptor_contract": {
                "required_fields": list(PROMOTION.REQUIRED_FIELDS),
                "strict_tolerances": {name: 0.0 for name in PROMOTION.REQUIRED_FIELDS},
                "corpus": "fixture.smi",
            }
        }
        diagnostics = {
            "schema_version": 1,
            "profile": "test",
            "chematic_version": "1.0.13",
            "corpus": "fixture.smi",
            "rows": 2,
            "parsed": 2,
            "parse_failures": 0,
            "rdkit_version": "2025.09.3",
            "unsupported_values": 0,
            "fields": self.complete_fields(),
        }
        diagnostics["fields"]["hba"]["mae"] = float("inf")
        errors = ACCURACY.validate(manifest, diagnostics, self.complete_holdout())
        self.assertTrue(any("mae is absent or non-finite" in error for error in errors))

        diagnostics["schema_version"] = 2
        diagnostics["raw_rows"] = [{"index": 0, "smiles": "C", "status": "ok"}]
        errors = ACCURACY.validate(manifest, diagnostics, self.complete_holdout())
        self.assertTrue(any("raw_rows must cover every diagnostics row" in error for error in errors))

    def test_diagnostics_rejects_more_than_fifty_strict_mismatches(self):
        rows = 60
        manifest = {
            "descriptor_contract": {
                "required_fields": list(PROMOTION.REQUIRED_FIELDS),
                "strict_tolerances": {name: 0.0 for name in PROMOTION.REQUIRED_FIELDS},
                "corpus": "fixture.smi",
            }
        }
        diagnostics = {
            "schema_version": 1,
            "profile": "test",
            "chematic_version": "1.0.13",
            "corpus": "fixture.smi",
            "rows": rows,
            "parsed": rows,
            "parse_failures": 0,
            "rdkit_version": "2025.09.3",
            "unsupported_values": 0,
            "fields": self.complete_fields(rows),
        }
        diagnostics["fields"]["hba"].update(
            {"matches": rows - 51, "strict_matches": rows - 51, "mismatches": 51}
        )
        errors = ACCURACY.validate(manifest, diagnostics, self.complete_holdout(rows))
        self.assertTrue(
            any("strict mismatch or incomplete coverage" in error for error in errors),
            errors,
        )

    def test_baseline_candidate_rejects_reused_artifacts_and_accepts_distinct_arms(self):
        def arm(role: str, suffix: str) -> dict:
            return {
                "schema_version": 1,
                "role": role,
                "contract": "rdkit_descriptor_semantics_v1",
                "profile": "rdkit_compat_v2",
                "provenance": {
                    "source_commit": f"commit-{suffix}",
                    "source_tree_sha256": f"tree-{suffix}",
                    "artifact_sha256": f"artifact-{suffix}",
                    "rdkit_version": "2025.09.3",
                    "corpus_sha256": "corpus-fixed",
                    "corpus_rows": 5000,
                },
                "fields": self.complete_fields(5000),
            }

        baseline = arm("baseline", "b")
        candidate = arm("candidate", "c")
        self.assertEqual(BASELINE_CANDIDATE.validate_pair(baseline, candidate), [])
        baseline["fields"]["molecular_weight"]["parsed"] = 4999
        self.assertEqual(BASELINE_CANDIDATE.validate_pair(baseline, candidate), [])
        candidate["provenance"]["artifact_sha256"] = baseline["provenance"]["artifact_sha256"]
        errors = BASELINE_CANDIDATE.validate_pair(baseline, candidate)
        self.assertTrue(any("same artifact" in error for error in errors))

        candidate["provenance"]["artifact_sha256"] = "artifact-c"
        candidate["provenance"]["rdkit_version"] = "2025.09.2"
        errors = BASELINE_CANDIDATE.validate_pair(baseline, candidate)
        self.assertTrue(any("rdkit_version differs" in error for error in errors))

        candidate["provenance"]["rdkit_version"] = "2025.09.3"
        candidate["provenance"]["corpus_sha256"] = "different-corpus"
        errors = BASELINE_CANDIDATE.validate_pair(baseline, candidate)
        self.assertTrue(any("corpus_sha256 differs" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
