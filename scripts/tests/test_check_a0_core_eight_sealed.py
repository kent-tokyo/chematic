import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/check_a0_core_eight_sealed.py"
SPEC = importlib.util.spec_from_file_location("check_a0_core_eight_sealed", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def checked_in_summary():
    return json.loads(MODULE.DEFAULT_SUMMARY.read_text(encoding="utf-8"))


def test_checked_in_a0_summary_is_accepted():
    assert MODULE.validate(checked_in_summary()) == []


def test_gate_rejects_one_field_mismatch():
    summary = checked_in_summary()
    summary["fields"]["tpsa"]["strict_matches"] = 7_999
    summary["fields"]["tpsa"]["mismatches"] = 1
    assert "tpsa is not strict-green on all rows" in MODULE.validate(summary)


def test_gate_rejects_incomplete_accounting_and_missing_provenance():
    summary = checked_in_summary()
    summary["accounting"]["candidate_parse_failures"] = 1
    del summary["provenance"]["candidate_binary_sha256"]
    errors = MODULE.validate(summary)
    assert "row accounting is incomplete" in errors
    assert "provenance candidate_binary_sha256 is missing" in errors
