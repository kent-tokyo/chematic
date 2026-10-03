"""Unit checks for the 265-row A6 source-wheel accounting diagnostic."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "summarize_a6_mmff94_termination.py"
SPEC = importlib.util.spec_from_file_location("a6_source_termination", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def rows() -> list[dict]:
    return [
        {
            "row_index": index,
            "arm": module.ARM,
            "status": "success",
            "force_field": {
                "converged": index == 0,
                "mmff94_termination": (
                    "gradient_converged" if index == 0 else "iteration_limit"
                ),
            },
            "final_validation": {"sound": True},
        }
        for index in range(265)
    ]


def test_complete_success_accounting_is_only_a_retention_gate():
    report = module.summarize(rows())
    assert report["success_retention_gate"] is True
    assert report["full_a6_quality_gate"] == "not_run"
    assert report["success_terminations"] == {
        "gradient_converged": 1,
        "iteration_limit": 264,
    }


def test_typed_failure_is_counted_and_blocks_retention():
    sample = rows()
    sample[246] = {
        "row_index": 246,
        "arm": module.ARM,
        "status": "typed_failure",
        "failure_cause": "excessive_residual_force",
    }
    report = module.summarize(sample)
    assert report["success_retention_gate"] is False
    assert report["failures"] == [
        {
            "row_index": 246,
            "status": "typed_failure",
            "cause": "excessive_residual_force",
        }
    ]


def test_missing_or_inconsistent_termination_is_rejected():
    sample = rows()
    sample[0]["force_field"]["mmff94_termination"] = "iteration_limit"
    with pytest.raises(ValueError, match="disagree"):
        module.summarize(sample)
    sample = rows()
    sample[0]["force_field"]["mmff94_termination"] = None
    with pytest.raises(ValueError, match="unknown"):
        module.summarize(sample)
