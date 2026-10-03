"""Fixed-row A6 source/published outcome comparison."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "compare_a6_source_published_rows.py"
SPEC = importlib.util.spec_from_file_location("a6_compare", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def rows() -> list[dict]:
    return [
        {
            "row_index": i,
            "name": f"row_{i}",
            "smiles": "CC",
            "arm": "chematic_pipeline_v2_mmff94_strict_stereo_safe",
            "status": "success",
        }
        for i in range(265)
    ]


def test_retention_is_separate_from_full_quality():
    source = rows()
    published = rows()
    source[53]["status"] = "typed_failure"
    source[53]["failure_cause"] = "stereo violation"
    result = module.compare(source, published)
    assert result["source_preserves_published_successes"] is False
    assert result["lost_published_successes"][0]["row_index"] == 53
    assert result["full_a6_quality_gate"] == "not_run"


def test_published_failure_is_not_counted_as_source_regression():
    source = rows()
    published = rows()
    published[53]["status"] = "typed_failure"
    source[53]["status"] = "typed_failure"
    assert module.compare(source, published)["source_preserves_published_successes"] is True
    source[53]["status"] = "success"
    assert module.compare(source, published)["gained_successes"] == [53]


def test_corpus_mismatch_is_rejected():
    source = rows()
    published = rows()
    source[0]["smiles"] = "CO"
    with pytest.raises(ValueError, match="corpus mismatch"):
        module.compare(source, published)


def test_row_count_mismatch_is_rejected_before_zipping():
    with pytest.raises(ValueError, match="different row counts"):
        module.compare(rows()[:-1], rows())
