"""A6 per-platform quality gate on the independent scorer (#739)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "a6_platform_gate.py"
SPEC = importlib.util.spec_from_file_location("a6_platform_gate", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

ARM = module.ARM


def rows(n: int = 265) -> list[dict]:
    return [
        {"tier": "A", "row_index": i, "name": f"row_{i}", "arm": ARM, "status": "success",
         "coords": [[float(i), 0.0, 0.0]]}
        for i in range(n)
    ]


def scored(n: int = 265) -> list[dict]:
    return [
        {"tier": "A", "name": f"row_{i}", "arm": ARM, "engine": "chematic", "status": "scored",
         "independently_sound": True, "gross_clash_count": 0,
         "stereo": {"declared": 1, "satisfied": 1, "violated": 0, "unevaluable": 0}}
        for i in range(n)
    ]


def test_all_sound_rows_pass():
    result = module.gate(rows(), scored(), 265)
    assert result["passed"] is True
    assert (result["successes"], result["independently_sound"], result["stereo_clean"]) == (265, 265, 265)


def test_final_stereo_violation_fails_the_gate():
    # Published v1.0.31 on Linux: rows 53 and 246 ended in FinalStereoViolation.
    data = rows()
    for i in (53, 246):
        data[i] = {**data[i], "status": "typed_failure", "coords": None,
                   "failure_cause": "PipelineV2Error: FinalStereoViolation"}
    result = module.gate(data, scored(), 265)
    assert result["passed"] is False
    assert result["successes"] == 263
    assert [f["row_index"] for f in result["failures"]] == [53, 246]


def test_scorer_violations_clashes_and_missing_rows_fail():
    bad = scored()
    bad[3]["stereo"] = {"violated": 1, "unevaluable": 0}
    bad[4]["gross_clash_count"] = 2
    bad[5]["independently_sound"] = False
    del bad[6]
    result = module.gate(rows(), bad, 265)
    assert result["passed"] is False
    assert result["stereo_violated_or_unevaluable"] == ["row_3"]
    assert result["clashing"] == ["row_4"]
    assert result["unsound"] == ["row_5"]
    assert result["unscored"] == ["row_6"]
    assert module.gate(rows(264), scored(264), 265)["passed"] is False


def test_coordinate_comparison_is_row_by_row():
    other = rows()
    other[7] = {**other[7], "coords": [[7.000000000000001, 0.0, 0.0]]}
    report = module.compare(rows(), other)
    assert report["rows_identical"] == 264
    assert report["rows_different"] == ["row_7"]
    assert module.coordinates_digest(rows()) != module.coordinates_digest(other)
