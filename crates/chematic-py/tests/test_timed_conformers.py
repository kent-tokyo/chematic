import math

import pytest

import chematic


def test_timed_conformer_generation_reports_mmff_unavailable():
    result = chematic.generate_conformers_timed("CB(C)C", 2, timeout_seconds=15)
    assert result["status"] == "completed"
    assert result["embedded"] > 0
    assert result["force_field_error"]
    assert result["conformers"]
    assert all(row[0] == -1 and row[1] is None for row in result["conformers"])


def test_timed_conformer_generation_terminates_and_reports_timeout():
    result = chematic.generate_conformers_timed("C" * 80, 20, timeout_seconds=0.1)
    assert result["status"] == "timeout"
    assert result["elapsed_seconds"] < 10


def test_timed_conformer_generation_reports_input_failure():
    result = chematic.generate_conformers_timed("not smiles", 1, timeout_seconds=15)
    assert result["status"] == "error"
    assert result["error_type"]


def test_timed_conformer_generation_rejects_a_different_worker_package(monkeypatch):
    monkeypatch.setattr(chematic, "__file__", "/a/different/chematic/__init__.py")
    result = chematic.generate_conformers_timed("CCO", 1, timeout_seconds=15)
    assert result["status"] == "error"
    assert "different chematic package" in result["error"]


def test_timed_conformer_generation_reports_unevaluable_stereo():
    result = chematic.generate_conformers_timed(
        "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F", 2, timeout_seconds=15
    )
    assert result["status"] == "completed"
    assert result["embedded"] > 0
    assert result["conformers"] == []
    assert result["stereo_unevaluable"] == result["embedded"]


@pytest.mark.parametrize("timeout", [0, -1, math.inf, math.nan])
def test_timed_conformer_generation_rejects_invalid_deadlines(timeout):
    with pytest.raises(ValueError):
        chematic.generate_conformers_timed("CCO", 1, timeout_seconds=timeout)
