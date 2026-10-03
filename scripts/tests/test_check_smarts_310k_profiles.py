"""Small contract tests for the pinned SMARTS profile accounting helper."""

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "check_smarts_310k_profiles.py"
SPEC = importlib.util.spec_from_file_location("check_smarts_310k_profiles", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_match_sets_ignore_embedding_and_match_order():
    assert MODULE.normalized_match_sets([[2, 1], [4]]) == MODULE.normalized_match_sets(
        [[4], [1, 2]]
    )


@pytest.mark.parametrize("bad", [None, {}, [[1, "2"]], [1]])
def test_match_sets_reject_malformed_rows(bad):
    with pytest.raises(ValueError, match="match sets"):
        MODULE.normalized_match_sets(bad)


def test_opt_in_typed_refusal_is_not_counted_as_false():
    assert MODULE.profile_matches({"parity_error": "RingModelAmbiguous", "parity": []}, "parity") is None
    assert MODULE.profile_matches({"parity_budget_exhausted": True, "parity": []}, "parity") is None
    assert MODULE.profile_matches({"shared_symmetrized": {"error": "ambiguous"}}, "shared_symmetrized") is None
    assert MODULE.profile_matches({"parity": []}, "parity") == []


def test_missing_profile_is_not_silently_typed_as_refusal():
    with pytest.raises(ValueError, match="neither matches nor a refusal"):
        MODULE.profile_matches({}, "parity")
    with pytest.raises(ValueError, match="profile is missing"):
        MODULE.profile_matches({}, "shared_symmetrized")
