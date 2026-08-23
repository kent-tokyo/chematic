#!/usr/bin/env python3
"""
Regression tests for `gen_pipeline_v2_vs_rdkit_report.py`'s
`GATE_WIDENING_EXPLANATIONS` -- the classification logic behind the
`newly_passing_unexplained` assertion that blocked report generation during
A1's fresh 265-molecule re-run (`docs/rfcs/a1_conformer_benchmark_failure_ledger.md`,
Finding 3: `chembl_tier_b_0030` flipped timeout->success under a gate
widening that turned out to be a mechanical no-op for that molecule, a
genuine explanation the assertion's original single check didn't recognize).

Three fixtures, not one:

1. POSITIVE (`uff_fallback_rescue`, pre-existing category, unchanged logic
   minus a fragile substring-match removal): a stricter policy times out,
   the wider gate's policy falls back to UFF because the newly-gated
   coverage dimension is missing. Synthetic data (this project's own fresh
   265-corpus re-run happened not to contain a live example of this
   specific category at the time this test was written -- the ONLY
   newly-passing case in the current committed data is chembl_tier_b_0030,
   which is the OTHER category below) -- still exercises the exact typed
   fields the check reads.
2. POSITIVE (`identical_coverage_timing_variance`, the new category): a
   frozen snapshot of the real row pair for `chembl_tier_b_0030`
   (`chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated` ->
   `..._complete_bonded_term_gated`) that originally tripped the assertion
   during the A1/A2 closeout round (2026-08-23) -- captured once, not
   re-read from the live rows file. This is a regression test of the
   classification LOGIC, not of that molecule's current benchmark outcome:
   `chembl_tier_b_0030` sits right at the 20s wall-clock timeout boundary
   under these two (otherwise unmodified) arms, so its real status flips
   between runs/machines on pure timing noise -- exactly the phenomenon
   this category exists to recognize, and exactly why re-reading it live
   would make this test assert on benchmark timing variance instead of on
   `_verify_identical_coverage_timing_variance`'s own logic. Found flipping
   for real during the best-of-N benchmark round (2026-08-24): timed out at
   21946ms in one run, succeeded cleanly at 6643ms in the next, both on
   identical unmodified code.
3. NEGATIVE (must keep failing both checks -- the one that matters, per
   this project's established two-control convention): a case that
   superficially resembles a rescue (earlier row times out, later row
   succeeds) but matches NEITHER recognized explanation -- later row still
   shows missing coverage in the gated dimension (not a fallback, not a
   verified no-op). A version of this check that silently widened to accept
   "any timeout -> success flip" would incorrectly pass this fixture.

Usage: python scripts/tests/test_gate_stage_delta_explanations.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from gen_pipeline_v2_vs_rdkit_report import (  # noqa: E402
    GATE_WIDENING_EXPLANATIONS,
    TOTAL_TIMEOUT_MS,
    _verify_identical_coverage_timing_variance,
    _verify_uff_fallback_rescue,
)


def positive_uff_fallback_rescue():
    earlier_row = {
        "status": "timeout",
        "failure_cause": "Timeout",
        "failure_stage": "ForceFieldMinimization",
    }
    later_row = {
        "status": "success",
        "force_field_actual": "UffOnly",
        "force_field_fallback": True,
        "stretch_bend_missing_count": 3,
    }
    ok, reason = _verify_uff_fallback_rescue(earlier_row, later_row, ["stretch_bend_missing_count"])
    assert ok, f"expected uff_fallback_rescue to match, got reason: {reason}"

    matched = [cat for cat, fn in GATE_WIDENING_EXPLANATIONS if fn(earlier_row, later_row, ["stretch_bend_missing_count"])[0]]
    assert matched == ["uff_fallback_rescue"], f"expected exactly ['uff_fallback_rescue'], got {matched}"


def positive_identical_coverage_timing_variance_chembl_tier_b_0030():
    """Frozen snapshot (captured 2026-08-23, A1/A2 closeout round) of the
    real row pair that originally tripped the assertion -- deliberately
    NEVER re-read from the live, regenerated-per-benchmark-round rows file.
    This function tests `_verify_identical_coverage_timing_variance`'s
    classification logic in isolation; the molecule's CURRENT benchmark
    outcome is irrelevant to that and does genuinely change from run to
    run (it sits right at the 20s wall-clock timeout boundary under these
    two otherwise-unmodified arms -- confirmed flipping timeout<->success
    on identical code between the A1/A2 closeout round and the best-of-N
    round, 2026-08-24). Loading live data here would make this test assert
    on benchmark timing noise instead of on the classifier."""
    earlier_row = {
        "arm": "chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated",
        "elapsed_ms": 21946,
        "failure_cause": "Timeout",
        "failure_stage": "ForceFieldMinimization",
        "name": "chembl_tier_b_0030",
        "status": "timeout",
        "tier": "B",
    }
    later_row = {
        "arm": "chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated",
        "elapsed_ms": 18914,
        "force_field_actual": "Mmff94BondAngleStrict",
        "force_field_fallback": False,
        "force_field_fallback_reason": None,
        "name": "chembl_tier_b_0030",
        "oop_missing_count": 0,
        "status": "success",
        "stretch_bend_missing_count": 0,
        "tier": "B",
        "torsion_missing_count": 0,
    }

    ok, reason = _verify_identical_coverage_timing_variance(
        earlier_row, later_row, ["torsion_missing_count", "oop_missing_count"]
    )
    assert ok, f"expected chembl_tier_b_0030 to match identical_coverage_timing_variance, got reason: {reason}"
    # Must NOT also match uff_fallback_rescue -- the two categories'
    # preconditions on force_field_fallback are mutually exclusive.
    ok2, _ = _verify_uff_fallback_rescue(earlier_row, later_row, ["torsion_missing_count", "oop_missing_count"])
    assert not ok2, "chembl_tier_b_0030's real row pair must not ALSO match uff_fallback_rescue"

    matched = [
        cat
        for cat, fn in GATE_WIDENING_EXPLANATIONS
        if fn(earlier_row, later_row, ["torsion_missing_count", "oop_missing_count"])[0]
    ]
    assert matched == ["identical_coverage_timing_variance"], f"expected exactly one match, got {matched}"


def negative_unexplained_case_still_fails():
    """A case that looks superficially like a rescue but is neither a
    verified UFF fallback NOR a verified no-op gate: the later row succeeded
    without falling back, but STILL shows missing coverage in the exact
    dimension this stage gates on -- meaning the gate was NOT a mechanical
    no-op here, so the "identical computation, just timing noise" story
    does not hold, and there is no fallback to invoke the other category
    either. Must be rejected by both checks."""
    earlier_row = {
        "status": "timeout",
        "failure_cause": "Timeout",
        "failure_stage": "ForceFieldMinimization",
    }
    later_row = {
        "status": "success",
        "force_field_actual": "Mmff94BondAngleStrict",
        "force_field_fallback": False,
        "torsion_missing_count": 2,  # gate did NOT verify as a no-op
        "oop_missing_count": 0,
        "elapsed_ms": 19000,
    }
    for category, check_fn in GATE_WIDENING_EXPLANATIONS:
        ok, reason = check_fn(earlier_row, later_row, ["torsion_missing_count", "oop_missing_count"])
        assert not ok, f"expected {category} to reject this unexplained case, but it matched"


def negative_fast_success_is_not_corroborated_as_timing_variance():
    """Even with fully-covered gate dimensions and no fallback, a later row
    that finishes trivially fast (nowhere near the timeout boundary) must
    NOT be accepted as "the same computation, just noisy" -- that would be
    accepting an unrelated, unexplained speed difference under a name that
    implies verification."""
    earlier_row = {
        "status": "timeout",
        "failure_cause": "Timeout",
        "failure_stage": "ForceFieldMinimization",
    }
    later_row = {
        "status": "success",
        "force_field_actual": "Mmff94BondAngleStrict",
        "force_field_fallback": False,
        "torsion_missing_count": 0,
        "oop_missing_count": 0,
        "elapsed_ms": 50,  # nowhere near TOTAL_TIMEOUT_MS
    }
    ok, reason = _verify_identical_coverage_timing_variance(
        earlier_row, later_row, ["torsion_missing_count", "oop_missing_count"]
    )
    assert not ok, f"a fast, uncorroborated success must not be accepted, got ok=True (TOTAL_TIMEOUT_MS={TOTAL_TIMEOUT_MS})"


def main():
    positive_uff_fallback_rescue()
    positive_identical_coverage_timing_variance_chembl_tier_b_0030()
    negative_unexplained_case_still_fails()
    negative_fast_success_is_not_corroborated_as_timing_variance()
    print("OK: 2 positive fixtures matched their expected category exactly; "
          "2 negative fixtures correctly rejected by every recognized category.")


if __name__ == "__main__":
    sys.exit(main())
