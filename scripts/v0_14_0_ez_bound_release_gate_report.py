#!/usr/bin/env python3
"""v0.14.0 release gate: E/Z bound constraint end-to-end effect on the 265-molecule
Tier A/B corpus, measured through the production `embed_pipeline_v2` entry point
(the same one `Mol.embed_pipeline_v2()`/`embed_pipeline_v2_json` call).

Two independent checks, per the release-gate plan:

1. Regression control: every pre-existing arm (the 11 that existed before this
   diagnostic arm was added) must be byte-identical in outcome to the last recorded
   v0.13.0 run -- `enforce_chirality: false` is untouched by the E/Z bound fix, so
   this run should show zero change.
2. Feature-effect measurement: `chematic_pipeline_v2_mmff94_strict_enforce_chirality`
   (new arm, `enforce_chirality: true`) vs `chematic_pipeline_v2_mmff94_strict`
   (existing baseline, `enforce_chirality: false`) -- identical in every other
   config field, so any delta between them isolates exactly what the E/Z bound
   constraint changes, end to end, through the real production pipeline.

Not a comparison against the `*_repair` (StereoPolicy::RepairAndVerify) arms --
those use a different, mutually-exclusive stereo mechanism (see pipeline_v2.rs's
InvalidConfiguration gate). This script only ever diffs same-run arm pairs or the
same arm across two runs -- never cross-arm-and-cross-run at once.

Usage:
    python3 scripts/v0_14_0_ez_bound_release_gate_report.py \\
        --new validation/results/pipeline_v2_vs_rdkit_v0_14_0_chematic_rows.jsonl \\
        --historical validation/results/pipeline_v2_vs_rdkit_v0_13_0_chematic_rows.jsonl
"""

import argparse
import json
import sys
from collections import defaultdict

BASELINE_ARM = "chematic_pipeline_v2_mmff94_strict"
NEW_ARM = "chematic_pipeline_v2_mmff94_strict_enforce_chirality"

# The 11 arms that existed before this release-gate diagnostic arm was added.
# `chematic_pipeline_v2_ring_torsion_failclosed_probe` and `chematic_legacy_etkdg`
# are also pre-existing but are separate one-off diagnostics, not part of
# PIPELINE_ARMS proper -- included here too since they're equally a regression
# control (enforce_chirality never touches them either).
PRE_EXISTING_ARMS = {
    "chematic_pipeline_v2_no_ff",
    "chematic_pipeline_v2_dreiding",
    "chematic_pipeline_v2_uff_only",
    "chematic_pipeline_v2_mmff94_strict",
    "chematic_pipeline_v2_mmff94_with_uff_fallback",
    "chematic_pipeline_v2_mmff94_strict_repair",
    "chematic_pipeline_v2_mmff94_with_uff_fallback_repair",
    "chematic_pipeline_v2_mmff94_strict_stretch_bend_gated",
    "chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated",
    "chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated",
    "chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated",
    "chematic_pipeline_v2_ring_torsion_failclosed_probe",
    "chematic_legacy_etkdg",
}

OUTCOME_FIELDS = (
    "status",
    "sound",
    "final_stereo_declared",
    "final_stereo_satisfied",
    "final_stereo_violations",
    "final_stereo_unevaluable",
)


def load_rows(path):
    rows = defaultdict(dict)  # (name, arm) -> row
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            rows[(row["name"], row["arm"])] = row
    return rows


def is_declared_ez(smiles):
    return "/" in smiles or "\\" in smiles


def summarize(rows, arm, predicate=None):
    """Per-molecule (name -> row) for one arm, optionally filtered by predicate(row)."""
    out = {}
    for (name, a), row in rows.items():
        if a != arm:
            continue
        if predicate is not None and not predicate(row):
            continue
        out[name] = row
    return out


def classify_transition(before_row, after_row):
    """One (before, after) row pair -> a short label describing what changed."""
    b_ok = before_row["status"] == "success"
    a_ok = after_row["status"] == "success"
    b_sound = bool(before_row.get("sound"))
    a_sound = bool(after_row.get("sound"))
    b_viol = before_row.get("final_stereo_violations", 0)
    a_viol = after_row.get("final_stereo_violations", 0)

    labels = []
    if b_ok and not a_ok:
        labels.append("success_to_failure")
    elif not b_ok and a_ok:
        labels.append("failure_to_success")

    if b_sound and not a_sound:
        labels.append("sound_to_unsound")
    elif not b_sound and a_sound:
        labels.append("unsound_to_sound")

    if b_ok and a_ok:
        if b_viol > 0 and a_viol == 0:
            labels.append("stereo_newly_fixed")
        elif b_viol == 0 and a_viol > 0:
            labels.append("stereo_newly_broken")

    return labels


def report_regression_control(new_rows, hist_rows):
    print("=" * 78)
    print("1. REGRESSION CONTROL — pre-existing arms, new run vs v0.13.0 historical")
    print("=" * 78)
    mismatches = []
    checked = 0
    for (name, arm), hist_row in hist_rows.items():
        if arm not in PRE_EXISTING_ARMS:
            continue
        new_row = new_rows.get((name, arm))
        if new_row is None:
            mismatches.append((name, arm, "MISSING in new run"))
            continue
        checked += 1
        diffs = []
        for field in OUTCOME_FIELDS:
            if hist_row.get(field) != new_row.get(field):
                diffs.append(f"{field}: {hist_row.get(field)!r} -> {new_row.get(field)!r}")
        if diffs:
            mismatches.append((name, arm, "; ".join(diffs)))

    print(f"Checked {checked} (molecule, arm) pairs across {len(PRE_EXISTING_ARMS)} arms.")
    if not mismatches:
        print("RESULT: zero mismatches. enforce_chirality: false path confirmed byte-identical "
              "in outcome to v0.13.0 for every pre-existing arm.")
    else:
        print(f"RESULT: {len(mismatches)} mismatch(es) found — investigate before declaring GO:")
        for name, arm, detail in mismatches:
            print(f"  [{arm}] {name}: {detail}")
    print()
    return mismatches


def report_feature_effect(new_rows, subset_name, predicate):
    baseline = summarize(new_rows, BASELINE_ARM, predicate)
    new = summarize(new_rows, NEW_ARM, predicate)
    names = sorted(set(baseline) & set(new))
    missing_baseline = sorted(set(new) - set(baseline))
    missing_new = sorted(set(baseline) - set(new))

    print("-" * 78)
    print(f"Subset: {subset_name} ({len(names)} molecules)")
    print("-" * 78)
    if missing_baseline or missing_new:
        print(f"  WARNING: {len(missing_baseline)} in new-arm-only, "
              f"{len(missing_new)} in baseline-only (should be 0/0)")

    # `final_stereo_*` fields are `null` (not 0, not absent) on a non-success row
    # (e.g. a timeout never reached the point of counting stereo elements) -- `or 0`
    # coerces that explicitly, `.get(field, 0)` alone does NOT (the key is present
    # with value None, so the default never kicks in).
    n_success_baseline = sum(1 for n in names if baseline[n]["status"] == "success")
    n_success_new = sum(1 for n in names if new[n]["status"] == "success")
    n_sound_baseline = sum(1 for n in names if baseline[n].get("sound"))
    n_sound_new = sum(1 for n in names if new[n].get("sound"))
    satisfied_baseline = sum(baseline[n].get("final_stereo_satisfied") or 0 for n in names)
    satisfied_new = sum(new[n].get("final_stereo_satisfied") or 0 for n in names)
    violated_baseline = sum(baseline[n].get("final_stereo_violations") or 0 for n in names)
    violated_new = sum(new[n].get("final_stereo_violations") or 0 for n in names)
    unevaluable_baseline = sum(baseline[n].get("final_stereo_unevaluable") or 0 for n in names)
    unevaluable_new = sum(new[n].get("final_stereo_unevaluable") or 0 for n in names)

    print(f"  pipeline success:      {n_success_baseline}/{len(names)} -> {n_success_new}/{len(names)}")
    print(f"  geometry sound:        {n_sound_baseline}/{len(names)} -> {n_sound_new}/{len(names)}")
    print(f"  stereo satisfied (sum): {satisfied_baseline} -> {satisfied_new}")
    print(f"  stereo violated (sum):  {violated_baseline} -> {violated_new}")
    print(f"  stereo unevaluable (sum): {unevaluable_baseline} -> {unevaluable_new}")

    transitions = defaultdict(list)
    for n in names:
        for label in classify_transition(baseline[n], new[n]):
            transitions[label].append(n)

    for label in (
        "success_to_failure",
        "failure_to_success",
        "sound_to_unsound",
        "unsound_to_sound",
        "stereo_newly_fixed",
        "stereo_newly_broken",
    ):
        mols = transitions.get(label, [])
        print(f"  {label}: {len(mols)}" + (f"  {mols}" if mols else ""))
    print()

    regression = bool(transitions.get("sound_to_unsound") or transitions.get("stereo_newly_broken")
                       or transitions.get("success_to_failure"))
    return regression


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True)
    ap.add_argument("--historical", required=True)
    args = ap.parse_args()

    new_rows = load_rows(args.new)
    hist_rows = load_rows(args.historical)

    mismatches = report_regression_control(new_rows, hist_rows)

    print("=" * 78)
    print("2. FEATURE EFFECT — chematic_pipeline_v2_mmff94_strict_enforce_chirality "
          "vs chematic_pipeline_v2_mmff94_strict (same run)")
    print("=" * 78)

    regression_all = report_feature_effect(new_rows, "ALL 265", lambda r: True)
    regression_declared = report_feature_effect(
        new_rows, "declared-stereo subset (stereo_before_declared > 0)",
        lambda r: (r.get("stereo_before_declared") or 0) > 0,
    )
    regression_ez = report_feature_effect(
        new_rows, "declared-E/Z subset (SMILES contains / or \\)",
        lambda r: is_declared_ez(r.get("smiles", "")),
    )

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    if mismatches:
        print("REGRESSION CONTROL: FAILED — pre-existing arms changed. STOP, do not declare GO.")
    else:
        print("REGRESSION CONTROL: PASS — pre-existing arms unchanged.")
    if regression_all or regression_declared or regression_ez:
        print("FEATURE EFFECT: REGRESSION DETECTED (newly-broken / sound-to-unsound / "
              "success-to-failure present) — STOP, do not declare GO.")
    else:
        print("FEATURE EFFECT: no newly-broken molecules, no soundness regression, "
              "no success-to-failure transitions.")

    if mismatches or regression_all or regression_declared or regression_ez:
        sys.exit(1)


if __name__ == "__main__":
    main()
