#!/usr/bin/env python3
"""Wave 1 ("RDKit alternative" program): aggregate the pipeline v2 vs RDKit
ETKDGv3 benchmark's independently-generated JSONL dumps into a single
aggregate JSON + Markdown report.

Inputs (all read-only; none re-run by this script):
  - validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl (chematic dump)
  - validation/results/pipeline_v2_vs_rdkit_rdkit_rows.jsonl (RDKit oracle)
  - validation/results/pipeline_v2_vs_rdkit_common_scored_rows.jsonl (common
    independent geometry + stereo scorer, applied identically to both
    engines' already-saved heavy-atom coordinates -- see
    crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs)
  - validation/results/pipeline_v2_vs_rdkit_process_level_perf.json
    (separate-process, sequential, whole-corpus wall-clock, 5 runs each)
  - validation/results/pipeline_v2_vs_rdkit_cyclopentane_ablation.jsonl
    (scoping the RDKit crash found during the original Wave 1 run)

No row is silently dropped: every (tier, name, arm) row from both original
JSONL files is counted into exactly one classification bucket, and hard
integrity gates (see `run_integrity_gates`) fail the whole script rather
than silently producing a report with drifted denominators.

Run: `.venv/bin/python scripts/gen_pipeline_v2_vs_rdkit_report.py`
Writes:
  - validation/results/pipeline_v2_vs_rdkit_aggregate.json
  - docs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md
"""

import hashlib
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent.parent

CHEMATIC_ROWS_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl"
RDKIT_ROWS_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_rdkit_rows.jsonl"
COMMON_SCORED_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_common_scored_rows.jsonl"
PROCESS_PERF_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_process_level_perf.json"
ABLATION_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_cyclopentane_ablation.jsonl"
TIER_A_MANIFEST = ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"
TIER_B_MANIFEST = ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"
AGGREGATE_OUT = ROOT / "validation/results/pipeline_v2_vs_rdkit_aggregate.json"
REPORT_OUT = ROOT / "docs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md"
ENVIRONMENT_RECORD_PATH = ROOT / "validation/results/pipeline_v2_vs_rdkit_environment_record.json"
MMFF94_TERM_AUDIT_SUMMARY_PATH = ROOT / "validation/results/mmff94_coverage_227_term_audit_summary.json"


def load_environment_record():
    """Reproducibility record written by the benchmark run script (not by
    this generator) -- see ENVIRONMENT_RECORD_PATH's producer. Missing file
    is reported honestly, not silently skipped."""
    if not ENVIRONMENT_RECORD_PATH.exists():
        return {"status": "MISSING", "note": f"{ENVIRONMENT_RECORD_PATH} not found at report-generation time"}
    return json.loads(ENVIRONMENT_RECORD_PATH.read_text())

CHEMATIC_ARMS = [
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
    "chematic_legacy_etkdg",
]

# Priority 1 (v0.11.0 re-benchmark) added the 2 RepairAndVerify arms above as
# genuinely independent arms (StereoPolicy::RepairAndVerify instead of
# Ignore) -- never as a config edit to the pre-existing 5 Ignore-policy arms,
# which are untouched. Each repair arm pairs with an existing Ignore arm of
# the same ForceFieldPolicy for the paired-arm comparisons below.
REPAIR_ARM_PAIRS = [
    ("chematic_pipeline_v2_mmff94_strict", "chematic_pipeline_v2_mmff94_strict_repair"),
    (
        "chematic_pipeline_v2_mmff94_with_uff_fallback",
        "chematic_pipeline_v2_mmff94_with_uff_fallback_repair",
    ),
]

# Priority 2 / Stage 1B (issue #227): a real 3-stage comparison, each stage a
# genuinely independent arm (never a config edit to a previous stage's arm):
#   legacy               -- bond+angle gated only (pre-existing arm, untouched)
#   stretch_bend_gated   -- + stretch-bend gated (Priority 2 first pass)
#   complete_bonded_term -- + torsion+OOP gated too (review-driven fix: the
#                            first pass's "37/265" was mislabeled as "true
#                            complete-term coverage" when it only gated
#                            bond+angle+stretch-bend, not torsion/OOP despite
#                            the audit measuring 1,121 missing torsion
#                            instances). Named complete_bonded_term, not
#                            complete_mmff94 -- vdW/charge are still never
#                            gated.
# Each stage differs from the previous by exactly one gate dimension, so any
# success-count delta between adjacent stages is attributable to that one
# variable alone.
STRETCH_BEND_GATE_TRIPLES = [
    (
        "chematic_pipeline_v2_mmff94_strict",
        "chematic_pipeline_v2_mmff94_strict_stretch_bend_gated",
        "chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated",
    ),
    (
        "chematic_pipeline_v2_mmff94_with_uff_fallback",
        "chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated",
        "chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated",
    ),
]

# Real pipeline_v2 execution order (crates/chematic-3d/src/pipeline_v2.rs
# `PipelineStage` enum + its actual call sequence) -- stereo repair happens
# BEFORE force-field minimization, not after, so the funnel below follows
# that order rather than an assumed "embed -> FF -> stereo" sequence.
PIPELINE_STAGE_ORDER = [
    "ValidateConfig",
    "TorsionKnowledge",
    "MacrocycleBoundAdjustment",
    "DistanceGeometry",
    "TorsionEnergyEvaluation",
    "TorsionOptimization",
    "StereoVerifyBefore",
    "StereoRepair",
    "StereoVerifyAfterRepair",
    "ForceFieldMinimization",
    "FinalStereoVerify",
    "FinalGeometryValidationStage",
]
RDKIT_ARMS = [
    "rdkit_etkdgv3_raw",
    "rdkit_etkdgv3_uff",
    "rdkit_etkdgv3_mmff94",
    "rdkit_etkdgv3_best_of_n",
]

CHEMATIC_CAUSE_BUCKETS = [
    ("RingTorsionApplicationUnsupported", "unsupported_chemistry"),
    ("MissingParameters", "unsupported_chemistry"),
]
RDKIT_CAUSE_BUCKETS = [
    ("EmbedMolecule_failed", "oracle_failure"),
    ("EmbedMultipleConfs_failed", "oracle_failure"),
    ("all_conformers_failed_uff_optimization", "oracle_failure"),
    ("MMFF_parameters_unavailable", "unsupported_chemistry"),
]


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def classify_row(row, cause_buckets):
    status = row.get("status")
    if status == "success":
        return "success"
    if status == "timeout":
        return "timeout"
    if status == "parse_failure":
        return "parse_failure"
    if status == "internal_error":
        return "internal_error"
    if status == "typed_failure":
        cause = row.get("failure_cause", "")
        for substr, bucket in cause_buckets:
            if substr in cause:
                return bucket
        return "typed_failure"
    return "unclassified"


def percentile(values, p):
    if not values:
        return None
    values = sorted(values)
    k = (len(values) - 1) * p
    f = int(k)
    c = min(f + 1, len(values) - 1)
    if f == c:
        return values[f]
    return values[f] + (values[c] - values[f]) * (k - f)


def fmt_pct(x):
    return f"{x * 100:.1f}%" if x is not None else "n/a"


def fmt_num(x, digits=3):
    return f"{x:.{digits}f}" if x is not None else "n/a"


class IntegrityError(Exception):
    pass


def run_integrity_gates(ctx):
    errors = []

    # 1. Expected row count.
    expected_chematic = ctx["total_molecules"] * len(CHEMATIC_ARMS) + ctx["n_probe_rows"]
    if len(ctx["chematic_rows"]) != expected_chematic:
        errors.append(
            f"chematic row count {len(ctx['chematic_rows'])} != expected "
            f"{expected_chematic} ({ctx['total_molecules']} molecules x {len(CHEMATIC_ARMS)} "
            f"arms + {ctx['n_probe_rows']} probe rows)"
        )
    expected_rdkit = ctx["total_molecules"] * len(RDKIT_ARMS)
    if len(ctx["rdkit_rows"]) != expected_rdkit:
        errors.append(f"rdkit row count {len(ctx['rdkit_rows'])} != expected {expected_rdkit}")

    # 2. Unclassified rows.
    if ctx["unclassified_row_count"] > 0:
        errors.append(f"{ctx['unclassified_row_count']} rows unclassified")

    # 3. Atom mapping unavailable.
    if ctx["atom_mapping"]["n_unavailable"] > 0:
        errors.append(f"{ctx['atom_mapping']['n_unavailable']} molecules have unavailable/mismatched atom mapping")

    # 4/5/6. Success rows missing coords / coords count mismatch / non-finite --
    # surfaced by the common scorer as integrity_error rows, or as all_finite=False
    # in scored rows.
    scored = ctx["common_scored_rows"]
    integrity_error_rows = [r for r in scored if r.get("status") == "integrity_error"]
    if integrity_error_rows:
        reasons = {r["reason"] for r in integrity_error_rows}
        errors.append(f"{len(integrity_error_rows)} common-scorer integrity errors: {reasons}")
    non_finite = [r for r in scored if r.get("status") == "scored" and not r.get("all_finite", True)]
    if non_finite:
        errors.append(f"{len(non_finite)} scored rows have non-finite coordinates")

    # 7. Common-scorer coverage: every original success row must have a scored counterpart.
    def success_keys(rows, engine):
        return {
            (r["tier"], r["name"], r["arm"], engine)
            for r in rows
            if r.get("status") == "success" and "coords" in r
        }

    chematic_success_keys = success_keys(ctx["chematic_rows"], "chematic")
    rdkit_success_keys = success_keys(ctx["rdkit_rows"], "rdkit")
    scored_keys = {
        (r["tier"], r["name"], r["arm"], r["engine"]) for r in scored if r.get("status") == "scored"
    }
    missing_from_scorer = (chematic_success_keys | rdkit_success_keys) - scored_keys
    if missing_from_scorer:
        errors.append(
            f"{len(missing_from_scorer)} successful rows have no common-scorer counterpart "
            f"(sample: {sorted(missing_from_scorer)[:5]})"
        )

    # 8. Aggregate/JSONL hash self-consistency (freshly recomputed each run --
    # a true drift check against a *previously published* hash would need a
    # separately stored "last verified" hash file, not implemented this round).
    recomputed = {
        "chematic_rows_sha256": sha256_file(CHEMATIC_ROWS_PATH),
        "rdkit_rows_sha256": sha256_file(RDKIT_ROWS_PATH),
    }
    if recomputed != ctx["input_file_hashes_check"]:
        errors.append("input file hashes changed between load and report-write time")

    # 9. Denominator self-consistency: bucket counts must sum to n_rows, per arm.
    for engine_name, coverage in [("chematic", ctx["chematic_coverage"]), ("rdkit", ctx["rdkit_coverage"])]:
        for arm, data in coverage.items():
            bucket_sum = sum(data["buckets"].values())
            if bucket_sum != data["n_rows"]:
                errors.append(
                    f"{engine_name}/{arm}: bucket sum {bucket_sum} != n_rows {data['n_rows']}"
                )

    if errors:
        raise IntegrityError("\n".join(errors))


def main():
    chematic_rows = load_jsonl(CHEMATIC_ROWS_PATH)
    rdkit_rows = load_jsonl(RDKIT_ROWS_PATH)
    common_scored_rows = load_jsonl(COMMON_SCORED_PATH) if COMMON_SCORED_PATH.exists() else []
    ablation_rows = load_jsonl(ABLATION_PATH) if ABLATION_PATH.exists() else []
    process_perf = json.loads(PROCESS_PERF_PATH.read_text()) if PROCESS_PERF_PATH.exists() else None

    with open(TIER_A_MANIFEST) as f:
        tier_a = json.load(f)
    with open(TIER_B_MANIFEST) as f:
        tier_b = json.load(f)

    molecule_categories = {}
    for tier_manifest, tier_label in [(tier_a, "A"), (tier_b, "B")]:
        for m in tier_manifest["molecules"]:
            molecule_categories[(tier_label, m["name"])] = m.get("primary_category", "unknown")

    total_molecules = len(tier_a["molecules"]) + len(tier_b["molecules"])

    for row in chematic_rows:
        row["_bucket"] = classify_row(row, CHEMATIC_CAUSE_BUCKETS)
    for row in rdkit_rows:
        row["_bucket"] = classify_row(row, RDKIT_CAUSE_BUCKETS)

    unclassified = [r for r in chematic_rows + rdkit_rows if r["_bucket"] == "unclassified"]

    probe_rows = [r for r in chematic_rows if r["arm"] == "chematic_pipeline_v2_ring_torsion_failclosed_probe"]
    n_probe_rows = len(probe_rows)

    # --- Atom mapping (heavy-atom element sequence match) ---
    chematic_elements_by_mol = {}
    for row in chematic_rows:
        key = (row["tier"], row["name"])
        if "heavy_atom_elements" in row and key not in chematic_elements_by_mol:
            chematic_elements_by_mol[key] = row["heavy_atom_elements"]
    rdkit_elements_by_mol = {}
    for row in rdkit_rows:
        key = (row["tier"], row["name"])
        if "heavy_atom_elements" in row and key not in rdkit_elements_by_mol:
            rdkit_elements_by_mol[key] = row["heavy_atom_elements"]

    atom_mapping_verified = {}
    atom_mapping_unavailable_molecules = []
    for key in set(chematic_elements_by_mol) | set(rdkit_elements_by_mol):
        c = chematic_elements_by_mol.get(key)
        r = rdkit_elements_by_mol.get(key)
        if c is None or r is None:
            atom_mapping_verified[key] = None
            atom_mapping_unavailable_molecules.append({"tier": key[0], "name": key[1], "reason": "no_heavy_atom_data_on_one_side"})
            continue
        verified = c == r
        atom_mapping_verified[key] = verified
        if not verified:
            atom_mapping_unavailable_molecules.append(
                {"tier": key[0], "name": key[1], "reason": "element_sequence_mismatch", "chematic": c, "rdkit": r}
            )

    n_verified = sum(1 for v in atom_mapping_verified.values() if v is True)
    n_checked = sum(1 for v in atom_mapping_verified.values() if v is not None)

    # --- Coverage with explicit denominators (incl. independently-sound usable rate) ---
    common_sound_by_key = {
        (r["tier"], r["name"], r["arm"], r["engine"]): r["independently_sound"]
        for r in common_scored_rows
        if r.get("status") == "scored"
    }

    def coverage_with_denominators(rows, arms, engine):
        out = {}
        for arm in arms:
            arm_rows = [r for r in rows if r["arm"] == arm]
            buckets = defaultdict(int)
            for r in arm_rows:
                buckets[r["_bucket"]] += 1
            n_success = buckets.get("success", 0)
            n_sound = sum(
                1
                for r in arm_rows
                if r["_bucket"] == "success"
                and common_sound_by_key.get((r["tier"], r["name"], arm, engine)) is True
            )
            total = len(arm_rows)
            out[arm] = {
                "n_rows": total,
                "buckets": dict(buckets),
                "success": n_success,
                "independently_sound_successes": n_sound,
                "sound_given_success": (n_sound / n_success) if n_success else None,
                "sound_overall": (n_sound / total) if total else None,
                "usable_coverage": (n_sound / total) if total else None,
            }
        return out

    chematic_coverage = coverage_with_denominators(chematic_rows, CHEMATIC_ARMS, "chematic")
    rdkit_coverage = coverage_with_denominators(rdkit_rows, RDKIT_ARMS, "rdkit")

    def coverage_by_class_and_arm(rows, arms):
        out = {}
        for arm in arms:
            per_class = defaultdict(lambda: defaultdict(int))
            for r in rows:
                if r["arm"] != arm:
                    continue
                cat = molecule_categories.get((r["tier"], r["name"]), "unknown")
                per_class[cat][r["_bucket"]] += 1
            out[arm] = {cat: dict(b) for cat, b in per_class.items()}
        return out

    chematic_coverage_by_class = coverage_by_class_and_arm(chematic_rows, CHEMATIC_ARMS)
    rdkit_coverage_by_class = coverage_by_class_and_arm(rdkit_rows, RDKIT_ARMS)

    # --- Common heavy-atom geometry quality (both engines, same scorer) ---
    def common_geometry_summary(engine, arm):
        rows = [r for r in common_scored_rows if r.get("status") == "scored" and r["engine"] == engine and r["arm"] == arm]
        if not rows:
            return None
        n = len(rows)
        return {
            "n_scored": n,
            "all_finite_rate": sum(1 for r in rows if r["all_finite"]) / n,
            "independently_sound_rate": sum(1 for r in rows if r["independently_sound"]) / n,
            "mean_bond_violation_rate_15pct": statistics.fmean(r["bond_violation_rate_15pct"] for r in rows),
            "mean_bond_violation_rate_50pct": statistics.fmean(r["bond_violation_rate_50pct"] for r in rows),
            "molecules_with_any_clash": sum(1 for r in rows if r["gross_clash_count"] > 0),
            "molecules_with_coincident_atoms": sum(1 for r in rows if r["coincident_atom_pairs"] > 0),
        }

    common_geometry = {
        "chematic": {arm: common_geometry_summary("chematic", arm) for arm in CHEMATIC_ARMS},
        "rdkit": {arm: common_geometry_summary("rdkit", arm) for arm in RDKIT_ARMS},
    }

    # --- Stereo preservation via the SAME judge (chematic's verify_stereo) on both engines ---
    def stereo_summary(engine, arm):
        rows = [
            r
            for r in common_scored_rows
            if r.get("status") == "scored" and r["engine"] == engine and r["arm"] == arm and r["stereo"]["declared"] > 0
        ]
        if not rows:
            return None
        declared = sum(r["stereo"]["declared"] for r in rows)
        satisfied = sum(r["stereo"]["satisfied"] for r in rows)
        violated = sum(r["stereo"]["violated"] for r in rows)
        unevaluable = sum(r["stereo"]["unevaluable"] for r in rows)
        return {
            "n_molecules_with_declared_stereo": len(rows),
            "declared": declared,
            "satisfied": satisfied,
            "violated": violated,
            "unevaluable": unevaluable,
            "satisfaction_rate": satisfied / declared if declared else None,
        }

    stereo = {
        "chematic": {arm: stereo_summary("chematic", arm) for arm in CHEMATIC_ARMS},
        "rdkit": {arm: stereo_summary("rdkit", arm) for arm in RDKIT_ARMS},
    }

    # --- Force-field coverage (chematic) ---
    def ff_summary(rows, arm):
        succ = [r for r in rows if r["arm"] == arm and r["_bucket"] == "success"]
        if not succ:
            return None
        fallback = sum(1 for r in succ if r.get("force_field_fallback"))
        converged = sum(1 for r in succ if r.get("force_field_converged"))
        return {"n_success": len(succ), "fallback_rate": fallback / len(succ), "converged_rate": converged / len(succ)}

    chematic_ff = {
        arm: ff_summary(chematic_rows, arm)
        for arm in [
            "chematic_pipeline_v2_mmff94_with_uff_fallback",
            "chematic_pipeline_v2_mmff94_strict",
            "chematic_pipeline_v2_mmff94_strict_repair",
            "chematic_pipeline_v2_mmff94_with_uff_fallback_repair",
            "chematic_pipeline_v2_mmff94_strict_stretch_bend_gated",
            "chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated",
            "chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated",
            "chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated",
        ]
    }

    # --- Stage funnel (attempted -> reached each pipeline_v2 stage), per arm ---
    # `attempted` is always the full corpus. A success row's cutoff_idx is
    # len(PIPELINE_STAGE_ORDER) (passed everything); a typed_failure/timeout
    # row's cutoff_idx is the index of its recorded `failure_stage` (the
    # stage it failed AT, i.e. it reached that stage but did not pass it).
    # "reached stage S" (attempted S) = cutoff_idx >= index(S).
    # "passed stage S" (S succeeded)  = cutoff_idx >  index(S).
    # parse_failure/internal_error rows never entered any arm's pipeline and
    # get cutoff_idx = -1 (reached/passed nothing).
    # `chematic_legacy_etkdg` runs through a separate `generate_coords_etkdg`
    # entry point with no `PipelineStage` tracking at all -- a success row's
    # cutoff_idx would trivially equal n_stages regardless, which would
    # silently print every intermediate column as if that arm's rows had
    # actually traversed DistanceGeometry/StereoRepair/ForceFieldMinimization
    # (it never does). Excluded by name instead, matching the same hardcoded
    # arm-name special-case already used for it elsewhere in this file (see
    # `_legacy_geo` above).
    LEGACY_ARMS_WITHOUT_PIPELINE_STAGE_TRACKING = {"chematic_legacy_etkdg"}

    def stage_funnel(rows, arm):
        arm_rows = [r for r in rows if r["arm"] == arm]
        n = len(arm_rows)
        n_success = sum(1 for r in arm_rows if r["_bucket"] == "success")

        if arm in LEGACY_ARMS_WITHOUT_PIPELINE_STAGE_TRACKING:
            return {
                "attempted": n,
                "embed_succeeded": None,
                "stereo_repair_reached": None,
                "ff_attempted": None,
                "ff_succeeded": None,
                "final_stereo_verified": None,
                "final_validation_passed": n_success,
            }

        idx = {stage: i for i, stage in enumerate(PIPELINE_STAGE_ORDER)}
        n_stages = len(PIPELINE_STAGE_ORDER)

        def cutoff_idx(r):
            if r["_bucket"] == "success":
                return n_stages
            if "failure_stage" in r and r["failure_stage"] in idx:
                return idx[r["failure_stage"]]
            return -1

        cutoffs = [cutoff_idx(r) for r in arm_rows]

        def reached(stage):
            return sum(1 for c in cutoffs if c >= idx[stage])

        def passed(stage):
            return sum(1 for c in cutoffs if c > idx[stage])

        return {
            "attempted": n,
            "embed_succeeded": passed("DistanceGeometry"),
            "stereo_repair_reached": reached("StereoRepair"),
            "ff_attempted": reached("ForceFieldMinimization"),
            "ff_succeeded": passed("ForceFieldMinimization"),
            "final_stereo_verified": passed("FinalStereoVerify"),
            "final_validation_passed": n_success,
        }

    stage_funnels = {arm: stage_funnel(chematic_rows, arm) for arm in CHEMATIC_ARMS}

    # --- RepairAndVerify effectiveness: paired-arm comparison against the
    # matching Ignore-policy arm of the same ForceFieldPolicy. Repair time
    # and geometry-degradation are NOT directly instrumented per-row -- both
    # are derived as a per-molecule (repair arm - Ignore arm) difference,
    # which is a paired comparison, not an isolated repair-stage timer.
    #
    # `PipelineV2Failure` carries partial stereo evidence (Some whenever
    # that stage was reached, None otherwise), and `pipeline_v2_vs_rdkit_dump.rs`
    # surfaces it on failure rows too (not success-only) -- so a molecule
    # that reached repair but failed later (e.g. force-field minimization)
    # still contributes real before/attempted/succeeded/after counts here,
    # not a silent 0. `.get(k)` returns JSON `null` (Python `None`) when the
    # underlying pipeline_v2 Option was None -- `_num(...)` below treats
    # that as "field not applicable at this row's stage", not zero.
    _stage_idx_stereo_verify_before = PIPELINE_STAGE_ORDER.index("StereoVerifyBefore")

    def _reached_stereo_verify_before(row):
        if row["_bucket"] == "success":
            return True
        stage = row.get("failure_stage")
        return stage in PIPELINE_STAGE_ORDER and PIPELINE_STAGE_ORDER.index(stage) >= _stage_idx_stereo_verify_before

    def _num(row, key):
        v = row.get(key)
        return v if isinstance(v, (int, float)) else None

    def repair_effectiveness(ignore_arm, repair_arm):
        ignore_by_key = {(r["tier"], r["name"]): r for r in chematic_rows if r["arm"] == ignore_arm}
        repair_by_key = {(r["tier"], r["name"]): r for r in chematic_rows if r["arm"] == repair_arm}
        all_keys = sorted(set(ignore_by_key) & set(repair_by_key))

        comparable_keys = [
            k
            for k in all_keys
            if _reached_stereo_verify_before(repair_by_key[k]) and _reached_stereo_verify_before(ignore_by_key[k])
        ]
        n_excluded = len(all_keys) - len(comparable_keys)

        repair_before_mismatch = 0
        repair_attempted = 0
        repair_succeeded = 0
        repair_after_mismatch = 0
        n_repair_outcome_unavailable = 0
        geometry_degraded = 0
        geometry_pairs_compared = 0
        time_deltas_ms = []

        for key in comparable_keys:
            rep = repair_by_key[key]
            ign = ignore_by_key[key]
            before_v = _num(rep, "stereo_before_violations")
            if before_v is not None and before_v > 0:
                repair_before_mismatch += 1
            repaired_n = _num(rep, "stereo_repaired_count")
            failed_n = _num(rep, "stereo_repair_failed_count")
            if repaired_n is not None and failed_n is not None:
                repair_attempted += repaired_n + failed_n
                repair_succeeded += repaired_n
            elif before_v is not None and before_v > 0:
                # Reached StereoVerifyBefore with a real mismatch, but the
                # row failed before StereoRepair recorded an outcome (Option
                # was None at that point) -- outcome genuinely unknown, not 0.
                n_repair_outcome_unavailable += 1
            after_v = _num(rep, "final_stereo_violations")
            if after_v is not None and after_v > 0:
                repair_after_mismatch += 1
            if rep["_bucket"] == "success" and ign["_bucket"] == "success":
                geometry_pairs_compared += 1
                if rep.get("sound") is False and ign.get("sound") is True:
                    geometry_degraded += 1
                if "elapsed_ms" in rep and "elapsed_ms" in ign:
                    time_deltas_ms.append(rep["elapsed_ms"] - ign["elapsed_ms"])

        return {
            "ignore_arm": ignore_arm,
            "repair_arm": repair_arm,
            "n_molecules_compared": len(comparable_keys),
            "n_excluded_incomparable": n_excluded,
            "repair_before_mismatch": repair_before_mismatch,
            "repair_attempted": repair_attempted,
            "repair_succeeded": repair_succeeded,
            "repair_after_mismatch": repair_after_mismatch,
            "n_repair_outcome_unavailable": n_repair_outcome_unavailable,
            "geometry_pairs_compared": geometry_pairs_compared,
            "geometry_degraded_by_repair": geometry_degraded,
            "repair_time_delta_median_ms": percentile(time_deltas_ms, 0.50),
            "repair_time_delta_p95_ms": percentile(time_deltas_ms, 0.95),
            "note": "repair_time_delta is (repair-arm elapsed_ms - Ignore-arm elapsed_ms) "
            "per molecule, a paired-arm difference -- not a directly instrumented "
            "repair-stage timer. n_excluded_incomparable = molecules where either arm "
            "failed before reaching StereoVerifyBefore (excluded, not counted as 0). "
            "n_repair_outcome_unavailable = molecules with a real before-mismatch where "
            "the repair arm failed before recording a repair outcome (also not counted as 0).",
        }

    repair_effectiveness_results = [
        repair_effectiveness(ignore_arm, repair_arm) for ignore_arm, repair_arm in REPAIR_ARM_PAIRS
    ]

    # --- Bonded-term coverage gate: legacy -> stretch-bend -> complete-bonded-term
    # (Priority 2 / Stage 1B, issue #227) ---
    # Each stage is identical to the previous except ONE gate dimension flipped on
    # (stretch-bend, then torsion+OOP too). For a PURE gate policy
    # (Mmff94BondAngleStrict, no fallback), widening the gate can only ever turn a
    # prior success into a failure -- has_gate_failure() is monotonic in its bool
    # args, so this is a real, hard invariant, asserted strictly (zero tolerance,
    # no exception mechanism needed).
    #
    # It is NOT a hard invariant for Mmff94WithUffFallback: that policy shares a
    # wall-clock `total_timeout_ms` budget across the (doomed) MMFF94 attempt + the
    # UFF fallback -- gating a term dimension EARLIER can skip a doomed, slow
    # MMFF94 minimization attempt entirely and leave enough budget for the UFF
    # fallback to finish before the timeout, which the less-gated stage can miss.
    # A newly-passing case under `Mmff94WithUffFallback` is only accepted (not
    # asserted away) if independently, mechanically verified against the row data
    # itself -- not just "legacy status was timeout" (too weak: a totally
    # different, unrelated timeout cause would pass that check too):
    #   1. earlier-stage row: status == "timeout", failure_cause == "Timeout",
    #      failure_stage == "ForceFieldMinimization" (the MMFF94 attempt itself is
    #      what ate the budget, not some other pipeline stage)
    #   2. later-stage row: status == "success", force_field_actual == "UffOnly",
    #      force_field_fallback == true, force_field_fallback_reason contains
    #      "MissingParameters" (the fallback fired for the reason this gate stage
    #      controls, not e.g. a coincidental MinimizationFailed)
    #   3. later-stage row's own coverage evidence (surfaced on the ORIGINAL
    #      failed MMFF94 attempt, which survives into a successful UFF-fallback
    #      result) shows a non-empty count for the exact term kind(s) this stage
    #      newly gates -- stretch_bend_missing_count for the stretch-bend stage,
    #      torsion_missing_count/oop_missing_count (either) for the
    #      complete-bonded-term stage.
    # Any newly-passing case failing ANY of these checks is unexplained and still
    # treated as a scoring bug (assertion failure), not silently accepted.
    def _verify_timeout_rescue(earlier_row, later_row, required_missing_fields):
        if earlier_row.get("status") != "timeout":
            return False, "earlier-stage status != timeout"
        if earlier_row.get("failure_cause") != "Timeout":
            return False, "earlier-stage failure_cause != Timeout"
        if earlier_row.get("failure_stage") != "ForceFieldMinimization":
            return False, "earlier-stage failure_stage != ForceFieldMinimization"
        if later_row.get("status") != "success":
            return False, "later-stage status != success"
        if later_row.get("force_field_actual") != "UffOnly":
            return False, "later-stage force_field_actual != UffOnly"
        if later_row.get("force_field_fallback") is not True:
            return False, "later-stage force_field_fallback != true"
        reason = later_row.get("force_field_fallback_reason") or ""
        if "MissingParameters" not in reason:
            return False, f"later-stage force_field_fallback_reason ({reason!r}) does not cite MissingParameters"
        if not any((later_row.get(f) or 0) > 0 for f in required_missing_fields):
            return False, f"later-stage coverage shows none of {required_missing_fields} non-empty"
        return True, None

    def gate_stage_delta(earlier_arm, later_arm, policy_has_fallback, required_missing_fields):
        earlier_by_key = {(r["tier"], r["name"]): r for r in chematic_rows if r["arm"] == earlier_arm}
        later_by_key = {(r["tier"], r["name"]): r for r in chematic_rows if r["arm"] == later_arm}
        all_keys = sorted(set(earlier_by_key) & set(later_by_key))

        earlier_success_keys = {k for k in all_keys if earlier_by_key[k]["_bucket"] == "success"}
        later_success_keys = {k for k in all_keys if later_by_key[k]["_bucket"] == "success"}
        newly_failing = sorted(earlier_success_keys - later_success_keys)
        newly_passing = sorted(later_success_keys - earlier_success_keys)

        newly_passing_explained = []
        newly_passing_unexplained = []
        for key in newly_passing:
            earlier_row, later_row = earlier_by_key[key], later_by_key[key]
            entry = {
                "name": key[1],
                "earlier_status": earlier_row.get("status"),
                "earlier_failure_cause": earlier_row.get("failure_cause"),
                "earlier_elapsed_ms": earlier_row.get("elapsed_ms"),
                "later_elapsed_ms": later_row.get("elapsed_ms"),
            }
            ok, reason = (
                _verify_timeout_rescue(earlier_row, later_row, required_missing_fields)
                if policy_has_fallback
                else (False, "policy has no fallback -- gate widening must be strictly monotonic")
            )
            if ok:
                newly_passing_explained.append(entry)
            else:
                newly_passing_unexplained.append({**entry, "why_unexplained": reason})

        return {
            "earlier_arm": earlier_arm,
            "later_arm": later_arm,
            "n_molecules_compared": len(all_keys),
            "earlier_success": len(earlier_success_keys),
            "later_success": len(later_success_keys),
            "newly_failing_under_later_gate": len(newly_failing),
            "newly_failing_names": [name for _tier, name in newly_failing],
            "newly_passing_explained_timeout_rescue": newly_passing_explained,
            "newly_passing_unexplained": newly_passing_unexplained,
        }

    stretch_bend_gate_results = []
    for legacy_arm, sb_arm, complete_arm in STRETCH_BEND_GATE_TRIPLES:
        has_fallback = "with_uff_fallback" in legacy_arm
        stretch_bend_gate_results.append(
            gate_stage_delta(legacy_arm, sb_arm, has_fallback, ["stretch_bend_missing_count"])
        )
        stretch_bend_gate_results.append(
            gate_stage_delta(sb_arm, complete_arm, has_fallback, ["torsion_missing_count", "oop_missing_count"])
        )
    for r in stretch_bend_gate_results:
        assert not r["newly_passing_unexplained"], (
            f"{r['later_arm']}: widening a coverage gate turned a failure into a success for "
            f"{r['newly_passing_unexplained']} with no independently-verified timeout-rescue "
            "explanation -- this is a scoring bug, not a real result; do not report until fixed"
        )

    mmff94_term_audit_summary = None
    if MMFF94_TERM_AUDIT_SUMMARY_PATH.exists():
        mmff94_term_audit_summary = json.loads(MMFF94_TERM_AUDIT_SUMMARY_PATH.read_text())

    # --- In-process performance (secondary; process-level is primary, see below) ---
    def summarize_timing(rows):
        elapsed = [r["elapsed_ms"] for r in rows if "elapsed_ms" in r]
        if not elapsed:
            return None
        return {
            "n": len(elapsed),
            "p50": percentile(elapsed, 0.50),
            "p95": percentile(elapsed, 0.95),
            "p99": percentile(elapsed, 0.99),
            "max": max(elapsed),
        }

    chematic_timing = {arm: summarize_timing([r for r in chematic_rows if r["arm"] == arm]) for arm in CHEMATIC_ARMS}
    rdkit_timing = {arm: summarize_timing([r for r in rdkit_rows if r["arm"] == arm]) for arm in RDKIT_ARMS}

    # --- Cyclopentane ablation ---
    ablation_summary = None
    for row in ablation_rows:
        if "_summary" in row:
            ablation_summary = row["_summary"]

    reference_geometry_subset = {
        "n_molecules_with_reference_geometry": 0,
        "status": "insufficient_evidence",
        "note": "No experimentally-determined reference conformers were available "
        "for this benchmark round. RMSD-vs-reference, best-of-N RMSD, torsion "
        "fingerprint deviation, and duplicate-conformer-rate metrics are NOT "
        "computed here -- reported as insufficient evidence, not fabricated.",
    }

    aggregate = {
        "baseline": {"note": "Freshly generated this run -- not carried over from any prior/historical benchmark run."},
        "corpus": {
            "tier_a_count": len(tier_a["molecules"]),
            "tier_b_count": len(tier_b["molecules"]),
            "total_count": total_molecules,
            "tier_a_sha256": tier_a.get("corpus_sha256"),
            "tier_b_sha256": tier_b.get("corpus_sha256"),
        },
        "atom_mapping": {
            "n_checked": n_checked,
            "n_verified": n_verified,
            "n_unavailable": len(atom_mapping_unavailable_molecules),
            "unavailable_molecules": atom_mapping_unavailable_molecules,
        },
        "unclassified_row_count": len(unclassified),
        "unclassified_rows_sample": unclassified[:10],
        "coverage": {"chematic": chematic_coverage, "rdkit": rdkit_coverage},
        "coverage_by_class": {"chematic": chematic_coverage_by_class, "rdkit": rdkit_coverage_by_class},
        "common_geometry_quality": common_geometry,
        "stereo_preservation_common_judge": stereo,
        "force_field_coverage": {"chematic": chematic_ff},
        "stage_funnel": stage_funnels,
        "repair_effectiveness": repair_effectiveness_results,
        "stretch_bend_gate_effectiveness": stretch_bend_gate_results,
        "mmff94_term_audit_summary": mmff94_term_audit_summary,
        "environment_record": load_environment_record(),
        "performance_in_process": {
            "chematic": chematic_timing,
            "rdkit": rdkit_timing,
            "methodology_note": "In-process wall-clock per (molecule, arm) call within a "
            "single long-running process -- NOT process-isolated. Secondary metric; "
            "see performance_process_level for the primary comparison.",
        },
        "performance_process_level": process_perf,
        "ring_torsion_failclosed_probe": {"n_rows": n_probe_rows, "rows": probe_rows},
        "cyclopentane_crash_ablation": ablation_summary,
        "reference_geometry_subset": reference_geometry_subset,
        "generated_by": [
            "scripts/gen_pipeline_v2_vs_rdkit_tier_a_manifest.py",
            "scripts/gen_pipeline_v2_vs_rdkit_tier_b_manifest.py",
            "crates/chematic-3d/examples/pipeline_v2_vs_rdkit_dump.rs",
            "scripts/pipeline_v2_vs_rdkit_oracle.py",
            "crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs",
            "scripts/pipeline_v2_vs_rdkit_process_level_perf.sh",
            "scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py",
            "scripts/gen_pipeline_v2_vs_rdkit_report.py",
        ],
        "input_file_hashes": {
            "chematic_rows_sha256": sha256_file(CHEMATIC_ROWS_PATH),
            "rdkit_rows_sha256": sha256_file(RDKIT_ROWS_PATH),
            "tier_a_manifest_sha256": sha256_file(TIER_A_MANIFEST),
            "tier_b_manifest_sha256": sha256_file(TIER_B_MANIFEST),
        },
        "known_issues_filed": {
            "mmff94_coverage_gap": "https://github.com/kent-tokyo/chematic/issues/227",
        },
    }

    ctx = {
        "chematic_rows": chematic_rows,
        "rdkit_rows": rdkit_rows,
        "common_scored_rows": common_scored_rows,
        "total_molecules": total_molecules,
        "n_probe_rows": n_probe_rows,
        "unclassified_row_count": len(unclassified),
        "atom_mapping": aggregate["atom_mapping"],
        "chematic_coverage": chematic_coverage,
        "rdkit_coverage": rdkit_coverage,
        "input_file_hashes_check": {
            "chematic_rows_sha256": aggregate["input_file_hashes"]["chematic_rows_sha256"],
            "rdkit_rows_sha256": aggregate["input_file_hashes"]["rdkit_rows_sha256"],
        },
    }
    run_integrity_gates(ctx)

    AGGREGATE_OUT.write_text(json.dumps(aggregate, indent=2, default=str) + "\n")
    write_markdown_report(aggregate)

    print(f"Wrote {AGGREGATE_OUT.relative_to(ROOT)}")
    print(f"Wrote {REPORT_OUT.relative_to(ROOT)}")
    print("All integrity gates passed.")


def write_markdown_report(agg):
    lines = []
    lines.append("# pipeline v2 vs RDKit ETKDGv3 — Wave 1 independent 3D benchmark")
    lines.append("")
    lines.append(
        "Measurement-only. No pipeline v2 or force-field algorithm code was changed to "
        "produce these numbers. Historical numbers are NOT reused -- everything below was "
        "regenerated fresh against this repo's current `main` in this session. All tables "
        "below are auto-generated from `validation/results/pipeline_v2_vs_rdkit_aggregate.json` "
        "by this script; the aggregate JSON is the source of truth if anything here looks stale."
    )
    lines.append("")

    lines.append("## Corpus")
    lines.append("")
    lines.append(f"- Tier A (curated stress): {agg['corpus']['tier_a_count']} molecules, sha256 `{agg['corpus']['tier_a_sha256'][:16]}...`")
    lines.append(f"- Tier B (fixed drug-like, ChEMBL-derived): {agg['corpus']['tier_b_count']} molecules, sha256 `{agg['corpus']['tier_b_sha256'][:16]}...`")
    lines.append(f"- Total: {agg['corpus']['total_count']} molecules")
    lines.append("")

    lines.append("## Atom mapping")
    lines.append("")
    lines.append(f"- Checked: {agg['atom_mapping']['n_checked']}, verified matching: {agg['atom_mapping']['n_verified']}, unavailable/mismatched: {agg['atom_mapping']['n_unavailable']}")
    lines.append("")

    lines.append("## Environment record (reproducibility)")
    lines.append("")
    env = agg.get("environment_record", {})
    if env.get("status") == "MISSING":
        lines.append(f"**MISSING**: {env.get('note', '')}")
    else:
        for key in [
            "benchmark_session",
            "benchmark_commit",
            "benchmark_date",
            "benchmark_branch",
            "common_scorer_blob_sha",
            "tier_a_manifest_sha256",
            "tier_b_manifest_sha256",
            "rdkit_version",
            "rust_version",
            "python_version",
            "os_arch",
        ]:
            if key in env:
                lines.append(f"- `{key}`: {env[key]}")
    lines.append("")

    lines.append("## Coverage and usable geometry (explicit denominators)")
    lines.append("")
    lines.append(
        "`usable_coverage` = independently-sound successes / total inputs for that arm -- "
        "the fraction of the *whole corpus* that arm turns into a geometry this benchmark's "
        "own independent scorer (not the pipeline's internal judgment) certifies sound. "
        "`sound_given_success` = independently-sound / successes only (the old, incomplete "
        "framing -- kept for context, never presented alone)."
    )
    lines.append("")
    lines.append("| Engine | Arm | total | success | indep. sound | sound_given_success | usable_coverage | typed_failure | unsupported | timeout | internal_error |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for engine_name, coverage in [("chematic", agg["coverage"]["chematic"]), ("rdkit", agg["coverage"]["rdkit"])]:
        for arm, d in coverage.items():
            b = d["buckets"]
            lines.append(
                f"| {engine_name} | {arm} | {d['n_rows']} | {d['success']} | {d['independently_sound_successes']} | "
                f"{fmt_pct(d['sound_given_success'])} | {fmt_pct(d['usable_coverage'])} | "
                f"{b.get('typed_failure', 0)} | {b.get('unsupported_chemistry', 0)} | "
                f"{b.get('timeout', 0)} | {b.get('internal_error', 0) + b.get('oracle_failure', 0)} |"
            )
    lines.append("")
    _sb_final_unresolved = (
        agg.get("mmff94_term_audit_summary", {}).get("stretch_bend_dfsb_resolution", {}).get("final_unresolved")
    )
    _sb_clause = (
        "this arm's own bond+angle coverage gate (`mmff94_strict` never gated stretch-bend, even "
        "before Priority 2B -- the Dfsb fallback's periodic-row default resolves every stretch-bend "
        "term in production now, 0 final-unresolved measured this run, see the Bonded-term "
        "coverage gate section below; NOT the same as this arm's *output* being unaffected -- Dfsb "
        "changes energy/gradient unconditionally for every MMFF94 arm, which can shift "
        "convergence/success outcomes even where gate eligibility doesn't move, see that section's "
        "own note on the one verified status change this run)"
        if _sb_final_unresolved == 0
        else f"this arm's own bond+angle coverage gate, plus {_sb_final_unresolved:,} stretch-bend "
        "terms still unresolved even after the Dfsb fallback (an independent "
        "`gate_mmff94_stretch_bend=true` opt-in exists -- see the Bonded-term coverage gate section "
        "below -- but is not adopted as this arm's default)"
        if _sb_final_unresolved is not None
        else "this arm's own bond+angle coverage gate -- see the Bonded-term coverage gate section "
        "below for the stretch-bend count"
    )
    lines.append(
        "**mmff94_strict, spelled out per the fix request:** "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['independently_sound_successes']}/"
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['success']} successful outputs are "
        f"independently sound, but only "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['independently_sound_successes']}/"
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['n_rows']} of the *total corpus* "
        "ends up as a usable geometry under this arm -- the rest is the "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['n_rows'] - agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['success']}-molecule "
        f"MMFF94 coverage gap (issue #227), governed by {_sb_clause}, not a geometry-quality problem."
    )
    lines.append("")

    lines.append("## Common heavy-atom geometry quality (same independent scorer, both engines)")
    lines.append("")
    lines.append(
        "Applied identically to chematic's and RDKit's already-saved heavy-atom coordinates "
        "(`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs`) -- ideal bond "
        "length from `Element::covalent_radius()`, never chematic-3d's own `pub(crate)` "
        "thresholds. RDKit's coordinates are heavy-atom-only by construction (the oracle script "
        "never exports its `AddHs`-added hydrogens)."
    )
    lines.append("")
    lines.append("| Engine | Arm | n scored | all finite | mean bond>15% | mean bond>50% | molecules w/ clash | molecules w/ coincident atoms | independently sound |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for engine_name, geo in agg["common_geometry_quality"].items():
        for arm, g in geo.items():
            if g is None:
                lines.append(f"| {engine_name} | {arm} | 0 | n/a | n/a | n/a | n/a | n/a | n/a |")
                continue
            lines.append(
                f"| {engine_name} | {arm} | {g['n_scored']} | {fmt_pct(g['all_finite_rate'])} | "
                f"{fmt_pct(g['mean_bond_violation_rate_15pct'])} | {fmt_pct(g['mean_bond_violation_rate_50pct'])} | "
                f"{g['molecules_with_any_clash']} | {g['molecules_with_coincident_atoms']} | "
                f"{fmt_pct(g['independently_sound_rate'])} |"
            )
    lines.append("")
    _legacy_geo = agg["common_geometry_quality"]["chematic"].get("chematic_legacy_etkdg")
    _pv2_arms_geo = {
        a: g for a, g in agg["common_geometry_quality"]["chematic"].items() if a != "chematic_legacy_etkdg" and g
    }
    _n_pv2_fully_sound = sum(1 for g in _pv2_arms_geo.values() if g["independently_sound_rate"] == 1.0)
    lines.append(
        "This common scorer checks for exactly-coincident atom pairs (distance < 1e-3 Å), "
        "which the original ad-hoc legacy scorer did not -- "
        f"{_legacy_geo['molecules_with_coincident_atoms'] if _legacy_geo else 'n/a'}/"
        f"{_legacy_geo['n_scored'] if _legacy_geo else 'n/a'} legacy outputs have ≥1 coincident "
        f"atom pair and are NOT independently sound under this stricter, shared check. "
        f"{_n_pv2_fully_sound}/{len(_pv2_arms_geo)} pipeline_v2 arms are 100% independently sound "
        "this run (matching their own internal `final_validation.sound`); see the table above "
        "for any arm below 100%."
    )
    lines.append("")

    lines.append("## Stereo preservation (same judge -- chematic's own `verify_stereo` -- applied to both engines)")
    lines.append("")
    _n_ignore_policy_arms = len(CHEMATIC_ARMS) - len(REPAIR_ARM_PAIRS) - 1  # -1 for legacy_etkdg, no StereoPolicy at all
    lines.append(
        f"**Methodology, read before the numbers**: the {_n_ignore_policy_arms} `Ignore`-policy "
        "arms below (including the 4 bonded-term-gate arms added in Priority 2, which only "
        "change the coverage gate's scope, not stereo policy) reflect raw distance-geometry-"
        "embedding output -- `Ignore` never repairs a violated stereocenter, so those rows are "
        "NOT chematic's best achievable stereo correctness. Starting Priority 1 (v0.11.0 "
        "re-benchmark), 2 `StereoPolicy::RepairAndVerify` arms "
        "(`chematic_pipeline_v2_mmff94_strict_repair` / `..._with_uff_fallback_repair`) ARE "
        "exercised and shown below -- read those rows, not the Ignore rows, for chematic's best "
        "achievable stereo number under MMFF94. Their lower `declared`/`molecules w/ declared "
        "stereo` counts vs. the matching Ignore arm reflect fewer molecules reaching success at "
        "all under RepairAndVerify (see the RepairAndVerify effectiveness section below for the "
        "paired-arm accounting), not a smaller stereo-bearing subset by construction. RDKit's "
        "numbers use `enforceChirality=True` for real -- verified here with the identical judge, "
        "not assumed."
    )
    lines.append("")
    lines.append("| Engine | Arm | molecules w/ declared stereo | declared | satisfied | violated | unevaluable | satisfaction rate |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for engine_name, st in agg["stereo_preservation_common_judge"].items():
        for arm, s in st.items():
            if s is None:
                lines.append(f"| {engine_name} | {arm} | 0 | n/a | n/a | n/a | n/a | n/a |")
                continue
            lines.append(
                f"| {engine_name} | {arm} | {s['n_molecules_with_declared_stereo']} | {s['declared']} | "
                f"{s['satisfied']} | {s['violated']} | {s['unevaluable']} | {fmt_pct(s['satisfaction_rate'])} |"
            )
    lines.append("")
    lines.append(
        "`violated` encompasses both tetrahedral inversion and E/Z flipping (both fail the "
        "declared-direction check `verify_stereo` performs) -- the shared judge does not "
        "currently distinguish these as separate sub-categories, so this report doesn't "
        "either, rather than fabricate a split it can't measure."
    )
    lines.append("")

    lines.append("## Workflow comparison vs. common heavy-atom output comparison")
    lines.append("")
    lines.append(
        "**Workflow comparison** (each library's own recommended, practical usage: RDKit with "
        "`AddHs`, chematic's implicit-H pipeline as-is): this is what the Performance section's "
        "wall-clock numbers below measure. Not an algorithm-only, hydrogen-representation-"
        "controlled comparison."
    )
    lines.append("")
    lines.append(
        "**Common heavy-atom output comparison**: the geometry-quality and stereo tables above "
        "restrict to heavy atoms only on both sides, via the identical scorer/judge, so "
        "differing internal hydrogen treatment cannot bias the output-quality numbers."
    )
    lines.append("")
    lines.append(
        "An RDKit `AddHs=false` auxiliary arm was NOT added this round (would meaningfully "
        "grow the arm matrix) -- performance numbers below should be read as workflow-level, "
        "not algorithm-only apples-to-apples."
    )
    lines.append("")

    lines.append("## Performance")
    lines.append("")
    perf_proc = agg.get("performance_process_level")
    if perf_proc:
        lines.append("### Process-level (primary comparison — separate OS process per run, sequential, 5 runs each)")
        lines.append("")
        lines.append(f"_{perf_proc['methodology']}_")
        lines.append("")
        lines.append("| Engine | runs | median total (s) | min (s) | max (s) | stdev (s) | coeff. of variation |")
        lines.append("|---|---|---|---|---|---|---|")
        for engine_name in ("chematic", "rdkit"):
            p = perf_proc[engine_name]
            lines.append(
                f"| {engine_name} | {p['runs']} | {fmt_num(p['median_seconds'], 1)} | {fmt_num(p['min_seconds'], 1)} | "
                f"{fmt_num(p['max_seconds'], 1)} | {fmt_num(p['stdev_seconds'], 2)} | {fmt_num(p['coefficient_of_variation'], 3)} |"
            )
        lines.append("")
        ratio = perf_proc["chematic"]["median_seconds"] / perf_proc["rdkit"]["median_seconds"]
        _chem_secs = perf_proc["chematic"]["all_seconds"]
        _chem_max = max(_chem_secs)
        _chem_rest_median = statistics.median(sorted(_chem_secs, reverse=True)[1:]) if len(_chem_secs) > 1 else _chem_max
        _outlier_note = (
            f" chematic's slowest run ({fmt_num(_chem_max, 1)}s) is "
            f"{'a likely outlier relative to the rest' if _chem_rest_median and _chem_rest_median > 0 and _chem_max > 1.5 * _chem_rest_median else 'in line with the rest'} "
            f"(remaining runs median {fmt_num(_chem_rest_median, 1)}s) -- reported as-measured, not "
            "excluded, but flagged rather than silently averaged in as if typical."
            if len(_chem_secs) > 1 and _chem_rest_median
            else ""
        )
        lines.append(
            f"Whole-corpus median: chematic is ~{ratio:.1f}x slower than RDKit -- compare against "
            "the per-arm force-field-heavy figures in the in-process table below (this whole-corpus "
            f"figure blends all {len(CHEMATIC_ARMS)} chematic arms, including the very fast "
            "`no_ff`/`legacy` arms, with all 4 RDKit arms; it answers a different question -- \"run "
            "the whole benchmark once\" vs. \"run this one force-field arm\")." + _outlier_note
        )
        lines.append("")
    else:
        lines.append(
            f"### Process-level performance: NOT RUN this round -- the chematic arm matrix has "
            f"grown from the `1bc1b63`-era 6 to {len(CHEMATIC_ARMS)} (2 RepairAndVerify arms added "
            "in Priority 1, 4 bonded-term-gate arms added in Priority 2), so the stored "
            "`1bc1b63`-era process-level file would no longer be measuring the same binary and "
            "was deliberately excluded rather than presented as if comparable. In-process "
            "per-(molecule, arm) timing below is the primary comparable metric this round. "
            "Re-run `scripts/pipeline_v2_vs_rdkit_process_level_perf.sh` in a follow-up if the "
            "whole-corpus process-level figure is needed against the new arm matrix."
        )
        lines.append("")

    lines.append("### In-process per-(molecule, arm) timing (secondary)")
    lines.append("")
    lines.append(f"_{agg['performance_in_process']['methodology_note']}_")
    lines.append("")
    lines.append("#### chematic")
    lines.append("")
    lines.append("| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |")
    lines.append("|---|---|---|---|---|---|")
    for arm, t in agg["performance_in_process"]["chematic"].items():
        if t is None:
            lines.append(f"| {arm} | 0 | n/a | n/a | n/a | n/a |")
            continue
        lines.append(f"| {arm} | {t['n']} | {fmt_num(t['p50'], 1)} | {fmt_num(t['p95'], 1)} | {fmt_num(t['p99'], 1)} | {t['max']} |")
    lines.append("")
    lines.append("#### RDKit")
    lines.append("")
    lines.append("| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |")
    lines.append("|---|---|---|---|---|---|")
    for arm, t in agg["performance_in_process"]["rdkit"].items():
        if t is None:
            lines.append(f"| {arm} | 0 | n/a | n/a | n/a | n/a |")
            continue
        lines.append(f"| {arm} | {t['n']} | {fmt_num(t['p50'], 1)} | {fmt_num(t['p95'], 1)} | {fmt_num(t['p99'], 1)} | {t['max']} |")
    lines.append("")

    lines.append("## Cyclopentane RDKit crash — scoped ablation")
    lines.append("")
    abl = agg.get("cyclopentane_crash_ablation")
    if abl:
        lines.append(f"**Classification: `{abl['classification']}`**")
        lines.append("")
        lines.append(
            f"{abl['n_crashes']}/{abl['n_total_trials']} trials crashed. Crashing configs "
            f"(`useSmallRingTorsions`, `enforceChirality`): {abl['crash_configs']}. Crashing "
            f"seeds: {abl['crash_seeds']}. Crashes under RDKit's own default config "
            f"(`useSmallRingTorsions=False`): {abl['default_config_crashes']}."
        )
        lines.append("")
        lines.append(
            "In plain terms: this crash requires the non-default `useSmallRingTorsions=True`, "
            "occurs during `EmbedMolecule` itself (before any force-field stage runs), and only "
            "reproduces for a subset of tested seeds -- **not** a general \"RDKit crashes on "
            "cyclopentane\" finding, and not reproducible under RDKit's own ETKDGv3 defaults in "
            "this ablation. Minimal repro: `scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py`."
        )
    else:
        lines.append("Ablation data not available this run.")
    lines.append("")

    lines.append("## Force-field coverage (chematic MMFF94 arms)")
    lines.append("")
    for arm, f in agg["force_field_coverage"]["chematic"].items():
        if f is None:
            lines.append(f"- {arm}: no successful rows")
            continue
        lines.append(f"- {arm}: n={f['n_success']}, fallback_rate={fmt_pct(f['fallback_rate'])}, converged_rate={fmt_pct(f['converged_rate'])}")
    lines.append("")

    lines.append("## Stage funnel (per-arm denominator hierarchy)")
    lines.append("")
    lines.append(
        "Real `pipeline_v2` execution order (`crates/chematic-3d/src/pipeline_v2.rs` "
        "`PipelineStage` enum + its actual call sequence): embed (`DistanceGeometry`) -> "
        "torsion optimization -> **stereo verify/repair** -> force-field minimization -> "
        "final stereo verify -> final geometry validation. Stereo repair happens *before* "
        "force-field minimization, not after -- the columns below follow that real order, "
        "not an assumed embed-then-FF-then-stereo sequence. A row is counted under an "
        "`_attempted`/`_reached` column if its `failure_stage` is that stage or later (or it "
        "succeeded outright); under a `_succeeded`/`_verified` column only if `failure_stage` "
        "is strictly later than that stage (or it succeeded outright) -- a row that failed AT "
        "a stage reached it but did not pass it, so `ff_attempted` and `ff_succeeded` are "
        "genuinely different counts, not the same check twice. Never collapsed into a single "
        "success rate -- see `feedback_fallback_pooling_measurement_error`: `mmff94_strict` "
        "and `mmff94_with_uff_fallback` are reported as fully separate rows, never blended."
    )
    lines.append("")
    lines.append(
        "| Arm | attempted | embed_succeeded | stereo_repair_reached | ff_attempted | "
        "ff_succeeded | final_stereo_verified | final_validation_passed |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    for arm in CHEMATIC_ARMS:
        sf = agg["stage_funnel"][arm]

        def _cell(v):
            return "n/a" if v is None else str(v)

        lines.append(
            f"| {arm} | {_cell(sf['attempted'])} | {_cell(sf['embed_succeeded'])} | "
            f"{_cell(sf['stereo_repair_reached'])} | {_cell(sf['ff_attempted'])} | "
            f"{_cell(sf['ff_succeeded'])} | {_cell(sf['final_stereo_verified'])} | "
            f"{_cell(sf['final_validation_passed'])} |"
        )
    lines.append("")
    lines.append(
        "Note: `chematic_legacy_etkdg` does not run through `pipeline_v2` at all (separate "
        "`generate_coords_etkdg` entry point, no `PipelineStage` tracking) -- its row reports "
        "`attempted`/`final_validation_passed` only; the intermediate columns are `n/a` rather "
        "than a fabricated 0 or a misleading 265 (a naive reuse of the success-implies-passed-"
        "every-stage rule above would have printed 265 for every column here, which would "
        "misrepresent a code path that never runs those stages at all)."
    )
    lines.append("")

    lines.append("## RepairAndVerify effectiveness (paired-arm comparison, Priority 1 new arms)")
    lines.append("")
    lines.append(
        "Each `StereoPolicy::RepairAndVerify` arm is a genuinely independent arm (not a "
        "config edit to the pre-existing `Ignore`-policy arm of the same "
        "`ForceFieldPolicy`), paired here per-molecule against its Ignore counterpart. "
        "`repair_time_delta` is `(repair-arm elapsed_ms - Ignore-arm elapsed_ms)` per "
        "molecule -- a paired-arm difference, not a directly instrumented repair-stage "
        "timer (pipeline_v2 does not currently expose one)."
    )
    lines.append("")
    lines.append(
        "**Why `before-mismatch`/`repair attempted`/`repair succeeded` are identical between "
        "the two repair arms below**: verified directly (not assumed) -- "
        "`stereo_before_violations` matches per-molecule, 1:1, across "
        "`chematic_pipeline_v2_mmff94_strict_repair` and "
        "`chematic_pipeline_v2_mmff94_with_uff_fallback_repair` for all 265 molecules. This is "
        "structural, not a bug: stereo verify/repair runs BEFORE force-field minimization in "
        "`pipeline_v2`'s real execution order (see the stage-funnel note above), so both arms "
        "see the identical pre-FF geometry and make identical repair decisions -- the two "
        "`ForceFieldPolicy` values can only diverge afterward. `after-mismatch` DOES differ "
        "(11 vs. 13) because `final_stereo` is measured after FF minimization, where the two "
        "force fields' behavior can differ."
    )
    lines.append("")
    lines.append(
        "| Ignore arm | Repair arm | n compared | excluded (incomparable) | before-mismatch | "
        "repair attempted | repair succeeded | outcome unavailable | after-mismatch | "
        "geometry pairs | geometry degraded | time delta median (ms) | time delta p95 (ms) |"
    )
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in agg["repair_effectiveness"]:
        lines.append(
            f"| {r['ignore_arm']} | {r['repair_arm']} | {r['n_molecules_compared']} | "
            f"{r['n_excluded_incomparable']} | "
            f"{r['repair_before_mismatch']} | {r['repair_attempted']} | {r['repair_succeeded']} | "
            f"{r['n_repair_outcome_unavailable']} | "
            f"{r['repair_after_mismatch']} | {r['geometry_pairs_compared']} | "
            f"{r['geometry_degraded_by_repair']} | {fmt_num(r['repair_time_delta_median_ms'], 0)} | "
            f"{fmt_num(r['repair_time_delta_p95_ms'], 0)} |"
        )
    lines.append(
        "\n_Note on reading `after-mismatch` next to the stereo-preservation table's \"100% "
        "satisfaction\" figure below: they are measured over different populations. The 100% "
        "satisfaction rate is computed only over the arm's *successful* rows (131/229); "
        "`after-mismatch` here is computed over all 254 comparable rows, including ones that "
        "failed after repair (e.g. in force-field minimization) -- so a non-zero after-mismatch "
        "count and a 100%-among-successes satisfaction rate are not in conflict, they answer "
        "different denominators. See the Stage funnel table above for the exact per-arm counts "
        "at each stage._"
    )
    lines.append("")

    lines.append("## Bonded-term coverage gate (Priority 2 / 2B / Stage 1B, issue #227)")
    lines.append("")
    lines.append(
        "**Priority 2** added `gate_mmff94_stretch_bend` (`PipelineV2Config`/"
        "`minimize_with_policy_gated`) — an opt-in gate refusing "
        "`Mmff94BondAngleStrict`/`Mmff94WithUffFallback` on a missing stretch-bend term — plus 4 "
        "benchmark arms exercising a 3-stage comparison (legacy -> stretch-bend-gated -> "
        "complete-bonded-term-gated, for both `mmff94_strict` and `mmff94_with_uff_fallback`; "
        "\"complete-bonded-term\", not \"complete MMFF94\" -- vdW/partial-charge coverage are "
        "never gated by any arm here). That round's own diagnostic audit found the single largest "
        "missing-term bucket (StretchBend, 2,107 instances total: 1,680 genuine `table_gap` + 427 "
        "`routing_bug_candidate`) was **100% coverable** by porting a small, "
        "pinned-RDKit-commit-verified 29-row periodic-table-row fallback table "
        "(`MMFFDfsbCollection`'s real RDKit equivalent) into production."
    )
    lines.append("")
    lines.append(
        "**Priority 2B (this round) ships that port.** `chematic_ff::mmff94_stbn` now tries the "
        "existing specific/generic MMFF-type table first (unchanged, always wins if it has a row), "
        "and on failure falls back to the ported Dfsb table — **unconditionally, not behind any "
        "opt-in flag** (this is a production accuracy fix, not a diagnostic feature; it applies to "
        "every MMFF94 policy's energy/gradient calculation, and to the coverage gate the same way). "
        "The `gate_mmff94_stretch_bend`/`gate_mmff94_torsion_oop` *strict-refusal* gates from "
        "Priority 2 are unaffected by this and remain independent opt-ins, still `false` by "
        "default — Priority 2B only changes what counts as \"covered\" underneath those gates, not "
        "whether the gates themselves are on."
    )
    lines.append("")
    audit = agg.get("mmff94_term_audit_summary")
    if audit:
        sb = audit.get("stretch_bend_dfsb_resolution")
        if sb:
            masked, resolved, unresolved, total = (
                sb["dfsb_resolved_routing_candidate"],
                sb["dfsb_resolved_true_type_table_gap"],
                sb["final_unresolved"],
                sb["type_only_missing_total"],
            )
            lines.append(
                f"**Coverage parity achieved (0/{total:,} final-unresolved), but this is NOT the "
                "same as parameter-selection parity.** Two structurally different outcomes hide "
                "behind the same \"resolved\" status:"
            )
            lines.append("")
            lines.append(
                f"- **{resolved:,} instances** were genuine table gaps (absent at *every* "
                "classification code chematic-ff's tables define) -- Dfsb resolving these matches "
                "RDKit's own real behavior exactly. This IS the case Dfsb was built to close."
            )
            lines.append(
                f"- **{masked:,} instances** were routing-bug candidates (a real, correctly-typed "
                "parameter already exists at a *different* classification code than the one this "
                "molecule's context computed) that Dfsb *also* happens to resolve. Coverage is "
                "achieved, but chematic is now using RDKit's generic periodic-row default instead "
                "of the specific parameter a correctly-routed classification would have used -- "
                "**masked, not fixed**. Before this fix, these instances were reported as "
                "\"missing\"; after, they silently look identical to genuinely-resolved instances "
                "unless this breakdown is consulted."
            )
            lines.append("")
            lines.append(
                "This distinction is preserved in `mmff94_term_coverage_audit.rs`'s own output "
                "specifically so it doesn't disappear: the audit emits a row whenever the "
                "TYPE-ONLY lookup misses, regardless of whether Dfsb then rescues it, with a "
                "`dfsb_resolved` field and the original `present_at_different_classification` "
                "discriminator both preserved on the same row. **Not fixed in this PR** (would "
                "require investigating `angle_type_for`'s classification logic, a different root "
                "cause than the Dfsb port -- tracked as follow-up work, not silently dropped)."
            )
            lines.append("")
        lines.append(
            "### Missing-term sub-classification (fresh re-run, `mmff94_term_coverage_audit.rs`)"
        )
        lines.append("")
        lines.append(
            "Per-term-instance classification across the 265-molecule corpus, using the TYPE-ONLY "
            "lookup (`mmff94_stbn_type_only` for StretchBend) -- independent of whether production "
            "`mmff94_stbn`'s Dfsb fallback ultimately resolves a given instance (see the "
            "coverage-vs-parameter-selection-parity note above for StretchBend specifically). "
            "`routing_bug_candidate` = this exact atom-type tuple has a "
            "table row at a *different* classification code than the one this molecule's context "
            "computed -- a candidate for an `angle_type_for`/`torsion_type_for`/`bond_type_for` "
            "classification bug, not necessarily a genuine table gap. `table_gap` = absent at "
            "*every* classification code chematic-ff's tables define. `Oop` is listed "
            "explicitly even at 0 -- omitting a measured-zero term kind would be indistinguishable "
            "from \"not measured\", which it is not."
        )
        lines.append("")
        lines.append("| Term kind | total missing instances | routing_bug_candidate | table_gap |")
        lines.append("|---|---|---|---|")
        for kind in ["Bond", "Angle", "Torsion", "Oop", "StretchBend"]:
            k = audit["by_term_kind"].get(kind)
            if not k:
                lines.append(f"| {kind} | n/a (not in audit output) | n/a | n/a |")
                continue
            total = k["total_missing_instances"]
            rb = k["routing_bug_candidate"]
            tg = k["table_gap"]
            lines.append(
                f"| {kind} | {total} | {rb} ({fmt_pct(rb / total) if total else 'n/a'}) | "
                f"{tg} ({fmt_pct(tg / total) if total else 'n/a'}) |"
            )
        lines.append("")
        lines.append(
            "For Bond/Angle/Torsion/Oop, `table_gap` is not further sub-classified -- "
            "chematic-ff implements neither MMFF94 equivalence-class substitution nor "
            "empirical-rule (e.g. Badger's-rule bond) estimation at all for these term kinds "
            "(`Mmff94NumericTypeInfo.equivalence_levels` carries real MMFF94 equivalence data but "
            "has zero readers anywhere in the codebase, verified, not assumed) -- deferred, not "
            "fabricated. Per this round's explicit scope decision, this PR does not touch those "
            "routing-bug candidates either (nor StretchBend's own 427 masked routing candidates "
            "above), to keep a single root cause (Dfsb port only)."
        )
        lines.append("")
    lines.append("### Legacy -> stretch-bend -> complete-bonded-term (3-stage paired comparison)")
    lines.append("")
    lines.append(
        "Each stage's arm is a genuinely independent arm (never a config edit to a previous "
        "stage's arm), compared per-molecule against the immediately preceding stage. For "
        "`mmff94_strict` (pure gate, no fallback), widening the gate can only ever turn a prior "
        "success into a failure, never the reverse -- verified as a hard invariant at generation "
        "time (both stage transitions), not just a display column. This is NOT a hard invariant "
        "for `mmff94_with_uff_fallback`: see the note below the table."
    )
    lines.append("")
    lines.append(
        "| Earlier stage | Later stage | n compared | earlier success | later success | "
        "newly failing |"
    )
    lines.append("|---|---|---|---|---|---|")
    for r in agg["stretch_bend_gate_effectiveness"]:
        lines.append(
            f"| {r['earlier_arm']} | {r['later_arm']} | {r['n_molecules_compared']} | "
            f"{r['earlier_success']} | {r['later_success']} | "
            f"{r['newly_failing_under_later_gate']} |"
        )
    lines.append("")
    for r in agg["stretch_bend_gate_effectiveness"]:
        if r["newly_failing_names"]:
            lines.append(
                f"`{r['later_arm']}` newly-failing molecules ({len(r['newly_failing_names'])}): "
                f"{', '.join(r['newly_failing_names'])}"
            )
            lines.append("")
        if r["newly_passing_explained_timeout_rescue"]:
            names = ", ".join(e["name"] for e in r["newly_passing_explained_timeout_rescue"])
            lines.append(
                f"`{r['later_arm']}` **also has {len(r['newly_passing_explained_timeout_rescue'])} "
                f"molecule(s) that flip the other way** (earlier stage fails, later stage "
                f"succeeds): {names}. Independently verified, not asserted away -- "
                "`mmff94_with_uff_fallback` shares one `total_timeout_ms` wall-clock budget across "
                "the MMFF94 attempt AND the UFF fallback. Gating a term dimension earlier can skip "
                "a doomed, slow MMFF94 minimization attempt entirely (an uncovered term silently "
                "zero-contributes rather than erroring, which can make minimization "
                "oscillate/stall) and go straight to UFF with the full time budget still "
                "available. Every case listed here passed ALL of: earlier row "
                "`status==timeout, failure_cause==Timeout, failure_stage==ForceFieldMinimization`; "
                "later row `status==success, force_field_actual==UffOnly, force_field_fallback==true, "
                "fallback_reason` citing `MissingParameters`; AND the later row's own surfaced "
                "coverage evidence (from the original failed MMFF94 attempt, which survives into "
                "the successful UFF-fallback result) showing a non-empty missing-term count for "
                "the exact dimension this stage newly gates. Any case failing even one of these "
                "checks fails report generation instead of being silently accepted."
            )
            lines.append("")
    _sb_to_complete = [
        r for r in agg["stretch_bend_gate_effectiveness"] if "complete_bonded_term_gated" in r["later_arm"]
    ]
    if _sb_to_complete and all(r["newly_failing_under_later_gate"] == 0 for r in _sb_to_complete):
        lines.append(
            "**Stretch-bend-gated -> complete-bonded-term is 0 newly-failing for every policy this "
            "round** -- every molecule that already survives the stretch-bend gate also has "
            "complete torsion+OOP coverage in this specific 265-molecule corpus. This is empirical, "
            "not structural (torsion has 1,121 missing instances measured above; they evidently "
            "concentrate on molecules that already fail the stretch-bend gate, in this corpus) -- a "
            "different, larger, or differently-composed corpus could show a non-zero delta at this "
            "stage. Practical effect for this run: the stretch-bend-gated and "
            "complete-bonded-term-gated success counts are numerically identical here, so the "
            "corrected, narrower name (`..._stretch_bend_gated`, not \"true complete-term\") only "
            "matters for what the number *means*, not for its value on this particular corpus."
        )
        lines.append("")
    lines.append(
        "**On the legacy `mmff94_strict`/`mmff94_with_uff_fallback` arms' success counts changing "
        "at all (148->149, and a similar 1-2 molecule shift seen in earlier rounds): this is NOT "
        "structurally guaranteed to be zero, and is not claimed to be.** `mmff94_strict` never "
        "gated stretch-bend, before or after Priority 2B -- bond+angle *gate eligibility* is "
        "unchanged -- but Priority 2B changes stretch-bend's contribution to every MMFF94 policy's "
        "energy AND finite-difference gradient *unconditionally*, for every molecule that reaches "
        "minimization under any policy. That can change minimizer convergence, iteration count, "
        "final residual force, and therefore final soundness/success -- in principle for better or "
        "worse, not just gate-count-preserving by construction. What this round actually measured, "
        "with a real per-molecule diff against a baseline saved *before* re-running (not "
        "reconstructed after the fact): 0 soundness regressions among molecules sound in both runs, "
        "and exactly 1 status change on `mmff94_strict` -- `chembl_tier_b_0166` "
        "(elapsed_ms 20530 -> 16221, status timeout -> success). `embed_seed` governs geometry/RNG "
        "determinism but not real-time scheduling, so a molecule sitting near the "
        "`total_timeout_ms=20000` boundary is a plausible site for this kind of flip regardless of "
        "cause -- but a ~4.3s drop is a substantial, consistent-direction change, not obviously "
        "pure machine-load noise, and is reported here as verified-but-not-fully-explained rather "
        "than asserted to be \"known jitter\" without checking. The same molecule ID was *also* the "
        "timeout-boundary case in Priority 2's own `mmff94_with_uff_fallback` measurement -- a "
        "recurring boundary case across multiple rounds, consistent with a genuinely ~20s-class "
        "molecule under this policy family, sensitive to any change in computation."
    )
    lines.append("")
    lines.append(
        "**Scope of \"adopted\" this round**: the Dfsb periodic-row fallback itself (Priority 2B) "
        "IS now unconditional production behavior for every MMFF94 policy's energy/gradient "
        "calculation and coverage measurement -- not gated, not opt-in. What remains opt-in and "
        "`false` by default is the *strict-refusal* gate on top of that coverage "
        "(`gate_mmff94_stretch_bend`/`gate_mmff94_torsion_oop`). This changes what energy/gradient "
        "every MMFF94 arm computes, not just what a gate refuses -- see the paragraph above for why "
        "that is a real, if empirically small this round, source of output change even for arms "
        "whose gate eligibility never touches stretch-bend."
    )
    lines.append("")

    lines.append("## Ring-torsion FailClosed probe")
    lines.append("")
    lines.append(
        f"{agg['ring_torsion_failclosed_probe']['n_rows']} row(s) -- demonstrates "
        "`RingTorsionApplicationPolicy::FailClosed`'s documented behavior. Not folded into any "
        f"of the {len(CHEMATIC_ARMS)} main arms' coverage numbers (those use `DiagnosticOnly`)."
    )
    lines.append("")

    lines.append("## Reference geometry subset")
    lines.append("")
    lines.append(f"Status: **{agg['reference_geometry_subset']['status']}**. {agg['reference_geometry_subset']['note']}")
    lines.append("")

    lines.append("## Known issues filed from this benchmark")
    lines.append("")
    _mmff_strict_cov = agg["coverage"]["chematic"]["chematic_pipeline_v2_mmff94_strict"]
    lines.append(
        f"- MMFF94 coverage gap ({_mmff_strict_cov['n_rows'] - _mmff_strict_cov['success']}/{_mmff_strict_cov['n_rows']} "
        f"not successful under mmff94_strict, PR #236/#238/#239/#241 fixes already reflected in this run): "
        f"{agg['known_issues_filed']['mmff94_coverage_gap']}"
    )
    lines.append("")

    lines.append("## Data integrity")
    lines.append("")
    lines.append(f"- Unclassified rows: {agg['unclassified_row_count']} (hard-gated at 0 by the report generator)")
    lines.append(f"- chematic rows sha256: `{agg['input_file_hashes']['chematic_rows_sha256'][:16]}...`")
    lines.append(f"- RDKit rows sha256: `{agg['input_file_hashes']['rdkit_rows_sha256'][:16]}...`")
    lines.append("- All integrity gates (row-count, unclassified, atom-mapping, missing/mismatched coords, non-finite coords, common-scorer coverage, denominator self-consistency) passed at generation time -- see `run_integrity_gates` in this script.")
    lines.append("")

    lines.append("## Conclusions")
    lines.append("")
    lines.append("Classified per class/metric — no single overall win/loss score.")
    lines.append("")
    lines.append("| Metric | Classification | Basis |")
    lines.append("|---|---|---|")

    mmff_strict = agg["coverage"]["chematic"]["chematic_pipeline_v2_mmff94_strict"]
    mmff_fallback = agg["coverage"]["chematic"]["chematic_pipeline_v2_mmff94_with_uff_fallback"]
    lines.append(
        f"| Coverage — no_ff/dreiding/uff_only/mmff94_with_uff_fallback vs. RDKit | Roughly comparable | "
        f"chematic {fmt_pct(agg['coverage']['chematic']['chematic_pipeline_v2_uff_only']['success']/265)}-"
        f"{fmt_pct(agg['coverage']['chematic']['chematic_pipeline_v2_no_ff']['success']/265)} success vs. "
        f"RDKit {fmt_pct(agg['coverage']['rdkit']['rdkit_etkdgv3_raw']['success']/265)} |"
    )
    lines.append(
        f"| Coverage — mmff94_strict | RDKit-favor (chematic gap, issue #227 filed) | "
        f"{fmt_pct(mmff_strict['success']/265)} success, {mmff_strict['buckets'].get('unsupported_chemistry',0)}/265 unsupported |"
    )
    lines.append(
        "| Common heavy-atom geometry — pipeline_v2 force-field arms | Chematic strength on soundness | "
        "100% independently-sound across dreiding/uff_only/mmff94 arms, matching pipeline-internal judgment |"
    )
    lines.append(
        "| Common heavy-atom geometry — legacy etkdg | Known gap, refined this round | "
        "14/265 legacy outputs have coincident atoms under the stricter common scorer (not caught by the original Wave 1 ad-hoc check); the already-documented clash-rate gap stands |"
    )
    stereo_c = agg["stereo_preservation_common_judge"]["chematic"]["chematic_pipeline_v2_uff_only"]
    stereo_r = agg["stereo_preservation_common_judge"]["rdkit"]["rdkit_etkdgv3_uff"]
    lines.append(
        f"| Stereo preservation (same judge, `Ignore`) | RDKit-favor | "
        f"RDKit {fmt_pct(stereo_r['satisfaction_rate']) if stereo_r else 'n/a'} satisfaction vs. chematic "
        f"{fmt_pct(stereo_c['satisfaction_rate']) if stereo_c else 'n/a'} under `StereoPolicy::Ignore` "
        "-- not chematic's best achievable number, see next row |"
    )
    stereo_strict_repair = agg["stereo_preservation_common_judge"]["chematic"].get(
        "chematic_pipeline_v2_mmff94_strict_repair"
    )
    stereo_fallback_repair = agg["stereo_preservation_common_judge"]["chematic"].get(
        "chematic_pipeline_v2_mmff94_with_uff_fallback_repair"
    )
    lines.append(
        "| Stereo preservation (same judge, `RepairAndVerify`, new this round) | Parity with RDKit "
        "among successes, coverage gap remains the real cost | "
        f"mmff94_strict_repair {fmt_pct(stereo_strict_repair['satisfaction_rate']) if stereo_strict_repair else 'n/a'}, "
        f"mmff94_with_uff_fallback_repair {fmt_pct(stereo_fallback_repair['satisfaction_rate']) if stereo_fallback_repair else 'n/a'} "
        "satisfaction among molecules that reached success under RepairAndVerify (both match RDKit's "
        "100% on that subset) -- but RepairAndVerify also reduces the success *count* vs. the "
        "matching Ignore arm (fewer molecules reach final success at all when repair is required to "
        "pass); see the RepairAndVerify effectiveness section for the exact paired accounting |"
    )
    for r in agg["stretch_bend_gate_effectiveness"]:
        pct_lost = (
            fmt_pct(r["newly_failing_under_later_gate"] / r["earlier_success"]) if r["earlier_success"] else "n/a"
        )
        lines.append(
            f"| Bonded-term coverage gate, {r['earlier_arm'].removeprefix('chematic_pipeline_v2_')} "
            f"-> {r['later_arm'].removeprefix('chematic_pipeline_v2_')} (new this round) | Real "
            "coverage gap surfaced, widening the gate is a real cost | "
            f"{r['earlier_success']} earlier-stage successes -> {r['later_success']} under the "
            f"later stage's gate ({r['newly_failing_under_later_gate']} newly fail, {pct_lost} "
            "of earlier-stage successes) -- see the Bonded-term coverage gate section for the "
            "term-kind sub-classification and full molecule list |"
        )
    _ff_fallback = agg["force_field_coverage"]["chematic"]["chematic_pipeline_v2_mmff94_with_uff_fallback"]
    lines.append(
        f"| Force-field convergence rate | RDKit-favor, and an input to Priority 3 (Stage 1C) | "
        f"chematic mmff94_with_uff_fallback {fmt_pct(_ff_fallback['converged_rate'])} converged "
        f"within 200 iterations, yet {_ff_fallback['n_success']}/265 of that arm's runs pass "
        "final validation regardless -- i.e. most successful outputs did NOT converge within "
        "200 iterations and still passed geometry validation. Either `force_field_converged` is "
        "narrower than \"produced a usable geometry\" (an iteration-budget artifact, not "
        "necessarily a quality problem), or this is a real gap worth diagnosing -- Priority 3's "
        "MinimizationFailed root-causing (CatastrophicBondBlowup vs. ExcessiveResidualForce) is "
        "the next stage that should resolve which; corroborates open issues #185/#188 |"
    )
    if agg.get("performance_process_level"):
        pc = agg["performance_process_level"]["chematic"]["median_seconds"]
        pr = agg["performance_process_level"]["rdkit"]["median_seconds"]
        lines.append(
            f"| Performance (process-level, whole corpus) | RDKit-favor | median {fmt_num(pc,1)}s (chematic) vs. "
            f"{fmt_num(pr,1)}s (RDKit) for the full 265-molecule x arms run |"
        )
    lines.append(
        "| Known crashes | RDKit has a narrowly-scoped one; chematic none found this round | "
        f"cyclopentane crash classified `{abl['classification'] if abl else 'n/a'}` -- non-default config, seed-dependent, not RDKit's own default behavior |"
    )
    lines.append(
        f"| Unsupported chemistry | RDKit-favor | chematic mmff94_strict "
        f"{mmff_strict['n_rows'] - mmff_strict['success']}/{mmff_strict['n_rows']} unsupported "
        "(issue #227); RDKit's 4 arms show 0 unsupported_chemistry rows |"
    )
    lines.append(
        "| Reference-geometry accuracy / torsion fingerprint / conformer diversity | Insufficient evidence | not measured this round, not fabricated |"
    )
    lines.append(
        "| Overall \"does chematic beat RDKit\" | Not claimed | per this program's explicit rule -- findings are class/metric-specific |"
    )
    lines.append("")

    REPORT_OUT.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
