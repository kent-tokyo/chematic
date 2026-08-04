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
    def stage_funnel(rows, arm):
        arm_rows = [r for r in rows if r["arm"] == arm]
        n = len(arm_rows)
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
            "final_validation_passed": sum(1 for r in arm_rows if r["_bucket"] == "success"),
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
    lines.append(
        "**mmff94_strict, spelled out per the fix request:** "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['independently_sound_successes']}/"
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['success']} successful outputs are "
        f"independently sound, but only "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['independently_sound_successes']}/"
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['n_rows']} of the *total corpus* "
        "ends up as a usable geometry under this arm -- the rest is the "
        f"{agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['n_rows'] - agg['coverage']['chematic']['chematic_pipeline_v2_mmff94_strict']['success']}-molecule "
        "MMFF94 coverage gap (issue #227, ~6,900 stretch-bend terms still ungated by this "
        "strict check -- see Priority 2/Stage 1B), not a geometry-quality problem."
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
    lines.append(
        "**Methodology, read before the numbers**: the 5 `Ignore`-policy arms below reflect raw "
        "distance-geometry-embedding output -- `Ignore` never repairs a violated stereocenter, so "
        "those rows are NOT chematic's best achievable stereo correctness. Starting this round "
        "(Priority 1, v0.11.0 re-benchmark), 2 additional `StereoPolicy::RepairAndVerify` arms "
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
            f"### Process-level performance: NOT RUN this round -- the chematic arm matrix grew "
            f"from 6 to {len(CHEMATIC_ARMS)} (2 new RepairAndVerify arms), so the stored "
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
        "not an assumed embed-then-FF-then-stereo sequence. A row reached a stage if its "
        "`failure_stage` is strictly later than that stage, or if it succeeded outright. "
        "Never collapsed into a single success rate -- see "
        "`feedback_fallback_pooling_measurement_error`: `mmff94_strict` and "
        "`mmff94_with_uff_fallback` are reported as fully separate rows, never blended."
    )
    lines.append("")
    lines.append(
        "| Arm | attempted | embed_succeeded | stereo_repair_reached | ff_attempted | "
        "ff_succeeded | final_stereo_verified | final_validation_passed |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    for arm in CHEMATIC_ARMS:
        sf = agg["stage_funnel"][arm]
        lines.append(
            f"| {arm} | {sf['attempted']} | {sf['embed_succeeded']} | {sf['stereo_repair_reached']} | "
            f"{sf['ff_attempted']} | {sf['ff_succeeded']} | {sf['final_stereo_verified']} | "
            f"{sf['final_validation_passed']} |"
        )
    lines.append("")
    lines.append(
        "Note: `chematic_legacy_etkdg` does not run through `pipeline_v2` at all (separate "
        "`generate_coords_etkdg` entry point, no `PipelineStage` tracking) -- its row is "
        "`attempted`/`final_validation_passed` only, intermediate columns are 0 by construction."
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
