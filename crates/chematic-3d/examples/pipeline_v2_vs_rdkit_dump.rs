//! Wave 1 ("RDKit alternative" program): chematic-side dump for the pipeline v2
//! vs RDKit ETKDGv3 independent benchmark.
//!
//! Calls `chematic_3d::pipeline_v2::embed_pipeline_v2` and the legacy
//! `etkdg::generate_coords_etkdg` entry point **directly** -- no Python binding
//! is used, so binding overhead/exception-conversion never mixes into the
//! measured 3D-algorithm numbers. Reads the Tier A/B corpus manifests
//! (`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`) and
//! writes one JSON row per (molecule, arm) to stdout as JSONL. Every row is
//! written -- nothing is silently dropped, including parse failures and
//! panics (caught via `catch_unwind`).
//!
//! This executable makes zero changes to `pipeline_v2.rs`, `minimize.rs`,
//! `distance_geometry_v2.rs`, or any other production algorithm file -- it is
//! a pure external caller (like `crates/chematic-py`/`crates/chematic-wasm`
//! already are), consuming only their existing public API.
//!
//! Run: `cargo run --release -p chematic-3d --example pipeline_v2_vs_rdkit_dump
//!   > validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl`

use std::collections::HashMap;
use std::panic::{self, AssertUnwindSafe};
use std::time::Instant;

use chematic_3d::coords::Coords3D;
use chematic_3d::distance_geometry_v2::EmbedParameters;
use chematic_3d::etkdg::generate_coords_etkdg;
use chematic_3d::etkdg_knowledge::TorsionOptimizationConfig;
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2::{
    self as pv2, PipelineV2Config, PipelineV2FailureCause, RingTorsionApplicationPolicy,
    StereoPolicy,
};
use chematic_3d::{ConformerDisposition, EnsembleV2Config, embed_ensemble_v2};
use chematic_core::Molecule;
use serde_json::{Value, json};

const EMBED_SEED: u64 = 20260801; // fixed for reproducibility; not a cross-platform bit-exactness claim
const MAX_ATTEMPTS: usize = 8;
const BEST_OF_N_ARM_NAME: &str = "chematic_pipeline_v2_uff_best_of_10";
const BEST_OF_N_COUNT: usize = 10;

#[derive(Clone, Copy)]
struct Arm {
    name: &'static str,
    force_field: ForceFieldPolicy,
    stereo_policy: StereoPolicy,
    gate_stretch_bend: bool,
    gate_torsion_oop: bool,
    /// `EmbedParameters.enforce_chirality` -- `false` for every pre-v0.14.0 arm
    /// (unaffected by this field, byte-identical to before it existed). The one
    /// `true` arm below is a v0.14.0 release-gate diagnostic, see its own comment.
    enforce_chirality: bool,
}

const PIPELINE_ARMS: &[Arm] = &[
    // Existing arms (Ignore StereoPolicy for backward compatibility with v0.10.0 baseline)
    Arm {
        name: "chematic_pipeline_v2_no_ff",
        force_field: ForceFieldPolicy::None,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_dreiding",
        force_field: ForceFieldPolicy::Dreiding,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_uff_only",
        force_field: ForceFieldPolicy::UffOnly,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_mmff94_strict",
        force_field: ForceFieldPolicy::Mmff94BondAngleStrict,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_mmff94_with_uff_fallback",
        force_field: ForceFieldPolicy::Mmff94WithUffFallback,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    // New arms for Priority 1 Wave 1 re-benchmark: RepairAndVerify variants
    Arm {
        name: "chematic_pipeline_v2_mmff94_strict_repair",
        force_field: ForceFieldPolicy::Mmff94BondAngleStrict,
        stereo_policy: StereoPolicy::RepairAndVerify,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_mmff94_with_uff_fallback_repair",
        force_field: ForceFieldPolicy::Mmff94WithUffFallback,
        stereo_policy: StereoPolicy::RepairAndVerify,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    // New arms for Priority 2 / Stage 1B (issue #227): "complete_term_strict_gate"
    // side of the legacy-vs-complete-term comparison -- identical to
    // chematic_pipeline_v2_mmff94_strict/..._with_uff_fallback (same
    // ForceFieldPolicy, same StereoPolicy::Ignore) except stretch-bend
    // coverage is also gated. Genuinely new, independent arms -- the
    // existing "legacy_strict_gate" arms above are NOT edited, so the delta
    // between a pair is attributable to exactly one variable.
    Arm {
        name: "chematic_pipeline_v2_mmff94_strict_stretch_bend_gated",
        force_field: ForceFieldPolicy::Mmff94BondAngleStrict,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: true,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated",
        force_field: ForceFieldPolicy::Mmff94WithUffFallback,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: true,
        gate_torsion_oop: false,
        enforce_chirality: false,
    },
    // Review-driven fix (Priority 2 blocker 1): the arms above only gate
    // bond+angle+stretch-bend -- NOT torsion/OOP, even though the audit
    // measures 1,121 missing torsion instances. These 2 new arms gate
    // bond+angle+stretch-bend+torsion+OOP -- the actual "complete bonded
    // term" coverage (still excludes vdW/charge, hence the name, not
    // "complete_mmff94"). Forms a real 3-stage comparison with the arms
    // above: legacy (bond+angle) -> stretch_bend_gated (+stretch-bend) ->
    // complete_bonded_term_gated (+torsion+OOP too).
    Arm {
        name: "chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated",
        force_field: ForceFieldPolicy::Mmff94BondAngleStrict,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: true,
        gate_torsion_oop: true,
        enforce_chirality: false,
    },
    Arm {
        name: "chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated",
        force_field: ForceFieldPolicy::Mmff94WithUffFallback,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: true,
        gate_torsion_oop: true,
        enforce_chirality: false,
    },
    // v0.14.0 release-gate diagnostic (issue #285): identical to
    // chematic_pipeline_v2_mmff94_strict in every field EXCEPT
    // enforce_chirality, isolating exactly one variable -- v0.14.0's new
    // pre-embedding declared-E/Z 1-4 bound constraint
    // (distance_geometry_v2.rs's `apply_declared_ez_bounds`). NOT a
    // comparison against the *_repair (RepairAndVerify) arms above -- those
    // use pipeline_v2's own post-hoc stage 7-11 stereo gate, a different,
    // mutually-exclusive mechanism (enforce_chirality=true requires
    // StereoPolicy::Ignore, see pipeline_v2.rs's InvalidConfiguration gate).
    // This arm's sole purpose is measuring what the new constraint changes,
    // end to end, through the actual production embed_pipeline_v2 entry
    // point real callers (Mol.embed_pipeline_v2(), embed_pipeline_v2_json)
    // use -- not a "which policy is better" comparison.
    Arm {
        name: "chematic_pipeline_v2_mmff94_strict_enforce_chirality",
        force_field: ForceFieldPolicy::Mmff94BondAngleStrict,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: true,
    },
];

fn base_config(
    force_field: ForceFieldPolicy,
    stereo_policy: StereoPolicy,
    gate_stretch_bend: bool,
    gate_torsion_oop: bool,
    enforce_chirality: bool,
) -> PipelineV2Config {
    PipelineV2Config {
        embed: EmbedParameters {
            random_seed: EMBED_SEED,
            max_attempts: MAX_ATTEMPTS,
            use_exp_torsions: true,
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            use_macrocycle_14_bounds: true,
            track_failures: true,
            enforce_chirality,
            ..EmbedParameters::default()
        },
        torsion_optimization: TorsionOptimizationConfig::default(),
        include_legacy_torsion_heuristic: false,
        // StereoPolicy::Ignore never gates success/failure on stereo, but
        // `stereo_before`/`stereo_after_repair`/`final_stereo` are still real,
        // non-fabricated evidence (see StereoPolicy::Ignore's own doc comment)
        // -- keeps "coverage" (typed success/failure) and "stereo correctness"
        // as genuinely separate metrics rather than conflating them.
        // StereoPolicy::RepairAndVerify attempts to correct stereo mismatches
        // after force-field minimization, and gates overall success on repair outcome.
        stereo_policy,
        fail_on_unevaluable_stereo: false,
        force_field_policy: force_field,
        force_field_max_iterations: 200,
        gate_mmff94_torsion_oop: gate_torsion_oop,
        gate_mmff94_stretch_bend: gate_stretch_bend,
        // DiagnosticOnly, not FailClosed: with use_small_ring_torsions/
        // use_macrocycle_torsions on, FailClosed rejects the whole pipeline
        // for nearly any ring-containing molecule (confirmed via a smoke
        // test before the full run -- cholesterol and a dedicated
        // ring-torsion fixture both failed on every one of the 5 arms).
        // That's a configuration artifact, not pipeline_v2's real achievable
        // geometry quality: DiagnosticOnly still scores those potentials
        // (ring_torsion_evidence.potentials[].applied_to_geometry stays
        // honestly `false`) without refusing to embed. FailClosed's
        // behavior on the dedicated known_fail_closed_case fixture is
        // reported separately in the aggregate report, not folded into the
        // 5 main arms' coverage numbers.
        ring_torsion_policy: RingTorsionApplicationPolicy::DiagnosticOnly,
        // A bounded safety net, not a scoring choice: an earlier full-corpus
        // run (no timeout) hung on chembl_tier_b_0124
        // ("COc1cc(C(=O)N2CCN(C(=O)c3cc(OC)c(OC)c(OC)c3)C(COC(=O)C3CCCCC3)C2)cc(OC)c1OC")
        // for several minutes on at least one arm -- a genuine slow/possibly
        // non-converging case worth reporting as a Timeout row, not a reason
        // to let the whole benchmark hang indefinitely (which would silently
        // omit that row forever, worse than reporting it as a timeout).
        total_timeout_ms: Some(20_000),
        expand_implicit_h_through_pipeline: false,
    }
}

fn coords_to_json(coords: &Coords3D) -> Value {
    let mut out = Vec::with_capacity(coords.atom_count());
    for i in 0..coords.atom_count() {
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        out.push(json!([p.x, p.y, p.z]));
    }
    Value::Array(out)
}

fn heavy_atom_elements(mol: &Molecule) -> Vec<&'static str> {
    (0..mol.atom_count())
        .map(|i| mol.atom(chematic_core::AtomIdx(i as u32)).element.symbol())
        .collect()
}

fn failure_cause_str(cause: &PipelineV2FailureCause) -> String {
    format!("{cause:?}")
}

fn stage_str(stage: &pv2::PipelineStage) -> String {
    format!("{stage:?}")
}

/// Independent (not reusing pipeline_v2's private `compute_final_validation`
/// or its `pub(crate)` constants) geometry-validity check, for the legacy
/// arm only, so pipeline_v2's own internal thresholds are never silently
/// duplicated into a second, possibly-drifting copy. Documented explicitly
/// as NOT bit-identical to `FinalGeometryValidation`'s formulas.
struct LegacyGeometryCheck {
    all_finite: bool,
    atom_count_unchanged: bool,
    worst_bond_length_ratio: f64,
    bond_violation_rate_15pct: f64,
    bond_violation_rate_50pct: f64,
    gross_clash_count: usize,
    sound: bool,
}

const LEGACY_CLASH_THRESHOLD_ANGSTROM: f64 = 1.2;
const LEGACY_MAX_SANE_BOND_LENGTH_RATIO: f64 = 3.0;

fn legacy_geometry_check(mol: &Molecule, coords: &Coords3D) -> LegacyGeometryCheck {
    let n = mol.atom_count();
    let all_finite = coords.is_finite();
    let atom_count_unchanged = coords.atom_count() == n;

    let mut worst_ratio = 0.0f64;
    let mut violations_15 = 0usize;
    let mut violations_50 = 0usize;
    let bond_count = mol.bond_count().max(1);
    for (_bidx, bond) in mol.bonds() {
        let a = bond.atom1;
        let b = bond.atom2;
        let actual = coords.get(a).distance(&coords.get(b));
        let ideal = (mol.atom(a).element.covalent_radius() as f64
            + mol.atom(b).element.covalent_radius() as f64)
            .max(0.3);
        let ratio = actual / ideal;
        let rel_error = (ratio - 1.0).abs();
        if rel_error > 0.15 {
            violations_15 += 1;
        }
        if rel_error > 0.50 {
            violations_50 += 1;
        }
        worst_ratio = worst_ratio.max(rel_error);
    }

    let mut clashes = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = coords
                .get(chematic_core::AtomIdx(i as u32))
                .distance(&coords.get(chematic_core::AtomIdx(j as u32)));
            if d < LEGACY_CLASH_THRESHOLD_ANGSTROM {
                clashes += 1;
            }
        }
    }

    let sound =
        all_finite && atom_count_unchanged && worst_ratio <= LEGACY_MAX_SANE_BOND_LENGTH_RATIO;

    LegacyGeometryCheck {
        all_finite,
        atom_count_unchanged,
        worst_bond_length_ratio: worst_ratio,
        bond_violation_rate_15pct: violations_15 as f64 / bond_count as f64,
        bond_violation_rate_50pct: violations_50 as f64 / bond_count as f64,
        gross_clash_count: clashes,
        sound,
    }
}

/// Dedicated probe demonstrating `RingTorsionApplicationPolicy::FailClosed`'s
/// documented behavior (rejects mechanically-unapplyable ring/macrocycle
/// torsion potentials rather than silently ignoring them) -- run only for
/// corpus rows tagged `known_fail_closed_case`, never as part of the 5 main
/// arms' coverage numbers (see `base_config`'s own comment for why the main
/// arms use `DiagnosticOnly` instead).
fn run_fail_closed_probe(mol: &Molecule) -> Value {
    let mut config = base_config(
        ForceFieldPolicy::Dreiding,
        StereoPolicy::Ignore,
        false,
        false,
        false,
    );
    config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
    let arm = Arm {
        name: "chematic_pipeline_v2_ring_torsion_failclosed_probe",
        force_field: ForceFieldPolicy::Dreiding,
        stereo_policy: StereoPolicy::Ignore,
        gate_stretch_bend: false,
        gate_torsion_oop: false,
        enforce_chirality: false,
    };
    run_pipeline_arm_with_config(mol, &arm, &config)
}

fn run_pipeline_arm(mol: &Molecule, arm: &Arm) -> Value {
    let config = base_config(
        arm.force_field,
        arm.stereo_policy,
        arm.gate_stretch_bend,
        arm.gate_torsion_oop,
        arm.enforce_chirality,
    );
    run_pipeline_arm_with_config(mol, arm, &config)
}

// `PipelineV2Failure` (production type, not touched here) is 1112 bytes;
// clippy's result_large_err fires on this benchmark-only call site because
// catch_unwind moves the whole Result across a generic boundary. Not a
// production defect to fix -- clippy's own suggested override.
#[allow(clippy::result_large_err)]
fn run_pipeline_arm_with_config(mol: &Molecule, arm: &Arm, config: &PipelineV2Config) -> Value {
    let start = Instant::now();
    let result = panic::catch_unwind(AssertUnwindSafe(|| pv2::embed_pipeline_v2(mol, config)));
    let elapsed_ms = start.elapsed().as_millis() as u64;

    match result {
        Err(_panic) => json!({
            "arm": arm.name,
            "status": "internal_error",
            "elapsed_ms": elapsed_ms,
        }),
        Ok(Ok(r)) => {
            let fv = &r.final_validation;
            json!({
                "arm": arm.name,
                "status": "success",
                "elapsed_ms": elapsed_ms,
                "atom_count": r.coords.atom_count(),
                "coords": coords_to_json(&r.coords),
                "embed_attempts_used": r.embed_stats.attempts_used,
                "all_finite": fv.all_finite,
                "atom_count_unchanged": fv.atom_count_unchanged,
                "worst_bond_length": fv.worst_bond_length,
                "bond_violation_rate_15pct": fv.bond_violation_rate_15pct,
                "bond_violation_rate_50pct": fv.bond_violation_rate_50pct,
                "gross_clash_count": fv.gross_clash_count,
                "ring_closure_delta": fv.ring_closure_delta,
                "sound": fv.sound,
                "stereo_before_declared": r.stereo_before.n_declared(),
                "stereo_before_satisfied": r.stereo_before.n_satisfied(),
                "stereo_before_violations": r.stereo_before.n_violations(),
                "stereo_before_unevaluable": r.stereo_before.n_unevaluable(),
                "final_stereo_declared": r.final_stereo.n_declared(),
                "final_stereo_satisfied": r.final_stereo.n_satisfied(),
                "final_stereo_violations": r.final_stereo.n_violations(),
                "final_stereo_unevaluable": r.final_stereo.n_unevaluable(),
                "stereo_repaired_count": r.stereo_repair.as_ref().map(|s| s.repaired.len()).unwrap_or(0),
                "stereo_repair_failed_count": r.stereo_repair.as_ref().map(|s| s.failures.len()).unwrap_or(0),
                "force_field_requested": format!("{:?}", r.force_field.requested_force_field),
                "force_field_actual": format!("{:?}", r.force_field.actual_force_field_used),
                "force_field_fallback": r.force_field.fallback_reason.is_some(),
                "force_field_fallback_reason": r.force_field.fallback_reason.as_ref().map(|e| format!("{e}")),
                "force_field_converged": r.force_field.converged,
                "force_field_iterations": r.force_field.iterations,
                // Only ever Some on a Mmff94WithUffFallback success-via-UFF
                // (the original failed MMFF94 attempt's coverage report
                // survives into the successful result specifically so a
                // paired-arm comparison can verify *why* the fallback fired,
                // not just that it did -- see the Priority 2 report's
                // timeout-rescue integrity check).
                "stretch_bend_missing_count": r.force_field.coverage.as_ref().map(|c| c.stretch_bend_missing.len()),
                "torsion_missing_count": r.force_field.coverage.as_ref().map(|c| c.torsions_missing.len()),
                "oop_missing_count": r.force_field.coverage.as_ref().map(|c| c.oop_missing.len()),
                "ring_torsion_potentials_total": r.ring_torsion_evidence.potentials.len(),
                "ring_torsion_potentials_applied": r
                    .ring_torsion_evidence
                    .potentials
                    .iter()
                    .filter(|p| p.applied_to_geometry)
                    .count(),
                "ring_torsion_diagnostic_only": r.ring_torsion_evidence.diagnostic_only,
            })
        }
        Ok(Err(f)) => {
            let cause_kind = failure_cause_str(&f.cause);
            let status = if matches!(f.cause, PipelineV2FailureCause::Timeout) {
                "timeout"
            } else {
                "typed_failure"
            };
            json!({
                "arm": arm.name,
                "status": status,
                "elapsed_ms": elapsed_ms,
                "failure_cause": cause_kind,
                "failure_stage": stage_str(&f.stage),
                "has_last_known_coords": f.last_known_coords.is_some(),
                // PipelineV2Failure carries partial stereo evidence computed
                // before the failure point (Some whenever that stage was
                // reached, regardless of eventual pipeline outcome) --
                // surfaced here so paired-arm RepairAndVerify analysis isn't
                // silently blind to molecules that reached repair but failed
                // at a later stage (e.g. force-field minimization).
                "stereo_before_declared": f.stereo_before.as_ref().map(|s| s.n_declared()),
                "stereo_before_satisfied": f.stereo_before.as_ref().map(|s| s.n_satisfied()),
                "stereo_before_violations": f.stereo_before.as_ref().map(|s| s.n_violations()),
                "stereo_before_unevaluable": f.stereo_before.as_ref().map(|s| s.n_unevaluable()),
                "stereo_repaired_count": f.stereo_repair.as_ref().map(|s| s.repaired.len()),
                "stereo_repair_failed_count": f.stereo_repair.as_ref().map(|s| s.failures.len()),
                "final_stereo_declared": f.final_stereo.as_ref().map(|s| s.n_declared()),
                "final_stereo_satisfied": f.final_stereo.as_ref().map(|s| s.n_satisfied()),
                "final_stereo_violations": f.final_stereo.as_ref().map(|s| s.n_violations()),
                "final_stereo_unevaluable": f.final_stereo.as_ref().map(|s| s.n_unevaluable()),
            })
        }
    }
}

fn run_legacy_arm(mol: &Molecule) -> Value {
    let start = Instant::now();
    let result = panic::catch_unwind(AssertUnwindSafe(|| generate_coords_etkdg(mol)));
    let elapsed_ms = start.elapsed().as_millis() as u64;

    match result {
        Err(_panic) => json!({
            "arm": "chematic_legacy_etkdg",
            "status": "internal_error",
            "elapsed_ms": elapsed_ms,
        }),
        Ok(coords) => {
            let check = legacy_geometry_check(mol, &coords);
            json!({
                "arm": "chematic_legacy_etkdg",
                "status": "success",
                "elapsed_ms": elapsed_ms,
                "atom_count": coords.atom_count(),
                "coords": coords_to_json(&coords),
                "all_finite": check.all_finite,
                "atom_count_unchanged": check.atom_count_unchanged,
                "worst_bond_length_ratio": check.worst_bond_length_ratio,
                "bond_violation_rate_15pct": check.bond_violation_rate_15pct,
                "bond_violation_rate_50pct": check.bond_violation_rate_50pct,
                "gross_clash_count": check.gross_clash_count,
                "sound": check.sound,
                "geometry_check_methodology": "independent_lightweight_not_pipeline_v2_internal",
            })
        }
    }
}

/// A2 best-of-N arm: `embed_ensemble_v2` with `count: BEST_OF_N_COUNT`,
/// force-field and seed matched to RDKit's own `rdkit_etkdgv3_best_of_n`
/// arm (`scripts/pipeline_v2_vs_rdkit_oracle.py`'s `run_best_of_n`:
/// `EmbedMultipleConfs(numConfs=10, randomSeed=20260801)` + per-conformer
/// UFF optimization + lowest-energy selection), so the "10 attempts, best
/// by energy" comparison means the same quantity on both sides.
#[allow(clippy::result_large_err)]
fn run_best_of_n_arm(mol: &Molecule) -> Value {
    let mut per_conformer = base_config(
        ForceFieldPolicy::UffOnly,
        StereoPolicy::Ignore,
        false,
        false,
        false,
    );
    // RDKit's EmbedMultipleConfs has no per-conformer retry loop; matching
    // its "10 draws" semantics means one attempt per conformer here too,
    // not this file's usual MAX_ATTEMPTS=8 (which retries a single slot on
    // a bad draw) -- leaving this at 8 would silently give chematic up to
    // 80 embedding attempts against RDKit's 10.
    per_conformer.embed.max_attempts = 1;
    let config = EnsembleV2Config {
        per_conformer,
        count: BEST_OF_N_COUNT,
        base_seed: EMBED_SEED,
        // Disables RMSD-dedup pruning: EmbedMultipleConfs doesn't RMSD-dedupe
        // before optimizing either, so "10 attempts" only means the same
        // thing on both sides with chematic's own pruning also off. RMSD
        // pruning is A2's own value-add, out of scope for this comparison.
        rmsd_threshold: 0.0,
        use_symmetric_rmsd_pruning: true,
        ensemble_timeout_ms: Some(20_000 * BEST_OF_N_COUNT as u64),
    };

    let start = Instant::now();
    let result = panic::catch_unwind(AssertUnwindSafe(|| embed_ensemble_v2(mol, &config)));
    let elapsed_ms = start.elapsed().as_millis() as u64;

    let r = match result {
        Err(_panic) => {
            return json!({
                "arm": BEST_OF_N_ARM_NAME,
                "status": "internal_error",
                "elapsed_ms": elapsed_ms,
            });
        }
        Ok(Err(cfg_err)) => {
            panic!("unexpected EnsembleV2ConfigError with a fixed valid config: {cfg_err}");
        }
        Ok(Ok(r)) => r,
    };

    let attempts_failed = r.attempts.iter().filter(|a| a.outcome.is_err()).count();
    let attempts_pruned = r
        .attempts
        .iter()
        .filter(|a| {
            matches!(
                &a.outcome,
                Ok(s) if matches!(s.disposition, ConformerDisposition::PrunedAsDuplicate { .. })
            )
        })
        .count();

    let Some(best_coords) = r.ensemble.get_conformer(0) else {
        // Not a new status vocabulary word: `gen_pipeline_v2_vs_rdkit_report.py`'s
        // `classify_row` exhaustively dispatches on `status`, and anything it
        // doesn't recognize becomes "unclassified" -- a hard integrity-gate
        // failure. "typed_failure" + a `failure_cause` substring keeps this
        // arm inside the report generator's existing closed vocabulary.
        return json!({
            "arm": BEST_OF_N_ARM_NAME,
            "status": "typed_failure",
            "failure_cause": "NoConformersKept",
            "elapsed_ms": elapsed_ms,
            "attempts_requested": r.requested_count,
            "conformers_kept": 0,
            "attempts_failed": attempts_failed,
            "attempts_pruned": attempts_pruned,
            "termination": format!("{:?}", r.termination),
        });
    };

    // `ensemble` is ordered group-then-ascending-energy; under this
    // single-policy (UffOnly) config there is exactly one force-field group
    // and UFF energy is always present on success, so index 0 is
    // unambiguously "the lowest-energy kept conformer." A future config
    // change to this arm (e.g. a mixed-policy ensemble) would break that
    // assumption -- reverse-mapping back to the originating attempt below
    // (rather than trusting index 0 alone) keeps `best_energy` honest.
    let best = r
        .attempts
        .iter()
        .find_map(|a| match &a.outcome {
            Ok(s)
                if matches!(
                    s.disposition,
                    ConformerDisposition::Kept { conformer_index: 0 }
                ) =>
            {
                Some(s)
            }
            _ => None,
        })
        .expect("ensemble.get_conformer(0) implies some attempt has Kept{conformer_index: 0}");

    let check = legacy_geometry_check(mol, best_coords);
    json!({
        "arm": BEST_OF_N_ARM_NAME,
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "atom_count": best_coords.atom_count(),
        "coords": coords_to_json(best_coords),
        "attempts_requested": r.requested_count,
        "conformers_kept": r.ensemble.conformer_count(),
        "attempts_failed": attempts_failed,
        "attempts_pruned": attempts_pruned,
        "termination": format!("{:?}", r.termination),
        "mixed_force_field": r.mixed_force_field,
        "best_energy": best.energy,
        "all_finite": check.all_finite,
        "atom_count_unchanged": check.atom_count_unchanged,
        "worst_bond_length_ratio": check.worst_bond_length_ratio,
        "bond_violation_rate_15pct": check.bond_violation_rate_15pct,
        "bond_violation_rate_50pct": check.bond_violation_rate_50pct,
        "gross_clash_count": check.gross_clash_count,
        "sound": check.sound,
        "geometry_check_methodology": "independent_lightweight_not_pipeline_v2_internal",
    })
}

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

fn main() {
    let mut manifests: Vec<(String, Value)> = Vec::new();
    for (tier, path) in [
        (
            "A",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
        ),
        (
            "B",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
        ),
    ] {
        manifests.push((tier.to_string(), load_manifest(path)));
    }

    let config_snapshot: HashMap<&str, Value> = PIPELINE_ARMS
        .iter()
        .map(|arm| {
            (
                arm.name,
                json!(format!(
                    "ff={:?} stereo={:?} gate_stretch_bend={} gate_torsion_oop={}",
                    arm.force_field, arm.stereo_policy, arm.gate_stretch_bend, arm.gate_torsion_oop
                )),
            )
        })
        .collect();
    eprintln!(
        "config_snapshot embed_seed={EMBED_SEED} max_attempts={MAX_ATTEMPTS} arms={config_snapshot:?}"
    );
    eprintln!(
        "config_snapshot {BEST_OF_N_ARM_NAME}: ff=UffOnly count={BEST_OF_N_COUNT} \
         base_seed={EMBED_SEED} max_attempts=1 (NOT {MAX_ATTEMPTS} -- one draw per \
         conformer, matching RDKit's EmbedMultipleConfs) rmsd_threshold=0.0 (pruning \
         disabled, for parity with RDKit's no-dedup best-of-N selection)"
    );

    for (tier, manifest) in &manifests {
        let molecules = manifest["molecules"].as_array().expect("molecules array");
        for m in molecules {
            let name = m["name"].as_str().unwrap();
            let smiles = m["smiles"].as_str().unwrap();
            let primary_category = m["primary_category"].as_str().unwrap_or("unknown");

            let parsed = panic::catch_unwind(AssertUnwindSafe(|| chematic_smiles::parse(smiles)));
            let mol = match parsed {
                Err(_) | Ok(Err(_)) => {
                    let row = json!({
                        "tier": tier,
                        "name": name,
                        "smiles": smiles,
                        "primary_category": primary_category,
                        "arm": "all",
                        "status": "parse_failure",
                    });
                    println!("{row}");
                    continue;
                }
                Ok(Ok(mol)) => mol,
            };

            let elements = heavy_atom_elements(&mol);

            for arm in PIPELINE_ARMS {
                let mut row = run_pipeline_arm(&mol, arm);
                row["tier"] = json!(tier);
                row["name"] = json!(name);
                row["smiles"] = json!(smiles);
                row["primary_category"] = json!(primary_category);
                row["heavy_atom_elements"] = json!(elements);
                println!("{row}");
            }

            let mut legacy_row = run_legacy_arm(&mol);
            legacy_row["tier"] = json!(tier);
            legacy_row["name"] = json!(name);
            legacy_row["smiles"] = json!(smiles);
            legacy_row["primary_category"] = json!(primary_category);
            legacy_row["heavy_atom_elements"] = json!(elements);
            println!("{legacy_row}");

            let mut best_of_n_row = run_best_of_n_arm(&mol);
            best_of_n_row["tier"] = json!(tier);
            best_of_n_row["name"] = json!(name);
            best_of_n_row["smiles"] = json!(smiles);
            best_of_n_row["primary_category"] = json!(primary_category);
            best_of_n_row["heavy_atom_elements"] = json!(elements);
            println!("{best_of_n_row}");

            if primary_category == "known_fail_closed_case" {
                let mut probe_row = run_fail_closed_probe(&mol);
                probe_row["tier"] = json!(tier);
                probe_row["name"] = json!(name);
                probe_row["smiles"] = json!(smiles);
                probe_row["primary_category"] = json!(primary_category);
                probe_row["heavy_atom_elements"] = json!(elements);
                println!("{probe_row}");
            }
        }
    }
}
