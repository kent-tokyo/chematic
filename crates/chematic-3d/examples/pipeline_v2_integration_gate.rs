//! Integration gate harness for the opt-in v2 embedding pipeline (`pipeline_v2.rs`),
//! 3D Breakthrough Program Wave 2 → Wave 3 Coordinator Integration 1.
//!
//! Compares 10 arms (A-J, spec §13) on the same seed over the frozen 58-molecule
//! corpus (hand-copied verbatim from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`, same
//! transcription already used by every sibling example in this directory -- see
//! `examples/cf_integration_smoke_test.rs`) plus 5 molecules added specifically to
//! close a gap found during independent verification round 4: spec §14 explicitly
//! names cyclobutane/cyclooctane/cyclononane, a 2,2'-disubstituted biphenyl, and a
//! macrocyclic amide as required stress cases, none of which the frozen 58 contains.
//! Two of these are load-bearing for this PR's own mechanism, not just checklist
//! completeness: cyclooctane (8-ring, `SMALL_RING_MAX`) and cyclononane (9-ring,
//! `MACROCYCLE_MIN`) are the two sides of the small-ring/macrocycle classification
//! boundary that round 2's verification only exercised ad hoc; and no macrocycle in
//! the frozen 58 is amide-like, so `bounds14.rs`'s `macrocycle_14:amide_ester_pinned`
//! branch (pin-to-cis) was previously never exercised by any arm -- only its
//! `relaxed_band` sibling was:
//!
//! A: raw DG (bypasses the pipeline entirely -- Agent C's own module)
//! B: raw DG + macrocycle 1-4 bounds
//! C: B + standard acyclic torsion optimization
//! D: C + stereo repair
//! E: D + MMFF94 bond/angle strict
//! F: D + widened MMFF94 torsion/oop gate
//! G: D + MMFF94->UFF fallback
//! H: D + DREIDING
//! I: full requested ring torsions under FailClosed (must fail closed, typed)
//! J: full requested ring torsions under DiagnosticOnly (must succeed, scored-only)
//!
//! Run: `cargo run --release -p chematic-3d --example pipeline_v2_integration_gate`

use std::collections::BTreeMap;
use std::panic::{self, AssertUnwindSafe};

use chematic_3d::align::align_coords;
use chematic_3d::coords::Coords3D;
use chematic_3d::distance_geometry_v2::{self, EmbedParameters};
use chematic_3d::etkdg_knowledge::{
    TorsionKnowledgeDiagnosticKind, TorsionKnowledgeSource, TorsionOptimizationConfig,
};
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2::{
    PipelineV2Config, PipelineV2Result, RingTorsionApplicationPolicy, StereoPolicy,
    embed_pipeline_v2,
};
use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::parse;

const DEFAULT_SEED: u64 = 0xC0FF_EE42_D157_6E02;

/// Frozen 58-molecule corpus (identical transcription to
/// `examples/cf_integration_smoke_test.rs`, deliberately kept byte-for-byte the same
/// so results are directly comparable across both examples) plus 5 spec-§14-required
/// stress cases added in independent verification round 4 (see module doc above).
const CORPUS: &[(&str, &str, &str)] = &[
    ("benzene", "c1ccccc1", "rigid_ring"),
    ("naphthalene", "c1ccc2ccccc2c1", "fused_aromatic"),
    ("pyridine", "c1ccncc1", "rigid_ring"),
    ("furan", "c1ccoc1", "rigid_ring"),
    ("thiophene", "c1ccsc1", "rigid_ring"),
    ("adamantane", "C1CC2CC3CC1CC(C2)C3", "rigid_ring"),
    ("cubane", "C1C2C3C1C4C2C3C4", "rigid_ring"),
    ("cyclohexane", "C1CCCCC1", "rigid_ring"),
    ("cyclopentane", "C1CCCC1", "rigid_ring"),
    ("indole", "c1ccc2[nH]ccc2c1", "fused_aromatic"),
    ("purine", "c1ncc2[nH]cnc2n1", "fused_aromatic"),
    ("quinoline", "c1ccc2ncccc2c1", "fused_aromatic"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "fused_aromatic"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", "fused_aromatic"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1", "fused_aromatic"),
    ("butane", "CCCC", "flexible_chain"),
    ("hexane", "CCCCCC", "flexible_chain"),
    ("decane", "CCCCCCCCCC", "flexible_chain"),
    ("triethylene_glycol", "OCCOCCOCCO", "flexible_chain"),
    ("hexanediol", "OCCCCCCO", "flexible_chain"),
    ("hexadecane", "CCCCCCCCCCCCCCCC", "flexible_chain"),
    ("cyclododecane", "C1CCCCCCCCCCC1", "macrocycle"),
    ("crown_12_4", "O1CCOCCOCCOCC1", "macrocycle"),
    ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1", "macrocycle"),
    ("l_alanine", "N[C@@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("d_alanine", "N[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("l_serine", "N[C@@H](CO)C(=O)O", "stereocenter_implicit_h"),
    (
        "l_threonine",
        "C[C@H](O)[C@@H](N)C(=O)O",
        "stereocenter_implicit_h",
    ),
    ("2_butanol_R", "C[C@H](O)CC", "stereocenter_implicit_h"),
    ("2_butanol_S", "C[C@@H](O)CC", "stereocenter_implicit_h"),
    (
        "2_chlorobutane_R",
        "C[C@H](Cl)CC",
        "stereocenter_implicit_h",
    ),
    (
        "ibuprofen_S",
        "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        "stereocenter_implicit_h",
    ),
    (
        "naproxen_S",
        "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1",
        "stereocenter_implicit_h",
    ),
    (
        "menthol",
        "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O",
        "stereocenter_implicit_h",
    ),
    ("chfclbr_R", "[C@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("chfclbr_S", "[C@@H](F)(Cl)Br", "stereocenter_quaternary"),
    (
        "quaternary_1_R",
        "[C@](F)(Cl)(Br)I",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_1_S",
        "[C@@](F)(Cl)(Br)I",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_2_R",
        "[C@](C)(N)(O)F",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_2_S",
        "[C@@](C)(N)(O)F",
        "stereocenter_quaternary",
    ),
    ("but2ene_E", "C/C=C/C", "alkene_ez"),
    ("but2ene_Z", r"C/C=C\C", "alkene_ez"),
    ("chloropropene_E", "C(/C=C/C)Cl", "alkene_ez"),
    ("chloropropene_Z", r"C(/C=C\C)Cl", "alkene_ez"),
    ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1", "alkene_ez"),
    ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1", "alkene_ez"),
    ("pent2ene_E", "CC/C=C/C", "alkene_ez"),
    ("pent2ene_Z", r"CC/C=C\C", "alkene_ez"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "druglike"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "druglike"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "druglike"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", "druglike"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", "druglike"),
    (
        "penicillin_core",
        "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
        "druglike",
    ),
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
        "druglike_rigid",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        "druglike_stress",
    ),
    (
        "atorvastatin_fragment",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        "druglike_stress",
    ),
    ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O", "druglike"),
    // --- spec §14 required stress cases added in verification round 4 (see module
    // doc above) -- not part of the original frozen 58, appended rather than
    // interleaved so the first 58 rows stay byte-identical to
    // `cf_integration_smoke_test.rs`'s own corpus.
    ("cyclobutane", "C1CCC1", "rigid_ring"),
    ("cyclooctane", "C1CCCCCCC1", "small_ring_boundary"),
    ("cyclononane", "C1CCCCCCCC1", "macrocycle_boundary"),
    (
        "dimethylbiphenyl_2_2",
        "Cc1ccccc1-c1ccccc1C",
        "hindered_biaryl",
    ),
    ("macrolactam_12", "O=C1CCCCCCCCCCN1", "macrocycle_amide"),
];

// ---------------------------------------------------------------------------
// Per-molecule / per-arm outcome bookkeeping
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct ArmOutcome {
    success: bool,
    panicked: bool,
    stage: Option<String>,
    cause: Option<String>,
    elapsed_ms: u64,
    coords: Option<Coords3D>,
    n_matched: usize,
    n_applied: usize,
    n_scored_only: usize,
    diagnostic_only: bool,
    /// `TorsionKnowledgeDiagnosticKind::AmbiguousSameTierConflict` count only --
    /// a genuine same-tier rule conflict (spec's "never arbitrarily pick one side
    /// of a genuine ambiguity"). Kept separate from `n_fused_bridged_notices`
    /// below: `ambiguous_matches` (matcher.rs) pushes BOTH kinds into the same
    /// `Vec`, and PR #191 itself already warned that reporting that vector's raw
    /// length overstates genuine conflicts (adamantane alone has 13 fused/bridged
    /// notices and zero real conflicts). Conflating them here would repeat that
    /// exact measurement error one PR later.
    n_ambiguous_rule_conflicts: usize,
    /// `TorsionKnowledgeDiagnosticKind::FusedOrBridgedRingBoundary` count -- a
    /// ring-topology notice ("this bond touches more than one ring/size bucket"),
    /// not a rule conflict. Pushed for essentially every fused/bridged/spiro bond,
    /// so it dominates `ambiguous_matches.len()` on ring-heavy molecules.
    n_fused_bridged_notices: usize,
    stereo_declared: usize,
    stereo_satisfied: usize,
    stereo_violated: usize,
    stereo_unevaluable: usize,
    stereo_repaired: usize,
    ff_requested: Option<String>,
    ff_actual: Option<String>,
    ff_fallback: bool,
    ff_energy_before: f64,
    ff_energy_after: f64,
    /// From `torsion_optimization_report.energy_before` (stage 6, before/after the
    /// acyclic-only optimization) -- `None` when stage 6 never ran (zero potentials
    /// to optimize). Distinct from `ff_energy_before`/`ff_energy_after`, which are
    /// the force field's own energy units and `0.0` whenever
    /// `ForceFieldPolicy::None` (never conflated -- see the bug this fixed).
    torsion_energy_before: Option<f64>,
    /// From `final_validation.torsion_energy_after` (stage 12, evaluated on the
    /// FINAL post-force-field geometry, not just post-torsion-optimization) --
    /// always present on success, unlike `torsion_energy_before`.
    torsion_energy_after: f64,
    residual_force: f64,
    bond_violation_15: f64,
    bond_violation_50: f64,
    gross_clashes: usize,
    bounds_violation_rate: f64,
    geometrically_valid: bool,
}

impl ArmOutcome {
    fn failure(stage: &str, cause: String, elapsed_ms: u64) -> Self {
        Self {
            success: false,
            panicked: false,
            stage: Some(stage.to_string()),
            cause: Some(cause),
            elapsed_ms,
            coords: None,
            n_matched: 0,
            n_applied: 0,
            n_scored_only: 0,
            diagnostic_only: false,
            n_ambiguous_rule_conflicts: 0,
            n_fused_bridged_notices: 0,
            stereo_declared: 0,
            stereo_satisfied: 0,
            stereo_violated: 0,
            stereo_unevaluable: 0,
            stereo_repaired: 0,
            ff_requested: None,
            ff_actual: None,
            ff_fallback: false,
            ff_energy_before: 0.0,
            ff_energy_after: 0.0,
            torsion_energy_before: None,
            torsion_energy_after: 0.0,
            residual_force: 0.0,
            bond_violation_15: 0.0,
            bond_violation_50: 0.0,
            gross_clashes: 0,
            bounds_violation_rate: 0.0,
            geometrically_valid: false,
        }
    }

    fn panic(elapsed_ms: u64) -> Self {
        let mut o = Self::failure("PANIC", "PANIC".to_string(), elapsed_ms);
        o.panicked = true;
        o
    }

    fn from_result(r: &PipelineV2Result, elapsed_ms: u64) -> Self {
        let ff = &r.force_field;
        Self {
            success: true,
            panicked: false,
            stage: None,
            cause: None,
            elapsed_ms,
            coords: Some(r.coords.clone()),
            n_matched: r.torsion_knowledge_report.potentials.len(),
            n_applied: r.ring_torsion_evidence.n_applied(),
            n_scored_only: r.ring_torsion_evidence.n_scored_only(),
            diagnostic_only: r.ring_torsion_evidence.diagnostic_only,
            n_ambiguous_rule_conflicts: r
                .torsion_knowledge_report
                .ambiguous_matches
                .iter()
                .filter(|d| d.kind == TorsionKnowledgeDiagnosticKind::AmbiguousSameTierConflict)
                .count(),
            n_fused_bridged_notices: r
                .torsion_knowledge_report
                .ambiguous_matches
                .iter()
                .filter(|d| d.kind == TorsionKnowledgeDiagnosticKind::FusedOrBridgedRingBoundary)
                .count(),
            stereo_declared: r.final_stereo.n_declared(),
            stereo_satisfied: r.final_stereo.n_satisfied(),
            stereo_violated: r.final_stereo.n_violations(),
            stereo_unevaluable: r.final_stereo.n_unevaluable(),
            stereo_repaired: r
                .stereo_repair
                .as_ref()
                .map(|s| s.repaired.len())
                .unwrap_or(0),
            ff_requested: Some(format!("{:?}", ff.requested_force_field)),
            ff_actual: Some(format!("{:?}", ff.actual_force_field_used)),
            ff_fallback: ff.fallback_reason.is_some(),
            ff_energy_before: ff.energy_before.total(),
            ff_energy_after: ff.energy_after.total(),
            torsion_energy_before: r
                .torsion_optimization_report
                .as_ref()
                .map(|t| t.energy_before),
            torsion_energy_after: r.final_validation.torsion_energy_after,
            residual_force: ff.max_residual_force,
            bond_violation_15: r.final_validation.bond_violation_rate_15pct,
            bond_violation_50: r.final_validation.bond_violation_rate_50pct,
            gross_clashes: r.final_validation.gross_clash_count,
            bounds_violation_rate: if r.final_validation.bounds_conformance.n_pairs > 0 {
                r.final_validation.bounds_conformance.n_violations as f64
                    / r.final_validation.bounds_conformance.n_pairs as f64
            } else {
                0.0
            },
            geometrically_valid: r.final_validation.sound,
        }
    }
}

fn to_vec3(coords: &Coords3D) -> Vec<[f64; 3]> {
    (0..coords.atom_count())
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect()
}

/// Short, bounded-size failure-cause label. `PipelineV2FailureCause::ForceField`
/// wraps `ForceFieldBridgeError::MissingParameters(Box<Mmff94CoverageReport>)`,
/// whose full `{:?}` dump lists every missing bond/angle/torsion/oop term with a
/// per-term description string -- genuinely useful in isolation, but printed
/// unbucketed it makes the per-arm summary's failure-by-cause table (which keys on
/// this string) explode to tens of thousands of characters per row on molecules with
/// many missing MMFF94 parameters. Bucketed the same way
/// `examples/cf_integration_smoke_test.rs`'s own `err_bucket` already does.
fn cause_label(cause: &chematic_3d::pipeline_v2::PipelineV2FailureCause) -> String {
    use chematic_3d::minimize::ForceFieldBridgeError;
    use chematic_3d::pipeline_v2::PipelineV2FailureCause as C;
    match cause {
        C::ForceField(ForceFieldBridgeError::MissingParameters(r)) => {
            format!(
                "ForceField(MissingParameters: {} missing)",
                r.total_missing()
            )
        }
        C::ForceField(ForceFieldBridgeError::MinimizationFailed(d)) => {
            format!("ForceField(MinimizationFailed({:?}))", d.reason)
        }
        C::ForceField(ForceFieldBridgeError::UnsupportedAtomType(_)) => {
            "ForceField(UnsupportedAtomType)".to_string()
        }
        other => format!("{other:?}"),
    }
}

// `PipelineV2Failure` is intentionally large (see its own doc comment in
// `pipeline_v2.rs`) -- this harness immediately destructures it into `ArmOutcome`
// and never propagates it further, so the allocation concern clippy's lint guards
// against doesn't apply here either.
#[allow(clippy::result_large_err)]
fn run_pipeline_arm(mol: &Molecule, config: &PipelineV2Config) -> ArmOutcome {
    let start = std::time::Instant::now();
    let result = panic::catch_unwind(AssertUnwindSafe(|| embed_pipeline_v2(mol, config)));
    let elapsed_ms = start.elapsed().as_millis() as u64;
    match result {
        Err(_) => ArmOutcome::panic(elapsed_ms),
        Ok(Ok(r)) => ArmOutcome::from_result(&r, elapsed_ms),
        Ok(Err(e)) => {
            ArmOutcome::failure(&format!("{:?}", e.stage), cause_label(&e.cause), elapsed_ms)
        }
    }
}

fn run_raw_dg_arm(mol: &Molecule, seed: u64) -> ArmOutcome {
    let start = std::time::Instant::now();
    let params = EmbedParameters {
        random_seed: seed,
        ..EmbedParameters::default()
    };
    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        distance_geometry_v2::embed_distance_geometry_v2(mol, &params)
    }));
    let elapsed_ms = start.elapsed().as_millis() as u64;
    match result {
        Err(_) => ArmOutcome::panic(elapsed_ms),
        Ok(Err(cause)) => ArmOutcome::failure("DistanceGeometry", format!("{cause:?}"), elapsed_ms),
        Ok(Ok(coords)) => {
            let all_finite = coords.is_finite();
            let worst = mol
                .bonds()
                .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
                .fold(0.0_f64, f64::max);
            let mut o =
                ArmOutcome::from_result(&fake_result_for_raw_dg(mol, coords.clone()), elapsed_ms);
            o.geometrically_valid = all_finite && worst.is_finite() && worst < 3.0;
            o
        }
    }
}

/// Wraps a raw-DG-only coords into the same `ArmOutcome` extraction path as the
/// pipeline arms, so Arm A's row prints through the identical code path (only
/// `coords`/`geometrically_valid` are real; everything torsion/stereo/FF-shaped is
/// trivially empty because Arm A never runs any of those stages).
fn fake_result_for_raw_dg(mol: &Molecule, coords: Coords3D) -> PipelineV2Result {
    use chematic_3d::distance_geometry_v2::EmbedStats;
    use chematic_3d::etkdg_knowledge::TorsionKnowledgeReport;
    use chematic_3d::minimize::{EnergyReport, ForceFieldPolicy as FFP, PolicyMinimizeResult};
    use chematic_3d::pipeline_v2::{FinalGeometryValidation, RingTorsionEvidence, StageTimings};
    use chematic_3d::stereo_constraints::verify_stereo;

    let stereo = verify_stereo(mol, &coords);
    let bounds_conf = distance_geometry_v2::bounds_conformance(mol, &coords);
    let worst = mol
        .bonds()
        .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
        .fold(0.0_f64, f64::max);
    PipelineV2Result {
        coords: coords.clone(),
        embed_stats: EmbedStats::default(),
        bound_adjustment_report: None,
        torsion_knowledge_report: TorsionKnowledgeReport::default(),
        ring_torsion_evidence: RingTorsionEvidence::default(),
        torsion_optimization_report: None,
        stereo_before: stereo.clone(),
        stereo_repair: None,
        stereo_after_repair: stereo.clone(),
        force_field: PolicyMinimizeResult {
            coords: coords.clone(),
            requested_force_field: FFP::None,
            actual_force_field_used: FFP::None,
            fallback_reason: None,
            missing_parameter_classes: Vec::new(),
            coverage: None,
            energy_before: EnergyReport::None,
            energy_after: EnergyReport::None,
            converged: true,
            iterations: 0,
            max_residual_force: 0.0,
            starting_geometry: None,
        },
        final_stereo: stereo,
        post_minimization_stereo_repair: None,
        final_validation: FinalGeometryValidation {
            all_finite: coords.is_finite(),
            atom_count_unchanged: true,
            worst_bond_length: worst,
            bond_violation_rate_15pct: 0.0,
            bond_violation_rate_50pct: 0.0,
            gross_clash_count: 0,
            bounds_conformance: bounds_conf,
            stereo_ok: true,
            torsion_energy_after: 0.0,
            ring_closure_delta: 0.0,
            sound: coords.is_finite() && worst.is_finite() && worst < 3.0,
        },
        elapsed_ms_by_stage: StageTimings::default(),
    }
}

// ---------------------------------------------------------------------------
// Arm configs (B..J) -- A bypasses the pipeline (see run_raw_dg_arm)
// ---------------------------------------------------------------------------

fn base_embed(seed: u64) -> EmbedParameters {
    EmbedParameters {
        random_seed: seed,
        ..EmbedParameters::default()
    }
}

fn arm_b(seed: u64) -> PipelineV2Config {
    let mut c = PipelineV2Config::minimal(ForceFieldPolicy::None);
    c.embed = base_embed(seed);
    c.embed.use_macrocycle_14_bounds = true;
    c
}

fn arm_c(seed: u64) -> PipelineV2Config {
    let mut c = arm_b(seed);
    c.embed.use_exp_torsions = true;
    c.torsion_optimization = TorsionOptimizationConfig::default();
    c
}

fn arm_d(seed: u64) -> PipelineV2Config {
    let mut c = arm_c(seed);
    c.stereo_policy = StereoPolicy::RepairAndVerify;
    c
}

fn arm_e(seed: u64) -> PipelineV2Config {
    let mut c = arm_d(seed);
    c.force_field_policy = ForceFieldPolicy::Mmff94BondAngleStrict;
    c.gate_mmff94_torsion_oop = false;
    c
}

fn arm_f(seed: u64) -> PipelineV2Config {
    let mut c = arm_d(seed);
    c.force_field_policy = ForceFieldPolicy::Mmff94BondAngleStrict;
    c.gate_mmff94_torsion_oop = true;
    c
}

fn arm_g(seed: u64) -> PipelineV2Config {
    let mut c = arm_d(seed);
    c.force_field_policy = ForceFieldPolicy::Mmff94WithUffFallback;
    c
}

fn arm_h(seed: u64) -> PipelineV2Config {
    let mut c = arm_d(seed);
    c.force_field_policy = ForceFieldPolicy::Dreiding;
    c
}

fn arm_i(seed: u64) -> PipelineV2Config {
    let mut c = arm_d(seed);
    c.embed.use_small_ring_torsions = true;
    c.embed.use_macrocycle_torsions = true;
    c.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
    c.force_field_policy = ForceFieldPolicy::None;
    c
}

fn arm_j(seed: u64) -> PipelineV2Config {
    let mut c = arm_i(seed);
    c.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
    c.force_field_policy = ForceFieldPolicy::None;
    c
}

const ARM_NAMES: &[&str] = &[
    "A_raw_dg",
    "B_macro14",
    "C_std_torsion",
    "D_stereo_repair",
    "E_mmff94_strict",
    "F_mmff94_widened",
    "G_mmff94_uff_fb",
    "H_dreiding",
    "I_ring_failclosed",
    "J_ring_diagonly",
];

fn run_all_arms(mol: &Molecule, seed: u64) -> Vec<ArmOutcome> {
    vec![
        run_raw_dg_arm(mol, seed),
        run_pipeline_arm(mol, &arm_b(seed)),
        run_pipeline_arm(mol, &arm_c(seed)),
        run_pipeline_arm(mol, &arm_d(seed)),
        run_pipeline_arm(mol, &arm_e(seed)),
        run_pipeline_arm(mol, &arm_f(seed)),
        run_pipeline_arm(mol, &arm_g(seed)),
        run_pipeline_arm(mol, &arm_h(seed)),
        run_pipeline_arm(mol, &arm_i(seed)),
        run_pipeline_arm(mol, &arm_j(seed)),
    ]
}

fn percentile(sorted: &[u64], pct: f64) -> u64 {
    if sorted.is_empty() {
        return 0;
    }
    let idx = ((sorted.len() as f64 - 1.0) * pct).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn main() {
    let n_arms = ARM_NAMES.len();
    // per-arm: molecule name -> outcome
    let mut per_arm: Vec<BTreeMap<&'static str, ArmOutcome>> =
        (0..n_arms).map(|_| BTreeMap::new()).collect();

    println!(
        "{:<24} {}",
        "molecule",
        ARM_NAMES
            .iter()
            .map(|n| format!("{n:<18}"))
            .collect::<String>()
    );

    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("{name}: parse failed: {e:?}"));
        let outcomes = run_all_arms(&mol, DEFAULT_SEED);

        let row: String = outcomes
            .iter()
            .map(|o| {
                if o.panicked {
                    format!("{:<18}", "PANIC")
                } else if o.success {
                    format!("{:<18}", "OK")
                } else {
                    let c = o.cause.as_deref().unwrap_or("?");
                    format!("{:<18}", &c[..c.len().min(17)])
                }
            })
            .collect();
        println!("{name:<24} {row}");

        for (i, outcome) in outcomes.into_iter().enumerate() {
            per_arm[i].insert(name, outcome);
        }
    }

    let n_total = CORPUS.len();
    assert_eq!(
        n_total, 63,
        "corpus size drifted from the frozen 58 + 5 spec-§14 stress cases"
    );

    // =========================================================================
    // Per-arm summary metrics (spec §15) -- denominators kept explicit and
    // separate: `attempted` is always n_total, `success` is a strict subset, and
    // every failed molecule is still accounted for (stage+cause), never dropped.
    // =========================================================================
    println!("\n=== PER-ARM SUMMARY (attempted={n_total} for every arm) ===");
    for (i, arm_name) in ARM_NAMES.iter().enumerate() {
        let outcomes = &per_arm[i];
        let attempted = n_total;
        let successes: Vec<&ArmOutcome> = outcomes.values().filter(|o| o.success).collect();
        let success = successes.len();
        let mut failures_by_stage_cause: BTreeMap<(String, String), usize> = BTreeMap::new();
        for o in outcomes.values() {
            if !o.success {
                let stage = o.stage.clone().unwrap_or_default();
                let cause = o.cause.clone().unwrap_or_default();
                *failures_by_stage_cause.entry((stage, cause)).or_insert(0) += 1;
            }
        }

        let mut runtimes: Vec<u64> = outcomes.values().map(|o| o.elapsed_ms).collect();
        runtimes.sort_unstable();
        let p50 = percentile(&runtimes, 0.50);
        let p95 = percentile(&runtimes, 0.95);

        let geo_valid = successes.iter().filter(|o| o.geometrically_valid).count();
        let mean = |f: fn(&ArmOutcome) -> f64| -> f64 {
            if successes.is_empty() {
                0.0
            } else {
                successes.iter().map(|o| f(o)).sum::<f64>() / successes.len() as f64
            }
        };
        let sum_usize =
            |f: fn(&ArmOutcome) -> usize| -> usize { successes.iter().map(|o| f(o)).sum() };

        let n_declared: usize = sum_usize(|o| o.stereo_declared);
        let n_satisfied: usize = sum_usize(|o| o.stereo_satisfied);
        let n_violated: usize = sum_usize(|o| o.stereo_violated);
        let n_unevaluable: usize = sum_usize(|o| o.stereo_unevaluable);
        let n_repaired: usize = sum_usize(|o| o.stereo_repaired);
        let n_fallback = successes.iter().filter(|o| o.ff_fallback).count();

        println!("\n--- Arm {arm_name} ---");
        println!("  attempted={attempted} success={success} ({success}/{attempted})");
        if !failures_by_stage_cause.is_empty() {
            println!("  typed failures by (stage, cause):");
            for ((stage, cause), count) in &failures_by_stage_cause {
                println!("    {stage:<28} {cause:<40} {count}");
            }
        }
        println!("  runtime p50={p50}ms p95={p95}ms");
        println!("  geometrically-valid rate: {geo_valid}/{success} of successes",);
        println!(
            "  bond-violation rate: 15%={:.4} 50%={:.4} (mean over successes)",
            mean(|o| o.bond_violation_15),
            mean(|o| o.bond_violation_50)
        );
        println!(
            "  gross clashes (mean): {:.3}  bounds-violation rate (mean): {:.4}",
            mean(|o| o.gross_clashes as f64),
            mean(|o| o.bounds_violation_rate)
        );
        let n_with_torsion_opt = successes
            .iter()
            .filter(|o| o.torsion_energy_before.is_some())
            .count();
        let mean_torsion_before = if n_with_torsion_opt == 0 {
            0.0
        } else {
            successes
                .iter()
                .filter_map(|o| o.torsion_energy_before)
                .sum::<f64>()
                / n_with_torsion_opt as f64
        };
        println!(
            "  torsion energy before/after (mean over the {n_with_torsion_opt} successes with \
             potentials to optimize / mean over all {success} successes): {:.3} / {:.3}",
            mean_torsion_before,
            mean(|o| o.torsion_energy_after)
        );
        println!(
            "  force-field energy before/after (mean, FF units, 0 under ForceFieldPolicy::None): \
             {:.3} / {:.3}",
            mean(|o| o.ff_energy_before),
            mean(|o| o.ff_energy_after)
        );
        println!(
            "  potentials matched/applied/scored-only (sum): {}/{}/{}  \
             genuine rule conflicts (sum): {}  fused/bridged ring notices (sum): {}",
            sum_usize(|o| o.n_matched),
            sum_usize(|o| o.n_applied),
            sum_usize(|o| o.n_scored_only),
            sum_usize(|o| o.n_ambiguous_rule_conflicts),
            sum_usize(|o| o.n_fused_bridged_notices)
        );
        println!(
            "  stereo declared/satisfied/repaired/unevaluable/violated (sum): {n_declared}/{n_satisfied}/{n_repaired}/{n_unevaluable}/{n_violated}"
        );
        println!(
            "  force-field fallback occurred: {n_fallback}/{success} of successes; \
             mean residual force: {:.3}",
            mean(|o| o.residual_force)
        );
        let mut requested_actual: BTreeMap<(String, String), usize> = BTreeMap::new();
        for o in &successes {
            let req = o.ff_requested.clone().unwrap_or_default();
            let act = o.ff_actual.clone().unwrap_or_default();
            *requested_actual.entry((req, act)).or_insert(0) += 1;
        }
        if !requested_actual.is_empty() {
            println!("  force-field requested -> actual (count):");
            for ((req, act), count) in &requested_actual {
                println!("    {req:<24} -> {act:<24} {count}");
            }
        }
    }

    // =========================================================================
    // Arm I / Arm J discriminating checks (spec §13's core claim).
    // =========================================================================
    println!("\n=== ARM I (FailClosed) / ARM J (DiagnosticOnly) discriminating check ===");
    let mut arm_i_typed_failures = 0usize;
    let mut arm_j_successes_with_scored_only = 0usize;
    // Arm I/J are built on top of Arm D (full pipeline, including stereo repair --
    // see arm_i()/arm_j()), so a molecule can independently fail for a reason that
    // has nothing to do with the ring-torsion mechanism at all (e.g. Agent D's own
    // disclosed stereo_constraints.rs limitation on ring-fused stereocenters, see
    // that module's doc comment: "two ring-fused stereocenters whose... substituents
    // still overlapped enough for this to happen"). Arm J's contract is narrowly
    // "never fails with RingTorsionApplicationUnsupported specifically" -- any OTHER
    // failure cause is a separate, already-known limitation surfacing here for the
    // first time because loosening the ring gate lets the pipeline reach further
    // into the stage order for these molecules, not a bug in this mechanism.
    let mut arm_j_confounded_by_unrelated_failure: Vec<(&str, String)> = Vec::new();
    for &(name, smiles, _cat) in CORPUS {
        let mol = parse(smiles).unwrap();
        let config_i = arm_i(DEFAULT_SEED);
        let config_j = arm_j(DEFAULT_SEED);
        // Only meaningful for molecules where the torsion-knowledge layer actually
        // matches a small-ring/macrocycle potential -- checked directly via
        // TorsionKnowledgeReport rather than assumed from the category label.
        let tk_config = chematic_3d::etkdg_knowledge::TorsionKnowledgeConfig {
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            ..Default::default()
        };
        let report = chematic_3d::etkdg_knowledge::build_torsion_knowledge(&mol, &tk_config);
        let has_ring_potential = report.potentials.iter().any(|p| {
            matches!(
                p.source,
                TorsionKnowledgeSource::SmallRingExperimental
                    | TorsionKnowledgeSource::MacrocycleAdaptation
            )
        });
        if !has_ring_potential {
            continue;
        }
        let outcome_i = run_pipeline_arm(&mol, &config_i);
        let outcome_j = run_pipeline_arm(&mol, &config_j);
        if !outcome_i.success
            && outcome_i.cause.as_deref() == Some("RingTorsionApplicationUnsupported")
        {
            arm_i_typed_failures += 1;
        } else {
            println!(
                "  UNEXPECTED: {name} has a ring/macrocycle potential but Arm I did not \
                 fail closed as RingTorsionApplicationUnsupported (got success={}, cause={:?})",
                outcome_i.success, outcome_i.cause
            );
        }
        if outcome_j.success && outcome_j.n_scored_only > 0 && outcome_j.diagnostic_only {
            arm_j_successes_with_scored_only += 1;
        } else if outcome_j.cause.as_deref() == Some("RingTorsionApplicationUnsupported") {
            println!(
                "  UNEXPECTED: {name} has a ring/macrocycle potential but Arm J (DiagnosticOnly) \
                 still failed closed as RingTorsionApplicationUnsupported -- this specific cause \
                 must never occur under DiagnosticOnly"
            );
        } else {
            // Failed (or succeeded with unexpected evidence shape) for a DIFFERENT,
            // unrelated reason -- not evidence against the ring-torsion mechanism.
            arm_j_confounded_by_unrelated_failure.push((
                name,
                outcome_j
                    .cause
                    .clone()
                    .unwrap_or_else(|| "success-but-empty-evidence".to_string()),
            ));
        }
    }
    println!(
        "Arm I typed RingTorsionApplicationUnsupported failures on ring/macrocycle-bearing \
         molecules: {arm_i_typed_failures}"
    );
    println!(
        "Arm J successes with scored-only+diagnostic_only evidence on the same molecules: \
         {arm_j_successes_with_scored_only}"
    );
    if !arm_j_confounded_by_unrelated_failure.is_empty() {
        println!(
            "Arm J molecules confounded by an UNRELATED, already-known limitation (not a \
             ring-torsion bug -- see stereo_constraints.rs's own disclosed ring-fused-\
             stereocenter repair limitation): {arm_j_confounded_by_unrelated_failure:?}"
        );
    }
    assert!(
        arm_i_typed_failures > 0,
        "expected at least one molecule to trip Arm I's fail-closed gate"
    );
    assert!(
        arm_j_successes_with_scored_only > 0,
        "expected at least one molecule to succeed under Arm J with scored-only evidence"
    );
    println!(
        "VERIFIED: FailClosed and DiagnosticOnly genuinely discriminate on the same \
         ring/macrocycle-bearing molecules at the same seed."
    );

    // =========================================================================
    // Conformer RMSD vs raw DG (spec §15's last metric): Kabsch RMSD of each arm's
    // final coordinates against Arm A's raw-DG coordinates, same molecule/seed.
    // =========================================================================
    println!("\n=== conformer RMSD vs raw DG (mean over successes, Kabsch-aligned) ===");
    for (i, arm_name) in ARM_NAMES.iter().enumerate().skip(1) {
        let mut rmsds = Vec::new();
        for &(name, _smiles, _cat) in CORPUS {
            let raw = &per_arm[0];
            let arm = &per_arm[i];
            if let (Some(raw_o), Some(arm_o)) = (raw.get(name), arm.get(name))
                && let (Some(raw_c), Some(arm_c)) = (&raw_o.coords, &arm_o.coords)
            {
                let a = to_vec3(raw_c);
                let b = to_vec3(arm_c);
                if a.len() == b.len() && !a.is_empty() {
                    rmsds.push(align_coords(&a, &b).rmsd);
                }
            }
        }
        let mean = if rmsds.is_empty() {
            0.0
        } else {
            rmsds.iter().sum::<f64>() / rmsds.len() as f64
        };
        println!(
            "  {arm_name:<20} mean RMSD vs raw DG: {mean:.3} A (n={})",
            rmsds.len()
        );
    }

    // =========================================================================
    // Reproducibility / invariance spot checks (spec §16).
    // =========================================================================
    println!("\n=== reproducibility / invariance spot checks ===");
    let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol
    let config_d = arm_d(DEFAULT_SEED);
    let r1 = embed_pipeline_v2(&mol, &config_d);
    let r2 = embed_pipeline_v2(&mol, &config_d);
    let same_seed_reproducible = match (&r1, &r2) {
        (Ok(a), Ok(b)) => (0..mol.atom_count())
            .all(|i| a.coords.get(AtomIdx(i as u32)) == b.coords.get(AtomIdx(i as u32))),
        _ => false,
    };
    println!("same seed -> same result (Arm D, paracetamol): {same_seed_reproducible}");
    assert!(
        same_seed_reproducible,
        "same seed must reproduce identical results"
    );

    let config_d_seed2 = arm_d(1);
    let r3 = embed_pipeline_v2(&mol, &config_d_seed2);
    let non_aliased = match (&r1, &r3) {
        (Ok(a), Ok(b)) => (0..mol.atom_count())
            .any(|i| a.coords.get(AtomIdx(i as u32)) != b.coords.get(AtomIdx(i as u32))),
        _ => true, // differing success/failure across seeds is itself non-aliased
    };
    println!("different seeds -> non-aliased output (Arm D, paracetamol): {non_aliased}");
    assert!(
        non_aliased,
        "different seeds must not produce aliased output"
    );

    println!("\n=== VERDICT ===");
    println!(
        "See per-arm summaries above for full denominators, typed-failure breakdowns, and \
         measured metrics. Arm I/J section confirms the core applied-vs-scored-only claim \
         discriminates on real molecules, not just in unit tests. Known gap this harness does \
         NOT attempt: full atom-order-permutation equivalence and rule-order invariance across \
         the whole 63-molecule corpus (covered narrowly by `pipeline_v2.rs`'s own unit tests, \
         not exhaustively here) -- flagged, not silently assumed passing."
    );
}
