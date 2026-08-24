//! Issue #227 Phase 2 follow-up: precise characterization of the one new
//! stereo violation the BCI fix introduces (`chembl_tier_b_0082`, State 2 ->
//! State 3). Measurement-only -- calls the new, purely additive
//! `chematic_3d::stereo_constraints::{debug_double_bond, debug_all_double_bonds}`
//! diagnostic (production `verify_double_bond`/`verify_stereo` untouched) and
//! `embed_pipeline_v2` with `ForceFieldPolicy::None` (the exact isolation
//! technique `pipeline_v2.rs`'s own module doc already used for the prior
//! `chembl_tier_b_0076`/`chembl_tier_b_0083` cases) to get the
//! pre-minimization ("stereo_before"-equivalent) geometry.
//!
//! Reads State 2's and State 3's already-saved final coordinates directly
//! from the committed dump JSONL (both states' `chembl_tier_b_0082` /
//! `chematic_pipeline_v2_mmff94_strict` rows) -- no re-run of the old,
//! pre-BCI-fix binary needed, since SMILES parsing/atom ordering is
//! unaffected by the BCI fix and is deterministic.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_bci_stereo_drift_diagnostic_227`

use chematic_3d::distance_geometry_v2::EmbedParameters;
use chematic_3d::etkdg_knowledge::TorsionOptimizationConfig;
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2::{
    self as pv2, PipelineV2Config, RingTorsionApplicationPolicy, StereoPolicy,
};
use chematic_3d::stereo_constraints::debug_all_double_bonds;
use chematic_core::AtomIdx;
use serde_json::Value;

const SMILES: &str = "COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC";
const NAME: &str = "chembl_tier_b_0082";
const ARM: &str = "chematic_pipeline_v2_mmff94_strict";
const EMBED_SEED: u64 = 20260801; // matches pipeline_v2_vs_rdkit_dump.rs
const MAX_ATTEMPTS: usize = 8;

fn coords_from_dump(path: &str, name: &str, arm: &str) -> Option<Vec<[f64; 3]>> {
    let text = std::fs::read_to_string(path).ok()?;
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(line).ok()?;
        if row["name"].as_str() == Some(name) && row["arm"].as_str() == Some(arm) {
            let arr = row["coords"].as_array()?;
            return Some(
                arr.iter()
                    .map(|p| {
                        let a = p.as_array().unwrap();
                        [
                            a[0].as_f64().unwrap(),
                            a[1].as_f64().unwrap(),
                            a[2].as_f64().unwrap(),
                        ]
                    })
                    .collect(),
            );
        }
    }
    None
}

fn to_coords3d(raw: &[[f64; 3]]) -> chematic_3d::coords::Coords3D {
    let mut c = chematic_3d::coords::Coords3D::new_zeroed(raw.len());
    for (i, p) in raw.iter().enumerate() {
        c.set(
            AtomIdx(i as u32),
            chematic_3d::coords::Point3::new(p[0], p[1], p[2]),
        );
    }
    c
}

fn report(label: &str, mol: &chematic_core::Molecule, coords: &chematic_3d::coords::Coords3D) {
    println!("--- {label} ---");
    for info in debug_all_double_bonds(mol, coords) {
        println!(
            "  bond={:?} end1={} end2={} sub1={} sub2={} declared_same_side={} \
             actual_angle_deg={:.3} actual_same_side={} margin_from_90deg={:.3} status={:?}",
            info.bond,
            info.end1.0,
            info.end2.0,
            info.sub1.0,
            info.sub2.0,
            info.declared_same_side,
            info.actual_angle_deg,
            info.actual_same_side,
            info.margin_from_boundary_deg,
            info.status
        );
    }
}

fn main() {
    let mol = chematic_smiles::parse(SMILES).expect("chembl_tier_b_0082 must parse");

    // Pre-minimization geometry: ForceFieldPolicy::None isolates embedding
    // (+ stage 6 torsion optimization, no MMFF94 at all), same technique
    // pipeline_v2.rs's own module doc used for chembl_tier_b_0076/0083.
    let config = PipelineV2Config {
        embed: EmbedParameters {
            random_seed: EMBED_SEED,
            max_attempts: MAX_ATTEMPTS,
            use_exp_torsions: true,
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            use_macrocycle_14_bounds: true,
            track_failures: true,
            ..EmbedParameters::default()
        },
        torsion_optimization: TorsionOptimizationConfig::default(),
        include_legacy_torsion_heuristic: false,
        stereo_policy: StereoPolicy::Ignore,
        fail_on_unevaluable_stereo: false,
        force_field_policy: ForceFieldPolicy::None,
        force_field_max_iterations: 200,
        gate_mmff94_torsion_oop: false,
        gate_mmff94_stretch_bend: false,
        ring_torsion_policy: RingTorsionApplicationPolicy::DiagnosticOnly,
        total_timeout_ms: Some(20_000),
        expand_implicit_h_through_pipeline: false,
    };
    let pre_min = pv2::embed_pipeline_v2(&mol, &config).expect("embedding must succeed");
    report(
        "pre-minimization (ForceFieldPolicy::None, embedding-only)",
        &mol,
        &pre_min.coords,
    );

    let state2_raw = coords_from_dump(
        "validation/results/mmff94_bci_gap_227_state2_post_torsion_fix_chematic_rows.jsonl",
        NAME,
        ARM,
    )
    .expect("state2 coords must be present");
    report(
        "State 2 final (post-torsion-fix, pre-BCI-fix)",
        &mol,
        &to_coords3d(&state2_raw),
    );

    let state3_raw = coords_from_dump(
        "validation/results/mmff94_bci_gap_227_state3_post_bci_fix_chematic_rows.jsonl",
        NAME,
        ARM,
    )
    .expect("state3 coords must be present");
    let state3_coords = to_coords3d(&state3_raw);
    report(
        "State 3 final (post-BCI-fix, this PR)",
        &mol,
        &state3_coords,
    );

    // Tier-3 fix feasibility check: does repair_stereo actually succeed on
    // this exact violated post-minimization geometry?
    println!("--- repair_stereo(mol, State3 final coords) ---");
    println!("  pre-repair:");
    soundness_check(&mol, &state3_coords);
    match chematic_3d::stereo_constraints::repair_stereo(&mol, &state3_coords) {
        Ok(outcome) => {
            println!("  repair_stereo: Ok, repaired={:?}", outcome.repaired);
            println!("  post-repair:");
            soundness_check(&mol, &outcome.coords);
            report("State 3 final, AFTER repair_stereo", &mol, &outcome.coords);
        }
        Err(e) => {
            println!("  repair_stereo: Err, failures={:?}", e.failures);
        }
    }
}

fn soundness_check(mol: &chematic_core::Molecule, coords: &chematic_3d::coords::Coords3D) {
    let n = mol.atom_count();
    let mut worst_ratio = 0.0f64;
    for (_bidx, bond) in mol.bonds() {
        let a = bond.atom1;
        let b = bond.atom2;
        let actual = coords.get(a).distance(&coords.get(b));
        let ideal = (mol.atom(a).element.covalent_radius() as f64
            + mol.atom(b).element.covalent_radius() as f64)
            .max(0.3);
        let rel_error = (actual / ideal - 1.0).abs();
        worst_ratio = worst_ratio.max(rel_error);
    }
    let mut clashes = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);
            if mol.bond_between(a, b).is_some() {
                continue;
            }
            let d = coords.get(a).distance(&coords.get(b));
            if d < 1.2 {
                clashes += 1;
            }
        }
    }
    println!(
        "  soundness: worst_bond_length_ratio={worst_ratio:.4} (sane<=3.0) gross_clash_count={clashes}"
    );
}
