//! Issue #291 Phase 0.5: does threading an H-expanded geometry through the full
//! pipeline stage sequence (DG embed -> torsion optimization -> stereo
//! verify/repair -> force-field minimization -> final verify) actually fix the
//! HEAVY-ATOM-ONLY coordinates `embed_pipeline_v2` returns for
//! testosterone/cholesterol -- or does it only fix the pipeline's own internal
//! (H-aware) bookkeeping while the publicly visible geometry stays exactly as
//! before?
//!
//! Motivating concern (raised in review, not yet checked before this file):
//! `repair_tetrahedral_center` always picks the SMALLEST bridge-eligible
//! component to move. A materialized implicit H is a 1-atom terminal
//! component, so it is essentially guaranteed to be the substituent
//! `repair_stereo` moves -- meaning the repair step, by itself, may leave every
//! HEAVY atom exactly where it was. If nothing downstream (torsion
//! optimization, force-field minimization) subsequently moves the heavy atoms
//! in response to the corrected H position, then "the H-expanded molecule's
//! own internal verify says Satisfied" (already established, see this
//! session's PR #380 and the reverted phantom-free-verify investigation) would
//! NOT imply "the coordinates `embed_pipeline_v2` actually returns are fixed" --
//! it would just mean the bookkeeping is right about a part of the structure
//! (the H position) that gets discarded before the caller ever sees it.
//!
//! This harness is NOT wired into `pipeline_v2.rs` or any public API, and
//! changes no default or exposed behavior. It manually re-runs
//! `pipeline_v2::embed_pipeline_v2`'s own stage 4/5/6/7/8/9/10/11 sequence
//! (reusing exactly the same functions that module calls, imported directly)
//! on an `add_hydrogens`-expanded copy of the molecule, and records the metrics
//! needed to answer the question above: heavy-atom RMSD and H-atom RMSD across
//! the repair step and the force-field step, which specific atoms moved,
//! expanded-side vs. original-molecule-plus-truncated-coords stereo verdicts,
//! agreement with an already-fully-explicit-H SMILES run through the
//! completely standard (unmodified) pipeline, and timing/atom-count overhead.
//!
//! Branch A (heavy-only output ends up independently correct): the existing
//! `embed_pipeline_v2` API can be fixed by threading the expanded state through
//! consistently. Branch B (heavy-only stays wrong even though the expanded
//! state is internally "Satisfied"): `embed_pipeline_v2`'s heavy-atom-only
//! contract cannot be fixed this way at all -- would need either a new API
//! that returns the explicit-H molecule+coords, or an approach that actually
//! moves heavy atoms (a real chiral-volume-penalty/local-constraint
//! optimization, not H-position bookkeeping).

use chematic_3d::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2_detail};
use chematic_3d::etkdg_knowledge::{
    TorsionKnowledgeConfig, TorsionOptimizationConfig, build_torsion_knowledge, optimize_torsions,
};
use chematic_3d::minimize::{ForceFieldPolicy, MinimizeConfig, minimize_with_policy_gated};
use chematic_3d::pipeline_v2::{PipelineV2Config, StereoPolicy, embed_pipeline_v2};
use chematic_3d::stereo_constraints::{repair_stereo, verify_stereo};
use chematic_3d::stereo3d::assign_stereo_from_3d;
use chematic_chem::assign_cip;
use chematic_core::AtomIdx;
use chematic_smiles::parse;

/// Mirrors `chematic_3d::minimize::MAX_SANE_BOND_LENGTH` (`pub(crate)`, not
/// reachable from an example binary outside the crate).
const MAX_SANE_BOND_LENGTH: f64 = 3.0;

const CASES: &[(&str, &str)] = &[
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    ),
];

const BASE_SEEDS: [u64; 5] = [0, 1, 2, 3, 4];
const MOVED_EPS: f64 = 1e-9;

fn rmsd_range(
    a: &chematic_3d::coords::Coords3D,
    b: &chematic_3d::coords::Coords3D,
    lo: usize,
    hi: usize,
) -> f64 {
    if hi <= lo {
        return 0.0;
    }
    let sum_sq: f64 = (lo..hi)
        .map(|i| {
            let d = a.get(AtomIdx(i as u32)).distance(&b.get(AtomIdx(i as u32)));
            d * d
        })
        .sum();
    (sum_sq / (hi - lo) as f64).sqrt()
}

fn moved_atoms_in_range(
    a: &chematic_3d::coords::Coords3D,
    b: &chematic_3d::coords::Coords3D,
    lo: usize,
    hi: usize,
) -> Vec<usize> {
    (lo..hi)
        .filter(|&i| a.get(AtomIdx(i as u32)).distance(&b.get(AtomIdx(i as u32))) > MOVED_EPS)
        .collect()
}

fn worst_bond_length(mol: &chematic_core::Molecule, coords: &chematic_3d::coords::Coords3D) -> f64 {
    mol.bonds()
        .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
        .fold(0.0_f64, f64::max)
}

fn truncate(coords: &chematic_3d::coords::Coords3D, n: usize) -> chematic_3d::coords::Coords3D {
    let mut out = chematic_3d::coords::Coords3D::new_zeroed(n);
    for i in 0..n {
        out.set(AtomIdx(i as u32), coords.get(AtomIdx(i as u32)));
    }
    out
}

fn main() {
    for &(name, smiles) in CASES {
        println!("\n================ {name} ================");
        let orig_mol = parse(smiles).unwrap();
        let n_heavy = orig_mol.atom_count();
        let expanded_mol = chematic_chem::add_hydrogens(&orig_mol);
        let n_total = expanded_mol.atom_count();
        let declared = assign_cip(&orig_mol);
        println!(
            "n_heavy={n_heavy} n_total_expanded={n_total} n_declared_stereocenters={}",
            declared.assignments.len()
        );

        // Baseline: fully-explicit-H SMILES (constructed from add_hydrogens) run
        // through the completely standard, unmodified pipeline_v2 -- no sentinel,
        // no repair-of-an-implicit-atom, no truncation. If this converges cleanly,
        // it's an independent ceiling on "what a correct answer looks like" that
        // doesn't depend on anything this harness or PR #380 built.
        {
            let mut baseline_ok = 0usize;
            for &seed in &BASE_SEEDS {
                let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
                config.embed.random_seed = seed;
                config.embed.max_attempts = 8;
                config.embed.enforce_chirality = true;
                config.stereo_policy = StereoPolicy::RepairAndVerify;
                if let Ok(r) = embed_pipeline_v2(&expanded_mol, &config)
                    && r.final_stereo.is_fully_satisfied()
                {
                    baseline_ok += 1;
                }
            }
            println!(
                "[baseline] fully-explicit-H molecule through STANDARD pipeline_v2: {baseline_ok}/5 seeds fully satisfied"
            );
        }

        let mut heavy_fixed_count = 0usize;
        let mut heavy_unchanged_count = 0usize;
        for &seed in &BASE_SEEDS {
            let t_start = std::time::Instant::now();
            // ---- Stage 4: DG embed on the expanded molecule ----
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                ..EmbedParameters::default()
            };
            let coords4 = match embed_distance_geometry_v2_detail(&expanded_mol, &params) {
                Ok((c, _)) => c,
                Err((cause, _)) => {
                    println!("seed={seed}: stage4 DG FAILED cause={cause:?}");
                    continue;
                }
            };
            let stage4_verify = verify_stereo(&expanded_mol, &coords4);

            // ---- Stage 5/6: torsion optimization (matches EmbedParameters::default(),
            // all knowledge flags off -- same no-op state used throughout this
            // investigation and by every other test in this file's own PR chain) ----
            let tk_config = TorsionKnowledgeConfig::default();
            let tk_report = build_torsion_knowledge(&expanded_mol, &tk_config);
            let coords6 = if tk_report.potentials.is_empty() {
                coords4.clone()
            } else {
                match optimize_torsions(
                    &expanded_mol,
                    &coords4,
                    &tk_report.potentials,
                    &TorsionOptimizationConfig::default(),
                ) {
                    Ok((c, _)) => c,
                    Err(e) => {
                        println!("seed={seed}: stage6 torsion opt FAILED {e:?}");
                        continue;
                    }
                }
            };
            let heavy_moved_by_torsion = moved_atoms_in_range(&coords4, &coords6, 0, n_heavy).len();

            // ---- Stage 7: verify (before repair) ----
            let stage7_verify = verify_stereo(&expanded_mol, &coords6);

            // ---- Stage 8: repair ----
            let coords8 = match repair_stereo(&expanded_mol, &coords6) {
                Ok(outcome) => outcome.coords,
                Err(failure) => failure.partial_coords,
            };
            let heavy_rmsd_repair = rmsd_range(&coords6, &coords8, 0, n_heavy);
            let h_rmsd_repair = rmsd_range(&coords6, &coords8, n_heavy, n_total);
            let heavy_moved_by_repair = moved_atoms_in_range(&coords6, &coords8, 0, n_heavy);
            let h_moved_by_repair =
                moved_atoms_in_range(&coords6, &coords8, n_heavy, n_total).len();

            // ---- Stage 9: verify (after repair) ----
            let stage9_verify = verify_stereo(&expanded_mol, &coords8);

            // ---- Stage 10: force-field minimization on the EXPANDED molecule ----
            let ff_config = MinimizeConfig::default();
            let ff_result = minimize_with_policy_gated(
                &expanded_mol,
                coords8.clone(),
                ForceFieldPolicy::Mmff94WithUffFallback,
                &ff_config,
                false,
                false,
            );
            let coords10 = match ff_result {
                Ok(r) => r.coords,
                Err(e) => {
                    println!("seed={seed}: stage10 FF FAILED {e:?}");
                    continue;
                }
            };
            let heavy_rmsd_ff = rmsd_range(&coords8, &coords10, 0, n_heavy);
            let heavy_moved_by_ff = moved_atoms_in_range(&coords8, &coords10, 0, n_heavy).len();

            // ---- Stage 11: final verify on the expanded state, plus the SAME
            // post-minimization repair-and-reverify mechanism pipeline_v2.rs already
            // runs in production for RepairAndVerify (mmff94/uff minimization has no
            // notion of declared stereo and can walk a satisfied geometry back across
            // a boundary -- this exact step is what already rescues that case for the
            // heavy-only path today; not yet tried here for the expanded path) ----
            let mut final_verify_expanded = verify_stereo(&expanded_mol, &coords10);
            let mut coords_final = coords10.clone();
            let mut post_min_repair_used = false;
            if final_verify_expanded.n_violations() > 0
                && let Ok(outcome) = repair_stereo(&expanded_mol, &coords10)
            {
                let reverified = verify_stereo(&expanded_mol, &outcome.coords);
                let sound_after = outcome.coords.is_finite()
                    && worst_bond_length(&expanded_mol, &outcome.coords) <= MAX_SANE_BOND_LENGTH;
                if reverified.n_violations() == 0 && sound_after {
                    coords_final = outcome.coords;
                    final_verify_expanded = reverified;
                    post_min_repair_used = true;
                }
            }
            let coords10 = coords_final;

            // ---- Diagnostic-only follow-up: if the post-min repair fired, does ONE
            // more short FF relaxation pass on the repaired (expanded) geometry
            // restore agreement between the expanded-side verify and the
            // production-equivalent (original-mol + truncated coords) verify? This
            // distinguishes "the repaired H position just needs relaxing, same as
            // the pre-minimization repair-then-refine pattern PR #380 already fixed"
            // from "the repaired geometry is genuinely still wrong." ----
            let mut coords10 = coords10;
            let mut extra_relax_note: Option<String> = None;
            if post_min_repair_used {
                let re_min = minimize_with_policy_gated(
                    &expanded_mol,
                    coords10.clone(),
                    ForceFieldPolicy::Mmff94WithUffFallback,
                    &ff_config,
                    false,
                    false,
                );
                if let Ok(r) = re_min {
                    let re_verify_expanded = verify_stereo(&expanded_mol, &r.coords);
                    let re_truncated = truncate(&r.coords, n_heavy);
                    let re_verify_original = verify_stereo(&orig_mol, &re_truncated);
                    extra_relax_note = Some(format!(
                        "expanded n_violations={} | production-equivalent is_fully_satisfied={} n_violations={}",
                        re_verify_expanded.n_violations(),
                        re_verify_original.is_fully_satisfied(),
                        re_verify_original.n_violations()
                    ));
                    if re_verify_expanded.n_violations() == 0 {
                        coords10 = r.coords;
                    }
                }
            }

            // ---- What production `embed_pipeline_v2` would actually return: the
            // heavy-only truncated coords, verified against the ORIGINAL molecule
            // (current, unmodified `verify_tetrahedral_center`/`phantom_neighbor_position`
            // code path -- exactly what stage 11 sees today in the real pipeline) ----
            let truncated_final = truncate(&coords10, n_heavy);
            let verify_original_truncated = verify_stereo(&orig_mol, &truncated_final);

            // ---- Did the HEAVY-ONLY truncated coordinates actually change end to end,
            // from the raw stage-4 embed through everything, or did repair+optimization
            // only ever move the H? ----
            let truncated_stage4 = truncate(&coords4, n_heavy);
            let heavy_rmsd_end_to_end = rmsd_range(&truncated_stage4, &truncated_final, 0, n_heavy);
            if heavy_rmsd_end_to_end > 1e-6 {
                heavy_fixed_count += 1;
            } else {
                heavy_unchanged_count += 1;
            }

            // ---- Independent oracle cross-check on the final truncated coords
            // (only checks the subset of declared centers with 4 real heavy
            // neighbors -- see this session's PR #380/#291 notes on that limit) ----
            let perceived = assign_stereo_from_3d(&orig_mol, &truncated_final);
            let mut oracle_mismatches = Vec::new();
            for &(idx, code) in &declared.assignments {
                if let Some(p) = perceived.get(idx)
                    && p != code
                {
                    oracle_mismatches.push((idx, code, p));
                }
            }

            let elapsed_ms = t_start.elapsed().as_millis();

            println!(
                "seed={seed} elapsed_ms={elapsed_ms} n_total={n_total} (x{:.2} atoms)",
                n_total as f64 / n_heavy as f64
            );
            println!(
                "  stage4(expanded,raw) n_violations={} | stage6(post-torsion) heavy_moved_by_torsion={heavy_moved_by_torsion} n_violations={}",
                stage4_verify.n_violations(),
                stage7_verify.n_violations()
            );
            println!(
                "  stage8(repair) heavy_rmsd={heavy_rmsd_repair:.6} h_rmsd={h_rmsd_repair:.6} heavy_atoms_moved={heavy_moved_by_repair:?} h_atoms_moved={h_moved_by_repair} | stage9 n_violations={}",
                stage9_verify.n_violations()
            );
            println!(
                "  stage10(force-field) heavy_rmsd={heavy_rmsd_ff:.6} heavy_atoms_moved={heavy_moved_by_ff} | stage11(expanded) n_violations={} post_min_repair_used={post_min_repair_used}",
                final_verify_expanded.n_violations()
            );
            if let Some(note) = &extra_relax_note {
                println!("  [extra-relax-after-post-min-repair] {note}");
            }
            println!(
                "  HEAVY-ONLY end-to-end rmsd (stage4 truncated -> final truncated) = {heavy_rmsd_end_to_end:.6}"
            );
            println!(
                "  production-equivalent verify (original mol + truncated coords): is_fully_satisfied={} n_violations={}",
                verify_original_truncated.is_fully_satisfied(),
                verify_original_truncated.n_violations()
            );
            println!("  oracle_mismatches={oracle_mismatches:?}");
        }
        println!(
            "{name}: heavy-only coords CHANGED end-to-end on {heavy_fixed_count}/5 seeds, UNCHANGED on {heavy_unchanged_count}/5 seeds"
        );
    }
}
