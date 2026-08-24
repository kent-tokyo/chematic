//! Issue #291 Step A: does `StereoPolicy::RepairAndVerify` already close the gap?
//!
//! Issue #291 measured that production `embed_pipeline_v2` under
//! `ForceFieldPolicy::UffOnly` + `StereoPolicy::Ignore` silently ships wrong
//! tetrahedral/E-Z stereochemistry on 18/29 (62%) of the declared-stereocenter
//! molecules in the 58-molecule corpus (`scripts/etkdg_vs_rdkit_gap.py::CORPUS`,
//! same 29-molecule subset `stereo_constraints_gap_check.rs`'s `stereo_subset`
//! filters to -- copied verbatim below). That measurement used
//! `Result::is_ok()` as its correctness signal, which under `Ignore` is
//! meaningless for stereo (violations never gate success).
//!
//! `pipeline_v2_vs_rdkit_dump.rs`'s existing `chematic_pipeline_v2_mmff94_strict_repair`
//! arm (`StereoPolicy::RepairAndVerify`, `enforce_chirality: false` -- already
//! shipped, no code change) already shows 100% stereo satisfaction on the
//! 265-molecule benchmark corpus, but only among rows that *succeed* -- and
//! that arm uses `Mmff94BondAngleStrict`, never `UffOnly`, and a different
//! corpus than #291's own 29-molecule population. This example fills that
//! specific gap, with 3 arms, against #291's *exact* population, across
//! multiple base seeds (#291's own caveat: n=29 single-seed is "one draw, not
//! a fully converged rate"), reading `final_stereo.n_violations()` directly
//! rather than `Result::is_ok()`:
//!
//! - **`ignore`**: `UffOnly` + `StereoPolicy::Ignore` -- reproduces #291's
//!   original measurement, to confirm the number is still current.
//! - **`repair_and_verify`**: `UffOnly` + `StereoPolicy::RepairAndVerify`,
//!   `enforce_chirality: false` -- zero code change, already shipped.
//! - **`enforce_chirality_repair_and_verify`**: same, plus
//!   `enforce_chirality: true` -- requires lifting `pipeline_v2.rs`'s
//!   `enforce_chirality && RepairAndVerify` guard (`InvalidConfiguration`),
//!   which existed because "composing the two repair mechanisms is a
//!   separate, not-yet-validated question (deliberately deferred, not
//!   decided by omission)" -- this arm is exactly that validation.
//!
//! **Result (2026-08-24, this corpus, 5 seeds)**: `ignore` reproduces #291
//! (58.6% silently wrong here vs. 62.1% in #291's single-seed number --
//! consistent). `repair_and_verify` alone already converts this to 0% silently
//! wrong, 86.2% correct *and* successful, 13.8% loud (honest) failure -- with
//! zero code changes. Lifting the guard for `enforce_chirality_repair_and_verify`
//! and re-running just the residual failures fully recovers naproxen_S,
//! ibuprofen_S, and penicillin_core (15/15 seed-runs each now correct); it
//! does NOT recover testosterone or cholesterol (5/5 seed-runs each still
//! fail, now at the embedding stage itself -- `EmbedFailureCause::
//! StereoConstraintFailed` -- rather than at post-repair verification). Those
//! two are exactly the ring-fused-stereocenter case `docs/rfcs/
//! etkdg_3d_gap_rfc.md` already diagnosed as unfixable by substituent-
//! reflection repair (no non-ring bridge-eligible substituent exists to
//! reflect) -- confirms that population needs the separately-scoped
//! chiral-volume-penalty-in-`refine_coords` work, not this fix. The two
//! remaining failures outside that pair (menthol, atorvastatin_fragment) are
//! an unrelated, pre-existing UFF `CatastrophicBondBlowup` issue, not stereo.
//!
//! Per-molecule, per-seed outcome is one of:
//! - `silently_wrong`: Ok, but `final_stereo.n_violations() > 0` (only possible
//!   under `Ignore` -- included here to reproduce #291's number before anything
//!   else, as a sanity check that it's still current).
//! - `correct_and_ok`: Ok, zero violations.
//! - `loud_failure_stereo`: Err, and the failure's own `final_stereo` (or, if
//!   unavailable, `stereo_after_repair`) shows a violation -- i.e. `RepairAndVerify`
//!   correctly refused to silently ship the wrong structure.
//! - `loud_failure_other`: Err for an unrelated reason (embedding never
//!   converged, timeout, etc.) -- not a stereo-policy outcome at all, kept
//!   separate so it isn't miscounted either way.
//!
//! Also reports, for every `correct_and_ok` row: whether the returned geometry
//! is still bond-length-sound (this project's established >50%-covalent-radius
//! blowup convention, matching `stereo_constraints_gap_check.rs` and
//! `distance_geometry_v2_gap_check.rs`) -- a stereo fix that leaves a strained
//! structure is not a real success.
//!
//! **4th arm added (issue #291 real implementation, following the "Phase 0.5"
//! feasibility measurement in `crates/chematic-3d/examples/
//! issue291_expanded_geometry_feasibility.rs`)**: `expand_implicit_h_repair_and_verify`
//! -- same as `enforce_chirality_repair_and_verify` plus
//! `PipelineV2Config::expand_implicit_h_through_pipeline: true`, threading an
//! `add_hydrogens`-expanded molecule through the whole pipeline instead of
//! truncating right after embed. Targets exactly testosterone/cholesterol,
//! the two molecules the prior arms could never fix (ring-fused declared
//! stereocenters with no non-ring substituent for `repair_stereo` to
//! reflect).
//!
//! **Result (2026-08-24, this corpus, 5 seeds)**: 144/145 (99.3%) correct_and_ok,
//! 0 silently_wrong, 0 loud_failure_stereo -- testosterone and cholesterol both
//! now succeed on every declared stereocenter, every seed (5/5 each), fully
//! closing the residual `enforce_chirality_repair_and_verify` could not. The
//! single remaining failure (cholesterol, one seed) is `loud_failure_other`
//! (`ForceField(MinimizationFailed(CatastrophicBondBlowup))` under `UffOnly`) --
//! the same pre-existing, unrelated UFF numerical issue this file's docs
//! already named for menthol/atorvastatin_fragment, not a stereo failure.
//!
//! Pure external caller via the existing public API, same convention as
//! `pipeline_v2_vs_rdkit_dump.rs` -- the only production change this needed
//! was lifting `pipeline_v2.rs`'s config-validation guard (see above), which
//! this file's own result justifies keeping.
//!
//! Run: `cargo run --release -p chematic-3d --example issue291_repair_policy_measurement`

use std::collections::BTreeMap;

use chematic_3d::distance_geometry_v2::EmbedParameters;
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2::{PipelineV2Config, StereoPolicy, embed_pipeline_v2};
use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};
use chematic_smiles::parse;

// Copied verbatim from `stereo_constraints_gap_check.rs`'s `CORPUS`, filtered to
// exactly the 29 declared-stereo molecules -- the same population issue #291
// measured (10 stereocenter_implicit_h, 6 stereocenter_quaternary, 8 alkene_ez,
// 5 stereo-bearing druglike*).
const STEREO_CORPUS: &[(&str, &str)] = &[
    ("l_alanine", "N[C@@H](C)C(=O)O"),
    ("d_alanine", "N[C@H](C)C(=O)O"),
    ("l_serine", "N[C@@H](CO)C(=O)O"),
    ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O"),
    ("2_butanol_R", "C[C@H](O)CC"),
    ("2_butanol_S", "C[C@@H](O)CC"),
    ("2_chlorobutane_R", "C[C@H](Cl)CC"),
    ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
    ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
    ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O"),
    ("chfclbr_R", "[C@H](F)(Cl)Br"),
    ("chfclbr_S", "[C@@H](F)(Cl)Br"),
    ("quaternary_1_R", "[C@](F)(Cl)(Br)I"),
    ("quaternary_1_S", "[C@@](F)(Cl)(Br)I"),
    ("quaternary_2_R", "[C@](C)(N)(O)F"),
    ("quaternary_2_S", "[C@@](C)(N)(O)F"),
    ("but2ene_E", "C/C=C/C"),
    ("but2ene_Z", r"C/C=C\C"),
    ("chloropropene_E", "C(/C=C/C)Cl"),
    ("chloropropene_Z", r"C(/C=C\C)Cl"),
    ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1"),
    ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
    ("pent2ene_E", "CC/C=C/C"),
    ("pent2ene_Z", r"CC/C=C\C"),
    (
        "penicillin_core",
        "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
    ),
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    ),
    (
        "atorvastatin_fragment",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    ),
    ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O"),
];

// 5 independent base seeds -- addresses #291's own "one draw, not a fully
// converged rate" caveat. Arbitrary, fixed for reproducibility.
const BASE_SEEDS: &[u64] = &[
    0xC0FF_EE42_D157_6E02,
    0x1234_5678_9ABC_DEF0,
    0x5EED_0001_0000_0001,
    0xA5A5_5A5A_1111_2222,
    0x9E37_79B9_7F4A_7C15,
];

const RDKIT_COVALENT_RADIUS: &[(u8, f64)] = &[
    (1, 0.31),
    (6, 0.76),
    (7, 0.71),
    (8, 0.66),
    (9, 0.57),
    (15, 1.07),
    (16, 1.05),
    (17, 1.02),
    (35, 1.20),
    (53, 1.39),
];

fn covalent_radius(z: u8) -> Option<f64> {
    RDKIT_COVALENT_RADIUS
        .iter()
        .find(|(n, _)| *n == z)
        .map(|(_, r)| *r)
}

fn bond_order_scale(order: BondOrder) -> f64 {
    match order {
        BondOrder::Double => 0.87,
        BondOrder::Triple => 0.78,
        BondOrder::Aromatic => 0.93,
        _ => 1.00,
    }
}

const BOND_BLOWUP_REL_ERROR: f64 = 0.5;

/// Max relative bond-length error across every bond with a known reference,
/// or `None` if the coords are missing/non-finite (caller treats as unsound).
fn max_bond_rel_error(mol: &Molecule, coords: &chematic_3d::coords::Coords3D) -> Option<f64> {
    let mut max_rel = 0.0_f64;
    for (_, bond) in mol.bonds() {
        let za = mol.atom(bond.atom1).element.atomic_number();
        let zb = mol.atom(bond.atom2).element.atomic_number();
        let (Some(ra), Some(rb)) = (covalent_radius(za), covalent_radius(zb)) else {
            continue;
        };
        let r0 = (ra + rb) * bond_order_scale(bond.order);
        let p1 = coords.get(bond.atom1);
        let p2 = coords.get(bond.atom2);
        if !p1.x.is_finite() || !p2.x.is_finite() {
            return None;
        }
        let r = p1.distance(&p2);
        let rel = (r - r0).abs() / r0;
        if rel > max_rel {
            max_rel = rel;
        }
    }
    Some(max_rel)
}

fn mol_has_declared_stereo(mol: &Molecule) -> bool {
    if (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None) {
        return true;
    }
    mol.bonds()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Up | BondOrder::Down))
}

fn base_config(
    stereo_policy: StereoPolicy,
    seed: u64,
    enforce_chirality: bool,
    expand_implicit_h_through_pipeline: bool,
) -> PipelineV2Config {
    PipelineV2Config {
        embed: EmbedParameters {
            random_seed: seed,
            // Same budget for every arm, including the new expanded-geometry one
            // -- matches what the Phase 0.5 measurement harness actually used
            // (8) and validated testosterone/cholesterol against; not bumped
            // speculatively.
            max_attempts: 8,
            enforce_chirality,
            ..EmbedParameters::default()
        },
        stereo_policy,
        expand_implicit_h_through_pipeline,
        ..PipelineV2Config::minimal(ForceFieldPolicy::UffOnly)
    }
}

#[derive(Default, Clone, Copy)]
struct Counts {
    silently_wrong: usize,
    correct_and_ok: usize,
    correct_and_ok_unsound_geometry: usize,
    loud_failure_stereo: usize,
    loud_failure_other: usize,
}

const ARMS: &[(&str, StereoPolicy, bool, bool)] = &[
    ("ignore", StereoPolicy::Ignore, false, false),
    (
        "repair_and_verify",
        StereoPolicy::RepairAndVerify,
        false,
        false,
    ),
    (
        "enforce_chirality_repair_and_verify",
        StereoPolicy::RepairAndVerify,
        true,
        false,
    ),
    // New arm (issue #291 real implementation): threads an add_hydrogens-expanded
    // molecule through the whole pipeline instead of just the embed. See
    // ROADMAP.md's #291 entry ("Phase 0.5") for the measurement this is based on.
    (
        "expand_implicit_h_repair_and_verify",
        StereoPolicy::RepairAndVerify,
        true,
        true,
    ),
];

fn main() {
    let mut counts: BTreeMap<&'static str, Counts> = BTreeMap::new();
    for &(arm_key, _, _, _) in ARMS {
        counts.insert(arm_key, Counts::default());
    }

    println!("molecule\tarm\tseed\toutcome\tviolations\tfailure_cause\tmax_bond_rel_error");

    for (name, smiles) in STEREO_CORPUS {
        let mol = parse(smiles).expect("corpus SMILES must parse");
        assert!(
            mol_has_declared_stereo(&mol),
            "{name}: expected to be in the declared-stereo subset"
        );

        for &(arm_key, policy, enforce_chirality, expand_implicit_h) in ARMS {
            for &seed in BASE_SEEDS {
                let config = base_config(policy, seed, enforce_chirality, expand_implicit_h);
                let acc = counts.get_mut(arm_key).unwrap();
                match embed_pipeline_v2(&mol, &config) {
                    Ok(result) => {
                        let violations = result.final_stereo.n_violations();
                        if violations == 0 {
                            let sound = max_bond_rel_error(&mol, &result.coords)
                                .map(|e| e <= BOND_BLOWUP_REL_ERROR)
                                .unwrap_or(false);
                            if sound {
                                acc.correct_and_ok += 1;
                            } else {
                                acc.correct_and_ok_unsound_geometry += 1;
                            }
                            println!(
                                "{name}\t{arm_key}\t{seed:#x}\tcorrect_and_ok\t{violations}\t-\t{:.4}",
                                max_bond_rel_error(&mol, &result.coords).unwrap_or(f64::NAN)
                            );
                        } else {
                            acc.silently_wrong += 1;
                            println!(
                                "{name}\t{arm_key}\t{seed:#x}\tsilently_wrong\t{violations}\t-\t-"
                            );
                        }
                    }
                    Err(failure) => {
                        let violations = failure
                            .final_stereo
                            .as_ref()
                            .or(failure.stereo_after_repair.as_ref())
                            .map(|s| s.n_violations());
                        let is_stereo_cause = matches!(
                            failure.cause,
                            chematic_3d::pipeline_v2::PipelineV2FailureCause::FinalStereoViolation
                                | chematic_3d::pipeline_v2::PipelineV2FailureCause::StereoRepairFailed
                                | chematic_3d::pipeline_v2::PipelineV2FailureCause::DistanceGeometry(
                                    chematic_3d::distance_geometry_v2::EmbedFailureCause::StereoConstraintFailed
                                )
                        );
                        if is_stereo_cause {
                            acc.loud_failure_stereo += 1;
                        } else {
                            acc.loud_failure_other += 1;
                        }
                        println!(
                            "{name}\t{arm_key}\t{seed:#x}\t{}\t{}\t{:?}\t-",
                            if is_stereo_cause {
                                "loud_failure_stereo"
                            } else {
                                "loud_failure_other"
                            },
                            violations.map(|v| v.to_string()).unwrap_or("?".into()),
                            failure.cause
                        );
                    }
                }
            }
        }
    }

    eprintln!("\n=== Summary (29 molecules x 5 seeds = 145 runs per arm) ===");
    for (arm, c) in &counts {
        eprintln!(
            "{arm}: silently_wrong={} correct_and_ok={} correct_and_ok_unsound_geometry={} loud_failure_stereo={} loud_failure_other={}",
            c.silently_wrong,
            c.correct_and_ok,
            c.correct_and_ok_unsound_geometry,
            c.loud_failure_stereo,
            c.loud_failure_other
        );
    }
}
