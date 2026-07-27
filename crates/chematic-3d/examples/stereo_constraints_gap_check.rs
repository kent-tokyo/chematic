//! Acceptance-gate measurement for `stereo_constraints::{verify_stereo, repair_stereo}`
//! (3D Breakthrough Program, Wave 2, Agent D). See `docs/3d_breakthrough_master_plan.md`
//! §0/§3 and the module doc comment in `crates/chematic-3d/src/stereo_constraints.rs`.
//!
//! Mirrors the pattern Agent C's `distance_geometry_v2_gap_check.rs` established: a
//! standalone example binary against the frozen 58-molecule corpus (copied verbatim
//! from that file, which itself copies `scripts/etkdg_vs_rdkit_gap.py::CORPUS`), with
//! a differential, multi-arm gate rather than one absolute pass/fail number.
//!
//! # 3-arm gate
//!
//! - **Arm (a) -- baseline**: embed via `embed_distance_geometry_v2` with
//!   `enforce_chirality: false` (today's default -- the embedder never looks at
//!   declared stereo), then run `verify_stereo`. This must show real, nonzero
//!   violations -- if it doesn't, the verifier isn't discriminating anything (the
//!   exact failure mode the master plan calls out from Wave 1).
//! - **Arm (b) -- verifier-only correctness, with controls**: hand-built
//!   known-correct/known-wrong tetrahedral (quaternary AND implicit-H) and E/Z
//!   fixtures (positive AND negative controls, this project's established
//!   convention). An independent cross-check against this crate's own
//!   `assign_stereo_from_3d` + `chematic_chem::assign_cip` was attempted first and
//!   deliberately dropped -- see "Discovered (not fixed): a pre-existing
//!   `assign_stereo_from_3d`/`assign_cip` inconsistency" below for why that oracle
//!   turned out to be unusable as a control, independent of this module.
//! - **Arm (c) -- full enforce**: same embed, then `repair_stereo`, then re-verify.
//!   Reports the post-repair violation rate, with genuinely-unrepairable cases (typed
//!   `RepairRejectionReason`) counted and named separately, never folded into
//!   "success". Also reports geometry disturbance: atoms moved, max displacement, and
//!   whether repair introduces any NEW bond-length blowup or clash elsewhere in the
//!   molecule (reusing this codebase's >50%-covalent-radius-blowup and 0.5 Å
//!   gross-clash conventions from Agent C's harness).
//!
//! Tetrahedral and E/Z are reported separately throughout, never pooled into one
//! number.
//!
//! No RDKit/Python dependency: pure Rust, runs in any `cargo test`/CI environment.
//!
//! # Discovered (not fixed): a pre-existing `assign_stereo_from_3d`/`assign_cip`
//! inconsistency
//!
//! While building arm (b)'s cross-check, a hand-built geometry for
//! `[C@](F)(Cl)(Br)I` was constructed and independently confirmed (two ways: this
//! module's own sign-convention tests, and a by-hand trace of
//! `chematic-chem/src/cip.rs`'s `assign_tetrahedral` parity-counting algorithm) to
//! correctly realize the declared `@` configuration. On that geometry,
//! `chematic_chem::assign_cip` declares `S` but `chematic_3d::stereo3d::
//! assign_stereo_from_3d` perceives `R` -- a disagreement, not a tie or a degenerate
//! case. The mirrored check (`[C@@](F)(Cl)(Br)I` with the mirrored geometry) shows
//! the same disagreement in the opposite direction (declared `R`, perceived `S`),
//! ruling out "one geometry built wrong" as the explanation. This is only observable
//! where `assign_stereo_from_3d` produces an answer at all, i.e. 4-heavy-neighbor
//! (no implicit H) centers -- confirmed separately: `assign_stereo_from_3d` returns
//! `None` for a 3-heavy-neighbor (implicit-H) center, consistent with its documented
//! `nb_count == 4` gate. This is the same disagreement pattern this file's own
//! cross-check (removed) found on 7/7 corpus 4-heavy-neighbor centers.
//!
//! Which of `assign_cip` / `assign_stereo_from_3d` (or the `cip_priority`/`rank4`
//! ranking either depends on) is actually wrong is **not determined here** -- this is
//! reported as a discovered inconsistency in already-shipped code, not attributed to
//! a specific root cause, and **not fixed by this PR**: neither `stereo3d.rs` nor
//! `chematic-perception::cip_priority` nor `chematic-chem::cip` is a file this PR
//! owns or touches. `stereo_constraints::verify_stereo` does not call either
//! function internally (see that module's doc comment on why it avoids CIP ranking
//! entirely), so this inconsistency does not affect this module's own correctness --
//! but it does mean an oracle built on top of those two functions is not a reliable
//! arm-(b) control, so this file uses only self-contained hand-built fixtures
//! instead.
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example stereo_constraints_gap_check
//! ```

use std::collections::BTreeMap;

use chematic_3d::coords::Coords3D;
use chematic_3d::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};
use chematic_3d::stereo_constraints::{
    RepairRejectionReason, StereoRejectionReason, StereoStatus, repair_stereo, verify_stereo,
};
use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};
use chematic_smiles::parse;

// ---------------------------------------------------------------------------
// Frozen 58-molecule corpus -- verbatim transcription of
// crates/chematic-3d/examples/distance_geometry_v2_gap_check.rs::CORPUS, which is
// itself a verbatim transcription of scripts/etkdg_vs_rdkit_gap.py::CORPUS.
// ---------------------------------------------------------------------------

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
];

// ---------------------------------------------------------------------------
// External geometry reference (same as distance_geometry_v2_gap_check.rs -- RDKit
// GetPeriodicTable().GetRcovalent(), not chematic's own table).
// ---------------------------------------------------------------------------

fn rdkit_covalent_radius(atomic_number: u8) -> Option<f64> {
    match atomic_number {
        1 => Some(0.31),
        6 => Some(0.76),
        7 => Some(0.71),
        8 => Some(0.66),
        9 => Some(0.57),
        15 => Some(1.07),
        16 => Some(1.05),
        17 => Some(1.02),
        35 => Some(1.20),
        53 => Some(1.39),
        _ => None,
    }
}

fn bond_order_scale(order: BondOrder) -> f64 {
    match order {
        BondOrder::Double => 0.87,
        BondOrder::Triple => 0.78,
        BondOrder::Aromatic => 0.93,
        _ => 1.00,
    }
}

fn ref_bond_length(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> Option<f64> {
    let za = mol.atom(a).element.atomic_number();
    let zb = mol.atom(b).element.atomic_number();
    let ra = rdkit_covalent_radius(za)?;
    let rb = rdkit_covalent_radius(zb)?;
    let order = mol
        .bond_between(a, b)
        .map(|(_, bond)| bond.order)
        .unwrap_or(BondOrder::Single);
    Some((ra + rb) * bond_order_scale(order))
}

const BOND_BLOWUP_REL_ERROR: f64 = 0.5; // matches distance_geometry_v2_gap_check.rs
const GROSS_CLASH_DIST: f64 = 0.5;

/// Max relative bond-length error across every bond with a known external reference.
fn max_bond_rel_error(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut max_rel = 0.0_f64;
    for (_, bond) in mol.bonds() {
        let Some(r0) = ref_bond_length(mol, bond.atom1, bond.atom2) else {
            continue;
        };
        let r = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        let rel = (r - r0).abs() / r0;
        if rel > max_rel {
            max_rel = rel;
        }
    }
    max_rel
}

fn has_gross_clash(coords: &Coords3D) -> bool {
    let n = coords.atom_count();
    for i in 0..n {
        for j in (i + 1)..n {
            if coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)))
                < GROSS_CLASH_DIST
            {
                return true;
            }
        }
    }
    false
}

fn mol_has_declared_stereo(mol: &Molecule) -> bool {
    if (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None) {
        return true;
    }
    mol.bonds()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Up | BondOrder::Down))
}

fn rate(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        "n/a".to_string()
    } else {
        format!("{:.4}", numerator as f64 / denominator as f64)
    }
}

#[derive(Default)]
struct Accum {
    declared: usize,
    satisfied: usize,
    violated: usize,
    unevaluable: usize,
    unevaluable_reasons: BTreeMap<&'static str, usize>,
}

fn reason_name(r: StereoRejectionReason) -> &'static str {
    match r {
        StereoRejectionReason::UnsupportedCoordination => "UnsupportedCoordination",
        StereoRejectionReason::DegenerateGeometry => "DegenerateGeometry",
        StereoRejectionReason::DegenerateImplicitHDirection => "DegenerateImplicitHDirection",
        StereoRejectionReason::TerminalOrCumulatedAlkene => "TerminalOrCumulatedAlkene",
        StereoRejectionReason::AmbiguousDirection => "AmbiguousDirection",
    }
}

fn repair_reason_name(r: RepairRejectionReason) -> &'static str {
    match r {
        RepairRejectionReason::NoBridgeEligibleSubstituent => "NoBridgeEligibleSubstituent",
        RepairRejectionReason::DegenerateReflectionPlane => "DegenerateReflectionPlane",
        RepairRejectionReason::DegenerateReflectionDirection => "DegenerateReflectionDirection",
        RepairRejectionReason::StillViolatedAfterRepair => "StillViolatedAfterRepair",
        RepairRejectionReason::RepairCausedNewViolation => "RepairCausedNewViolation",
    }
}

fn main() {
    let stereo_subset: Vec<&(&str, &str, &str)> = CORPUS
        .iter()
        .filter(|(_, smiles, _)| mol_has_declared_stereo(&parse(smiles).unwrap()))
        .collect();

    println!("=== Corpus ===");
    println!(
        "{} / {} molecules have declared stereo (@/@@ or /\\\\).",
        stereo_subset.len(),
        CORPUS.len()
    );

    // =========================================================================
    // ARM (a): BASELINE -- enforce_chirality: false (today's default), verify only.
    // =========================================================================
    let mut tet_a = Accum::default();
    let mut ez_a = Accum::default();
    let mut n_embed_failed = 0usize;

    // Also carries data forward into arm (c) so embedding is only done once per
    // molecule.
    let mut embedded: Vec<(&str, Molecule, Coords3D)> = Vec::new();

    for &&(name, smiles, _category) in &stereo_subset {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("corpus SMILES failed ({name}): {e:?}"));
        let params = EmbedParameters::default(); // enforce_chirality: false
        match embed_distance_geometry_v2(&mol, &params) {
            Ok(coords) => {
                let report = verify_stereo(&mol, &coords);
                for t in &report.tetrahedral {
                    match t.status {
                        StereoStatus::Satisfied => tet_a.satisfied += 1,
                        StereoStatus::Violated => tet_a.violated += 1,
                        StereoStatus::Unevaluable(r) => {
                            tet_a.unevaluable += 1;
                            *tet_a.unevaluable_reasons.entry(reason_name(r)).or_insert(0) += 1;
                        }
                    }
                    tet_a.declared += 1;
                }
                for d in &report.double_bond {
                    match d.status {
                        StereoStatus::Satisfied => ez_a.satisfied += 1,
                        StereoStatus::Violated => ez_a.violated += 1,
                        StereoStatus::Unevaluable(r) => {
                            ez_a.unevaluable += 1;
                            *ez_a.unevaluable_reasons.entry(reason_name(r)).or_insert(0) += 1;
                        }
                    }
                    ez_a.declared += 1;
                }
                embedded.push((name, mol, coords));
            }
            Err(cause) => {
                n_embed_failed += 1;
                eprintln!("embed failed for {name}: {cause:?} (excluded from stereo arms)");
            }
        }
    }

    println!("\n=== ARM (a): BASELINE (enforce_chirality: false, verify only) ===");
    println!("embed failures (excluded from all arms below): {n_embed_failed}");
    println!(
        "tetrahedral: declared={} satisfied={} violated={} unevaluable={} | violation_rate(among evaluable)={}",
        tet_a.declared,
        tet_a.satisfied,
        tet_a.violated,
        tet_a.unevaluable,
        rate(tet_a.violated, tet_a.satisfied + tet_a.violated)
    );
    println!(
        "  tetrahedral unevaluable reasons: {:?}",
        tet_a.unevaluable_reasons
    );
    println!(
        "E/Z:         declared={} satisfied={} violated={} unevaluable={} | violation_rate(among evaluable)={}",
        ez_a.declared,
        ez_a.satisfied,
        ez_a.violated,
        ez_a.unevaluable,
        rate(ez_a.violated, ez_a.satisfied + ez_a.violated)
    );
    println!("  E/Z unevaluable reasons: {:?}", ez_a.unevaluable_reasons);
    if tet_a.violated == 0 && ez_a.violated == 0 {
        eprintln!(
            "\nWARNING: arm (a) shows ZERO violations. Since the embedder is stereo-blind \
             (enforce_chirality: false), this would mean the verifier is not discriminating \
             anything -- the exact Wave-1 gate-design failure this program's own docs warn \
             against. Treat this as a red flag on the verifier, not evidence the embedder is \
             secretly stereo-aware."
        );
    } else {
        println!(
            "(nonzero violations confirm the baseline embedder is genuinely stereo-blind and \
             the verifier discriminates real cases -- not a dead gate.)"
        );
    }

    // =========================================================================
    // ARM (b): VERIFIER-ONLY CORRECTNESS -- positive/negative controls.
    //
    // An independent cross-check against assign_stereo_from_3d + assign_cip was
    // attempted here and removed -- see this file's module doc comment ("Discovered
    // (not fixed): a pre-existing assign_stereo_from_3d/assign_cip inconsistency")
    // for why that oracle pair disagrees with itself on a geometry independently
    // confirmed correct two other ways, making it unusable as a control.
    // =========================================================================
    println!("\n=== ARM (b): VERIFIER-ONLY CORRECTNESS (controls) ===");
    run_fixture_controls();

    // =========================================================================
    // ARM (c): FULL ENFORCE -- embed, repair, re-verify.
    // =========================================================================
    println!("\n=== ARM (c): FULL ENFORCE (embed -> repair_stereo -> re-verify) ===");
    let mut tet_c = Accum::default();
    let mut ez_c = Accum::default();
    let mut n_repair_ok = 0usize;
    let mut n_repair_failed = 0usize;
    let mut repair_failure_reasons: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut total_atoms_moved = 0usize;
    let mut max_displacement_seen = 0.0_f64;
    let mut n_new_blowup_introduced = 0usize;
    let mut n_new_clash_introduced = 0usize;
    let mut n_molecules_with_violations = 0usize;
    // Reconciliation counters: track the fate of the EXACT SAME originally-violated
    // elements arm (a) counted (23 tetrahedral + 4 E/Z), not just "declared elements
    // in fully-repaired molecules" (a different, smaller denominator that a reviewer
    // re-deriving the numbers would otherwise have to guess at).
    let mut baseline_violated_tet_total = 0usize;
    let mut baseline_violated_tet_fixed = 0usize;
    let mut baseline_violated_ez_total = 0usize;
    let mut baseline_violated_ez_fixed = 0usize;

    for (name, mol, coords) in &embedded {
        let baseline = verify_stereo(mol, coords);
        if baseline.n_violations() == 0 {
            continue; // nothing to repair for this molecule
        }
        n_molecules_with_violations += 1;
        let baseline_violated_tet: Vec<_> = baseline
            .tetrahedral
            .iter()
            .filter(|t| t.status.is_violated())
            .collect();
        let baseline_violated_ez: Vec<_> = baseline
            .double_bond
            .iter()
            .filter(|d| d.status.is_violated())
            .collect();
        baseline_violated_tet_total += baseline_violated_tet.len();
        baseline_violated_ez_total += baseline_violated_ez.len();

        let pre_max_rel = max_bond_rel_error(mol, coords);
        let pre_clash = has_gross_clash(coords);

        match repair_stereo(mol, coords) {
            Ok(outcome) => {
                n_repair_ok += 1;
                for r in &outcome.repaired {
                    total_atoms_moved += r.atoms_moved;
                    if r.max_displacement > max_displacement_seen {
                        max_displacement_seen = r.max_displacement;
                    }
                }
                let post_max_rel = max_bond_rel_error(mol, &outcome.coords);
                let post_clash = has_gross_clash(&outcome.coords);
                if post_max_rel > BOND_BLOWUP_REL_ERROR && pre_max_rel <= BOND_BLOWUP_REL_ERROR {
                    n_new_blowup_introduced += 1;
                    eprintln!(
                        "{name}: repair introduced a NEW bond-length blowup (pre={pre_max_rel:.3} post={post_max_rel:.3})"
                    );
                }
                if post_clash && !pre_clash {
                    n_new_clash_introduced += 1;
                    eprintln!("{name}: repair introduced a NEW gross clash");
                }

                let post = verify_stereo(mol, &outcome.coords);
                for t in &post.tetrahedral {
                    tet_c.declared += 1;
                    match t.status {
                        StereoStatus::Satisfied => tet_c.satisfied += 1,
                        StereoStatus::Violated => tet_c.violated += 1,
                        StereoStatus::Unevaluable(r) => {
                            tet_c.unevaluable += 1;
                            *tet_c.unevaluable_reasons.entry(reason_name(r)).or_insert(0) += 1;
                        }
                    }
                }
                for d in &post.double_bond {
                    ez_c.declared += 1;
                    match d.status {
                        StereoStatus::Satisfied => ez_c.satisfied += 1,
                        StereoStatus::Violated => ez_c.violated += 1,
                        StereoStatus::Unevaluable(r) => {
                            ez_c.unevaluable += 1;
                            *ez_c.unevaluable_reasons.entry(reason_name(r)).or_insert(0) += 1;
                        }
                    }
                }
                // repair_stereo's Ok contract guarantees every originally-violated
                // element in this molecule is now Satisfied -- confirm directly
                // rather than just trusting the contract.
                for t in &baseline_violated_tet {
                    let now = post
                        .tetrahedral
                        .iter()
                        .find(|x| x.atom == t.atom)
                        .map(|x| x.status);
                    assert_eq!(
                        now,
                        Some(StereoStatus::Satisfied),
                        "{name}: Ok(outcome) must fix every originally-violated tetrahedral center"
                    );
                    baseline_violated_tet_fixed += 1;
                }
                for d in &baseline_violated_ez {
                    let now = post
                        .double_bond
                        .iter()
                        .find(|x| x.bond == d.bond)
                        .map(|x| x.status);
                    assert_eq!(
                        now,
                        Some(StereoStatus::Satisfied),
                        "{name}: Ok(outcome) must fix every originally-violated E/Z bond"
                    );
                    baseline_violated_ez_fixed += 1;
                }
            }
            Err(failure) => {
                n_repair_failed += 1;
                for (elem, reason) in &failure.failures {
                    *repair_failure_reasons
                        .entry(repair_reason_name(*reason))
                        .or_insert(0) += 1;
                    eprintln!("{name}: repair FAILED for {elem:?}: {reason:?}");
                }
                // Some elements in this molecule may still have succeeded
                // individually even though the overall call returned Err (a
                // DIFFERENT element failed) -- measure directly against
                // partial_coords rather than assume, for the same reconciliation
                // reason as the Ok branch above.
                let partial = verify_stereo(mol, &failure.partial_coords);
                for t in &baseline_violated_tet {
                    let now = partial
                        .tetrahedral
                        .iter()
                        .find(|x| x.atom == t.atom)
                        .map(|x| x.status);
                    if now == Some(StereoStatus::Satisfied) {
                        baseline_violated_tet_fixed += 1;
                    }
                }
                for d in &baseline_violated_ez {
                    let now = partial
                        .double_bond
                        .iter()
                        .find(|x| x.bond == d.bond)
                        .map(|x| x.status);
                    if now == Some(StereoStatus::Satisfied) {
                        baseline_violated_ez_fixed += 1;
                    }
                }
            }
        }
    }

    println!(
        "molecules with >=1 baseline violation: {n_molecules_with_violations} \
         (repair_stereo fully succeeded: {n_repair_ok}, failed: {n_repair_failed})"
    );
    println!("repair failure reasons: {repair_failure_reasons:?}");
    println!(
        "reconciliation vs arm (a): of {baseline_violated_tet_total} baseline tetrahedral \
         violations, {baseline_violated_tet_fixed} are now Satisfied; of \
         {baseline_violated_ez_total} baseline E/Z violations, {baseline_violated_ez_fixed} are \
         now Satisfied (remainder = the typed failures listed above)"
    );
    println!(
        "NOTE: the two tables below count every DECLARED element (not just the originally- \
         violated ones) in the {n_repair_ok} molecules whose repair fully succeeded -- a smaller, \
         different denominator than arm (a)'s corpus-wide 39/8 (this table's `declared` also \
         includes elements that were already Satisfied at baseline and molecules whose repair \
         failed are excluded entirely). See the reconciliation line above for the apples-to-apples \
         originally-violated-element comparison against arm (a)."
    );
    println!(
        "tetrahedral (post-repair, successful repairs only): declared={} satisfied={} violated={} \
         unevaluable={} | violation_rate(among evaluable)={}",
        tet_c.declared,
        tet_c.satisfied,
        tet_c.violated,
        tet_c.unevaluable,
        rate(tet_c.violated, tet_c.satisfied + tet_c.violated)
    );
    println!(
        "E/Z (post-repair, successful repairs only):         declared={} satisfied={} violated={} \
         unevaluable={} | violation_rate(among evaluable)={}",
        ez_c.declared,
        ez_c.satisfied,
        ez_c.violated,
        ez_c.unevaluable,
        rate(ez_c.violated, ez_c.satisfied + ez_c.violated)
    );
    println!(
        "\n--- geometry disturbance (successful repairs only) ---\n\
         total atoms moved across all repairs: {total_atoms_moved}\n\
         max single-atom displacement seen: {max_displacement_seen:.4} \u{c5}\n\
         NOTE: unaffected-atom RMSD is 0.0000 by construction (repair only translates/rotates \
         the targeted subtree, never touches the rest) -- that is not evidence of \"free\" repair \
         by itself; atoms-moved / max-displacement / new-blowup / new-clash above are the metrics \
         that actually discriminate a bad repair.\n\
         molecules where repair introduced a NEW bond-length blowup (>{BOND_BLOWUP_REL_ERROR:.1} \
         rel. error) that wasn't already present: {n_new_blowup_introduced}\n\
         molecules where repair introduced a NEW gross clash (<{GROSS_CLASH_DIST} \u{c5}) that \
         wasn't already present: {n_new_clash_introduced}"
    );

    println!(
        "\n=== SUMMARY ===\n\
         stereo-declared subset: {}/{} corpus molecules\n\
         arm (a) baseline violation rate: tetrahedral {} ({}/{}), E/Z {} ({}/{})\n\
         arm (b) fixture controls: see above (all must pass)\n\
         arm (c) post-repair violation rate: tetrahedral {} ({}/{}), E/Z {} ({}/{}), \
         {n_repair_failed} molecule(s) with an unrepairable (typed) case",
        stereo_subset.len(),
        CORPUS.len(),
        rate(tet_a.violated, tet_a.satisfied + tet_a.violated),
        tet_a.violated,
        tet_a.satisfied + tet_a.violated,
        rate(ez_a.violated, ez_a.satisfied + ez_a.violated),
        ez_a.violated,
        ez_a.satisfied + ez_a.violated,
        rate(tet_c.violated, tet_c.satisfied + tet_c.violated),
        tet_c.violated,
        tet_c.satisfied + tet_c.violated,
        rate(ez_c.violated, ez_c.satisfied + ez_c.violated),
        ez_c.violated,
        ez_c.satisfied + ez_c.violated,
    );

    assert_eq!(CORPUS.len(), 58, "corpus size drifted from the frozen 58");
}

// ---------------------------------------------------------------------------
// Arm (b): hand-built positive/negative control fixtures.
// ---------------------------------------------------------------------------

fn run_fixture_controls() {
    use chematic_3d::coords::Point3;

    let mut n_pass = 0usize;
    let mut n_total = 0usize;
    let mut check = |label: &str,
                     mol: &Molecule,
                     coords: &Coords3D,
                     want_tet: Option<StereoStatus>,
                     want_ez: Option<StereoStatus>| {
        let report = verify_stereo(mol, coords);
        if let Some(want) = want_tet {
            n_total += 1;
            let got = report.tetrahedral.first().map(|r| r.status);
            let ok = got == Some(want);
            if ok {
                n_pass += 1;
            }
            println!(
                "  [{}] tetrahedral: want={want:?} got={got:?} -> {}",
                label,
                if ok { "PASS" } else { "FAIL" }
            );
        }
        if let Some(want) = want_ez {
            n_total += 1;
            let got = report.double_bond.first().map(|r| r.status);
            let ok = got == Some(want);
            if ok {
                n_pass += 1;
            }
            println!(
                "  [{}] E/Z: want={want:?} got={got:?} -> {}",
                label,
                if ok { "PASS" } else { "FAIL" }
            );
        }
    };

    // Positive control: chfclbr_R declared `[C@H](F)(Cl)Br`, geometry independently
    // verified to be Satisfied -- NOT captured live from `verify_stereo` (that would
    // make this check unable to fail regardless of whether the geometry is actually
    // correct, the same gate flaw the "worked at A, same shape at B" incident
    // elsewhere in this program warns about). The expected value below was derived
    // via a standalone hand computation (phantom position + signed volume, written
    // independently of this module, i.e. not calling `verify_stereo`) and
    // cross-checked against `verify_stereo`'s own output -- both agree. `order[0]`
    // is the implicit-H sentinel for `[C@H](...)` (H inserted at position 0: no
    // preceding atom in the SMILES); phantom position is derived from the 3 real
    // neighbors, so it cannot be placed directly -- the real neighbors' *relative*
    // arrangement is what determines the sign, and CCW/CW intuition from staring at
    // xyz coordinates is exactly what got this wrong the first time, hence the
    // by-hand recomputation rather than another guess.
    // Negative control: same declared molecule, two substituents swapped -> inverts
    // parity -> must read as Violated (this is the ARM the sign was originally,
    // incorrectly, anchored on).
    let m = parse("[C@H](F)(Cl)Br").unwrap();
    let order = m.stereo_neighbor_order(AtomIdx(0)).unwrap().to_vec();
    let real: Vec<u32> = order
        .iter()
        .copied()
        .filter(|&n| n != chematic_core::STEREO_H_SENTINEL)
        .collect();
    let mut correct = Coords3D::new_zeroed(m.atom_count());
    correct.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
    correct.set(AtomIdx(real[0]), Point3::new(1.0, 0.0, -0.3));
    correct.set(AtomIdx(real[1]), Point3::new(-0.5, -0.87, -0.3));
    correct.set(AtomIdx(real[2]), Point3::new(-0.5, 0.87, -0.3));
    let mut wrong = correct.clone();
    let p0 = wrong.get(AtomIdx(real[0]));
    let p1 = wrong.get(AtomIdx(real[1]));
    wrong.set(AtomIdx(real[0]), p1);
    wrong.set(AtomIdx(real[1]), p0);

    check(
        "chfclbr_R positive control",
        &m,
        &correct,
        Some(StereoStatus::Satisfied),
        None,
    );
    check(
        "chfclbr_R negative control (swapped substituents)",
        &m,
        &wrong,
        Some(StereoStatus::Violated),
        None,
    );

    // Second implicit-H positive/negative pair, on a real drug-relevant amino acid
    // shape (l_alanine: N[C@@H](C)C(=O)O) rather than the halomethane fixture above,
    // so this control isn't specific to one molecule shape. Same independent
    // derivation + cross-check as chfclbr_R above -- hardcoded, not self-referential.
    let m_ala = parse("N[C@@H](C)C(=O)O").unwrap();
    let ala_idx = AtomIdx(1);
    let ala_order = m_ala.stereo_neighbor_order(ala_idx).unwrap().to_vec();
    let ala_real: Vec<u32> = ala_order
        .iter()
        .copied()
        .filter(|&n| n != chematic_core::STEREO_H_SENTINEL)
        .collect();
    let mut ala_correct = Coords3D::new_zeroed(m_ala.atom_count());
    ala_correct.set(ala_idx, Point3::new(0.0, 0.0, 0.0));
    ala_correct.set(AtomIdx(ala_real[0]), Point3::new(1.0, 0.0, -0.3));
    ala_correct.set(AtomIdx(ala_real[1]), Point3::new(-0.5, -0.87, -0.3));
    ala_correct.set(AtomIdx(ala_real[2]), Point3::new(-0.5, 0.87, -0.3));
    let mut ala_wrong = ala_correct.clone();
    let ap0 = ala_wrong.get(AtomIdx(ala_real[0]));
    let ap1 = ala_wrong.get(AtomIdx(ala_real[1]));
    ala_wrong.set(AtomIdx(ala_real[0]), ap1);
    ala_wrong.set(AtomIdx(ala_real[1]), ap0);
    check(
        "l_alanine positive control",
        &m_ala,
        &ala_correct,
        Some(StereoStatus::Satisfied),
        None,
    );
    check(
        "l_alanine negative control (swapped substituents)",
        &m_ala,
        &ala_wrong,
        Some(StereoStatus::Violated),
        None,
    );

    // E/Z positive/negative controls: but2ene_E with trans (correct) and cis (wrong) geometry.
    let ez_mol = parse("C/C=C/C").unwrap();
    let mut trans = Coords3D::new_zeroed(ez_mol.atom_count());
    trans.set(AtomIdx(0), Point3::new(-1.5, 0.5, 0.0));
    trans.set(AtomIdx(1), Point3::new(-0.67, 0.0, 0.0));
    trans.set(AtomIdx(2), Point3::new(0.67, 0.0, 0.0));
    trans.set(AtomIdx(3), Point3::new(1.5, -0.5, 0.0));
    let mut cis = trans.clone();
    cis.set(AtomIdx(3), Point3::new(1.5, 0.5, 0.0));
    check(
        "but2ene_E positive control (trans geometry)",
        &ez_mol,
        &trans,
        None,
        Some(StereoStatus::Satisfied),
    );
    check(
        "but2ene_E negative control (cis geometry)",
        &ez_mol,
        &cis,
        None,
        Some(StereoStatus::Violated),
    );

    // Second E/Z pair with UNSYMMETRIC substituents on both alkene ends (mirroring
    // cinnamic_acid_E's shape: OC(=O)/C=C/c1ccccc1 -- a carboxyl-bearing carbon on
    // one end, an aryl-bearing carbon on the other), so the E/Z control isn't
    // specific to a simple, symmetric but-2-ene-style alkene.
    let asym_mol = parse(r"OC(=O)/C=C/c1ccccc1").unwrap();
    // Atoms: 0=O(H), 1=C(=O), 2=O(=), 3=C(alkene, bears the carboxyl branch via atom1),
    // 4=C(alkene, bears the aryl ring via atom5), 5..=ring.
    let mut asym_trans = Coords3D::new_zeroed(asym_mol.atom_count());
    asym_trans.set(AtomIdx(1), Point3::new(-1.5, 0.5, 0.0)); // carboxyl carbon substituent
    asym_trans.set(AtomIdx(3), Point3::new(-0.67, 0.0, 0.0)); // alkene C
    asym_trans.set(AtomIdx(4), Point3::new(0.67, 0.0, 0.0)); // alkene C
    asym_trans.set(AtomIdx(5), Point3::new(1.5, -0.5, 0.0)); // aryl ipso carbon substituent
    let mut asym_cis = asym_trans.clone();
    asym_cis.set(AtomIdx(5), Point3::new(1.5, 0.5, 0.0));
    check(
        "cinnamic_acid_E-shaped positive control (trans geometry, unsymmetric substituents)",
        &asym_mol,
        &asym_trans,
        None,
        Some(StereoStatus::Satisfied),
    );
    check(
        "cinnamic_acid_E-shaped negative control (cis geometry, unsymmetric substituents)",
        &asym_mol,
        &asym_cis,
        None,
        Some(StereoStatus::Violated),
    );

    // Degenerate-geometry control: coplanar tetrahedral substituents must be
    // Unevaluable, not silently Satisfied or Violated.
    let m2 = parse("[C@](F)(Cl)(Br)I").unwrap();
    let order2 = m2.stereo_neighbor_order(AtomIdx(0)).unwrap().to_vec();
    let mut coplanar = Coords3D::new_zeroed(m2.atom_count());
    coplanar.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
    coplanar.set(AtomIdx(order2[0]), Point3::new(1.0, 0.0, 0.0));
    coplanar.set(AtomIdx(order2[1]), Point3::new(-1.0, 0.0, 0.0));
    coplanar.set(AtomIdx(order2[2]), Point3::new(0.0, 1.0, 0.0));
    coplanar.set(AtomIdx(order2[3]), Point3::new(0.0, -1.0, 0.0));
    check(
        "quaternary coplanar (degenerate) control",
        &m2,
        &coplanar,
        Some(StereoStatus::Unevaluable(
            StereoRejectionReason::DegenerateGeometry,
        )),
        None,
    );

    println!("fixture controls: {n_pass}/{n_total} passed");
    if n_pass != n_total {
        eprintln!("GATE FAILED: one or more fixture controls did not match expectation.");
        std::process::exit(1);
    }
}
