//! Acceptance-gate measurement for `distance_geometry_v2::embed_distance_geometry_v2`
//! (3D Breakthrough Program, Wave 1, Agent C — see `docs/rfcs/3d_breakthrough_master_plan.md`
//! §4 for the exact gate this answers).
//!
//! Measures the RAW embedder output (bounds construction → smoothing → Gram/
//! eigendecomposition → bounds-force refinement) on the frozen 58-molecule corpus
//! from `scripts/etkdg_vs_rdkit_gap.py::CORPUS` (transcribed verbatim below — if that
//! list ever changes, update this one too), with **no** MMFF94/DREIDING minimization
//! pass applied afterward, per the master plan's explicit gate scoping. This isolates
//! Agent C's own deliverable from Agent F's (force-field minimization) and Agent E's
//! (torsion knowledge, Wave 2).
//!
//! Deliberately calls `embed_distance_geometry_v2` directly, never through
//! `Mol.conformer_ensemble()` / the live `etkdg.rs` path (which is Coordinator-only
//! and not touched by this PR).
//!
//! # Non-circularity of the validity check
//!
//! `distance_geometry_v2`'s own bounds construction (`dg_fft::ideal_bond_length`) uses
//! `chematic_core::Element::covalent_radius()`. This example deliberately does **not**
//! reuse that table for the bond-length pass/fail check -- that would make the gate
//! close to tautological (measuring whether the bounds hit the targets they were built
//! from). Instead it hardcodes RDKit's own `GetPeriodicTable().GetRcovalent()` values
//! for the elements this corpus needs, dumped read-only from an ISOLATED venv (never
//! the shared repo `.venv` -- see this crate's PR body) against installed RDKit
//! 2025.09.2, matching `scripts/etkdg_vs_rdkit_gap.py::ref_bond_length`'s own
//! external-reference methodology and its `_BOND_ORDER_SCALE` factors exactly.
//!
//! The **angle** check below is different: it reuses `dg_fft::ideal_bond_angle`
//! (the same generic ~109.5°/120° model the embedder's own bounds are built from),
//! so it is explicitly an *internal-consistency* check (does the final geometry
//! agree with the model that constrained it?), not an external-oracle check like
//! the bond-length one. Reported separately, never mixed into the same table.
//!
//! # Formerly a known limitation: 3-membered rings
//!
//! Every 3-membered ring (cyclopropane, epoxide, aziridine, thiirane, ...) used to
//! fail closed with `EmbedFailureCause::BoundsConstructionFailed`, root-caused to
//! `dg_fft::build_bond_angle_bounds`'s angle-constraint loop treating a ring-closing
//! bonded pair as a generic 1-3 (through-center) relationship. Fixed by skipping the
//! angle-derived bound for any neighbor pair that is itself directly bonded -- see
//! `distance_geometry_v2::tests::three_membered_rings_embed_successfully`.
//! **The frozen 58-molecule corpus below still contains zero 3-membered rings**, so
//! the "58/58" gate result never exercised this path either way -- stated explicitly
//! so the gate doesn't read as broader coverage than it has.
//!
//! # 4-layer acceptance gate
//!
//! Redesigned (Coordinator-mandated) from a single "58/58 not-torn" pass/fail into
//! four separately-reported layers, all against the frozen 58-molecule corpus:
//!
//! 1. **Safety gate** (hard-fail on any violation): panics, NaN/Inf coordinates,
//!    all-zero coordinates, catastrophic bond blow-up (defined as bond length >2x
//!    the external RDKit-covalent-radius reference, i.e. relative error >1.0 --
//!    distinct from and much coarser than the 0.5 "torn" threshold `classify()`
//!    uses for the ok/not-ok bucket below).
//! 2. **Precision gate** (reported, not gated): bond-length violation rate @ 10%
//!    and @ 15% tolerance (external RDKit reference), angle violation rate
//!    (internal-consistency reference, see above), gross-clash rate,
//!    bounds-conformance residual (`distance_geometry_v2::bounds_conformance`),
//!    ring-planarity deviation for aromatic rings (Newell's-method best-fit plane).
//! 3. **Novelty gate**: this PR's seeded-stochastic-metrization engine vs. (a)
//!    `dg_fft::generate_coords_dg` (the *unmodified*, deterministic-midpoint
//!    engine that already sits in this same file post-bounds-fix -- see the PR
//!    body for why this comparison matters: it isolates the bounds fix from the
//!    stochastic-metrization novelty) and (b) `dg::generate_coords` (the existing,
//!    separate rule-based DFS placer this module is meant to obsolete). Same
//!    corpus, same tolerance definitions, all three engines.
//! 4. **Stochastic gate**: same-seed reproducibility (bit-identical), adjacent-seed
//!    non-aliasing through the *full* embedding pipeline (not just at the `Prng`
//!    level -- see `prng::from_seed`'s SplitMix64 fix), per-attempt derived-seed
//!    uniqueness (checked at the lib-test level, see
//!    `distance_geometry_v2::tests::derive_attempt_seed_is_distinct_across_attempts`),
//!    and nonzero geometric diversity across a seed ensemble on floppy molecules
//!    (RMSD spread via `chematic_3d::align::align_coords`).
//!
//! No RDKit/Python dependency: this binary is pure Rust so it can run in any `cargo
//! test`/CI environment without a Python venv. A separate, hand-run Python spot check
//! against RDKit (RMSD, chirality coverage) is reported in the PR body directly,
//! using an isolated venv created for this PR (never the shared repo `.venv`).
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example distance_geometry_v2_gap_check
//! ```

use std::collections::BTreeMap;
use std::panic::{self, AssertUnwindSafe};
use std::time::Instant;

use chematic_3d::align::align_coords;
use chematic_3d::coords::{Coords3D, Point3};
use chematic_3d::dg;
use chematic_3d::dg_fft;
use chematic_3d::distance_geometry_v2::{
    EmbedParameters, bounds_conformance, embed_distance_geometry_v2_detail,
};
use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;
use chematic_smiles::parse;

// ---------------------------------------------------------------------------
// Frozen 58-molecule corpus -- verbatim transcription of
// scripts/etkdg_vs_rdkit_gap.py::CORPUS (name, SMILES, category).
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
// External reference (RDKit-sourced, NOT chematic's own covalent_radius table --
// see module docs for why this must stay independent of the embedder's own bounds).
// ---------------------------------------------------------------------------

/// RDKit `GetPeriodicTable().GetRcovalent()` values (Å) for every element this
/// corpus uses. Dumped read-only via an isolated venv (`pip install rdkit` into a
/// throwaway venv, never the shared repo `.venv`) against installed RDKit 2025.09.2:
/// `Chem.GetPeriodicTable().GetRcovalent(Chem.GetPeriodicTable().GetAtomicNumber(sym))`
/// for sym in H,C,N,O,F,P,S,Cl,Br,I. See PR body for the exact dump script/output.
fn rdkit_covalent_radius(atomic_number: u8) -> Option<f64> {
    match atomic_number {
        1 => Some(0.31),  // H
        6 => Some(0.76),  // C
        7 => Some(0.71),  // N
        8 => Some(0.66),  // O
        9 => Some(0.57),  // F
        15 => Some(1.07), // P
        16 => Some(1.05), // S
        17 => Some(1.02), // Cl
        35 => Some(1.20), // Br
        53 => Some(1.39), // I
        _ => None,
    }
}

/// Same bond-order length-scale factors `scripts/etkdg_vs_rdkit_gap.py::_BOND_ORDER_SCALE`
/// uses against RDKit's own `BondType`.
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

// Distinct thresholds for distinct purposes -- kept separate, never collapsed:
//   - BOND_LEN_TOL_FRAC_10 / _15 (0.10 / 0.15): per-bond "violation" counts in
//     bond_violations() (precision-gate signal -- how many individual bonds are
//     off by more than X%, external RDKit-radius reference).
//   - BOND_BLOWUP_REL_ERROR (0.5): the "is this molecule torn" status-bucket
//     decision (checked against `max_rel_error`, not `n_violations`), used by
//     `classify()` below.
//   - CATASTROPHIC_BLOWUP_REL_ERROR (1.0, i.e. >2x the reference length): the
//     safety-gate's "catastrophic bond blow-up" definition -- much coarser than
//     the 0.5 "torn" threshold, on purpose (safety gate should basically never
//     fire; 0.5 already does via classify()'s bond_length_blowup bucket).
const BOND_LEN_TOL_FRAC_10: f64 = 0.10;
const BOND_LEN_TOL_FRAC_15: f64 = 0.15; // matches scripts/etkdg_vs_rdkit_gap.py
const BOND_BLOWUP_REL_ERROR: f64 = 0.5; // matches scripts/etkdg_vs_rdkit_gap.py
const CATASTROPHIC_BLOWUP_REL_ERROR: f64 = 1.0; // >2x reference bond length
const GROSS_CLASH_DIST: f64 = 0.5; // matches scripts/etkdg_vs_rdkit_gap.py
/// Internal-consistency tolerance (degrees) around `dg_fft::ideal_bond_angle`'s
/// generic 109.5°/120° model -- NOT an external-oracle tolerance (see module docs).
const ANGLE_TOL_DEG: f64 = 15.0;

struct BondCheck {
    n_bonds: usize,
    n_violations: usize,
    max_rel_error: f64,
}

/// Identical check to `scripts/etkdg_vs_rdkit_gap.py::bond_violations`, reimplemented
/// in Rust against the external RDKit-radius reference table above. `tol_frac`
/// selects which precision-gate column this call is for (0.10, 0.15) or the
/// safety-gate's catastrophic-blowup threshold (1.0) -- one implementation, three
/// call sites, never three copies of the same loop.
fn bond_violations(mol: &Molecule, coords: &Coords3D, tol_frac: f64) -> BondCheck {
    let mut n_bonds = 0;
    let mut n_violations = 0;
    let mut max_rel_error = 0.0_f64;
    for (_, bond) in mol.bonds() {
        let Some(r0) = ref_bond_length(mol, bond.atom1, bond.atom2) else {
            continue; // element outside this corpus's reference table
        };
        let r = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        let frac = (r - r0).abs() / r0;
        n_bonds += 1;
        if frac > max_rel_error {
            max_rel_error = frac;
        }
        if frac > tol_frac {
            n_violations += 1;
        }
    }
    BondCheck {
        n_bonds,
        n_violations,
        max_rel_error,
    }
}

struct AngleCheck {
    n_angles: usize,
    n_violations: usize,
    max_abs_error_deg: f64,
}

/// Internal-consistency angle check: for every atom with >=2 neighbors, compare the
/// *actual* neighbor-center-neighbor angle in `coords` against
/// `dg_fft::ideal_bond_angle`'s generic model (the same one `build_bound_matrix`
/// used to build this geometry's own angle bounds). Not an external-oracle check --
/// see module docs.
fn angle_violations(mol: &Molecule, coords: &Coords3D) -> AngleCheck {
    let mut n_angles = 0;
    let mut n_violations = 0;
    let mut max_abs_error_deg = 0.0_f64;
    for center_idx in 0..mol.atom_count() {
        let center = AtomIdx(center_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();
        if neighbors.len() < 2 {
            continue;
        }
        let ideal_deg = dg_fft::ideal_bond_angle(mol, center).to_degrees();
        let pc = coords.get(center);
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let vi = coords.get(neighbors[i]).sub(&pc);
                let vj = coords.get(neighbors[j]).sub(&pc);
                let (ni, nj) = (vi.norm(), vj.norm());
                if ni < 1e-9 || nj < 1e-9 {
                    continue; // degenerate/coincident placement, not a real angle
                }
                let cos_theta = (vi.dot(&vj) / (ni * nj)).clamp(-1.0, 1.0);
                let actual_deg = cos_theta.acos().to_degrees();
                let err = (actual_deg - ideal_deg).abs();
                n_angles += 1;
                if err > max_abs_error_deg {
                    max_abs_error_deg = err;
                }
                if err > ANGLE_TOL_DEG {
                    n_violations += 1;
                }
            }
        }
    }
    AngleCheck {
        n_angles,
        n_violations,
        max_abs_error_deg,
    }
}

/// Best-fit-plane deviation for one ring: Newell's method for the plane normal
/// (robust to non-convex/near-degenerate point sets, no SVD needed), then the max
/// perpendicular distance of any ring atom from the centroid plane, in Å.
fn ring_planarity_deviation(ring: &[AtomIdx], coords: &Coords3D) -> f64 {
    let pts: Vec<Point3> = ring.iter().map(|&idx| coords.get(idx)).collect();
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut normal = Point3::new(0.0, 0.0, 0.0);
    let mut centroid = Point3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let p = pts[i];
        let q = pts[(i + 1) % n];
        normal = Point3::new(
            normal.x + (p.y - q.y) * (p.z + q.z),
            normal.y + (p.z - q.z) * (p.x + q.x),
            normal.z + (p.x - q.x) * (p.y + q.y),
        );
        centroid = centroid.add(&p);
    }
    centroid = centroid.scale(1.0 / n as f64);
    let Some(normal) = normal.try_normalize() else {
        return 0.0; // degenerate (collinear points) -- shouldn't happen for a real ring
    };
    pts.iter()
        .map(|p| p.sub(&centroid).dot(&normal).abs())
        .fold(0.0_f64, f64::max)
}

/// (count, max deviation) across every all-aromatic SSSR ring in `mol` -- count is
/// the number of *rings*, not molecules, so a fused system (e.g. naphthalene, 2
/// aromatic rings) contributes 2, not 1. Corpus uses lowercase aromatic SMILES
/// throughout, so `atom.aromatic` is already set by the parser -- no
/// `apply_aromaticity` call needed.
fn aromatic_ring_deviations(mol: &Molecule, coords: &Coords3D) -> (usize, f64) {
    let sssr = find_sssr(mol);
    let mut n_rings = 0usize;
    let mut max_dev = 0.0_f64;
    for ring in sssr.rings() {
        if ring.iter().all(|&idx| mol.atom(idx).aromatic) {
            n_rings += 1;
            max_dev = max_dev.max(ring_planarity_deviation(ring, coords));
        }
    }
    (n_rings, max_dev)
}

fn gross_clash(coords: &Coords3D) -> bool {
    let n = coords.atom_count();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)));
            if d < GROSS_CLASH_DIST {
                return true;
            }
        }
    }
    false
}

fn all_finite(coords: &Coords3D) -> bool {
    (0..coords.atom_count()).all(|i| {
        let p = coords.get(AtomIdx(i as u32));
        p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
    })
}

/// Safety-gate check: every coordinate exactly (0,0,0) for a molecule with >1 atom
/// is a real bug signature (the only legitimate all-zero case is the n==1 trivial
/// molecule, not present in this corpus).
fn all_zero_coords(coords: &Coords3D) -> bool {
    coords.atom_count() > 1
        && (0..coords.atom_count()).all(|i| {
            let p = coords.get(AtomIdx(i as u32));
            p.x == 0.0 && p.y == 0.0 && p.z == 0.0
        })
}

/// Status bucket for one engine's output on one molecule, mirroring
/// scripts/etkdg_vs_rdkit_gap.py's named-bucket status strings (no silent drops).
fn classify(mol: &Molecule, coords: Option<&Coords3D>, embed_err: Option<String>) -> String {
    if let Some(cause) = embed_err {
        return format!("embed_failed:{cause}");
    }
    let coords = coords.expect("coords must be present when embed_err is None");
    if !all_finite(coords) {
        return "nonfinite_coords".to_string();
    }
    if gross_clash(coords) {
        return "gross_clash".to_string();
    }
    let check = bond_violations(mol, coords, BOND_LEN_TOL_FRAC_15);
    if check.max_rel_error > BOND_BLOWUP_REL_ERROR {
        "bond_length_blowup".to_string()
    } else {
        "ok".to_string()
    }
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn rate_str(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        "null".to_string()
    } else {
        format!("{:.4}", numerator as f64 / denominator as f64)
    }
}

/// Re-measure `geometrically_valid_rate` for the new embedder only, at a given seed --
/// used by the seed-robustness sweep below so the headline 100% isn't reported from a
/// single cherry-picked seed. Deliberately duplicates only the (cheap) status
/// classification, not the full row/JSON reporting in `main`.
fn measure_validity_at_seed(seed: u64) -> (usize, usize) {
    let mut n_ok = 0usize;
    let mut n_total = 0usize;
    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse ({name}): {e:?}"));
        let params = EmbedParameters {
            random_seed: seed,
            ..EmbedParameters::default()
        };
        n_total += 1;
        let status = match embed_distance_geometry_v2_detail(&mol, &params) {
            Ok((coords, _stats)) => classify(&mol, Some(&coords), None),
            Err((cause, _stats)) => classify(&mol, None, Some(format!("{cause:?}"))),
        };
        if status == "ok" {
            n_ok += 1;
        }
    }
    (n_ok, n_total)
}

fn main() {
    let mut rows = Vec::new();
    let mut new_status_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut dg_status_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut dgfft_status_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut n_total = 0usize;
    let mut n_new_ok = 0usize;
    let mut n_dg_ok = 0usize;
    let mut n_dgfft_ok = 0usize;
    let mut n_both_ok = 0usize;
    let mut n_regressions = 0usize; // dg (legacy rule-based) ok, new NOT ok -- must stay 0
    let mut regression_names: Vec<&str> = Vec::new();
    let mut new_bonds_checked_10 = 0usize;
    let mut new_bond_violations_10 = 0usize;
    let mut new_bonds_checked_15 = 0usize;
    let mut new_bond_violations_15 = 0usize;
    let mut dg_bonds_checked = 0usize;
    let mut dg_bond_violations = 0usize;
    let mut dgfft_bonds_checked = 0usize;
    let mut dgfft_bond_violations = 0usize;
    // Diagnostic (NOT the acceptance gate): how well does the final geometry
    // satisfy the smoothed distance bounds it was built from, over ALL atom pairs
    // (not just bonded ones)? This is a different question from
    // `smoothing_preserves_invariants` (bounds-vs-bounds, before any geometry
    // exists) -- this is geometry-vs-bounds, after refinement.
    let mut bc_pairs_total = 0usize;
    let mut bc_violations_total = 0usize;
    let mut bc_max_rel_violation = 0.0_f64;
    let mut bc_worst_molecule = "";
    // Precision gate: angle + ring-planarity accumulators (new engine only).
    let mut angle_total = 0usize;
    let mut angle_violations_total = 0usize;
    let mut angle_max_abs_error_deg = 0.0_f64;
    let mut ring_planarity_n_rings = 0usize;
    let mut ring_planarity_max_dev = 0.0_f64;
    let mut ring_planarity_worst_molecule = "";
    // Safety gate accumulators (new engine only -- this is the gate this PR must
    // never fail).
    let mut n_panics = 0usize;
    let mut n_nonfinite = 0usize;
    let mut n_all_zero = 0usize;
    let mut n_catastrophic_blowup = 0usize;

    for &(name, smiles, category) in CORPUS {
        n_total += 1;
        let mol = parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse ({name}): {e:?}"));

        // --- safety gate: does the call itself panic? ---
        let embed_result = panic::catch_unwind(AssertUnwindSafe(|| {
            embed_distance_geometry_v2_detail(&mol, &EmbedParameters::default())
        }));
        let embed_result = match embed_result {
            Ok(r) => r,
            Err(_) => {
                n_panics += 1;
                eprintln!("SAFETY GATE VIOLATION: {name} panicked during embedding");
                continue; // no coords to check further for this molecule
            }
        };

        // --- new embedder (this PR's deliverable), raw output, no minimization ---
        let (new_coords, new_err, new_stats_line) = match embed_result {
            Ok((coords, stats)) => (
                Some(coords),
                None,
                format!(
                    "attempts_used={} negative_eigs={} max_neg_mag={:.4} used_random_coords={}",
                    stats.attempts_used,
                    stats.negative_eigenvalues_beyond_embedding_dim,
                    stats.max_negative_eigenvalue_magnitude,
                    stats.used_random_coords
                ),
            ),
            Err((cause, stats)) => (
                None,
                Some(format!("{cause:?}")),
                format!("attempts_used={}", stats.attempts_used),
            ),
        };
        let new_status = classify(&mol, new_coords.as_ref(), new_err.clone());
        *new_status_counts.entry(new_status.clone()).or_insert(0) += 1;
        let new_ok = new_status == "ok";
        if new_ok {
            n_new_ok += 1;
        }

        // --- remaining safety-gate checks on the returned coords ---
        if let Some(c) = new_coords.as_ref() {
            if !all_finite(c) {
                n_nonfinite += 1;
                eprintln!("SAFETY GATE VIOLATION: {name} produced non-finite coordinates");
            }
            if all_zero_coords(c) {
                n_all_zero += 1;
                eprintln!("SAFETY GATE VIOLATION: {name} produced all-zero coordinates");
            }
            if bond_violations(&mol, c, CATASTROPHIC_BLOWUP_REL_ERROR).n_violations > 0 {
                n_catastrophic_blowup += 1;
                eprintln!("SAFETY GATE VIOLATION: {name} has a bond >2x the reference length");
            }
        }

        // --- novelty gate, arm "dgfft_fixed_midpoint": dg_fft::generate_coords_dg
        // (unmodified, deterministic-midpoint, but on TODAY's bounds-fixed
        // build_bound_matrix -- see PR body for the historical "dgfft_unfixed"
        // baseline measured separately) ---
        let dgfft_coords = dg_fft::generate_coords_dg(&mol);
        let dgfft_status = classify(&mol, Some(&dgfft_coords), None);
        *dgfft_status_counts.entry(dgfft_status.clone()).or_insert(0) += 1;
        if dgfft_status == "ok" {
            n_dgfft_ok += 1;
        }
        let dgfft_bv = bond_violations(&mol, &dgfft_coords, BOND_LEN_TOL_FRAC_15);
        dgfft_bonds_checked += dgfft_bv.n_bonds;
        dgfft_bond_violations += dgfft_bv.n_violations;

        // --- novelty gate, arm "dg_legacy_dfs": dg::generate_coords (existing
        // rule-based DFS placer this module is meant to obsolete) ---
        let dg_coords = dg::generate_coords(&mol);
        let dg_status = classify(&mol, Some(&dg_coords), None);
        *dg_status_counts.entry(dg_status.clone()).or_insert(0) += 1;
        let dg_ok = dg_status == "ok";
        if dg_ok {
            n_dg_ok += 1;
        }

        if new_ok && dg_ok {
            n_both_ok += 1;
        }
        if dg_ok && !new_ok {
            n_regressions += 1;
            regression_names.push(name);
        }

        let new_max_rel = new_coords.as_ref().map(|c| {
            let bv10 = bond_violations(&mol, c, BOND_LEN_TOL_FRAC_10);
            new_bonds_checked_10 += bv10.n_bonds;
            new_bond_violations_10 += bv10.n_violations;
            let bv15 = bond_violations(&mol, c, BOND_LEN_TOL_FRAC_15);
            new_bonds_checked_15 += bv15.n_bonds;
            new_bond_violations_15 += bv15.n_violations;

            let bc = bounds_conformance(&mol, c);
            bc_pairs_total += bc.n_pairs;
            bc_violations_total += bc.n_violations;
            if bc.max_rel_violation > bc_max_rel_violation {
                bc_max_rel_violation = bc.max_rel_violation;
                bc_worst_molecule = name;
            }

            let ac = angle_violations(&mol, c);
            angle_total += ac.n_angles;
            angle_violations_total += ac.n_violations;
            if ac.max_abs_error_deg > angle_max_abs_error_deg {
                angle_max_abs_error_deg = ac.max_abs_error_deg;
            }

            let (n_rings, dev) = aromatic_ring_deviations(&mol, c);
            ring_planarity_n_rings += n_rings;
            if n_rings > 0 && dev > ring_planarity_max_dev {
                ring_planarity_max_dev = dev;
                ring_planarity_worst_molecule = name;
            }

            bv15.max_rel_error
        });
        let dg_bv = bond_violations(&mol, &dg_coords, BOND_LEN_TOL_FRAC_15);
        dg_bonds_checked += dg_bv.n_bonds;
        dg_bond_violations += dg_bv.n_violations;
        let dg_max_rel = dg_bv.max_rel_error;

        rows.push(format!(
            "{{\"name\":\"{}\",\"category\":\"{}\",\"n_atoms\":{},\"new_status\":\"{}\",\"new_max_rel_error\":{},\"new_stats\":\"{}\",\"dgfft_status\":\"{}\",\"dg_status\":\"{}\",\"dg_max_rel_error\":{:.4}}}",
            json_escape(name),
            json_escape(category),
            mol.atom_count(),
            json_escape(&new_status),
            new_max_rel.map(|v| format!("{v:.4}")).unwrap_or_else(|| "null".to_string()),
            json_escape(&new_stats_line),
            json_escape(&dgfft_status),
            json_escape(&dg_status),
            dg_max_rel,
        ));

        println!("{new_status:<28} (dgfft: {dgfft_status:<12} dg: {dg_status:<12}) {name}");
    }

    // =========================================================================
    // LAYER 1: SAFETY GATE (hard-fail on any nonzero count)
    // =========================================================================
    println!("\n=== LAYER 1: SAFETY GATE (new embedder, {n_total} molecules) ===");
    println!("panics:                        {n_panics} (must be 0)");
    println!("NaN/Inf coordinates:           {n_nonfinite} (must be 0)");
    println!("all-zero coordinates:          {n_all_zero} (must be 0)");
    println!(
        "catastrophic bond blow-up:     {n_catastrophic_blowup} (must be 0; defined as any bond \
         >2x the external RDKit-covalent-radius reference length, i.e. relative error >{CATASTROPHIC_BLOWUP_REL_ERROR:.1} \
         -- distinct from and much coarser than the {BOND_BLOWUP_REL_ERROR:.1} \"torn\" threshold below)"
    );
    let safety_gate_violations = n_panics + n_nonfinite + n_all_zero + n_catastrophic_blowup;
    if safety_gate_violations > 0 {
        eprintln!(
            "\nSAFETY GATE FAILED: {safety_gate_violations} total violation(s) across the checks above."
        );
        std::process::exit(1);
    }
    println!("SAFETY GATE PASSED: 0 panics, 0 NaN/Inf, 0 all-zero, 0 catastrophic blow-up.");

    // =========================================================================
    // LAYER 2: PRECISION GATE (reported, not a pass/fail threshold)
    // =========================================================================
    let geometrically_valid_rate = n_new_ok as f64 / n_total as f64;
    println!("\n=== LAYER 2: PRECISION GATE (new embedder, {n_total} molecules) ===");
    println!(
        "not-torn rate (classify()==\"ok\", {BOND_BLOWUP_REL_ERROR:.1} max-rel-error threshold): \
         {geometrically_valid_rate:.4} ({n_new_ok}/{n_total})"
    );
    println!(
        "bond-length violation rate @ 10% tol (external RDKit reference): {} ({new_bond_violations_10}/{new_bonds_checked_10} bonds)",
        rate_str(new_bond_violations_10, new_bonds_checked_10)
    );
    println!(
        "bond-length violation rate @ 15% tol (external RDKit reference): {} ({new_bond_violations_15}/{new_bonds_checked_15} bonds)",
        rate_str(new_bond_violations_15, new_bonds_checked_15)
    );
    println!(
        "angle violation rate @ {ANGLE_TOL_DEG:.0} deg tol (INTERNAL-CONSISTENCY reference, \
         dg_fft::ideal_bond_angle's own generic 109.5/120 deg model -- not RDKit-referenced): \
         {} ({angle_violations_total}/{angle_total} angles, max abs error {angle_max_abs_error_deg:.1} deg)",
        rate_str(angle_violations_total, angle_total)
    );
    println!(
        "gross-clash rate (any non-bonded pair < {GROSS_CLASH_DIST} Angstrom): {}",
        rate_str(*new_status_counts.get("gross_clash").unwrap_or(&0), n_total)
    );
    println!(
        "bounds-conformance residual (all pairs vs. the smoothed bounds they were built from): \
         {} ({bc_violations_total}/{bc_pairs_total} pairs), max relative violation {bc_max_rel_violation:.4} ({})",
        rate_str(bc_violations_total, bc_pairs_total),
        json_escape(bc_worst_molecule)
    );
    println!(
        "ring-planarity: {ring_planarity_n_rings} all-aromatic SSSR ring(s) checked across the \
         corpus, max out-of-plane deviation {ring_planarity_max_dev:.4} Angstrom ({})",
        json_escape(ring_planarity_worst_molecule)
    );

    // =========================================================================
    // LAYER 3: NOVELTY GATE -- 3-arm comparison, same corpus, same tolerances
    //
    // Arm names are deliberately descriptive, not lettered (a)/(b)/(c) -- this is
    // the exact same 3-arm comparison the PR body's attribution table reports, and
    // a lettering scheme that doesn't match 1:1 between this binary's stdout and
    // the PR body's prose is a real, easy-to-introduce mismatch (a 4th, purely
    // historical "unfixed dg_fft" arm is reported ONLY in the PR body, measured via
    // a one-time uncommitted local revert -- see the note printed below).
    // =========================================================================
    let dgfft_valid_rate = n_dgfft_ok as f64 / n_total as f64;
    let dg_valid_rate = n_dg_ok as f64 / n_total as f64;
    println!("\n=== LAYER 3: NOVELTY GATE (3-arm comparison, {n_total} molecules) ===");
    println!(
        "this_pr_stochastic   embed_distance_geometry_v2 (seeded stochastic metrization): \
         not-torn {geometrically_valid_rate:.4} ({n_new_ok}/{n_total}), bond-violation-rate@15% {}",
        rate_str(new_bond_violations_15, new_bonds_checked_15)
    );
    println!(
        "dgfft_fixed_midpoint dg_fft::generate_coords_dg (deterministic-midpoint, SAME bounds-fixed \
         build_bound_matrix as this PR): not-torn {dgfft_valid_rate:.4} ({n_dgfft_ok}/{n_total}), \
         bond-violation-rate@15% {}",
        rate_str(dgfft_bond_violations, dgfft_bonds_checked)
    );
    println!(
        "dg_legacy_dfs        dg::generate_coords (existing, separate rule-based DFS placer): \
         not-torn {dg_valid_rate:.4} ({n_dg_ok}/{n_total}), bond-violation-rate@15% {}",
        rate_str(dg_bond_violations, dg_bonds_checked)
    );
    println!(
        "n_both_ok (this_pr_stochastic AND dg_legacy_dfs both not-torn): {n_both_ok}; \
         regressions (dg_legacy_dfs not-torn, this_pr_stochastic torn, must be 0): {n_regressions} {regression_names:?}"
    );
    println!(
        "See PR body for the separately-measured historical arm, dgfft_unfixed: dg_fft with the PRE-this-PR \
         (buggy, 9-entry-table) ideal_bond_length/vdw_radius -- not reproducible from this \
         binary alone since that table no longer exists in this file."
    );

    let status_counts_json = |m: &BTreeMap<String, usize>| {
        m.iter()
            .map(|(k, v)| format!("\"{}\":{}", json_escape(k), v))
            .collect::<Vec<_>>()
            .join(",")
    };

    let summary = format!(
        "{{\n  \"n_molecules\": {n_total},\n  \"new_status_counts\": {{{}}},\n  \"dgfft_status_counts\": {{{}}},\n  \"dg_status_counts\": {{{}}},\n  \"geometrically_valid_rate_new\": {geometrically_valid_rate:.4},\n  \"geometrically_valid_rate_dgfft\": {dgfft_valid_rate:.4},\n  \"geometrically_valid_rate_dg\": {dg_valid_rate:.4},\n  \"n_both_ok\": {n_both_ok},\n  \"n_regressions_dg_ok_new_not_ok\": {n_regressions},\n  \"regression_names\": {:?},\n  \"new_bond_violation_rate_at_10pct_tol\": {},\n  \"new_bond_violation_rate_at_15pct_tol\": {},\n  \"dgfft_bond_violation_rate_at_15pct_tol\": {},\n  \"dg_bond_violation_rate_at_15pct_tol\": {},\n  \"angle_violation_rate\": {},\n  \"bounds_conformance_all_pairs_violation_rate\": {},\n  \"bounds_conformance_max_rel_violation\": {:.4},\n  \"bounds_conformance_worst_molecule\": \"{}\",\n  \"ring_planarity_n_rings_checked\": {ring_planarity_n_rings},\n  \"ring_planarity_max_deviation\": {ring_planarity_max_dev:.4}\n}}",
        status_counts_json(&new_status_counts),
        status_counts_json(&dgfft_status_counts),
        status_counts_json(&dg_status_counts),
        regression_names,
        rate_str(new_bond_violations_10, new_bonds_checked_10),
        rate_str(new_bond_violations_15, new_bonds_checked_15),
        rate_str(dgfft_bond_violations, dgfft_bonds_checked),
        rate_str(dg_bond_violations, dg_bonds_checked),
        rate_str(angle_violations_total, angle_total),
        rate_str(bc_violations_total, bc_pairs_total),
        bc_max_rel_violation,
        json_escape(bc_worst_molecule),
    );

    println!("\n--- summary (JSON) ---");
    println!("{summary}");

    println!("\n--- rows (JSONL) ---");
    for row in &rows {
        println!("{row}");
    }

    assert_eq!(
        n_total, 58,
        "corpus size drifted from the frozen 58 -- re-sync with scripts/etkdg_vs_rdkit_gap.py::CORPUS"
    );

    if n_regressions > 0 {
        eprintln!(
            "\nGATE FAILED: {n_regressions} molecule(s) where the existing dg::generate_coords \
             was not-torn but this PR's embedder was not. See regression_names above."
        );
        std::process::exit(1);
    }
    println!(
        "\nnot-torn rate (this PR): {geometrically_valid_rate:.4} on the frozen 58 (raw \
         embedder, pre-minimization; see the precision/novelty gates above for what \
         \"not-torn\" does and doesn't tell you -- this single number is not, by itself, \
         evidence of anything beyond the {BOND_BLOWUP_REL_ERROR:.1}-threshold torn/not-torn split)."
    );

    // --- seed-robustness sweep: don't report the not-torn rate from one cherry-picked seed ---
    println!("\n--- seed-robustness sweep (new embedder only, not-torn rate per seed) ---");
    let mut all_seeds_perfect = true;
    for seed in [0u64, 1, 2, 42, 999, 0xDEAD_BEEF, u64::MAX] {
        let (n_ok, n_tot) = measure_validity_at_seed(seed);
        let rate = n_ok as f64 / n_tot as f64;
        println!("seed={seed:#x}  not_torn_rate={rate:.4}  ({n_ok}/{n_tot})");
        if n_ok != n_tot {
            all_seeds_perfect = false;
        }
    }
    if all_seeds_perfect {
        println!("All swept seeds reach 1.0 -- the not-torn rate is not a single-seed artifact.");
    } else {
        println!(
            "NOTE: not every swept seed reaches 1.0 (retries/max_attempts absorb per-attempt \
             draws via EmbedParameters::default()'s max_attempts=8 -- see stats.attempts_used \
             per molecule in the JSONL rows above for the default seed)."
        );
    }

    // =========================================================================
    // LAYER 4: STOCHASTIC GATE
    // =========================================================================
    println!("\n=== LAYER 4: STOCHASTIC GATE ===");

    // --- same-seed reproducibility: bit-identical, across the whole corpus ---
    let default_seed = EmbedParameters::default().random_seed;
    let (n_mismatch, n_checked) = same_seed_reproducible_all(default_seed);
    println!(
        "same-seed reproducibility (seed={default_seed:#x}, two independent calls, bit-identical \
         coordinates required): {}/{n_checked} molecules reproduced exactly \
         (mismatches: {n_mismatch})",
        n_checked - n_mismatch
    );

    // --- adjacent-seed non-aliasing through the FULL embedding pipeline, not just
    // the Prng level (that's covered by prng::tests::from_seed_adjacent_seeds_do_not_collide) ---
    let (n_differ, n_checked_adj) = adjacent_seed_differs_all(0);
    println!(
        "adjacent-seed non-aliasing (seed 0 vs seed 1, full pipeline): {n_differ}/{n_checked_adj} \
         molecules produced DIFFERENT coordinates (expected: all of them -- a collision here \
         would mean the pre-SplitMix64 `seed | 1` bug, or an equivalent, survived somewhere \
         downstream of Prng::from_seed itself)"
    );

    // --- per-attempt derived-seed uniqueness: checked at the lib-test level, see
    // distance_geometry_v2::tests::derive_attempt_seed_is_distinct_across_attempts ---
    println!(
        "per-attempt derived-seed uniqueness: checked at the lib-test level \
         (distance_geometry_v2::tests::derive_attempt_seed_is_distinct_across_attempts) -- \
         32 attempts, 5 base seeds, all pairwise distinct."
    );

    // --- nonzero diversity across a seed ensemble on floppy molecules ---
    println!(
        "\nensemble diversity (pairwise RMSD in Angstrom across {} seeds, after Kabsch alignment, \
         on floppy corpus molecules with real rotatable bonds):",
        DIVERSITY_SEEDS.len()
    );
    for &mol_name in DIVERSITY_MOLECULES {
        let smiles = CORPUS
            .iter()
            .find(|&&(name, _, _)| name == mol_name)
            .map(|&(_, smiles, _)| smiles)
            .expect("diversity molecule must be in CORPUS");
        let (min, mean, max) = measure_ensemble_diversity(mol_name, smiles);
        println!("  {mol_name:<20} min={min:.4}  mean={mean:.4}  max={max:.4}");
        assert!(
            mean > 0.0,
            "{mol_name}: expected nonzero geometric diversity across different seeds"
        );
    }

    // --- runtime: honest per-engine comparison, not just the new embedder's own cost ---
    // Separate dedicated timing pass (not mixed into the JSON-building loop above) so
    // JSON/string-formatting overhead doesn't pollute the measurement. Single-threaded,
    // wall-clock, same machine/run as everything else in this binary -- reported as
    // p50/p95/max, not a single mean (this repo's own issue #70 lesson: don't reduce
    // noisy wall-clock measurements to one number).
    println!(
        "\n--- runtime (new embedder vs. dg::generate_coords, single-threaded, wall-clock) ---"
    );
    let mut new_times_ms: Vec<f64> = Vec::with_capacity(CORPUS.len());
    let mut dg_times_ms: Vec<f64> = Vec::with_capacity(CORPUS.len());
    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles).unwrap();
        let params = EmbedParameters::default();
        let t0 = Instant::now();
        let _ = embed_distance_geometry_v2_detail(&mol, &params);
        new_times_ms.push(t0.elapsed().as_secs_f64() * 1000.0);

        let t0 = Instant::now();
        let _ = dg::generate_coords(&mol);
        dg_times_ms.push(t0.elapsed().as_secs_f64() * 1000.0);
        let _ = name;
    }
    let (new_p50, new_p95, new_max) = percentiles(&mut new_times_ms);
    let (dg_p50, dg_p95, dg_max) = percentiles(&mut dg_times_ms);
    println!("new embedder:      p50={new_p50:.3}ms  p95={new_p95:.3}ms  max={new_max:.3}ms");
    println!("dg::generate_coords: p50={dg_p50:.3}ms  p95={dg_p95:.3}ms  max={dg_max:.3}ms");
    println!(
        "(new embedder is ~{:.0}x slower than the DFS placer at the median -- expected and \
         acceptable for real distance geometry, i.e. an actual bounds/Gram/eigendecomposition \
         solve vs. a template walk; stated plainly, not omitted.)",
        new_p50 / dg_p50.max(1e-9)
    );
}

/// (p50, p95, max) in milliseconds. Sorts `values` in place.
fn percentiles(values: &mut [f64]) -> (f64, f64, f64) {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = values.len();
    let p50 = values[n / 2];
    let p95_idx = (((n as f64) * 0.95) as usize).min(n - 1);
    let p95 = values[p95_idx];
    let max = values[n - 1];
    (p50, p95, max)
}

/// Stochastic gate: run `embed_distance_geometry_v2_detail` twice with the SAME seed
/// on every corpus molecule and check the returned coordinates are bit-identical
/// (`Point3` derives `PartialEq` on exact f64 equality -- both calls run the same
/// deterministic arithmetic given the same seed, so this is a real "bit-identical",
/// not a "close enough" check). Returns (n_mismatch, n_checked).
fn same_seed_reproducible_all(seed: u64) -> (usize, usize) {
    let mut n_mismatch = 0usize;
    let mut n_checked = 0usize;
    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles).unwrap();
        let params = EmbedParameters {
            random_seed: seed,
            ..EmbedParameters::default()
        };
        let r1 = embed_distance_geometry_v2_detail(&mol, &params);
        let r2 = embed_distance_geometry_v2_detail(&mol, &params);
        n_checked += 1;
        let identical = match (&r1, &r2) {
            (Ok((c1, _)), Ok((c2, _))) => {
                (0..c1.atom_count()).all(|i| c1.get(AtomIdx(i as u32)) == c2.get(AtomIdx(i as u32)))
            }
            (Err((e1, _)), Err((e2, _))) => e1 == e2,
            _ => false,
        };
        if !identical {
            n_mismatch += 1;
            eprintln!("STOCHASTIC GATE VIOLATION: {name} not reproducible at seed {seed:#x}");
        }
    }
    (n_mismatch, n_checked)
}

/// Stochastic gate: run the FULL embedding pipeline (not just `Prng::from_seed`
/// directly) at `base` and `base + 1` on every corpus molecule and check the
/// returned coordinates actually differ. This is the pipeline-level guarantee the
/// SplitMix64 fix in `prng::from_seed` is supposed to provide -- checked here, not
/// assumed to hold just because the `Prng`-level unit tests pass. Returns
/// (n_differ, n_checked).
fn adjacent_seed_differs_all(base: u64) -> (usize, usize) {
    let mut n_differ = 0usize;
    let mut n_checked = 0usize;
    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles).unwrap();
        let p1 = EmbedParameters {
            random_seed: base,
            ..EmbedParameters::default()
        };
        let p2 = EmbedParameters {
            random_seed: base.wrapping_add(1),
            ..EmbedParameters::default()
        };
        let r1 = embed_distance_geometry_v2_detail(&mol, &p1);
        let r2 = embed_distance_geometry_v2_detail(&mol, &p2);
        n_checked += 1;
        if let (Ok((c1, _)), Ok((c2, _))) = (r1, r2) {
            let differs = (0..c1.atom_count())
                .any(|i| c1.get(AtomIdx(i as u32)) != c2.get(AtomIdx(i as u32)));
            if differs {
                n_differ += 1;
            } else {
                eprintln!(
                    "STOCHASTIC GATE VIOLATION: {name} produced IDENTICAL coordinates for seeds \
                     {base:#x} and {:#x}",
                    base.wrapping_add(1)
                );
            }
        }
    }
    (n_differ, n_checked)
}

/// Floppy corpus molecules (real rotatable bonds) used for the ensemble-diversity
/// measurement -- picked from `flexible_chain`/`druglike` categories, mirroring
/// `scripts/etkdg_vs_rdkit_gap.py::ENSEMBLE_SUBSET`'s intent.
const DIVERSITY_MOLECULES: &[&str] = &["decane", "diphenhydramine", "hexadecane"];
const DIVERSITY_SEEDS: &[u64] = &[1, 2, 3, 4, 5, 6, 7, 8];

/// (min, mean, max) pairwise RMSD (Angstrom, after Kabsch alignment via
/// `chematic_3d::align::align_coords`) across `DIVERSITY_SEEDS`' worth of independent
/// embeddings of the named molecule. A translation/rotation-invariant diversity
/// signal, not a raw coordinate diff (which would conflate rigid-body placement
/// differences with real conformational diversity).
fn measure_ensemble_diversity(name: &str, smiles: &str) -> (f64, f64, f64) {
    let mol = parse(smiles).unwrap();
    let mut geoms: Vec<Vec<[f64; 3]>> = Vec::with_capacity(DIVERSITY_SEEDS.len());
    for &seed in DIVERSITY_SEEDS {
        let params = EmbedParameters {
            random_seed: seed,
            ..EmbedParameters::default()
        };
        let (coords, _stats) =
            embed_distance_geometry_v2_detail(&mol, &params).unwrap_or_else(|(cause, _)| {
                panic!("{name} failed to embed at seed {seed}: {cause:?}")
            });
        let pts: Vec<[f64; 3]> = (0..coords.atom_count())
            .map(|i| {
                let p = coords.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        geoms.push(pts);
    }
    let mut rmsds = Vec::new();
    for i in 0..geoms.len() {
        for j in (i + 1)..geoms.len() {
            rmsds.push(align_coords(&geoms[i], &geoms[j]).rmsd);
        }
    }
    let min = rmsds.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = rmsds.iter().cloned().fold(0.0_f64, f64::max);
    let mean = rmsds.iter().sum::<f64>() / rmsds.len() as f64;
    (min, mean, max)
}
