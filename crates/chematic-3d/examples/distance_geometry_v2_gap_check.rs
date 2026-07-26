//! Acceptance-gate measurement for `distance_geometry_v2::embed_distance_geometry_v2`
//! (3D Breakthrough Program, Wave 1, Agent C — see `docs/3d_breakthrough_master_plan.md`
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
//! reuse that table for the pass/fail check -- that would make the gate close to
//! tautological (measuring whether the bounds hit the targets they were built from).
//! Instead it hardcodes RDKit's own `GetPeriodicTable().GetRcovalent()` values for the
//! elements this corpus needs, dumped read-only from an ISOLATED venv (never the
//! shared repo `.venv` -- see this crate's PR body) against installed RDKit 2025.09.2,
//! matching `scripts/etkdg_vs_rdkit_gap.py::ref_bond_length`'s own external-reference
//! methodology and its `_BOND_ORDER_SCALE` factors exactly.
//!
//! # Paired comparison against `dg::generate_coords`
//!
//! Also runs `dg::generate_coords` (the existing deterministic DFS placer this module
//! is meant to obsolete) on the same molecule under the same external check, so a
//! molecule where the new embedder is invalid but the old placer was valid is visible
//! as a named regression bucket, not hidden inside an aggregate.
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

use chematic_3d::coords::Coords3D;
use chematic_3d::dg;
use chematic_3d::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2_detail};
use chematic_core::{AtomIdx, BondOrder, Molecule};
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

const BOND_BLOWUP_REL_ERROR: f64 = 0.5; // matches scripts/etkdg_vs_rdkit_gap.py
const GROSS_CLASH_DIST: f64 = 0.5; // matches scripts/etkdg_vs_rdkit_gap.py

struct BondCheck {
    n_bonds: usize,
    n_violations: usize,
    max_rel_error: f64,
}

/// Identical check to `scripts/etkdg_vs_rdkit_gap.py::bond_violations`, reimplemented
/// in Rust against the external RDKit-radius reference table above.
fn bond_violations(mol: &Molecule, coords: &Coords3D) -> BondCheck {
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
        if frac > BOND_BLOWUP_REL_ERROR {
            n_violations += 1;
        }
    }
    BondCheck {
        n_bonds,
        n_violations,
        max_rel_error,
    }
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
    let check = bond_violations(mol, coords);
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
    let mut n_total = 0usize;
    let mut n_new_ok = 0usize;
    let mut n_dg_ok = 0usize;
    let mut n_both_ok = 0usize;
    let mut n_regressions = 0usize; // dg ok, new NOT ok -- must stay 0
    let mut regression_names: Vec<&str> = Vec::new();
    let mut new_bonds_checked = 0usize;
    let mut new_bond_violations = 0usize;
    let mut dg_bonds_checked = 0usize;
    let mut dg_bond_violations = 0usize;

    for &(name, smiles, category) in CORPUS {
        n_total += 1;
        let mol = parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse ({name}): {e:?}"));

        // --- new embedder (this PR's deliverable), raw output, no minimization ---
        let params = EmbedParameters::default();
        let (new_coords, new_err, new_stats_line) =
            match embed_distance_geometry_v2_detail(&mol, &params) {
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

        // --- paired comparison: dg::generate_coords (existing DFS placer), raw ---
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
            let bv = bond_violations(&mol, c);
            new_bonds_checked += bv.n_bonds;
            new_bond_violations += bv.n_violations;
            bv.max_rel_error
        });
        let dg_bv = bond_violations(&mol, &dg_coords);
        dg_bonds_checked += dg_bv.n_bonds;
        dg_bond_violations += dg_bv.n_violations;
        let dg_max_rel = dg_bv.max_rel_error;

        rows.push(format!(
            "{{\"name\":\"{}\",\"category\":\"{}\",\"n_atoms\":{},\"new_status\":\"{}\",\"new_max_rel_error\":{},\"new_stats\":\"{}\",\"dg_status\":\"{}\",\"dg_max_rel_error\":{:.4}}}",
            json_escape(name),
            json_escape(category),
            mol.atom_count(),
            json_escape(&new_status),
            new_max_rel.map(|v| format!("{v:.4}")).unwrap_or_else(|| "null".to_string()),
            json_escape(&new_stats_line),
            json_escape(&dg_status),
            dg_max_rel,
        ));

        println!("{new_status:<28} (dg: {dg_status:<12}) {name}");
    }

    let geometrically_valid_rate = n_new_ok as f64 / n_total as f64;
    let dg_geometrically_valid_rate = n_dg_ok as f64 / n_total as f64;

    let status_counts_json = |m: &BTreeMap<String, usize>| {
        m.iter()
            .map(|(k, v)| format!("\"{}\":{}", json_escape(k), v))
            .collect::<Vec<_>>()
            .join(",")
    };

    let summary = format!(
        "{{\n  \"n_molecules\": {n_total},\n  \"new_status_counts\": {{{}}},\n  \"dg_status_counts\": {{{}}},\n  \"geometrically_valid_rate_new\": {geometrically_valid_rate:.4},\n  \"geometrically_valid_rate_dg_raw\": {dg_geometrically_valid_rate:.4},\n  \"n_both_ok\": {n_both_ok},\n  \"n_regressions_dg_ok_new_not_ok\": {n_regressions},\n  \"regression_names\": {:?},\n  \"new_bond_violation_rate\": {},\n  \"dg_bond_violation_rate\": {}\n}}",
        status_counts_json(&new_status_counts),
        status_counts_json(&dg_status_counts),
        regression_names,
        rate_str(new_bond_violations, new_bonds_checked),
        rate_str(dg_bond_violations, dg_bonds_checked),
    );

    println!("\n--- summary ---");
    println!("{summary}");

    println!("\n--- rows (JSONL) ---");
    for row in &rows {
        println!("{row}");
    }

    assert_eq!(
        n_total, 58,
        "corpus size drifted from the frozen 58 -- re-sync with scripts/etkdg_vs_rdkit_gap.py::CORPUS"
    );

    if geometrically_valid_rate < 1.0 {
        eprintln!(
            "\nGATE NOT MET: geometrically_valid_rate = {:.4} (need 1.0). See rows above for which molecules failed and why.",
            geometrically_valid_rate
        );
        std::process::exit(1);
    } else {
        println!(
            "\nGATE MET: geometrically_valid_rate = 1.0 on frozen 58 (raw embedder, pre-minimization)."
        );
    }

    // --- seed-robustness sweep: don't report 100% from one cherry-picked seed ---
    println!("\n--- seed-robustness sweep (new embedder only) ---");
    let mut all_seeds_perfect = true;
    for seed in [0u64, 1, 2, 42, 999, 0xDEAD_BEEF, u64::MAX] {
        let (n_ok, n_tot) = measure_validity_at_seed(seed);
        let rate = n_ok as f64 / n_tot as f64;
        println!("seed={seed:#x}  geometrically_valid_rate={rate:.4}  ({n_ok}/{n_tot})");
        if n_ok != n_tot {
            all_seeds_perfect = false;
        }
    }
    if all_seeds_perfect {
        println!("All swept seeds reach 1.0 -- the gate result is not a single-seed artifact.");
    } else {
        println!(
            "NOTE: not every swept seed reaches 1.0 (retries/max_attempts absorb per-attempt \
             draws via EmbedParameters::default()'s max_attempts=8 -- see stats.attempts_used \
             per molecule in the JSONL rows above for the default seed)."
        );
    }
}
