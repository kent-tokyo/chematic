//! Measures MMFF94 parameter coverage (bond/angle/torsion/out-of-plane) and
//! the mechanism-3 "silent zero-gradient" fix, over the frozen 58-molecule
//! corpus from `scripts/etkdg_vs_rdkit_gap.py::CORPUS` (`docs/rfcs/etkdg_3d_gap_rfc.md`).
//!
//! This is a self-contained validation entry point owned by this PR (Agent
//! F, `feat/3d-full-force-field-bridge`) — it does not write into
//! `validation/` (that directory is Agent H's territory per
//! `docs/rfcs/3d_breakthrough_master_plan.md` §3) and does not touch the shared
//! repo `.venv`; it is a pure-Rust `cargo run --example` tool. The SMILES
//! list below is hand-copied from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`
//! (58 entries) so this can run without Python at all.
//!
//! For each molecule this reports:
//! - old path: `minimize_mmff94` (existing, crippled, silently-zeroing energy
//!   function) worst-bond-length after minimization, starting from
//!   `dg::generate_coords`.
//! - new path: `ForceFieldPolicy::Mmff94BondAngleStrict` coverage (bond/angle/
//!   torsion/oop missing counts, cited by atom element pairs) and, via
//!   `Mmff94WithUffFallback`, the fixed worst-bond-length.
//!
//! Usage: `cargo run -p chematic-3d --release --example mmff94_bridge_coverage_report`

use chematic_3d::dg::generate_coords;
use chematic_3d::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_mmff94, minimize_with_policy,
    minimize_with_policy_gated,
};
use chematic_core::Molecule;
use chematic_smiles::parse;

/// Hand-copied verbatim from `scripts/etkdg_vs_rdkit_gap.py::CORPUS` (58
/// entries, name/SMILES/category — category dropped, not needed here).
const CORPUS: &[(&str, &str)] = &[
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("pyridine", "c1ccncc1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("adamantane", "C1CC2CC3CC1CC(C2)C3"),
    ("cubane", "C1C2C3C1C4C2C3C4"),
    ("cyclohexane", "C1CCCCC1"),
    ("cyclopentane", "C1CCCC1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("purine", "c1ncc2[nH]cnc2n1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("butane", "CCCC"),
    ("hexane", "CCCCCC"),
    ("decane", "CCCCCCCCCC"),
    ("triethylene_glycol", "OCCOCCOCCO"),
    ("hexanediol", "OCCCCCCO"),
    ("hexadecane", "CCCCCCCCCCCCCCCC"),
    ("cyclododecane", "C1CCCCCCCCCCC1"),
    ("crown_12_4", "O1CCOCCOCCOCC1"),
    ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1"),
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
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
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

fn worst_bond(mol: &Molecule, coords: &chematic_3d::Coords3D) -> f64 {
    let mut worst = 0.0_f64;
    for (_, bond) in mol.bonds() {
        let d = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        if d > worst {
            worst = d;
        }
    }
    worst
}

/// Aggregation key: `Mmff94MissingTerm::description` with the trailing
/// per-molecule "(atom indices [...])" suffix stripped -- NOT a hand-rolled
/// re-derivation from element symbols. `description` already carries each
/// atom's specific numeric MMFF94 type (not just its element symbol -- many
/// distinct MMFF94 types share an element) and is built from
/// `canonicalize_term_atoms`-ordered atoms, so physically-equivalent
/// citations (a bond read C->N vs N->C, etc.) already collapse to the same
/// string -- reusing it here means this aggregation can't silently drift
/// from the bridge's own canonicalization the way a second, independent
/// implementation could.
fn pattern_key(term: &chematic_3d::minimize::Mmff94MissingTerm) -> String {
    match term.description.find(" (atom indices") {
        Some(pos) => term.description[..pos].to_string(),
        None => term.description.clone(),
    }
}

fn main() {
    let config = MinimizeConfig::default();

    let mut n_total = 0usize;
    let mut n_parsed = 0usize;
    let mut n_strict_ok = 0usize;
    let mut n_strict_missing = 0usize;
    // Same gate, widened to also require torsion+oop coverage -- measures
    // whether `Mmff94BondAngleStrict`'s current bond+angle-only scope is still the
    // right call post-chematic-ff-#183, or whether torsion/oop coverage is
    // now good enough that widening the gate (or renaming this policy to be
    // honest about a narrower scope) should be reconsidered. Separate
    // denominator from `n_strict_ok`/`n_strict_missing` -- never pooled.
    let mut n_strict_gated_ok = 0usize;
    let mut n_strict_gated_missing = 0usize;
    let mut n_old_blowup = 0usize; // old minimize_mmff94 worst bond > 3 A
    let mut n_fixed_via_fallback = 0usize; // was blown up, now <= 3 A (Ok, sound)
    let mut n_still_blown_up = 0usize; // was blown up, still > 3 A after fallback (Ok>3A or typed Err)
    let mut n_new_regression = 0usize; // was fine (<=3 A), fallback made it > 3 A (Ok>3A or typed Err)
    let mut n_new_regression_typed_err = 0usize; // ...of which now a typed Err(MinimizationFailed)
    // Bug check: an Ok result must NEVER have a blown-up (>3 A) worst bond --
    // that would mean check_minimization_soundness has a gap and a blown-up
    // geometry is silently escaping as "success". MUST be 0.
    let mut n_ok_but_blown_up_bug = 0usize;
    // Disclosure, not a bug check: Ok results that are sound (no soundness-gate
    // trigger) but simply didn't converge within the default iteration budget
    // -- expected to be nonzero (see check_minimization_soundness's doc).
    let mut n_fallback_ok_not_converged = 0usize;

    // pattern -> occurrence count across the whole corpus.
    let mut missing_bond_patterns: std::collections::BTreeMap<String, usize> = Default::default();
    let mut missing_angle_patterns: std::collections::BTreeMap<String, usize> = Default::default();
    let mut missing_torsion_patterns: std::collections::BTreeMap<String, usize> =
        Default::default();
    let mut missing_oop_patterns: std::collections::BTreeMap<String, usize> = Default::default();

    println!(
        "{:<24} {:>4} {:>4} {:>4} {:>4}  {:>10} {:>10} {:>6} {:>12}  strict",
        "molecule", "bond", "angl", "tors", "oop", "old_worst", "new_worst", "conv?", "max_resid"
    );

    for &(name, smiles) in CORPUS {
        n_total += 1;
        let mol = match parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                println!("{name:<24} PARSE ERROR: {e}");
                continue;
            }
        };
        n_parsed += 1;

        let coords = generate_coords(&mol);

        // Old (crippled, silently-zeroing) path.
        let old_result = minimize_mmff94(&mol, coords.clone());
        let old_worst = worst_bond(&mol, &old_result);
        if old_worst > 3.0 {
            n_old_blowup += 1;
        }

        // New: strict coverage gate.
        let strict = minimize_with_policy(
            &mol,
            coords.clone(),
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        );
        let (n_bond_miss, n_angle_miss, n_tors_miss, n_oop_miss, strict_label);
        match &strict {
            Ok(r) => {
                n_strict_ok += 1;
                let c = r
                    .coverage
                    .as_ref()
                    .expect("Mmff94BondAngleStrict always reports coverage");
                n_bond_miss = c.bonds_missing.len();
                n_angle_miss = c.angles_missing.len();
                n_tors_miss = c.torsions_missing.len();
                n_oop_miss = c.oop_missing.len();
                strict_label = "OK".to_string();
            }
            Err(ForceFieldBridgeError::MissingParameters(c)) => {
                n_strict_missing += 1;
                n_bond_miss = c.bonds_missing.len();
                n_angle_miss = c.angles_missing.len();
                n_tors_miss = c.torsions_missing.len();
                n_oop_miss = c.oop_missing.len();
                strict_label = "MISSING".to_string();
                for t in &c.bonds_missing {
                    *missing_bond_patterns.entry(pattern_key(t)).or_insert(0) += 1;
                }
                for t in &c.angles_missing {
                    *missing_angle_patterns.entry(pattern_key(t)).or_insert(0) += 1;
                }
                for t in &c.torsions_missing {
                    *missing_torsion_patterns.entry(pattern_key(t)).or_insert(0) += 1;
                }
                for t in &c.oop_missing {
                    *missing_oop_patterns.entry(pattern_key(t)).or_insert(0) += 1;
                }
            }
            Err(ForceFieldBridgeError::UnsupportedAtomType(e)) => {
                n_strict_missing += 1;
                n_bond_miss = 0;
                n_angle_miss = 0;
                n_tors_miss = 0;
                n_oop_miss = 0;
                strict_label = format!("UNSUPPORTED({e})");
            }
            Err(ForceFieldBridgeError::MinimizationFailed(d)) => {
                // The strict MMFF94 attempt itself produced an unsound
                // geometry despite full bond+angle coverage -- distinct from
                // "missing parameters", still counted in the same
                // denominator since it's still "did not return Ok".
                n_strict_missing += 1;
                n_bond_miss = 0;
                n_angle_miss = 0;
                n_tors_miss = 0;
                n_oop_miss = 0;
                strict_label = format!("UNSOUND({:?})", d.reason);
            }
        }

        // Same gate widened to torsion+oop too -- feeds the Mmff94BondAngleStrict
        // naming/scope decision (see PR body). Separate denominator, not
        // pooled with the bond+angle-only n_strict_ok/n_strict_missing above.
        match minimize_with_policy_gated(
            &mol,
            coords.clone(),
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
            true,
            false,
        ) {
            Ok(_) => n_strict_gated_ok += 1,
            Err(_) => n_strict_gated_missing += 1,
        }

        // New: fallback path -- NOT infallible (see ForceFieldPolicy::Mmff94WithUffFallback's
        // doc): the UFF fallback itself can be unsound, in which case this
        // is now a typed Err(MinimizationFailed) instead of the
        // Ok(converged=false) it used to be. new_worst/new_conv/new_resid
        // are pulled from whichever side (Ok or the typed Err's detail)
        // actually ran, so a fixed-but-still-blown molecule still reports
        // its real geometry, not a blank.
        let fallback_result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
        );
        let (new_worst, new_conv, new_resid, fallback_label);
        match &fallback_result {
            Ok(r) => {
                new_worst = worst_bond(&mol, &r.coords);
                new_conv = r.converged;
                new_resid = r.max_residual_force;
                fallback_label = "OK".to_string();
                if new_worst > 3.0 {
                    // Must never happen: check_minimization_soundness should
                    // have caught this and returned Err instead.
                    n_ok_but_blown_up_bug += 1;
                }
                if !new_conv {
                    n_fallback_ok_not_converged += 1;
                }
            }
            Err(ForceFieldBridgeError::MinimizationFailed(d)) => {
                new_worst = d.worst_bond_length;
                new_conv = d.converged;
                new_resid = d.max_residual_force;
                fallback_label = format!("TYPED_FAIL({:?})", d.reason);
            }
            Err(other) => {
                // Structurally shouldn't happen (Mmff94WithUffFallback's
                // only own failure mode is MinimizationFailed), but don't
                // let an unanticipated variant silently vanish from the
                // blow-up counting below.
                new_worst = f64::INFINITY;
                new_conv = false;
                new_resid = f64::INFINITY;
                fallback_label = format!("UNEXPECTED_ERR({other})");
            }
        }
        if old_worst > 3.0 && new_worst <= 3.0 {
            n_fixed_via_fallback += 1;
        } else if old_worst > 3.0 && new_worst > 3.0 {
            n_still_blown_up += 1;
        } else if old_worst <= 3.0 && new_worst > 3.0 {
            n_new_regression += 1;
            if fallback_result.is_err() {
                n_new_regression_typed_err += 1;
            }
        }

        println!(
            "{name:<24} {n_bond_miss:>4} {n_angle_miss:>4} {n_tors_miss:>4} {n_oop_miss:>4}  {old_worst:>10.2} {new_worst:>10.2} {:>6} {:>12.2}  {strict_label} | fallback={fallback_label}",
            new_conv, new_resid
        );
    }

    println!();
    println!("=== summary ===");
    println!("total corpus molecules:                       {n_total}");
    println!("parsed:                                        {n_parsed}");
    println!("Mmff94BondAngleStrict (bond+angle gate) fully covered (OK):   {n_strict_ok}");
    println!("Mmff94BondAngleStrict (bond+angle gate) missing params:        {n_strict_missing}");
    println!(
        "Mmff94BondAngleStrict widened (+torsion+oop gate) fully covered (OK): {n_strict_gated_ok}"
    );
    println!(
        "Mmff94BondAngleStrict widened (+torsion+oop gate) missing params:      {n_strict_gated_missing}"
    );
    println!("old minimize_mmff94 worst>3A (blown up):        {n_old_blowup}");
    println!("  ...of which FIXED by Mmff94WithUffFallback:  {n_fixed_via_fallback}");
    println!("  ...of which STILL blown up (Ok>3A or typed Err):  {n_still_blown_up}");
    println!(
        "NEW regressions (old<=3A, fallback made it >3A -- Ok>3A or typed Err):    {n_new_regression}"
    );
    println!(
        "  ...of which now surfaced as a typed Err(MinimizationFailed):  {n_new_regression_typed_err} \
         (should equal {n_new_regression} -- an Ok>3A here would mean check_minimization_soundness \
         has a gap)"
    );
    println!(
        "Ok result with worst bond >3A (soundness-gate bug check, MUST be 0): {n_ok_but_blown_up_bug}"
    );
    println!(
        "Ok fallback results with converged=false (sound geometry, just didn't fully converge \
         within the default iteration budget -- expected, not a bug): {n_fallback_ok_not_converged}"
    );
    println!();
    println!(
        "=== distinct missing bond element-type patterns (pattern: occurrence count across corpus) ==="
    );
    for (pattern, count) in &missing_bond_patterns {
        println!("  {pattern}: {count}");
    }
    println!("=== distinct missing angle element-type patterns ===");
    for (pattern, count) in &missing_angle_patterns {
        println!("  {pattern}: {count}");
    }
    println!(
        "=== distinct missing torsion element-type patterns ({} distinct) ===",
        missing_torsion_patterns.len()
    );
    for (pattern, count) in &missing_torsion_patterns {
        println!("  {pattern}: {count}");
    }
    println!("=== distinct missing oop element-type patterns ===");
    for (pattern, count) in &missing_oop_patterns {
        println!("  {pattern}: {count}");
    }
}
