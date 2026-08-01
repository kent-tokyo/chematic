//! Phase 1A audit for issue #227 -- DIAGNOSTIC PROBE, KNOWN-WRONG BY DESIGN.
//!
//! This is NOT a proposed fix and `corrected_type` below is NOT a
//! recommended remapping -- do not lift it into production code. It exists
//! only to produce a *negative result* for the audit report: a crude,
//! context-blind 3-value swap (63<->37, 38->64) of chmatic's aromatic-carbon
//! numeric MMFF94 types, applied uniformly regardless of ring context,
//! measurably:
//!   - flips only 17/221 (7.7%) of currently-gate-failing molecules to pass
//!     (nowhere near the 80% bar for a same-session fix), and
//!   - REGRESSES furan (was passing, now fails) -- proof that the real fix
//!     needs full ring-context-aware, C+N+O+S-together re-derivation against
//!     a correct reference (see the diagnosis doc), not a blind constant
//!     swap.
//!
//! Does not modify production code. Run only to regenerate this audit's
//! reported numbers, not as a basis for any merge decision.

use chematic_core::{AtomIdx, Molecule};
use chematic_ff::{
    angle_type_for, assign_mmff94_numeric_types, bond_type_for, mmff94_angle_energy,
    mmff94_bond_energy,
};
use chematic_perception::find_sssr;
use serde_json::Value;

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

/// High-confidence carbon-type correction only (dominant single mismatch
/// cluster from the RDKit oracle diff, 2451/3177 = 77% of all atom-level
/// mismatches). Leaves nitrogen/other mismatches (messier, multi-target,
/// likely a distinct aromaticity-perception issue) untouched.
fn corrected_type(t: u8) -> u8 {
    match t {
        63 => 37, // was: "CB 6-ring" (wrong) -> real CB 6-ring is 37
        37 => 63, // was: "C5A 5-ring alpha" (wrong) -> real C5-ring alpha is 63
        38 => 64, // was: "C5B 5-ring beta" (wrong) -> real C 5-ring beta is 64
        other => other,
    }
}

fn bond_angle_gate_ok(mol: &Molecule, types: &[u8]) -> bool {
    for (_, bond) in mol.bonds() {
        let (a1, a2) = (bond.atom1, bond.atom2);
        let (t1, t2) = (types[a1.0 as usize], types[a2.0 as usize]);
        let bt = bond_type_for(t1, t2, bond.order);
        if mmff94_bond_energy(bt, t1, t2).is_none() {
            return false;
        }
    }
    let rings = find_sssr(mol);
    let rings = rings.rings();
    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let (a, c) = (neighbors[i], neighbors[j]);
                let (ta, tc) = (types[a.0 as usize], types[c.0 as usize]);
                let at = angle_type_for(mol, rings, a.0 as usize, b_idx, c.0 as usize, types);
                if mmff94_angle_energy(at, ta, types[b_idx], tc).is_none() {
                    return false;
                }
            }
        }
    }
    true
}

fn main() {
    let mut corpus: Vec<(String, String, String)> = Vec::new();
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
        let manifest = load_manifest(path);
        for m in manifest["molecules"].as_array().unwrap() {
            corpus.push((
                tier.to_string(),
                m["name"].as_str().unwrap().to_string(),
                m["smiles"].as_str().unwrap().to_string(),
            ));
        }
    }

    let mut n_total = 0;
    let mut n_before_fail = 0;
    let mut n_after_fail = 0;
    let mut n_flipped_to_pass = 0;
    let mut n_before_pass_after_fail_regression = 0; // MUST be 0
    let mut flipped_names: Vec<String> = Vec::new();
    let mut still_failing_names: Vec<String> = Vec::new();

    for (_tier, name, smiles) in &corpus {
        let mol = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(_) => continue,
        };
        let types = match assign_mmff94_numeric_types(&mol) {
            Ok(t) => t,
            Err(_) => continue,
        };
        n_total += 1;

        let before_ok = bond_angle_gate_ok(&mol, &types);
        let corrected: Vec<u8> = types.iter().map(|&t| corrected_type(t)).collect();
        let after_ok = bond_angle_gate_ok(&mol, &corrected);

        if !before_ok {
            n_before_fail += 1;
        }
        if !after_ok {
            n_after_fail += 1;
            if !before_ok {
                still_failing_names.push(name.clone());
            }
        }
        if !before_ok && after_ok {
            n_flipped_to_pass += 1;
            flipped_names.push(name.clone());
        }
        if before_ok && !after_ok {
            n_before_pass_after_fail_regression += 1;
            eprintln!("REGRESSION: {name} was passing, now fails after the C-type correction");
        }
    }

    eprintln!("=== typing-fix simulation (C-type 63<->37, 38->64 remap only) ===");
    eprintln!("total molecules typed: {n_total}");
    eprintln!("bond+angle gate failing BEFORE: {n_before_fail}");
    eprintln!("bond+angle gate failing AFTER:  {n_after_fail}");
    eprintln!(
        "flipped fail->pass: {n_flipped_to_pass} ({:.1}% of before-failures)",
        100.0 * n_flipped_to_pass as f64 / n_before_fail.max(1) as f64
    );
    eprintln!("REGRESSIONS (pass->fail, must be 0): {n_before_pass_after_fail_regression}");
    eprintln!();
    eprintln!(
        "still failing after correction ({}): {:?}",
        still_failing_names.len(),
        still_failing_names
    );
}
