//! Phase 1A audit for issue #227: faithful re-measurement of
//! `ForceFieldPolicy::Mmff94BondAngleStrict` over the Wave 1 265-molecule
//! corpus, using the exact same production entry point
//! (`minimize::minimize_with_policy`) and starting geometry
//! (`dg::generate_coords`) as PR #226's methodology -- NOT the simplified
//! bond/angle-only `Some`/`None` check `mmff94_term_coverage_audit.rs` uses,
//! so the "216 -> N" headline number is measured the same way the issue was
//! originally filed, on today's main.
//!
//! Per-molecule JSONL rows (stdout) additionally carry every field already on
//! [`chematic_3d::minimize::MinimizationFailureDetail`] for `MinimizationFailed`
//! rows (`reason` -- `CatastrophicBondBlowup`/`ExcessiveResidualForce`/
//! `NonFiniteCoordinates` -- `converged`, `iterations`, `max_residual_force`,
//! `worst_bond_length`, `distance_geometry_v2_retry_attempted`) so Priority 3's
//! diagnosis reads directly off this file; no separate classifier needed, the
//! production soundness check already computed and attached this evidence.
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example mmff94_strict_gate_remeasure_227 \
//!   > validation/results/mmff94_strict_gate_remeasure_227_rows.jsonl \
//!   2> validation/results/mmff94_strict_gate_remeasure_227_stderr.log
//! ```

use chematic_3d::dg::generate_coords;
use chematic_3d::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_with_policy,
};
use serde_json::{Value, json};

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

fn main() {
    let config = MinimizeConfig::default();
    let mut n_total = 0;
    let mut n_ok = 0;
    let mut n_missing_params = 0;
    let mut n_unsupported_atom_type = 0;
    let mut n_minimization_failed = 0;
    let mut n_parse_or_other = 0;

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
            let name = m["name"].as_str().unwrap_or("<unnamed>");
            let smiles = m["smiles"].as_str().unwrap();
            let mol = match chematic_smiles::parse(smiles) {
                Ok(m) => m,
                Err(e) => {
                    n_total += 1;
                    n_parse_or_other += 1;
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles,
                               "status": "ParseError", "detail": e.to_string()})
                    );
                    continue;
                }
            };
            n_total += 1;
            let coords = generate_coords(&mol);
            let row = match minimize_with_policy(
                &mol,
                coords,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                &config,
            ) {
                Ok(_) => {
                    n_ok += 1;
                    json!({"tier": tier, "name": name, "smiles": smiles, "status": "Ok"})
                }
                Err(ForceFieldBridgeError::MissingParameters(r)) => {
                    n_missing_params += 1;
                    json!({"tier": tier, "name": name, "smiles": smiles,
                           "status": "MissingParameters",
                           "bonds_missing": r.bonds_missing.len(),
                           "angles_missing": r.angles_missing.len(),
                           "torsions_missing": r.torsions_missing.len(),
                           "oop_missing": r.oop_missing.len(),
                           "stretch_bend_missing": r.stretch_bend_missing.len()})
                }
                Err(ForceFieldBridgeError::UnsupportedAtomType(e)) => {
                    n_unsupported_atom_type += 1;
                    json!({"tier": tier, "name": name, "smiles": smiles,
                           "status": "UnsupportedAtomType", "detail": e.to_string()})
                }
                Err(ForceFieldBridgeError::MinimizationFailed(d)) => {
                    n_minimization_failed += 1;
                    json!({"tier": tier, "name": name, "smiles": smiles,
                           "status": "MinimizationFailed",
                           "reason": format!("{:?}", d.reason),
                           "converged": d.converged,
                           "iterations": d.iterations,
                           "max_residual_force": d.max_residual_force,
                           "worst_bond_length": d.worst_bond_length,
                           "distance_geometry_v2_retry_attempted":
                               d.distance_geometry_v2_retry_attempted})
                }
            };
            println!("{row}");
        }
    }

    eprintln!("=== issue #227 faithful re-measurement (production minimize_with_policy) ===");
    eprintln!("total: {n_total}");
    eprintln!("Ok (success): {n_ok}");
    eprintln!("Err(MissingParameters) [\"unsupported\" in the original issue]: {n_missing_params}");
    eprintln!("Err(UnsupportedAtomType): {n_unsupported_atom_type}");
    eprintln!(
        "Err(MinimizationFailed) [\"typed_failure\" in the original issue]: {n_minimization_failed}"
    );
    eprintln!("parse/other: {n_parse_or_other}");
}
