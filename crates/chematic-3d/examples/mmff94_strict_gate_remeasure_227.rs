//! Phase 1A audit for issue #227: faithful re-measurement of
//! `ForceFieldPolicy::Mmff94BondAngleStrict` over the Wave 1 265-molecule
//! corpus, using the exact same production entry point
//! (`minimize::minimize_with_policy`) and starting geometry
//! (`dg::generate_coords`) as PR #226's methodology -- NOT the simplified
//! bond/angle-only `Some`/`None` check `mmff94_term_coverage_audit.rs` uses,
//! so the "216 -> N" headline number is measured the same way the issue was
//! originally filed, on today's main.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_strict_gate_remeasure_227`

use chematic_3d::dg::generate_coords;
use chematic_3d::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_with_policy,
};
use serde_json::Value;

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

    for (_tier, path) in [
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
            let smiles = m["smiles"].as_str().unwrap();
            let mol = match chematic_smiles::parse(smiles) {
                Ok(m) => m,
                Err(_) => {
                    n_total += 1;
                    n_parse_or_other += 1;
                    continue;
                }
            };
            n_total += 1;
            let coords = generate_coords(&mol);
            match minimize_with_policy(
                &mol,
                coords,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                &config,
            ) {
                Ok(_) => n_ok += 1,
                Err(ForceFieldBridgeError::MissingParameters(_)) => n_missing_params += 1,
                Err(ForceFieldBridgeError::UnsupportedAtomType(_)) => n_unsupported_atom_type += 1,
                Err(ForceFieldBridgeError::MinimizationFailed(_)) => n_minimization_failed += 1,
            }
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
