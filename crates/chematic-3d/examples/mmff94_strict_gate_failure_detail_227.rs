//! Phase 1B-0 audit for issue #227: per-molecule detail for
//! `ForceFieldPolicy::Mmff94BondAngleStrict`'s `UnsupportedAtomType` and
//! `MinimizationFailed` failure modes (new post-fix categories worth
//! investigating individually, not just counting).
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_strict_gate_failure_detail_227`

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
            let name = m["name"].as_str().unwrap();
            let smiles = m["smiles"].as_str().unwrap();
            let mol = match chematic_smiles::parse(smiles) {
                Ok(m) => m,
                Err(_) => continue,
            };
            let coords = generate_coords(&mol);
            match minimize_with_policy(
                &mol,
                coords,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                &config,
            ) {
                Err(ForceFieldBridgeError::UnsupportedAtomType(e)) => {
                    println!("[{tier}] {name} ({smiles}): UnsupportedAtomType: {e}");
                }
                Err(ForceFieldBridgeError::MinimizationFailed(d)) => {
                    println!(
                        "[{tier}] {name} ({smiles}): MinimizationFailed: reason={:?} converged={} worst_bond={:.3}",
                        d.reason, d.converged, d.worst_bond_length
                    );
                }
                _ => {}
            }
        }
    }
}
