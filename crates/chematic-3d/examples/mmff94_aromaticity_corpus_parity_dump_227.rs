//! Issue #227 Phase 0.4: dumps chematic's post-Phase-0.1/0.2/0.3
//! `compute_mmff94_aromatic_view` atom + bond aromaticity flags for every
//! heavy atom/bond of the 265-molecule Wave 1 corpus, for cross-checking
//! against `scripts/mmff94_aromaticity_corpus_parity_227.py`'s RDKit
//! oracle. Atom-index-aligned per PR #226's already-verified mapping.
//!
//! Run: `cargo run --release -p chematic-3d --example \
//!   mmff94_aromaticity_corpus_parity_dump_227 \
//!   > validation/results/mmff94_aromaticity_corpus_parity_227_chematic.jsonl`

use chematic_core::AtomIdx;
use chematic_ff::mmff94_numeric::compute_mmff94_aromatic_view;
use serde_json::{Value, json};

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

fn main() {
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
        for m in manifest["molecules"].as_array().expect("molecules array") {
            let name = m["name"].as_str().unwrap();
            let smiles = m["smiles"].as_str().unwrap();
            let mol = match chematic_smiles::parse(smiles) {
                Ok(mol) => mol,
                Err(e) => {
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles, "status": "parse_failure", "error": e.to_string()})
                    );
                    continue;
                }
            };

            let rings = chematic_perception::find_symmetrized_sssr(&mol)
                .rings()
                .to_vec();
            let view = match compute_mmff94_aromatic_view(&mol, &rings) {
                Ok(v) => v,
                Err(e) => {
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles, "status": "reperception_failure", "error": e.to_string()})
                    );
                    continue;
                }
            };

            let n = view.atom_count();
            let mut atom_aromatic = serde_json::Map::new();
            for i in 0..n {
                atom_aromatic.insert(
                    i.to_string(),
                    Value::Bool(view.atom(AtomIdx(i as u32)).aromatic),
                );
            }
            let mut bond_aromatic = serde_json::Map::new();
            for (_, b) in view.bonds() {
                let (a1, a2) = (b.atom1.0, b.atom2.0);
                let key = format!("{}-{}", a1.min(a2), a1.max(a2));
                bond_aromatic.insert(
                    key,
                    Value::Bool(b.order == chematic_core::BondOrder::Aromatic),
                );
            }

            println!(
                "{}",
                json!({
                    "tier": tier,
                    "name": name,
                    "smiles": smiles,
                    "status": "ok",
                    "n_heavy_atoms": n,
                    "atom_aromatic": atom_aromatic,
                    "bond_aromatic": bond_aromatic,
                })
            );
        }
    }
}
