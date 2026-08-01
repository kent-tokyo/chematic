//! Phase 1A audit for issue #227: dumps chematic's `assign_mmff94_numeric_types`
//! output for every heavy atom of the 265-molecule Wave 1 corpus, for
//! cross-checking against an RDKit oracle (`scripts/mmff94_rdkit_type_oracle.py`)
//! atom-index-aligned (both parsers preserve SMILES left-to-right atom order,
//! per PR #226's already-verified 265/265 atom mapping).
//!
//! Run: `cargo run --release -p chematic-ff --example mmff94_numeric_type_dump \
//!   > validation/results/mmff94_chematic_numeric_types.jsonl`

use chematic_ff::assign_mmff94_numeric_types;
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
            match assign_mmff94_numeric_types(&mol) {
                Ok(types) => {
                    let atoms: Vec<Value> = (0..mol.atom_count())
                        .map(|i| {
                            let a = mol.atom(chematic_core::AtomIdx(i as u32));
                            json!({
                                "index": i,
                                "element": a.element.symbol(),
                                "aromatic": a.aromatic,
                                "chematic_numeric_type": types[i],
                            })
                        })
                        .collect();
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles, "status": "ok", "atoms": atoms})
                    );
                }
                Err(e) => {
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles, "status": "typing_failure", "error": e.to_string()})
                    );
                }
            }
        }
    }
}
