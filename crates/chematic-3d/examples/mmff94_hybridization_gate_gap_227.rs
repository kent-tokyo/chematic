//! Issue #227 Phase 0.3: dumps chematic's `total_degree(atom) > 3` gate
//! decision for every ring C/N heavy atom of the 265-molecule Wave 1 corpus,
//! for cross-checking against RDKit's real `atom->getHybridization() !=
//! Atom::SP2` gate (`Code/GraphMol/Aromaticity.cpp` line 1023, pinned commit
//! -- see `scripts/mmff94_provenance/PROVENANCE.md`) via
//! `scripts/mmff94_hybridization_gate_gap_227.py`. Atom-index-aligned with
//! RDKit's own parser per PR #226's already-verified 265/265 mapping.
//!
//! Mirrors `compute_mmff94_aromatic_view`'s own Kekulization step exactly
//! (same `kekulize`/`apply_kekule` sequence) so `total_degree` is computed
//! on the same molecule state the real gate uses -- Kekulization does not
//! change bond *count* but can change implicit-H count via the valence
//! model, so this must not be measured on the pre-Kekulization molecule.
//!
//! Run: `cargo run --release -p chematic-ff --example \
//!   mmff94_hybridization_gate_gap_227 \
//!   > validation/results/mmff94_hybridization_gate_gap_227_chematic.jsonl`

use chematic_core::{AtomIdx, Element, implicit_hcount};
use serde_json::{Value, json};

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

fn total_degree(mol: &chematic_core::Molecule, idx: AtomIdx) -> usize {
    mol.neighbors(idx).count() + implicit_hcount(mol, idx) as usize
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

            let rings = chematic_perception::find_sssr(&mol).rings().to_vec();
            if rings.is_empty() {
                println!(
                    "{}",
                    json!({"tier": tier, "name": name, "smiles": smiles, "status": "ok", "atoms": []})
                );
                continue;
            }

            let kmol = match chematic_core::kekulize(&mol) {
                Ok(kek) if kek.is_empty() => mol.clone(),
                Ok(kek) => chematic_core::apply_kekule(&mol, &kek),
                Err(e) => {
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles, "status": "kekulize_failure", "error": e.detail})
                    );
                    continue;
                }
            };

            let ring_atoms: std::collections::HashSet<u32> =
                rings.iter().flatten().map(|a| a.0).collect();

            let atoms: Vec<Value> = ring_atoms
                .iter()
                .copied()
                .filter(|&i| {
                    let e = kmol.atom(AtomIdx(i)).element;
                    e == Element::C || e == Element::N
                })
                .map(|i| {
                    let idx = AtomIdx(i);
                    let td = total_degree(&kmol, idx);
                    json!({
                        "index": i,
                        "element": kmol.atom(idx).element.symbol(),
                        "total_degree": td,
                        "gate_fires_reject": td > 3,
                    })
                })
                .collect();

            println!(
                "{}",
                json!({"tier": tier, "name": name, "smiles": smiles, "status": "ok", "atoms": atoms})
            );
        }
    }
}
