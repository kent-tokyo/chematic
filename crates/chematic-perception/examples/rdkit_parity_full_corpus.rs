//! Aromaticity-A1-1b-0: run `rdkit_parity_aromaticity` against a full SMILES
//! corpus (one SMILES per line/first CSV column) and emit per-atom and
//! per-bond rows, for a real set-level join against real RDKit (not just a
//! count comparison, which is blind to same-count/different-atoms
//! mismatches) at the same scale as this project's existing 99.44%/98.82%
//! baseline (`docs/rdkit_compat.md`).
//!
//! Diagnostic only. Run (requires the `diagnostics` feature):
//! ```text
//! cargo run -p chematic-perception --release --features diagnostics \
//!     --example rdkit_parity_full_corpus \
//!     -- ~/Downloads/SMILES.csv \
//!     > validation/results/rdkit_parity_full_corpus_trace.jsonl
//! ```

use std::fs;

use chematic_perception::diagnostics::rdkit_parity_aromaticity;
use serde_json::json;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: rdkit_parity_full_corpus <smiles.csv>");
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut n_ok = 0usize;
    let mut n_fail = 0usize;
    let mut fail_smiles: Vec<String> = Vec::new();

    for line in content.lines() {
        let smi = line.split(',').next().unwrap_or("").trim();
        if smi.is_empty() || smi.eq_ignore_ascii_case("smiles") {
            continue;
        }
        let raw = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(_) => {
                n_fail += 1;
                fail_smiles.push(smi.to_string());
                continue;
            }
        };
        let mol = match chematic_core::kekulize(&raw) {
            Ok(k) => chematic_core::apply_kekule(&raw, &k),
            Err(_) => {
                n_fail += 1;
                fail_smiles.push(smi.to_string());
                continue;
            }
        };
        let (atoms, bonds) = rdkit_parity_aromaticity(&mol);

        for (idx, _atom) in mol.atoms() {
            let row = json!({
                "kind": "atom",
                "smiles": smi,
                "atom_idx": idx.0,
                "rdkit_parity_atom_aromatic": atoms.contains(&idx),
            });
            println!("{row}");
        }
        for (bidx, bond) in mol.bonds() {
            let (a1, a2) = (
                bond.atom1.0.min(bond.atom2.0),
                bond.atom1.0.max(bond.atom2.0),
            );
            let row = json!({
                "kind": "bond",
                "smiles": smi,
                "bond_atoms": [a1, a2],
                "rdkit_parity_bond_aromatic": bonds.contains(&bidx),
            });
            println!("{row}");
        }
        n_ok += 1;
    }

    eprintln!("processed {n_ok} molecules, {n_fail} parse/kekulize failures");
    for smi in &fail_smiles {
        eprintln!("FAILED: {smi}");
    }
}
