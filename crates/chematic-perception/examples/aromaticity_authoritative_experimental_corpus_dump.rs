//! Corpus-verification tool for the opt-in
//! `apply_aromaticity_authoritative_experimental`/
//! `assign_aromaticity_authoritative_experimental` engine
//! (crates/chematic-perception/src/aromaticity.rs) -- the fused-diazine
//! ring-fusion fix plus bidirectional (promote+demote) atom-flag
//! authoritative demotion, shipped opt-in per coordinator decision (the
//! default `apply_aromaticity`/`apply_aromaticity_ex` stay byte-identical
//! to pre-K2b behavior; see `aromaticity_flag_corpus_dump.rs`, the sibling
//! tool for the default engine).
//!
//! Same shape/self-consistency check as `aromaticity_flag_corpus_dump.rs`
//! (see its own doc comment for the invariant and the two calling
//! conventions) but pointed at the opt-in engine instead, so its corpus-wide
//! behavior stays independently verifiable and diffable across commits.
//!
//! Diagnostic only. Calls only existing public API.
//!
//! Run (default corpus: `scripts/descriptor_census_corpus.smi`):
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_authoritative_experimental_corpus_dump \
//!     -- scripts/descriptor_census_corpus.smi \
//!     > /tmp/authoritative_experimental_corpus_dump.jsonl
//! ```

use std::fs;

use chematic_core::BondOrder;
use chematic_perception::{
    apply_aromaticity_authoritative_experimental, assign_aromaticity_authoritative_experimental,
};
use serde_json::json;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/descriptor_census_corpus.smi".to_string());
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut n_total = 0usize;
    let mut n_parse_fail = 0usize;

    for line in text.lines() {
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        n_total += 1;

        let raw = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                n_parse_fail += 1;
                eprintln!("PARSE_FAIL: {smiles:?}: {e}");
                continue;
            }
        };

        // Same two calling conventions as `aromaticity_flag_corpus_dump.rs`:
        // "raw" (un-Kekulized, as parsed) and "kekulized" (explicit
        // Kekulize-first, matching RDKit's own sanitization order and the
        // 40-fixture oracle's methodology).
        let kekulize_result = chematic_core::kekulize(&raw);
        let kekulize_ok = kekulize_result.is_ok();
        let kek_mol = match &kekulize_result {
            Ok(k) => chematic_core::apply_kekule(&raw, k),
            Err(_) => raw.clone(),
        };

        let model_aromatic_atom_count =
            assign_aromaticity_authoritative_experimental(&raw).aromatic_atom_count();
        let applied_raw = apply_aromaticity_authoritative_experimental(&raw);
        let applied_kekulized = apply_aromaticity_authoritative_experimental(&kek_mol);

        let dump_one = |applied: &chematic_core::Molecule| {
            // Same self-consistency invariant as `aromaticity_flag_corpus_dump.rs`.
            let mut inconsistent_atoms: Vec<u32> = Vec::new();
            for (_, bond) in applied.bonds() {
                if bond.order != BondOrder::Aromatic {
                    continue;
                }
                for a in [bond.atom1, bond.atom2] {
                    if !applied.atom(a).aromatic {
                        inconsistent_atoms.push(a.0);
                    }
                }
            }
            for (idx, atom) in applied.atoms() {
                if !atom.aromatic {
                    continue;
                }
                let has_ring_bond = applied.neighbors(idx).next().is_some();
                let has_aromatic_bond = applied
                    .neighbors(idx)
                    .any(|(_, bidx)| applied.bond(bidx).order == BondOrder::Aromatic);
                if has_ring_bond && !has_aromatic_bond {
                    inconsistent_atoms.push(idx.0);
                }
            }
            inconsistent_atoms.sort_unstable();
            inconsistent_atoms.dedup();

            let atoms: Vec<_> = applied
                .atoms()
                .map(|(idx, atom)| json!({"idx": idx.0, "aromatic": atom.aromatic}))
                .collect();
            let bonds: Vec<_> = applied
                .bonds()
                .map(|(idx, bond)| {
                    json!({
                        "idx": idx.0,
                        "a1": bond.atom1.0,
                        "a2": bond.atom2.0,
                        "aromatic": bond.order == BondOrder::Aromatic,
                    })
                })
                .collect();
            json!({
                "consistent": inconsistent_atoms.is_empty(),
                "inconsistent_atom_idxs": inconsistent_atoms,
                "atoms": atoms,
                "bonds": bonds,
            })
        };

        let row = json!({
            "smiles": smiles,
            "kekulize_ok": kekulize_ok,
            "model_aromatic_atom_count": model_aromatic_atom_count,
            "raw": dump_one(&applied_raw),
            "kekulized": dump_one(&applied_kekulized),
        });
        println!("{row}");
    }

    eprintln!(
        "dumped {} molecules, {n_parse_fail} parse failures",
        n_total - n_parse_fail
    );
}
