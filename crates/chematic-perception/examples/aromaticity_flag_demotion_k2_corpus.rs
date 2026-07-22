//! fix/aromaticity-flag-demotion-k2: full-corpus dump of `apply_aromaticity`'s
//! per-atom/per-bond aromatic flags plus a self-consistency check ("no
//! `BondOrder::Aromatic` bond has a non-aromatic endpoint atom"), for a
//! before/after comparison of this PR's fix to `build_molecule_from_model`
//! (crates/chematic-perception/src/aromaticity.rs).
//!
//! Diagnostic only. Calls only existing public API
//! (`chematic_smiles::parse`, `chematic_core::kekulize`/`apply_kekule`,
//! `chematic_perception::{apply_aromaticity, assign_aromaticity}`) exactly as
//! an external caller would -- this file itself is unchanged by the fix (it
//! only calls the library), so it can be run byte-identically against the
//! pre-fix and post-fix crate to produce the two sides of the diff described
//! in the PR: run it once on `main` (or any commit before this fix), once on
//! this branch, and diff the two JSONL dumps by `smiles`.
//!
//! Run (default corpus: `scripts/descriptor_census_corpus.smi`, reused from
//! PR #137's descriptor census -- not rebuilt for this PR):
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_flag_demotion_k2_corpus \
//!     -- scripts/descriptor_census_corpus.smi \
//!     > /tmp/k2_corpus_after.jsonl
//! ```

use std::fs;

use chematic_core::BondOrder;
use chematic_perception::{apply_aromaticity, assign_aromaticity};
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

        // Two calling conventions, both documented as supported by
        // `apply_aromaticity`'s own doc comment ("may be kekulized... or may
        // retain Aromatic bond orders from the SMILES parser"):
        //   - "raw": call `apply_aromaticity` directly on the freshly-parsed
        //     molecule (bonds may still be `BondOrder::Aromatic`). When the
        //     model can't independently confirm a ring that's still in this
        //     representation, `build_molecule_from_model`'s bond-consistency
        //     fallback (see its doc comment) keeps that ring's atoms/bonds
        //     aromatic together rather than demoting the atom against a
        //     bond order it has no independently-verified replacement for.
        //   - "kekulized": Kekulize first (mirroring RDKit's own
        //     `Kekulize`-before-`setAromaticity` sanitization order), then
        //     call `apply_aromaticity`. Every bond is then a real
        //     Single/Double value, so the model's demotion is never masked
        //     by that fallback -- this is the pathway the 40-fixture oracle
        //     (`aromaticity_rdkit_parity_dump.rs`) itself uses.
        let kekulize_result = chematic_core::kekulize(&raw);
        let kekulize_ok = kekulize_result.is_ok();
        let kek_mol = match &kekulize_result {
            Ok(k) => chematic_core::apply_kekule(&raw, k),
            Err(_) => raw.clone(),
        };

        let model_aromatic_atom_count = assign_aromaticity(&raw).aromatic_atom_count();
        let applied_raw = apply_aromaticity(&raw);
        let applied_kekulized = apply_aromaticity(&kek_mol);

        let dump_one = |applied: &chematic_core::Molecule| {
            // Self-consistency invariant, checked in BOTH directions -- this
            // is the precise shape of the pre-K2 bug (RFC §1b): an
            // `Aromatic`-order bond may not have a non-aromatic endpoint
            // atom (direction the experimental engine's own
            // `validate_aromaticity_invariants` already checks), AND an
            // `aromatic: true` ring atom may not have every one of its ring
            // bonds be a non-aromatic order (that second direction is
            // exactly what selenophene/azulene looked like before this fix:
            // atom.aromatic=true sitting on all-Single/Double ring bonds).
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
