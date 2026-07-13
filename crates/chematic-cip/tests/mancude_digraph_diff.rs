//! Milestone 3B-0: for each of the 98 aromatic/MANCUDE-adjacent corpus cases, build the
//! CURRENT (pre-M3B) `CipDigraph` under multiple chemically-equivalent input
//! representations (aromatic notation, Kekulé respelling, atom renumbering, SMILES
//! respelling) and report the first point where they structurally diverge -- the
//! evidence artifact Milestone 3B-1 designs its production representation against. Zero
//! production behavior change this round.

mod common;

use serde_json::{Value, json};

use chematic_cip::CipBudget;
use chematic_cip::digraph_diff::{
    DigraphDivergence, find_atom_by_signature, first_divergence, renumber_molecule,
};
use chematic_core::AtomIdx;
use chematic_core::kekulization::{apply_kekule, kekulize};

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cip_label_corpus.jsonl"
));

fn divergence_json(result: Result<Option<DigraphDivergence>, chematic_cip::CipError>) -> Value {
    match result {
        Ok(None) => json!({"diverges": false}),
        Ok(Some(d)) => json!({
            "diverges": true,
            "depth": d.depth,
            "left": d.left,
            "right": d.right,
        }),
        Err(e) => json!({"error": e.to_string()}),
    }
}

#[test]
fn digraph_diff_report() {
    let mut reports = Vec::new();
    let mut cases_with_kekule_comparison = 0usize;
    let mut checked = 0usize;

    for line in CORPUS.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value: Value = serde_json::from_str(line).unwrap();
        let Some(smiles) = value.get("smiles").and_then(|v| v.as_str()) else {
            continue; // manifest line
        };
        let atom_idx = value.get("atom_idx").and_then(|v| v.as_u64()).unwrap() as u32;
        let bucket = value.get("bucket").and_then(|v| v.as_str());
        let modern = value.get("modern").and_then(|v| v.as_str()).unwrap();

        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse: {smiles}: {e:?}"));

        if !common::in_mancude_scope(&mol, atom_idx, bucket, modern) {
            continue;
        }
        checked += 1;

        let atom = AtomIdx(atom_idx);
        let budget = CipBudget::default_budget();
        let mut representations = serde_json::Map::new();

        // Kekulé-respelled: apply_kekule preserves atom indices, so no atom-mapping is
        // needed -- the same atom_idx names the same physical atom on both sides.
        if let Ok(kekule) = kekulize(&mol) {
            let kekule_mol = apply_kekule(&mol, &kekule);
            representations.insert(
                "kekule".to_string(),
                divergence_json(first_divergence(&mol, atom, &kekule_mol, atom, budget)),
            );
            cases_with_kekule_comparison += 1;
        } else {
            representations.insert(
                "kekule".to_string(),
                json!({"skipped": "kekulize() failed for this molecule"}),
            );
        }

        // Atom-renumbered: a fixed reversal permutation, tracked via the returned map.
        let n = mol.atom_count();
        let perm: Vec<usize> = (0..n).rev().collect();
        let (renumbered, old_to_new) = renumber_molecule(&mol, &perm);
        let new_atom = AtomIdx(old_to_new[atom_idx as usize]);
        representations.insert(
            "renumbered".to_string(),
            divergence_json(first_divergence(&mol, atom, &renumbered, new_atom, budget)),
        );

        // SMILES-respelled: re-serialize via chematic_smiles's own canonicalizer (a
        // genuinely different traversal order than the corpus's RDKit-authored input,
        // reusing existing infrastructure instead of writing a new "alternate root"
        // writer this milestone doesn't need) and locate the corresponding atom by its
        // local structural signature.
        let text = chematic_smiles::canonical_smiles(&mol);
        match chematic_smiles::parse(&text) {
            Ok(respelled) => match find_atom_by_signature(&mol, atom, &respelled) {
                Some(respelled_atom) => {
                    representations.insert(
                        "canonical_smiles".to_string(),
                        divergence_json(first_divergence(
                            &mol,
                            atom,
                            &respelled,
                            respelled_atom,
                            budget,
                        )),
                    );
                }
                None => {
                    representations.insert(
                        "canonical_smiles".to_string(),
                        json!({"skipped": "stereocenter atom could not be uniquely relocated"}),
                    );
                }
            },
            Err(_) => {
                representations.insert(
                    "canonical_smiles".to_string(),
                    json!({"skipped": "re-parse of the respelled SMILES failed"}),
                );
            }
        }

        reports.push(json!({
            "smiles": smiles,
            "atom_idx": atom_idx,
            "bucket": bucket.unwrap_or("uncharacterized (BucketMisclassified)"),
            "representations": representations,
        }));
    }

    println!("\n=== mancude digraph-diff report (Milestone 3B-0) ===");
    println!("{}", serde_json::to_string_pretty(&reports).unwrap());
    println!("\n=== {checked} cases, {cases_with_kekule_comparison} with a Kekulé comparison ===");

    assert_eq!(
        checked, 98,
        "expected the same 98-case scope as the classification test"
    );
}
