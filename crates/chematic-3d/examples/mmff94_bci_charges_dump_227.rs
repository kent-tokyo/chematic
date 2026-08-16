//! Issue #227 Phase 2 (BCI investigation): dumps chematic's CURRENT,
//! UNMODIFIED `chematic_ff::mmff94_charges_numeric` output per heavy atom
//! for every molecule in the 265-molecule Wave 1 corpus. Measurement-only --
//! no chematic-ff production code is touched by this file. Compared
//! downstream (`scripts/mmff94_bci_charges_compare_227.py`) against a live
//! RDKit oracle (`scripts/mmff94_bci_charges_oracle_227.py`) to establish
//! whether a real gap exists BEFORE any production fix is written, per the
//! issue #227 Phase 2 directive's falsify-before-fix requirement.
//!
//! Run this example again, unmodified, after any production fix to
//! `mmff94_charges_numeric` lands -- the "before" and "after" snapshots are
//! produced by literally the same tool, so a diff between the two runs is a
//! true before/after comparison, not two different measurement methods.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_bci_charges_dump_227
//!   > validation/results/mmff94_bci_charges_chematic_<label>.jsonl`

use chematic_core::AtomIdx;
use serde_json::{Value, json};
use std::collections::HashMap;

fn load_manifest(path: &str) -> Vec<(String, String)> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    let v: Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("bad json in {path}: {e}"));
    v["molecules"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            (
                m["name"].as_str().unwrap().to_string(),
                m["smiles"].as_str().unwrap().to_string(),
            )
        })
        .collect()
}

fn main() {
    let manifests: Vec<(&str, &str)> = vec![
        (
            "A",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
        ),
        (
            "B",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
        ),
    ];

    let mut n_ok = 0usize;
    let mut n_fail = 0usize;
    let mut fail_reasons: HashMap<String, usize> = HashMap::new();

    for (tier, path) in manifests {
        for (name, smiles) in load_manifest(path) {
            let mol = match chematic_smiles::parse(&smiles) {
                Ok(m) => m,
                Err(e) => {
                    n_fail += 1;
                    *fail_reasons.entry("parse_failure".to_string()).or_insert(0) += 1;
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles,
                               "status": "parse_failure", "reason": format!("{e}")})
                    );
                    continue;
                }
            };

            let n_heavy = mol.atom_count();
            match chematic_ff::mmff94_charges_numeric(&mol) {
                Ok(charges) => {
                    let atom_types = chematic_ff::assign_mmff94_numeric_types(&mol).ok();
                    let charges_json: Vec<Value> = (0..n_heavy)
                        .map(|i| {
                            let a = mol.atom(AtomIdx(i as u32));
                            json!({
                                "index": i,
                                "element": a.element.symbol(),
                                "mmff94_numeric_type": atom_types.as_ref().map(|t| t[i]),
                                "chematic_partial_charge": charges[i],
                            })
                        })
                        .collect();
                    n_ok += 1;
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles,
                               "status": "ok", "n_heavy": n_heavy, "charges": charges_json})
                    );
                }
                Err(e) => {
                    n_fail += 1;
                    *fail_reasons
                        .entry(
                            format!("{e:?}")
                                .split('(')
                                .next()
                                .unwrap_or("unknown")
                                .to_string(),
                        )
                        .or_insert(0) += 1;
                    println!(
                        "{}",
                        json!({"tier": tier, "name": name, "smiles": smiles,
                               "status": "charges_error", "reason": format!("{e:?}")})
                    );
                }
            }
        }
    }

    eprintln!("total_ok={n_ok} total_fail={n_fail} fail_reasons={fail_reasons:?}");
}
