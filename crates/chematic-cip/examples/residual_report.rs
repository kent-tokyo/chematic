//! Milestone 4A-0: for a frozen residual set (rows where `assign_cip_accurate_experimental`
//! disagrees with the modern RDKit `rdCIPLabeler` oracle), print all three engines'
//! outputs side by side -- `FastApproximate` (`chematic_chem::assign_cip`),
//! `AccurateExperimental` (this crate), and the oracle (`modern`, supplied by the input
//! file, produced by `scripts/cip_ground_truth.py`-style RDKit regeneration since this
//! crate has no RDKit binding of its own).
//!
//! Deliberately narrow, not a general corpus tool: `scripts/cip_accurate_full_corpus_report.py`
//! already does full-corpus accounting. This exists for the small, fixed residual set a
//! milestone's own classification work needs to inspect closely, row by row.
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example residual_report -- <residual.jsonl>
//!
//! Input: one JSON object per line, `{"smiles": "...", "atom_idx": N, "modern": "R"}`
//! (exactly what a residual-diff step, e.g. M4A-0's own analysis script, produces).
//! Output: one JSON object per line, adding `legacy` (FastApproximate) and `accurate`
//! (AccurateExperimental, or `skip:<reason>`) fields.

use std::env;
use std::fs;

use chematic_cip::{CipBudget, SkipReason, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, CipCode};
use serde_json::Value;

fn code_str(c: CipCode) -> &'static str {
    match c {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
        CipCode::LowerR => "r",
        CipCode::LowerS => "s",
    }
}

fn skip_str(r: SkipReason) -> &'static str {
    match r {
        SkipReason::NotFourSubstituents => "skip:not4",
        SkipReason::Tied => "skip:tied",
        SkipReason::BudgetExceeded => "skip:budget",
        SkipReason::OracleUnstable => "skip:oracle-unstable",
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: residual_report <residual.jsonl>");
        std::process::exit(64);
    };
    let content = fs::read_to_string(path).expect("read residual file");

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(line).expect("valid JSON line");
        let smi = row["smiles"].as_str().expect("smiles field").to_string();
        let atom_idx = row["atom_idx"].as_u64().expect("atom_idx field") as u32;
        let modern = row["modern"].as_str().expect("modern field").to_string();

        let mol = match chematic_smiles::parse(&smi) {
            Ok(m) => m,
            Err(e) => {
                println!(
                    r#"{{"smiles":{smi:?},"atom_idx":{atom_idx},"modern":{modern:?},"legacy":"PARSE_ERR","accurate":"PARSE_ERR: {e}"}}"#
                );
                continue;
            }
        };

        let legacy = chematic_chem::assign_cip(&mol)
            .get(AtomIdx(atom_idx))
            .map(code_str)
            .unwrap_or("none");

        let accurate = match assign_cip_accurate_experimental(&mol, CipBudget::default_budget()) {
            Ok(result) => result
                .assignments
                .iter()
                .find(|(idx, _)| idx.0 == atom_idx)
                .map(|(_, code)| code_str(*code).to_string())
                .or_else(|| {
                    result
                        .skipped
                        .iter()
                        .find(|(idx, _)| idx.0 == atom_idx)
                        .map(|(_, reason)| skip_str(*reason).to_string())
                })
                .unwrap_or_else(|| "MISSING".to_string()),
            Err(e) => format!("ERR: {e:?}"),
        };

        println!(
            r#"{{"smiles":{smi:?},"atom_idx":{atom_idx},"modern":{modern:?},"legacy":{legacy:?},"accurate":{accurate:?}}}"#
        );
    }
}
