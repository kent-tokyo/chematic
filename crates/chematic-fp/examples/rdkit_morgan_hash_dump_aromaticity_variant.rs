//! Morgan M4-A0 residual-mechanism confirmation: identical to
//! `rdkit_morgan_hash_dump` except it feeds the RDKit-parity aromaticity
//! reference engine (`apply_aromaticity_rdkit_parity_experimental`,
//! `crates/chematic-perception/src/rdkit_parity.rs`) instead of production
//! `apply_aromaticity`, to test whether `rdkit_morgan_hash_dump`'s residual
//! mismatches against real RDKit trace to Hueckel-vs-RDKit aromaticity
//! *perception* disagreement rather than a hash defect. Falls back to
//! `apply_aromaticity` on kekulization failure (the RDKit-parity engine is
//! fallible) -- `fallback_to_hueckel: true` is recorded on those rows so a
//! fallback-driven "match" can't be mistaken for a parity-engine-driven one.
//!
//! Load-bearing for the M4-A0 report's residual-localization claim -- run
//! on the FULL corpus, not just the known-mismatching subset, so it can
//! also catch a regression among molecules that already matched under
//! Hueckel aromaticity. Kept as a permanent diagnostic artifact (not
//! removed after use) so the claim stays reproducible from the PR alone.
//! Requires the `diagnostics` feature:
//! ```text
//! cargo run -p chematic-fp --release --features diagnostics \
//!     --example rdkit_morgan_hash_dump_aromaticity_variant -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::diagnostics::rdkit_morgan_raw_trace;
use chematic_perception::{apply_aromaticity, apply_aromaticity_rdkit_parity_experimental};
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args.get(1).cloned().unwrap_or_else(|| {
        panic!("usage: rdkit_morgan_hash_dump_aromaticity_variant <SMILES.csv> <out.jsonl>")
    });
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "rdkit_morgan_hash_aromaticity_variant.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let mut rows = 0usize;
    let mut fallback_count = 0usize;
    let mut parse_fail = 0usize;

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let row_id = rows;

        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                parse_fail += 1;
                writeln!(
                    f,
                    "{}",
                    json!({"row_id": row_id, "smiles": smi, "parse_ok": false, "error": e.to_string()})
                )
                .unwrap();
                rows += 1;
                continue;
            }
        };

        let (mol, fallback_to_hueckel) = match apply_aromaticity_rdkit_parity_experimental(&mol) {
            Ok(m) => (m, false),
            Err(_) => {
                fallback_count += 1;
                (apply_aromaticity(&mol), true)
            }
        };

        let mut entries: Vec<(u32, u32, Option<u32>, Option<u32>)> =
            rdkit_morgan_raw_trace(&mol, 2)
                .into_iter()
                .map(|e| {
                    (
                        e.atom_idx,
                        e.radius,
                        e.raw_identifier_full,
                        e.raw_identifier_default,
                    )
                })
                .collect();
        entries.sort_unstable();

        writeln!(
            f,
            "{}",
            json!({
                "row_id": row_id,
                "smiles": smi,
                "parse_ok": true,
                "atom_count": mol.atom_count(),
                "fallback_to_hueckel": fallback_to_hueckel,
                "entries": entries,
            })
        )
        .unwrap();
        rows += 1;
    }
    eprintln!(
        "rows={rows} parse_fail={parse_fail} fallback_to_hueckel={fallback_count} out={out_path}"
    );
}
