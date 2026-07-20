//! Morgan M4-A0 residual-mechanism confirmation: identical to
//! `rdkit_morgan_hash_dump` except it feeds the RDKit-parity aromaticity
//! reference engine (`apply_aromaticity_rdkit_parity_experimental`,
//! `crates/chematic-perception/src/rdkit_parity.rs`) instead of production
//! `apply_aromaticity`, to test whether `rdkit_morgan_hash_dump`'s residual
//! mismatches against real RDKit trace to Hueckel-vs-RDKit aromaticity
//! *perception* disagreement rather than a hash defect.
//!
//! **No silent fallback.** An earlier version of this example fell back to
//! production `apply_aromaticity` whenever the RDKit-parity engine returned
//! `Err`, and folded that Hueckel-derived row into the "RDKit-parity"
//! exact-match count -- exactly the accident the RDKit-parity engine's own
//! `Result` contract exists to prevent (fail explicitly, never silently
//! substitute a different algorithm's answer). Fixed: every row is tagged
//! with one of four `aromaticity_status` values (`rdkit_parity_success` /
//! `rdkit_parity_kekulization_failed` / `rdkit_parity_internal_error` /
//! `parse_failed`); a row whose RDKit-parity preprocessing fails carries
//! `entries: null` and an `aromaticity_error` message, never a Hueckel-
//! derived trace. `scripts/ecfp_rdkit_raw_identifier_parity.py`'s
//! exact-match denominator only ever counts `rdkit_parity_success` rows.
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
use chematic_perception::{AromaticityError, apply_aromaticity_rdkit_parity_experimental};
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
    let mut parse_fail = 0usize;
    let mut kekulization_failed = 0usize;
    let mut internal_error = 0usize;
    let mut success = 0usize;

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
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": false,
                        "aromaticity_status": "parse_failed",
                        "aromaticity_error": e.to_string(),
                        "entries": serde_json::Value::Null,
                    })
                )
                .unwrap();
                rows += 1;
                continue;
            }
        };

        match apply_aromaticity_rdkit_parity_experimental(&mol) {
            Ok(aromatized) => {
                success += 1;
                let mut entries: Vec<(u32, u32, Option<u32>, Option<u32>)> =
                    rdkit_morgan_raw_trace(&aromatized, 2)
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
                        "aromaticity_status": "rdkit_parity_success",
                        "atom_count": aromatized.atom_count(),
                        "entries": entries,
                    })
                )
                .unwrap();
            }
            Err(err) => {
                let status = match &err {
                    AromaticityError::KekulizationFailed { .. } => {
                        kekulization_failed += 1;
                        "rdkit_parity_kekulization_failed"
                    }
                    AromaticityError::InternalInvariantViolation { .. } => {
                        internal_error += 1;
                        "rdkit_parity_internal_error"
                    }
                };
                writeln!(
                    f,
                    "{}",
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": true,
                        "aromaticity_status": status,
                        "aromaticity_error": err.to_string(),
                        "entries": serde_json::Value::Null,
                    })
                )
                .unwrap();
            }
        }
        rows += 1;
    }
    eprintln!(
        "rows={rows} parse_fail={parse_fail} rdkit_parity_success={success} \
         rdkit_parity_kekulization_failed={kekulization_failed} \
         rdkit_parity_internal_error={internal_error} out={out_path}"
    );
}
