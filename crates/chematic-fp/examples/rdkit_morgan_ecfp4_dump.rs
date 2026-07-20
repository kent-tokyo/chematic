//! Phase B acceptance-gate dump: `rdkit_morgan_ecfp4_experimental`'s real
//! output on the M4-A0 5,048-input corpus, for comparison against the same
//! RDKit oracle rows M4-A0 itself validated against
//! (`scripts/gen_ecfp_rdkit_environment_oracle.py`'s `default`
//! (`includeRedundantEnvironments=false`) lifecycle -- the only lifecycle
//! this production API ever computes).
//!
//! Every row is tagged with an explicit `status`, mirroring
//! `rdkit_morgan_hash_dump_aromaticity_variant`'s established taxonomy --
//! no fallback, no pooling a different status into `"success"`. Denominator
//! discipline for the comparator: only `status == "success"` rows enter the
//! exact-match rate.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example rdkit_morgan_ecfp4_dump \
//!     -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::{RdkitMorganError, rdkit_morgan_ecfp4_experimental};
use chematic_perception::AromaticityError;
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: rdkit_morgan_ecfp4_dump <SMILES.csv> <out.jsonl>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "rdkit_morgan_ecfp4_dump.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));

    let mut rows = 0usize;
    let mut status_counts = std::collections::BTreeMap::<&'static str, usize>::new();

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let row_id = rows;

        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                *status_counts.entry("parse_failed").or_default() += 1;
                writeln!(
                    f,
                    "{}",
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": false,
                        "status": "parse_failed",
                        "error": e.to_string(),
                    })
                )
                .unwrap();
                rows += 1;
                continue;
            }
        };

        match rdkit_morgan_ecfp4_experimental(&mol) {
            Ok(r) => {
                *status_counts.entry("success").or_default() += 1;

                let mut default_pairs: Vec<(u32, u32, u32)> = r
                    .raw_bit_info
                    .iter()
                    .flat_map(|(&raw_id, envs)| envs.iter().map(move |&(a, rad)| (a, rad, raw_id)))
                    .collect();
                default_pairs.sort_unstable();

                let mut sparse_counts: Vec<(u32, u32)> = r.sparse_counts.into_iter().collect();
                sparse_counts.sort_unstable();

                let mut folded_on_bits: Vec<usize> =
                    (0..2048usize).filter(|&b| r.fingerprint.get(b)).collect();
                folded_on_bits.sort_unstable();

                let mut folded_bit_info: Vec<(usize, Vec<(u32, u32)>)> = r
                    .folded_bit_info
                    .into_iter()
                    .map(|(bit, mut envs)| {
                        envs.sort_unstable();
                        (bit, envs)
                    })
                    .collect();
                folded_bit_info.sort_unstable_by_key(|(bit, _)| *bit);

                writeln!(
                    f,
                    "{}",
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": true,
                        "status": "success",
                        "default_pairs": default_pairs,
                        "sparse_counts": sparse_counts,
                        "folded_on_bits": folded_on_bits,
                        "folded_bit_info": folded_bit_info,
                    })
                )
                .unwrap();
            }
            Err(err) => {
                let status = match &err {
                    RdkitMorganError::Aromaticity(AromaticityError::KekulizationFailed {
                        ..
                    }) => "rdkit_parity_kekulization_failed",
                    RdkitMorganError::Aromaticity(
                        AromaticityError::InternalInvariantViolation { .. },
                    ) => "rdkit_parity_internal_error",
                    RdkitMorganError::UnsupportedBondOrder { .. } => "unsupported_bond_order",
                    RdkitMorganError::InternalInvariantViolation { .. } => {
                        "internal_invariant_violation"
                    }
                };
                *status_counts.entry(status).or_default() += 1;
                writeln!(
                    f,
                    "{}",
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": true,
                        "status": status,
                        "error": err.to_string(),
                    })
                )
                .unwrap();
            }
        }
        rows += 1;
    }

    eprintln!("rows={rows} status_counts={status_counts:?} out={out_path}");
}
