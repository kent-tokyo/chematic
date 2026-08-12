//! Diagnosis-only dump for the ECFP4 bit-exact API-design workstream
//! (`diag/ecfp4-bitexact-api`). Not production code, not reused by any
//! production path -- see `docs/rfcs/ecfp4_bitexact_api_rfc.md`.
//!
//! For every fixture in `scripts/ecfp4_bitexact_matrix_fixtures.csv`, dumps:
//! - The production, single-config path
//!   (`chematic_fp::rdkit_morgan_ecfp4_experimental`, radius=2, fpSize=2048,
//!   `includeRedundantEnvironments=false`, `useChirality=false`,
//!   `useBondTypes=true`): status (success/error+reason), the unfolded
//!   `sparse_counts` (raw_id -> emission count, radius 0-2, already public --
//!   used by the diagnosis script to test the nBits and count-vs-binary axes
//!   by *folding this same data*, not by touching any production code), and
//!   the folded 2048 on-bits.
//! - A radius sweep (`chematic_fp::diagnostics::rdkit_morgan_raw_trace`,
//!   `diagnostics` feature) up to radius 3, `raw_identifier_default` only
//!   (RDKit's real `includeRedundantEnvironments=false` lifecycle -- the
//!   `full` lifecycle is a distinct diagnostic surface, not compared here).
//!   Run on the SAME aromaticity engine the production path itself uses
//!   (`apply_aromaticity_rdkit_parity_experimental`), so a radius-3 mismatch
//!   can't be mistaken for an aromaticity-engine difference.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --features diagnostics \
//!     --example rdkit_ecfp4_bitexact_matrix_dump -- \
//!     scripts/ecfp4_bitexact_matrix_fixtures.csv \
//!     validation/results/ecfp4_bitexact_matrix_dump.jsonl
//! ```

use chematic_fp::diagnostics::rdkit_morgan_raw_trace;
use chematic_fp::{RdkitMorganError, rdkit_morgan_ecfp4_experimental};
use chematic_perception::{AromaticityError, apply_aromaticity_rdkit_parity_experimental};
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

const MAX_RADIUS_PROBE: u32 = 3;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args.get(1).cloned().unwrap_or_else(|| {
        panic!("usage: rdkit_ecfp4_bitexact_matrix_dump <fixtures.csv> <out.jsonl>")
    });
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "ecfp4_bitexact_matrix_dump.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));

    let mut rows = 0usize;
    let mut status_counts = std::collections::BTreeMap::<&'static str, usize>::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("id|tags|smiles") {
            continue;
        }
        let mut parts = line.splitn(3, '|');
        let fixture_id = parts.next().unwrap_or_default();
        let tags = parts.next().unwrap_or_default();
        let smi = parts.next().unwrap_or_default().trim();
        if smi.is_empty() {
            continue;
        }

        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                *status_counts.entry("parse_failed").or_default() += 1;
                writeln!(
                    f,
                    "{}",
                    json!({
                        "fixture_id": fixture_id,
                        "tags": tags,
                        "smiles": smi,
                        "status": "parse_failed",
                        "error": e.to_string(),
                    })
                )
                .unwrap();
                rows += 1;
                continue;
            }
        };

        // ---- Production, single-config path ----
        let (status, error, sparse_counts, folded_on_bits) =
            match rdkit_morgan_ecfp4_experimental(&mol) {
                Ok(r) => {
                    *status_counts.entry("success").or_default() += 1;
                    let mut counts: Vec<(u32, u32)> =
                        r.sparse_counts.iter().map(|(&k, &v)| (k, v)).collect();
                    counts.sort_unstable();
                    let mut bits: Vec<usize> =
                        (0..2048).filter(|&b| r.fingerprint.get(b)).collect();
                    bits.sort_unstable();
                    ("success", None, counts, bits)
                }
                Err(e) => {
                    let bucket = match &e {
                        RdkitMorganError::Aromaticity(AromaticityError::KekulizationFailed {
                            ..
                        }) => "error_kekulization_failed",
                        RdkitMorganError::Aromaticity(_) => "error_aromaticity_other",
                        RdkitMorganError::UnsupportedBondOrder { .. } => {
                            "error_unsupported_bond_order"
                        }
                        RdkitMorganError::InternalInvariantViolation { .. } => {
                            "error_internal_invariant"
                        }
                    };
                    *status_counts.entry(bucket).or_default() += 1;
                    (bucket, Some(e.to_string()), Vec::new(), Vec::new())
                }
            };

        // ---- Radius sweep (diagnostics feature), same aromaticity engine ----
        let radius_default_pairs: Vec<(u32, u32, u32)> =
            match apply_aromaticity_rdkit_parity_experimental(&mol) {
                Ok(aromatized) => {
                    let trace = rdkit_morgan_raw_trace(&aromatized, MAX_RADIUS_PROBE);
                    let mut pairs: Vec<(u32, u32, u32)> = trace
                        .iter()
                        .filter_map(|e| {
                            e.raw_identifier_default
                                .map(|rid| (e.atom_idx, e.radius, rid))
                        })
                        .collect();
                    pairs.sort_unstable();
                    pairs
                }
                Err(_) => Vec::new(), // same aromaticity failure as the production path above
            };

        writeln!(
            f,
            "{}",
            json!({
                "fixture_id": fixture_id,
                "tags": tags,
                "smiles": smi,
                "status": status,
                "error": error,
                "sparse_counts": sparse_counts,
                "folded_on_bits_2048": folded_on_bits,
                "radius_default_pairs": radius_default_pairs,
            })
        )
        .unwrap();
        rows += 1;
    }

    eprintln!("wrote {rows} rows to {out_path}");
    eprintln!("status counts: {status_counts:?}");
}
