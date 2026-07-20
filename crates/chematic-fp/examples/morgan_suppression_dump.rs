//! Phase B (redundant-environment suppression) validation: dumps chematic's
//! *emitted* `(atom_idx, radius)` pairs under
//! `ecfp_with_bitinfo_rdkit_environment_experimental` — i.e. RDKit-invariant
//! atom typing plus RDKit-equivalent redundant-environment suppression — for
//! a SMILES corpus, one JSON object per molecule.
//!
//! Uses only the public API (no `diagnostics` feature needed). Compared
//! against `row["default"]["sparse_bit_info"]` (flattened) from
//! `scripts/gen_ecfp_rdkit_environment_oracle.py`'s oracle rows by
//! `scripts/ecfp_rdkit_suppression_parity.py` — RDKit's own default
//! (`includeRedundantEnvironments=False`) generator's real emitted-pair set,
//! already produced by PR #120's unmodified oracle script.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example morgan_suppression_dump \
//!     -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::{EcfpConfig, ecfp_with_bitinfo_rdkit_environment_experimental};
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: morgan_suppression_dump <SMILES.csv> <out.jsonl>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "suppression.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let config = EcfpConfig::default(); // radius 2, nbits 2048, chirality off

    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let mut input_lines = 0usize;
    let mut parse_fail = 0usize;
    let mut rows = 0usize;

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        input_lines += 1;
        // 0-based position in this file's row stream -- must match the RDKit
        // oracle's row_id exactly (same sync discipline as PR #120's
        // morgan_rdkit_environment_trace.rs).
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

        let (_fp, info) = ecfp_with_bitinfo_rdkit_environment_experimental(&mol, &config);
        let mut emitted: Vec<[u32; 2]> = info
            .values()
            .flat_map(|envs| envs.iter().map(|&(atom, radius)| [atom, radius]))
            .collect();
        emitted.sort_unstable();

        writeln!(
            f,
            "{}",
            json!({
                "row_id": row_id,
                "smiles": smi,
                "parse_ok": true,
                "atom_count": mol.atom_count(),
                "emitted": emitted,
            })
        )
        .unwrap();
        rows += 1;
    }

    eprintln!("input_lines={input_lines} parse_fail={parse_fail} rows={rows} out={out_path}");
}
