//! Phase B (redundant-environment suppression) validation: dumps chematic's
//! *emitted* `(atom_idx, radius, raw_environment_id)` triples under
//! RDKit-invariant atom typing plus RDKit-equivalent redundant-environment
//! suppression, for a SMILES corpus, one JSON object per molecule.
//!
//! Uses `chematic_fp::diagnostics::suppressed_environments_diagnostic` (the
//! `diagnostics` feature) rather than the public
//! `ecfp_with_bitinfo_rdkit_environment_experimental` API, since the raw
//! (unfolded) `raw_environment_id` is needed to measure real sparse-count
//! *multiplicity* (how many distinct environments hash to the same id) --
//! the public API only exposes `(atom, radius)` pairs keyed by *folded* bit,
//! which can't distinguish two different raw ids that happen to fold to the
//! same bit. Compared by `scripts/ecfp_rdkit_suppression_parity.py` against
//! `row["default"]["sparse_bit_info"]` (pair-set match) and
//! `row["default"]["sparse_counts"]` (count-multiset shape match) from
//! `scripts/gen_ecfp_rdkit_environment_oracle.py`'s oracle rows.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --features diagnostics \
//!     --example morgan_suppression_dump -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::EcfpConfig;
use chematic_fp::diagnostics::suppressed_environments_diagnostic;
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

        let mut emitted: Vec<[u64; 3]> = suppressed_environments_diagnostic(&mol, &config)
            .into_iter()
            .map(|(atom, radius, raw_id)| [atom as u64, radius as u64, raw_id])
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
