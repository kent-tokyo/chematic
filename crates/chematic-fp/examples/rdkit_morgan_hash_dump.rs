//! Morgan M4-A0 diagnostic: dumps chematic's RDKit-exact-hash reference
//! trace (`chematic_fp::diagnostics::rdkit_morgan_raw_trace`, see
//! `crates/chematic-fp/src/rdkit_morgan_hash.rs`) for a SMILES corpus, one
//! JSON object per molecule, `entries: [[atom_idx, radius,
//! raw_identifier_full_or_null, raw_identifier_default_or_null], ...]` — the
//! two identifiers are independent (see `RdkitMorganRawTraceEntry`'s doc
//! comment for why they can legitimately differ), not one shared value plus
//! emitted flags.
//!
//! Applies `chematic_perception::apply_aromaticity` before tracing (see
//! `rdkit_morgan_hash.rs`'s caller contract: un-aromatized Kekule input
//! diverges from RDKit at every radius >= 1, since RDKit's own
//! `Chem.MolFromSmiles` always sanitizes before hashing).
//!
//! Compared by `scripts/ecfp_rdkit_raw_identifier_parity.py` against
//! `row["full"]["sparse_bit_info"]` (radius 0/1/2 numeric parity) and
//! `row["default"]["sparse_bit_info"]`/`["sparse_counts"]`/
//! `["folded_on_bits"]`/`["folded_bit_info"]` (representative-selection,
//! sparse-count, folded-bit, bitInfo parity) from
//! `scripts/gen_ecfp_rdkit_environment_oracle.py`'s oracle rows.
//!
//! Diagnostic only. Requires the `diagnostics` feature:
//! ```text
//! cargo run -p chematic-fp --release --features diagnostics \
//!     --example rdkit_morgan_hash_dump -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::diagnostics::rdkit_morgan_raw_trace;
use chematic_perception::apply_aromaticity;
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: rdkit_morgan_hash_dump <SMILES.csv> <out.jsonl>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "rdkit_morgan_hash.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));

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
        let mol = apply_aromaticity(&mol);

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
                "entries": entries,
            })
        )
        .unwrap();
        rows += 1;
    }

    eprintln!("input_lines={input_lines} parse_fail={parse_fail} rows={rows} out={out_path}");
}
