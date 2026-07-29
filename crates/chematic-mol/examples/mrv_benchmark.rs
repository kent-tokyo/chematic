//! IO-3 performance record (not a gate, not a Criterion-registered
//! benchmark -- see `feedback_criterion_gate_pseudo_replication` and the
//! precedent of `tdt_benchmark.rs`/`smiles_table_benchmark.rs`). Performance
//! here is purely informational.
//!
//! Unlike SMILES-table/TDT (one file, many records), MRV is one document
//! per molecule -- this reports documents/sec parsing every `.mrv` file in
//! a directory (e.g. the 206-fixture pool from `gen_rdkit_mrv_oracle.py`).
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example mrv_benchmark -- <fixtures_dir>
//! ```

use chematic_mol::parse_mrv;
use std::time::Instant;

fn main() {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| panic!("usage: mrv_benchmark <fixtures_dir>"));

    let mut paths: Vec<_> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read_dir {dir}: {e}"))
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|ext| ext == "mrv"))
        .collect();
    paths.sort();

    if paths.is_empty() {
        panic!("no .mrv files found in {dir}");
    }

    let texts: Vec<String> = paths
        .iter()
        .map(|p| std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {p:?}: {e}")))
        .collect();

    let start = Instant::now();
    let mut success = 0usize;
    let mut errors = 0usize;
    for text in &texts {
        match parse_mrv(text) {
            Ok(_) => success += 1,
            Err(_) => errors += 1,
        }
    }
    let elapsed = start.elapsed();
    let total = success + errors;
    let docs_per_sec = total as f64 / elapsed.as_secs_f64();

    println!(
        "total={total} success={success} errors={errors} elapsed={elapsed:?} docs_per_sec={docs_per_sec:.0}"
    );
}
