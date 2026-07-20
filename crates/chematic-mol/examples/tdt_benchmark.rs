//! IO-2 performance record (not a gate, not a Criterion-registered
//! benchmark -- see `feedback_criterion_gate_pseudo_replication` and the
//! precedent of `smiles_table_benchmark.rs`). Performance here is purely
//! informational.
//!
//! Reports records/sec over a full streaming pass of a 10,000-record TDT
//! file (SMI + NAME + one property, ~2% deliberately malformed rows).
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example tdt_benchmark -- <10k.tdt>
//! ```

use chematic_mol::{TdtReaderOptions, TdtRecordReader};
use std::io::BufReader;
use std::time::Instant;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| panic!("usage: tdt_benchmark <10k.tdt>"));

    let file = std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let reader = TdtRecordReader::new(BufReader::new(file), TdtReaderOptions::default());

    let start = Instant::now();
    let mut success = 0usize;
    let mut errors = 0usize;
    for result in reader {
        match result {
            Ok(_) => success += 1,
            Err(_) => errors += 1,
        }
    }
    let elapsed = start.elapsed();
    let total = success + errors;
    let records_per_sec = total as f64 / elapsed.as_secs_f64();

    println!(
        "total={total} success={success} errors={errors} elapsed={elapsed:?} records_per_sec={records_per_sec:.0}"
    );
}
