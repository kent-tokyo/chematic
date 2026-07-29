//! IO-1 performance record (not a gate, not a Criterion-registered
//! benchmark -- see `feedback_criterion_gate_pseudo_replication` for why
//! this repo treats the automated Criterion CI gate's samples as
//! non-independent; this is a one-off, manually-reported measurement,
//! matching the precedent of `morgan_suppression_benchmark.rs`/
//! `rdkit_morgan_ecfp4_benchmark.rs`). Performance here is purely
//! informational -- parser correctness is never traded for speed.
//!
//! Reports records/sec over a full streaming pass of a 10,000-record
//! SMILES table file (title line + name column + 2 extra property columns,
//! ~2% deliberately malformed rows to also measure invalid-record recovery
//! throughput), using the streaming `SmilesRecordReader` (no
//! `read_to_string` of the whole file).
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example smiles_table_benchmark -- <10k.smi>
//! ```

use chematic_mol::{SmilesReaderOptions, SmilesRecordReader};
use std::io::BufReader;
use std::time::Instant;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| panic!("usage: smiles_table_benchmark <10k.smi>"));

    let file = std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let reader = SmilesRecordReader::new(
        BufReader::new(file),
        SmilesReaderOptions {
            title_line: true,
            ..Default::default()
        },
    );

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
