//! Same-process A/B benchmark for SD-field serialization.
//!
//! Usage:
//! `cargo run --release -p chematic-mol --example sdf_writer_ab -- path/to/input.sdf`

use std::collections::HashMap;
use std::fs::File;
use std::hint::black_box;
use std::io::BufReader;
use std::time::{Duration, Instant};

use chematic_core::Molecule;
use chematic_mol::mol2000::write_sdf_record_into;
use chematic_mol::{MolMetadata, SdfFileReader, write_mol_with_coords, write_sdf_record};

fn baseline_record(mol: &Molecule, meta: &MolMetadata, props: &HashMap<String, String>) -> String {
    let mut out = write_mol_with_coords(mol, meta, &[]);
    for (key, value) in props {
        if !key.starts_with('_') {
            out.push_str(&format!("> <{key}>\n{value}\n\n"));
        }
    }
    out.push_str("$$$$\n");
    out
}

fn measure(
    records: &[(Molecule, MolMetadata, HashMap<String, String>)],
    repeats: usize,
    optimized: bool,
) -> (Duration, usize) {
    let started = Instant::now();
    let mut bytes = 0usize;
    for _ in 0..repeats {
        for (mol, meta, props) in records {
            let record = if optimized {
                write_sdf_record(mol, meta, &[], props)
            } else {
                baseline_record(mol, meta, props)
            };
            bytes = bytes.wrapping_add(black_box(record.len()));
        }
    }
    (started.elapsed(), bytes)
}

fn measure_reused(
    records: &[(Molecule, MolMetadata, HashMap<String, String>)],
    repeats: usize,
) -> (Duration, usize) {
    let started = Instant::now();
    let mut bytes = 0usize;
    let mut record = String::new();
    for _ in 0..repeats {
        for (mol, meta, props) in records {
            write_sdf_record_into(&mut record, mol, meta, &[], props);
            bytes = bytes.wrapping_add(black_box(record.len()));
        }
    }
    (started.elapsed(), bytes)
}

fn median(values: &mut [Duration]) -> Duration {
    values.sort_unstable();
    values[values.len() / 2]
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("pass an SDF corpus path as the first argument");
    let file = File::open(&path).expect("open SDF corpus");
    let records: Vec<_> = SdfFileReader::fast(BufReader::new(file))
        .map(|record| {
            let record = record.expect("parse SDF record");
            (record.mol, record.meta, record.properties)
        })
        .collect();
    assert!(!records.is_empty(), "SDF corpus must contain records");

    let repeats = 50;
    let rounds = 9;
    for (mol, meta, props) in &records {
        assert_eq!(
            baseline_record(mol, meta, props),
            write_sdf_record(mol, meta, &[], props),
            "optimized serializer changed record bytes"
        );
    }
    let _ = measure(&records, 1, false);
    let _ = measure(&records, 1, true);

    let mut baseline = Vec::with_capacity(rounds);
    let mut optimized = Vec::with_capacity(rounds);
    let mut reused = Vec::with_capacity(rounds);
    let mut expected_bytes = None;
    for round in 0..rounds {
        let (first_optimized, second_optimized) = if round % 2 == 0 {
            (false, true)
        } else {
            (true, false)
        };
        for is_optimized in [first_optimized, second_optimized] {
            let (elapsed, bytes) = measure(&records, repeats, is_optimized);
            if let Some(expected) = expected_bytes {
                assert_eq!(
                    bytes, expected,
                    "serializers produced different byte counts"
                );
            } else {
                expected_bytes = Some(bytes);
            }
            if is_optimized {
                optimized.push(elapsed);
            } else {
                baseline.push(elapsed);
            }
        }
        let (elapsed, bytes) = measure_reused(&records, repeats);
        assert_eq!(bytes, expected_bytes.expect("one serializer measurement"));
        reused.push(elapsed);
    }

    let baseline = median(&mut baseline);
    let optimized = median(&mut optimized);
    let reused = median(&mut reused);
    let operations = records.len() * repeats;
    println!(
        "records={} repeats={} baseline_us_per_record={:.3} optimized_us_per_record={:.3} reused_us_per_record={:.3} existing_to_reused_speedup_x={:.3} historical_to_reused_speedup_x={:.3}",
        records.len(),
        repeats,
        baseline.as_secs_f64() * 1e6 / operations as f64,
        optimized.as_secs_f64() * 1e6 / operations as f64,
        reused.as_secs_f64() * 1e6 / operations as f64,
        optimized.as_secs_f64() / reused.as_secs_f64(),
        baseline.as_secs_f64() / reused.as_secs_f64(),
    );
}
