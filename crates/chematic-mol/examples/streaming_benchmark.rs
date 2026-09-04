//! Measure file-backed SDF/MOL/XYZ record streaming in the current workspace.
//!
//! Usage:
//! `cargo run -p chematic-mol --example streaming_benchmark -- --format sdf path`

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use chematic_mol::{SdfFileReader, XyzFileReader};

fn arg(name: &str, default: &str) -> String {
    let mut args = env::args().skip(1);
    while let Some(value) = args.next() {
        if value == name {
            return args.next().unwrap_or_else(|| default.to_string());
        }
    }
    default.to_string()
}

fn main() {
    let format = arg("--format", "sdf");
    let path = arg("--path", "benchmarks/fixtures/streaming.sdf");
    let repeats: usize = arg("--repeats", "20")
        .parse()
        .expect("--repeats must be an integer");
    assert!(repeats > 0, "--repeats must be positive");
    let bytes = std::fs::metadata(&path)
        .expect("benchmark input must exist")
        .len();
    let started = Instant::now();
    let mut records = 0usize;
    let mut failures = 0usize;

    for _ in 0..repeats {
        match format.as_str() {
            "sdf" | "mol" => {
                let input = File::open(&path).expect("open SDF/MOL input");
                for result in SdfFileReader::fast(BufReader::new(input)) {
                    match result {
                        Ok(_) => records += 1,
                        Err(_) => failures += 1,
                    }
                }
            }
            "xyz" => {
                let input = File::open(&path).expect("open XYZ input");
                for result in XyzFileReader::new(BufReader::new(input)) {
                    match result {
                        Ok(_) => records += 1,
                        Err(_) => failures += 1,
                    }
                }
            }
            other => panic!("unsupported format {other}; choose sdf, mol, or xyz"),
        }
    }

    let elapsed = started.elapsed().as_secs_f64();
    let input_bytes = bytes as usize * repeats;
    println!(
        "{{\"format\":\"{format}\",\"path\":\"{}\",\"repeats\":{repeats},\"records\":{records},\"failures\":{failures},\"input_bytes\":{input_bytes},\"seconds\":{elapsed:.6},\"records_per_second\":{:.2},\"bytes_per_second\":{:.2}}}",
        Path::new(&path).display(),
        records as f64 / elapsed,
        input_bytes as f64 / elapsed,
    );
}
