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

fn arg_usize(name: &str, default: usize) -> usize {
    arg(name, &default.to_string())
        .parse()
        .unwrap_or_else(|_| panic!("{name} must be an integer"))
}

fn main() {
    let format = arg("--format", "sdf");
    let path = arg("--path", "benchmarks/fixtures/streaming.sdf");
    let repeats: usize = arg("--repeats", "20")
        .parse()
        .expect("--repeats must be an integer");
    assert!(repeats > 0, "--repeats must be positive");
    let max_input_bytes = arg_usize("--max-input-bytes", 1 << 30);
    let max_record_bytes = arg_usize("--max-record-bytes", 16 << 20);
    let max_line_bytes = arg_usize("--max-line-bytes", 16 << 20);
    let max_records = arg_usize("--max-records", 100_000);
    let max_atoms = arg_usize("--max-atoms", 1_000_000);
    let materialized = matches!(
        format.as_str(),
        "v3000" | "mol2" | "cml" | "cdxml" | "mmcif"
    );
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
                let limits = chematic_mol::SdfParseLimits {
                    max_input_bytes,
                    max_record_bytes,
                    max_line_bytes,
                    max_records,
                };
                for result in SdfFileReader::with_limits(BufReader::new(input), limits) {
                    match result {
                        Ok(_) => records += 1,
                        Err(_) => failures += 1,
                    }
                }
            }
            "xyz" => {
                let input = File::open(&path).expect("open XYZ input");
                let limits = chematic_mol::XyzParseLimits {
                    max_input_bytes,
                    max_atoms_per_frame: max_atoms,
                    max_frames: max_records,
                    max_line_bytes,
                };
                for result in XyzFileReader::with_limits(BufReader::new(input), limits) {
                    match result {
                        Ok(_) => records += 1,
                        Err(_) => failures += 1,
                    }
                }
            }
            "v3000" | "mol2" | "cml" | "cdxml" | "mmcif" => {
                let text = std::fs::read_to_string(&path)
                    .unwrap_or_else(|error| panic!("read {format} input: {error}"));
                if text.len() > max_input_bytes || max_records == 0 {
                    failures += 1;
                    continue;
                }
                let parsed = match format.as_str() {
                    "v3000" => chematic_mol::parse_mol_v3000(&text).is_ok(),
                    "mol2" => chematic_mol::parse_mol2(&text).is_ok(),
                    "cml" => chematic_mol::parse_cml(&text).is_ok(),
                    "cdxml" => chematic_mol::parse_cdxml(&text).is_ok(),
                    "mmcif" => chematic_mol::parse_mmcif(&text).is_ok(),
                    _ => unreachable!(),
                };
                if parsed {
                    records += 1;
                } else {
                    failures += 1;
                }
            }
            other => panic!(
                "unsupported format {other}; choose sdf, mol, xyz, v3000, mol2, cml, cdxml, or mmcif"
            ),
        }
    }

    let elapsed = started.elapsed().as_secs_f64();
    let input_bytes = bytes as usize * repeats;
    println!(
        "{{\"format\":\"{format}\",\"execution_mode\":\"{}\",\"path\":\"{}\",\"repeats\":{repeats},\"records\":{records},\"failures\":{failures},\"input_bytes\":{input_bytes},\"limits\":{{\"max_input_bytes\":{max_input_bytes},\"max_record_bytes\":{max_record_bytes},\"max_line_bytes\":{max_line_bytes},\"max_records\":{max_records},\"max_atoms\":{max_atoms}}},\"seconds\":{elapsed:.6},\"records_per_second\":{:.2},\"bytes_per_second\":{:.2}}}",
        if materialized {
            "materialized_one_shot"
        } else {
            "file_backed_bufread"
        },
        Path::new(&path).display(),
        records as f64 / elapsed,
        input_bytes as f64 / elapsed,
    );
}
