//! Measure file-backed SDF/MOL/XYZ/Extended XYZ record streaming in the current workspace.
//!
//! Usage:
//! `cargo run -p chematic-mol --example streaming_benchmark -- --format sdf path`
//!
//! The JSON result includes aggregate `records`/`failures` plus deterministic
//! `failure_kinds` counts keyed by parser error variant.

use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::time::Instant;

use chematic_mol::{ExtxyzFileReader, SdfFileReader, XyzFileReader};
use flate2::read::GzDecoder;

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

fn has_flag(name: &str) -> bool {
    env::args().skip(1).any(|value| value == name)
}

fn read_text(path: &str, gzip: bool, max_input_bytes: usize) -> Option<String> {
    let file = File::open(path).unwrap_or_else(|error| panic!("open {path}: {error}"));
    let mut text = String::new();
    if gzip {
        let decoder = GzDecoder::new(file);
        decoder
            .take(max_input_bytes.saturating_add(1) as u64)
            .read_to_string(&mut text)
            .unwrap_or_else(|error| panic!("decompress {path}: {error}"));
    } else {
        file.take(max_input_bytes.saturating_add(1) as u64)
            .read_to_string(&mut text)
            .unwrap_or_else(|error| panic!("read {path}: {error}"));
    }
    (text.len() <= max_input_bytes).then_some(text)
}

fn count_sdf<R: BufRead>(
    reader: R,
    limits: chematic_mol::SdfParseLimits,
) -> (usize, usize, BTreeMap<String, usize>) {
    SdfFileReader::with_limits(reader, limits).fold(
        (0, 0, BTreeMap::new()),
        |(records, failures, mut kinds), result| {
            if result.is_ok() {
                (records + 1, failures, kinds)
            } else {
                let kind = match result.as_ref() {
                    Err(error) => format!("{error:?}")
                        .split(['{', '('])
                        .next()
                        .unwrap_or("unknown")
                        .trim()
                        .to_string(),
                    Ok(_) => unreachable!("result was checked as an error"),
                };
                *kinds.entry(kind).or_default() += 1;
                (records, failures + 1, kinds)
            }
        },
    )
}

fn count_xyz<R: BufRead>(
    reader: R,
    limits: chematic_mol::XyzParseLimits,
) -> (usize, usize, BTreeMap<String, usize>) {
    XyzFileReader::with_limits(reader, limits).fold(
        (0, 0, BTreeMap::new()),
        |(records, failures, mut kinds), result| {
            if result.is_ok() {
                (records + 1, failures, kinds)
            } else {
                let kind = match result.as_ref() {
                    Err(error) => format!("{error:?}")
                        .split(['{', '('])
                        .next()
                        .unwrap_or("unknown")
                        .trim()
                        .to_string(),
                    Ok(_) => unreachable!("result was checked as an error"),
                };
                *kinds.entry(kind).or_default() += 1;
                (records, failures + 1, kinds)
            }
        },
    )
}

fn count_extxyz<R: BufRead>(
    reader: R,
    limits: chematic_mol::XyzParseLimits,
) -> (usize, usize, BTreeMap<String, usize>) {
    ExtxyzFileReader::with_limits(reader, limits).fold(
        (0, 0, BTreeMap::new()),
        |(records, failures, mut kinds), result| {
            if result.is_ok() {
                (records + 1, failures, kinds)
            } else {
                let kind = match result.as_ref() {
                    Err(error) => format!("{error:?}")
                        .split(['{', '('])
                        .next()
                        .unwrap_or("unknown")
                        .trim()
                        .to_string(),
                    Ok(_) => unreachable!("result was checked as an error"),
                };
                *kinds.entry(kind).or_default() += 1;
                (records, failures + 1, kinds)
            }
        },
    )
}

fn merge_failure_kinds(into: &mut BTreeMap<String, usize>, from: BTreeMap<String, usize>) {
    for (kind, count) in from {
        *into.entry(kind).or_default() += count;
    }
}

fn record_result_failure<T, E: std::fmt::Debug>(
    result: &Result<T, E>,
    failures: &mut BTreeMap<String, usize>,
) {
    if let Err(error) = result {
        let kind = format!("{error:?}")
            .split(['{', '('])
            .next()
            .unwrap_or("unknown")
            .trim()
            .to_string();
        *failures.entry(kind).or_default() += 1;
    }
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
    let gzip = has_flag("--gzip");
    let materialized = matches!(
        format.as_str(),
        "v3000" | "mol2" | "cml" | "cdxml" | "mmcif" | "pdb"
    );
    let bytes = std::fs::metadata(&path)
        .expect("benchmark input must exist")
        .len();
    let started = Instant::now();
    let mut records = 0usize;
    let mut failures = 0usize;
    let mut failure_kinds = BTreeMap::<String, usize>::new();

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
                let (ok, bad, kinds) = if gzip {
                    count_sdf(BufReader::new(GzDecoder::new(input)), limits)
                } else {
                    count_sdf(BufReader::new(input), limits)
                };
                records += ok;
                failures += bad;
                merge_failure_kinds(&mut failure_kinds, kinds);
            }
            "xyz" => {
                let input = File::open(&path).expect("open XYZ input");
                let limits = chematic_mol::XyzParseLimits {
                    max_input_bytes,
                    max_atoms_per_frame: max_atoms,
                    max_frames: max_records,
                    max_line_bytes,
                };
                let (ok, bad, kinds) = if gzip {
                    count_xyz(BufReader::new(GzDecoder::new(input)), limits)
                } else {
                    count_xyz(BufReader::new(input), limits)
                };
                records += ok;
                failures += bad;
                merge_failure_kinds(&mut failure_kinds, kinds);
            }
            "extxyz" => {
                let input = File::open(&path).expect("open Extended XYZ input");
                let limits = chematic_mol::XyzParseLimits {
                    max_input_bytes,
                    max_atoms_per_frame: max_atoms,
                    max_frames: max_records,
                    max_line_bytes,
                };
                let (ok, bad, kinds) = if gzip {
                    count_extxyz(BufReader::new(GzDecoder::new(input)), limits)
                } else {
                    count_extxyz(BufReader::new(input), limits)
                };
                records += ok;
                failures += bad;
                merge_failure_kinds(&mut failure_kinds, kinds);
            }
            "v3000" | "mol2" | "cml" | "cdxml" | "mmcif" => {
                let Some(text) = read_text(&path, gzip, max_input_bytes) else {
                    failures += 1;
                    *failure_kinds
                        .entry("ResourceLimit".to_string())
                        .or_default() += 1;
                    continue;
                };
                if max_records == 0 {
                    failures += 1;
                    *failure_kinds
                        .entry("ResourceLimit".to_string())
                        .or_default() += 1;
                    continue;
                }
                let parsed = match format.as_str() {
                    "v3000" => {
                        let result = chematic_mol::parse_mol_v3000(&text);
                        record_result_failure(&result, &mut failure_kinds);
                        result.is_ok()
                    }
                    "mol2" => {
                        let result = chematic_mol::parse_mol2_with_limits(
                            &text,
                            &chematic_mol::Mol2ParseLimits {
                                max_input_bytes,
                                max_line_bytes,
                                max_lines: max_records,
                                max_atoms,
                                max_bonds: max_atoms.saturating_mul(2),
                                ..chematic_mol::Mol2ParseLimits::default()
                            },
                        );
                        record_result_failure(&result, &mut failure_kinds);
                        result.is_ok()
                    }
                    "cml" => {
                        let result = chematic_mol::parse_cml_with_limits(
                            &text,
                            &chematic_mol::CmlParseLimits {
                                max_input_bytes,
                                max_line_bytes,
                                max_lines: max_records,
                                max_atoms,
                                max_bonds: max_atoms.saturating_mul(2),
                                max_xml_elements: max_records.saturating_mul(4),
                            },
                        );
                        record_result_failure(&result, &mut failure_kinds);
                        result.is_ok()
                    }
                    "cdxml" => {
                        let result = chematic_mol::parse_cdxml_with_limits(
                            &text,
                            &chematic_mol::CdxmlParseLimits {
                                max_input_bytes,
                                max_line_bytes,
                                max_lines: max_records,
                                max_attribute_bytes: max_line_bytes,
                                max_atoms,
                                max_bonds: max_atoms.saturating_mul(2),
                                max_fragments: max_records,
                            },
                        );
                        record_result_failure(&result, &mut failure_kinds);
                        result.is_ok()
                    }
                    "mmcif" => {
                        let result = chematic_mol::parse_mmcif_with_limits(
                            &text,
                            &chematic_mol::MmcifParseLimits {
                                max_input_bytes,
                                max_atoms,
                                max_line_len: max_line_bytes,
                            },
                        );
                        record_result_failure(&result, &mut failure_kinds);
                        result.is_ok()
                    }
                    _ => unreachable!(),
                };
                if parsed {
                    records += 1;
                } else {
                    failures += 1;
                }
            }
            "pdb" => {
                let Some(text) = read_text(&path, gzip, max_input_bytes) else {
                    failures += 1;
                    *failure_kinds
                        .entry("ResourceLimit".to_string())
                        .or_default() += 1;
                    continue;
                };
                if max_records == 0 {
                    failures += 1;
                    *failure_kinds
                        .entry("ResourceLimit".to_string())
                        .or_default() += 1;
                    continue;
                }
                let parsed = chematic_3d::parse_pdb_atoms_with_limits(
                    &text,
                    &chematic_3d::PdbParseLimits {
                        max_input_bytes,
                        max_line_bytes,
                        max_atoms,
                        max_models: max_records,
                    },
                );
                record_result_failure(&parsed, &mut failure_kinds);
                if parsed.is_ok() {
                    records += 1;
                } else {
                    failures += 1;
                }
            }
            other => panic!(
                "unsupported format {other}; choose sdf, mol, xyz, extxyz, v3000, mol2, cml, cdxml, mmcif, or pdb"
            ),
        }
    }

    let elapsed = started.elapsed().as_secs_f64();
    let input_bytes = bytes as usize * repeats;
    println!(
        "{{\"format\":\"{format}\",\"compression\":\"{}\",\"execution_mode\":\"{}\",\"path\":\"{}\",\"repeats\":{repeats},\"records\":{records},\"failures\":{failures},\"failure_kinds\":{},\"input_bytes\":{input_bytes},\"limits\":{{\"max_input_bytes\":{max_input_bytes},\"max_record_bytes\":{max_record_bytes},\"max_line_bytes\":{max_line_bytes},\"max_records\":{max_records},\"max_atoms\":{max_atoms}}},\"seconds\":{elapsed:.6},\"records_per_second\":{:.2},\"bytes_per_second\":{:.2}}}",
        if gzip { "gzip" } else { "none" },
        if materialized {
            "materialized_one_shot"
        } else {
            "file_backed_bufread"
        },
        Path::new(&path).display(),
        serde_json::to_string(&failure_kinds).expect("failure kinds serialize"),
        records as f64 / elapsed,
        input_bytes as f64 / elapsed,
    );
}
