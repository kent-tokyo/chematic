//! One isolated bounded parser invocation for the parser-security gate.
//!
//! This intentionally reports malformed input as JSON instead of panicking or
//! using process failure as a parser result. The Python gate applies the
//! process-level wall-time and address-space limits around this executable.

use std::{env, fs};

#[cfg(target_os = "linux")]
fn peak_rss_kib() -> Option<u64> {
    fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

#[cfg(not(target_os = "linux"))]
fn peak_rss_kib() -> Option<u64> {
    None
}

fn peak_rss_json_value() -> String {
    peak_rss_kib()
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string())
}

fn argument(name: &str) -> String {
    let mut args = env::args().skip(1);
    while let Some(value) = args.next() {
        if value == name {
            return args
                .next()
                .unwrap_or_else(|| panic!("{name} requires a value"));
        }
    }
    panic!("missing {name}");
}

fn main() {
    let format = argument("--format");
    let input_path = argument("--input");
    let bytes = match fs::read(&input_path) {
        Ok(bytes) if bytes.len() <= 1 << 20 => bytes,
        Ok(_) => {
            println!(
                "{{\"format\":{format:?},\"status\":\"rejected\",\"kind\":\"input_bytes_limit\",\"peak_rss_kib\":{}}}",
                peak_rss_json_value()
            );
            return;
        }
        Err(error) => panic!("read {input_path}: {error}"),
    };
    let text = match String::from_utf8(bytes) {
        Ok(text) => text,
        Err(_) => {
            println!(
                "{{\"format\":{format:?},\"status\":\"rejected\",\"kind\":\"utf8\",\"peak_rss_kib\":{}}}",
                peak_rss_json_value()
            );
            return;
        }
    };
    let accepted = match format.as_str() {
        "smiles" => chematic_smiles::parse_with_limits(
            &text,
            &chematic_smiles::SmilesParseLimits {
                max_input_bytes: 1 << 20,
                max_atoms: 10_000,
                max_bonds: 20_000,
            },
        )
        .is_ok(),
        "smarts" => chematic_smarts::parse_smarts(&text).is_ok(),
        "mol_v2000" => chematic_mol::parse_mol(&text).is_ok(),
        "mol_v3000" => chematic_mol::parse_mol_v3000(&text).is_ok(),
        "sdf" => chematic_mol::parse_sdf_with_limits(
            &text,
            chematic_mol::SdfParseLimits {
                max_input_bytes: 1 << 20,
                max_record_bytes: 1 << 20,
                max_line_bytes: 1 << 20,
                max_records: 1,
            },
        )
        .is_ok(),
        other => panic!("unsupported format {other}"),
    };
    let status = if accepted { "accepted" } else { "rejected" };
    println!(
        "{{\"format\":{format:?},\"status\":{status:?},\"peak_rss_kib\":{}}}",
        peak_rss_json_value()
    );
}
