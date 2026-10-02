//! Timing-only worker compiled against one pinned published crates.io graph.
//! The parent script builds both versions with `--locked` and checks outputs.

use std::env;
use std::fs;
use std::hint::black_box;
use std::time::Instant;

use chematic::{chem, core, fp, smiles};
use serde_json::json;
use sha2::{Digest, Sha256};

// The fingerprint value is intentionally stored inline: boxing would add a
// per-molecule allocation that is not part of the published library call.
#[allow(clippy::large_enum_variant)]
enum Output {
    Hba(u32),
    Morgan(fp::BitVec2048),
}

fn evaluate(operation: &str, mol: &core::Molecule) -> Result<Output, Box<dyn std::error::Error>> {
    match operation {
        "hba" => Ok(Output::Hba(chem::rdkit_hba_count(mol) as u32)),
        "morgan" => Ok(Output::Morgan(fp::rdkit_morgan_ecfp4_bitvec(mol)?)),
        _ => Err("operation must be hba or morgan".into()),
    }
}

fn encode(value: Output) -> String {
    match value {
        Output::Hba(count) => count.to_string(),
        Output::Morgan(bits) => {
            let mut out = String::with_capacity(512);
            for byte_index in 0..256 {
                let mut byte = 0u8;
                for bit in 0..8 {
                    if bits.get(byte_index * 8 + bit) {
                        byte |= 1 << bit;
                    }
                }
                out.push_str(&format!("{byte:02x}"));
            }
            out
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 5 {
        return Err("usage: isolated_speed CORPUS 5000 hba|morgan parse_inclusive|prepared_first_use|precomputed".into());
    }
    let limit: usize = args[2].parse()?;
    if limit == 0 {
        return Err("limit must be positive".into());
    }
    let operation = args[3].as_str();
    let mode = args[4].as_str();
    if !matches!(
        mode,
        "parse_inclusive" | "prepared_first_use" | "precomputed"
    ) {
        return Err("unknown measurement mode".into());
    }
    let corpus = fs::read_to_string(&args[1])?;
    let rows: Vec<_> = corpus
        .lines()
        .take(limit)
        .map(|line| line.split_whitespace().next().unwrap_or("").to_owned())
        .collect();
    if rows.len() != limit || rows.iter().any(String::is_empty) {
        return Err("corpus has fewer than requested rows".into());
    }

    // Disjoint from the timed objects. Per-molecule caches are not prewarmed
    // unless the mode explicitly says so.
    let warm = smiles::parse("CCO")?;
    black_box(evaluate(operation, &warm)?);

    let setup_start = Instant::now();
    let prepared = if mode == "parse_inclusive" {
        None
    } else {
        Some(
            rows.iter()
                .map(|smiles| smiles::parse(smiles))
                .collect::<Result<Vec<_>, _>>()?,
        )
    };
    if mode == "precomputed" {
        for mol in prepared.as_ref().expect("precomputed molecules") {
            black_box(evaluate(operation, mol)?);
        }
    }
    let setup_ns = setup_start.elapsed().as_nanos();

    let started = Instant::now();
    let values: Vec<Output> = if let Some(molecules) = &prepared {
        molecules
            .iter()
            .map(|mol| evaluate(operation, mol))
            .collect::<Result<_, _>>()?
    } else {
        rows.iter()
            .map(|text| evaluate(operation, &smiles::parse(text)?))
            .collect::<Result<_, _>>()?
    };
    let operation_ns = started.elapsed().as_nanos();
    let mut bytes = Vec::new();
    for value in values {
        bytes.extend_from_slice(encode(value).as_bytes());
        bytes.push(b'\n');
    }
    println!(
        "{}",
        json!({
            "operation": operation,
            "mode": mode,
            "input_count": limit,
            "setup_ns": setup_ns,
            "operation_ns": operation_ns,
            "output_sha256": format!("{:x}", Sha256::digest(&bytes)),
        })
    );
    Ok(())
}
