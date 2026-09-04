//! Repeatable canonical-SMILES throughput benchmark for a line-based corpus.
//!
//! Parsing is performed once and excluded from the timed region. Each round
//! canonicalizes the complete corpus; the reported value is the median round.
//!
//! Usage:
//! `cargo run --release -p chematic-smiles --example canonical_throughput -- corpus.smi 9`

use std::fs;
use std::hint::black_box;
use std::time::{Duration, Instant};

use chematic_smiles::{canonical_smiles, parse};

fn output_digest(mols: &[chematic_core::Molecule]) -> u64 {
    // Stable FNV-1a over length-delimited canonical outputs. This is not a
    // cryptographic digest; it is a lightweight guard that lets repeated
    // benchmark records detect accidental output changes.
    let mut digest = 0xcbf29ce484222325_u64;
    for mol in mols {
        let output = canonical_smiles(mol);
        for byte in (output.len() as u64)
            .to_le_bytes()
            .into_iter()
            .chain(output.bytes())
        {
            digest ^= u64::from(byte);
            digest = digest.wrapping_mul(0x100000001b3);
        }
    }
    digest
}

fn median(values: &mut [Duration]) -> Duration {
    values.sort_unstable();
    values[values.len() / 2]
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("pass a line-based SMILES corpus");
    let rounds: usize = args
        .next()
        .as_deref()
        .unwrap_or("9")
        .parse()
        .expect("rounds");
    assert!(
        rounds > 0 && rounds % 2 == 1,
        "rounds must be positive and odd"
    );

    let input = fs::read_to_string(&path).expect("read corpus");
    let mols: Vec<_> = input
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .filter(|smiles| !smiles.is_empty() && !smiles.starts_with('#'))
        .map(|smiles| parse(smiles).unwrap_or_else(|error| panic!("{smiles}: {error}")))
        .collect();
    assert!(!mols.is_empty(), "corpus must contain molecules");
    let digest = output_digest(&mols);

    let mut durations = Vec::with_capacity(rounds);
    let mut expected_bytes = None;
    for _ in 0..rounds {
        let started = Instant::now();
        let mut bytes = 0usize;
        for mol in &mols {
            bytes = bytes.wrapping_add(black_box(canonical_smiles(mol)).len());
        }
        durations.push(started.elapsed());
        if let Some(expected) = expected_bytes {
            assert_eq!(
                bytes, expected,
                "canonical output length changed between rounds"
            );
        } else {
            expected_bytes = Some(bytes);
        }
    }

    let elapsed = median(&mut durations);
    println!(
        "molecules={} rounds={} total_ms={:.3} us_per_molecule={:.3} output_bytes={} output_fnv1a={digest:016x}",
        mols.len(),
        rounds,
        elapsed.as_secs_f64() * 1e3,
        elapsed.as_secs_f64() * 1e6 / mols.len() as f64,
        expected_bytes.unwrap_or(0),
    );
}
