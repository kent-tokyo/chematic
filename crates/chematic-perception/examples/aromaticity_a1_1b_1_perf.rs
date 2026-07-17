//! Aromaticity-A1-1b-1: record (not optimize) per-molecule wall-clock
//! timing for `apply_aromaticity_ex(.., RdkitLike)` (current production
//! default) vs `apply_aromaticity_rdkit_parity_experimental` (the new
//! opt-in engine), over the full SMILES corpus.
//!
//! This PR does not tune either engine for speed -- these numbers are a
//! baseline record for a future performance-focused round (A1-1c), not a
//! gate this PR is judged against.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_a1_1b_1_perf -- ~/Downloads/SMILES.csv
//! ```

use std::fs;
use std::time::Instant;

use chematic_perception::apply_aromaticity_rdkit_parity_experimental;
use chematic_perception::{AromaticityAlgorithm, apply_aromaticity_ex};

fn percentile(sorted_ns: &[u64], p: f64) -> u64 {
    if sorted_ns.is_empty() {
        return 0;
    }
    let idx = ((sorted_ns.len() as f64 - 1.0) * p).round() as usize;
    sorted_ns[idx]
}

fn report(label: &str, mut ns: Vec<u64>) {
    ns.sort_unstable();
    let n = ns.len() as u64;
    let sum: u64 = ns.iter().sum();
    let mean = sum as f64 / n as f64 / 1000.0;
    let p50 = percentile(&ns, 0.50) as f64 / 1000.0;
    let p95 = percentile(&ns, 0.95) as f64 / 1000.0;
    let max = *ns.last().unwrap_or(&0) as f64 / 1000.0;
    println!("{label}: n={n} mean={mean:.2}us p50={p50:.2}us p95={p95:.2}us max={max:.2}us");
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: aromaticity_a1_1b_1_perf <smiles.csv>");
    let content =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let smiles: Vec<&str> = content
        .lines()
        .map(|l| l.split(',').next().unwrap_or("").trim())
        .filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("smiles"))
        .collect();

    let mols: Vec<_> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    println!("parsed {}/{} molecules", mols.len(), smiles.len());

    // Warm-up pass (page faults, allocator warm-up) -- excluded from timing.
    for mol in &mols {
        let _ = apply_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike);
        let _ = apply_aromaticity_rdkit_parity_experimental(mol);
    }

    let mut default_ns = Vec::with_capacity(mols.len());
    for mol in &mols {
        let start = Instant::now();
        let _ = apply_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike);
        default_ns.push(start.elapsed().as_nanos() as u64);
    }
    report("RdkitLike (current production default)", default_ns);

    let mut experimental_ns = Vec::with_capacity(mols.len());
    let mut n_kekulize_failed = 0usize;
    for mol in &mols {
        let start = Instant::now();
        let result = apply_aromaticity_rdkit_parity_experimental(mol);
        experimental_ns.push(start.elapsed().as_nanos() as u64);
        if result.is_err() {
            n_kekulize_failed += 1;
        }
    }
    report("RdkitParityExperimental", experimental_ns);
    println!("(experimental KekulizationFailed count during timing: {n_kekulize_failed})");
}
