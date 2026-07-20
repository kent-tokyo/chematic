//! Phase B performance record (not a gate, not a Criterion-registered
//! benchmark -- see [[feedback_criterion_gate_pseudo_replication]] for why
//! this repo treats the automated Criterion CI gate's samples as
//! non-independent; this is a one-off, manually-reported measurement, not a
//! claim of statistically rigorous comparison). Compares the new suppression
//! path (`ecfp4_rdkit_environment_experimental`) against the existing
//! baseline (`ecfp4_rdkit_invariants`, same atom-invariant mode, no
//! suppression) on the same inputs: median wall time per call across N
//! repetitions (median of medians across a few independent batches), plus an
//! aggregate full-corpus pass.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example morgan_suppression_benchmark \
//!     -- <SMILES.csv>
//! ```
//! (`<SMILES.csv>` is optional -- omit to skip the full-corpus pass and only
//! run the named/synthetic per-molecule cases.)

use chematic_core::Molecule;
use chematic_fp::{ecfp4_rdkit_environment_experimental, ecfp4_rdkit_invariants};
use chematic_smiles::parse;
use std::time::{Duration, Instant};

fn median(mut durations: Vec<Duration>) -> Duration {
    durations.sort_unstable();
    durations[durations.len() / 2]
}

fn time_calls<F: FnMut()>(mut f: F, reps: usize) -> Duration {
    // Warmup, not counted.
    for _ in 0..(reps / 10).max(3) {
        f();
    }
    let mut samples = Vec::with_capacity(reps);
    for _ in 0..reps {
        let start = Instant::now();
        f();
        samples.push(start.elapsed());
    }
    median(samples)
}

fn bench_molecule(label: &str, mol: &Molecule, reps: usize) {
    let baseline = time_calls(
        || {
            ecfp4_rdkit_invariants(mol);
        },
        reps,
    );
    let suppression = time_calls(
        || {
            ecfp4_rdkit_environment_experimental(mol);
        },
        reps,
    );
    let ratio = suppression.as_secs_f64() / baseline.as_secs_f64();
    println!(
        "{label:32} atoms={:5} baseline={baseline:>12?} suppression={suppression:>12?} ratio={ratio:.3}x",
        mol.atom_count()
    );
}

fn linear_alkane(n: usize) -> Molecule {
    parse(&"C".repeat(n)).expect("linear alkane must parse")
}

/// CLI: `morgan_suppression_benchmark [<SMILES.csv>] [--corpus-only=baseline|suppression|both]`
/// `--corpus-only` restricts the full-corpus pass to just one path (for
/// isolated process-level wall-time / peak-RSS measurement via an external
/// tool like `/usr/bin/time -l`); omit it to run both in one process (for
/// the quick in-process ratio print).
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let corpus_only = args
        .iter()
        .find_map(|a| a.strip_prefix("--corpus-only=").map(str::to_string));
    let csv_path = args.iter().find(|a| !a.starts_with("--")).cloned();

    if corpus_only.is_none() {
        println!("=== per-molecule median wall time (Instant-based, not Criterion) ===");
        bench_molecule(
            "typical drug-like (aspirin)",
            &parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(),
            2000,
        );
        bench_molecule(
            "symmetric ring (benzene)",
            &parse("c1ccccc1").unwrap(),
            2000,
        );
        bench_molecule(
            "fused/polycyclic (steroid-like)",
            &parse("C[C@]12CCC3C(CCC4=CC(=O)CC[C@@]43C3CO3)C1CCC2=O").unwrap(),
            2000,
        );
        bench_molecule("~100-atom (linear alkane C100)", &linear_alkane(100), 500);
        bench_molecule("~200-atom (linear alkane C200)", &linear_alkane(200), 200);
    }

    if let Some(csv_path) = csv_path {
        let content =
            std::fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
        let mols: Vec<Molecule> = content
            .lines()
            .filter_map(|l| {
                let smi = l.trim();
                if smi.is_empty() {
                    None
                } else {
                    parse(smi).ok()
                }
            })
            .collect();

        match corpus_only.as_deref() {
            Some("baseline") => {
                let start = Instant::now();
                for m in &mols {
                    ecfp4_rdkit_invariants(m);
                }
                println!("baseline_total={:?} n={}", start.elapsed(), mols.len());
            }
            Some("suppression") => {
                let start = Instant::now();
                for m in &mols {
                    ecfp4_rdkit_environment_experimental(m);
                }
                println!("suppression_total={:?} n={}", start.elapsed(), mols.len());
            }
            _ => {
                println!(
                    "\n=== full corpus ({} molecules), single in-process pass ===",
                    mols.len()
                );
                let start = Instant::now();
                for m in &mols {
                    ecfp4_rdkit_invariants(m);
                }
                let baseline_total = start.elapsed();
                let start = Instant::now();
                for m in &mols {
                    ecfp4_rdkit_environment_experimental(m);
                }
                let suppression_total = start.elapsed();
                let ratio = suppression_total.as_secs_f64() / baseline_total.as_secs_f64();
                println!(
                    "baseline_total={baseline_total:?} suppression_total={suppression_total:?} ratio={ratio:.3}x"
                );
            }
        }
    }
}
