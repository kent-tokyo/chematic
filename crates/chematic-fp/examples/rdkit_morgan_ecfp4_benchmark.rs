//! Phase B performance record (not a gate, not a Criterion-registered
//! benchmark -- see `feedback_criterion_gate_pseudo_replication` for why this
//! repo treats the automated Criterion CI gate's samples as non-independent;
//! this is a one-off, manually-reported measurement, matching
//! `morgan_suppression_benchmark.rs`'s own precedent). Compares
//! `rdkit_morgan_ecfp4_experimental` (candidate: RDKit-bit-exact, does its
//! own kekulization + RDKit-parity aromaticity perception internally every
//! call) against `ecfp4_rdkit_environment_experimental` (baseline: reads
//! whatever aromatic flags are already on the input `Molecule`, no
//! aromaticity engine call of its own) on the same inputs.
//!
//! These two functions do genuinely different amounts of preprocessing work
//! by design (see the module docs on `rdkit_morgan_ecfp4.rs` -- the
//! bit-exactness guarantee *requires* the candidate's own aromaticity
//! perception step). A `--split` run reports how much of the candidate's
//! time is preprocessing (`apply_aromaticity_rdkit_parity_experimental`
//! alone) vs. the hash expansion itself, so a >2x ratio is attributable
//! rather than treated as an undifferentiated "regression."
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example rdkit_morgan_ecfp4_benchmark \
//!     -- <SMILES.csv> [--corpus-only=baseline|candidate|preprocessing]
//! ```
//! (`<SMILES.csv>` is optional -- omit to skip the full-corpus pass and only
//! run the named/synthetic per-molecule cases. `--corpus-only=...` runs one
//! side only, for external multi-process median/RSS measurement.)

use chematic_core::Molecule;
use chematic_fp::{ecfp4_rdkit_environment_experimental, rdkit_morgan_ecfp4_experimental};
use chematic_perception::apply_aromaticity_rdkit_parity_experimental;
use chematic_smiles::parse;
use std::time::{Duration, Instant};

fn median(mut durations: Vec<Duration>) -> Duration {
    durations.sort_unstable();
    durations[durations.len() / 2]
}

fn time_calls<F: FnMut()>(mut f: F, reps: usize) -> Duration {
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
            ecfp4_rdkit_environment_experimental(mol);
        },
        reps,
    );
    let preprocessing = time_calls(
        || {
            let _ = apply_aromaticity_rdkit_parity_experimental(mol);
        },
        reps,
    );
    let candidate = time_calls(
        || {
            let _ = rdkit_morgan_ecfp4_experimental(mol);
        },
        reps,
    );
    let ratio = candidate.as_secs_f64() / baseline.as_secs_f64();
    let preprocessing_share = preprocessing.as_secs_f64() / candidate.as_secs_f64();
    println!(
        "{label:32} atoms={:5} baseline={baseline:>12?} candidate={candidate:>12?} \
         ratio={ratio:.3}x preprocessing={preprocessing:>12?} ({:.0}% of candidate)",
        mol.atom_count(),
        preprocessing_share * 100.0,
    );
}

fn linear_alkane(n: usize) -> Molecule {
    parse(&"C".repeat(n)).unwrap()
}

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
                    ecfp4_rdkit_environment_experimental(m);
                }
                println!("baseline_total={:?} n={}", start.elapsed(), mols.len());
            }
            Some("candidate") => {
                let start = Instant::now();
                let mut errors = 0usize;
                for m in &mols {
                    if rdkit_morgan_ecfp4_experimental(m).is_err() {
                        errors += 1;
                    }
                }
                println!(
                    "candidate_total={:?} n={} errors={errors}",
                    start.elapsed(),
                    mols.len()
                );
            }
            Some("preprocessing") => {
                let start = Instant::now();
                for m in &mols {
                    let _ = apply_aromaticity_rdkit_parity_experimental(m);
                }
                println!("preprocessing_total={:?} n={}", start.elapsed(), mols.len());
            }
            _ => {
                println!(
                    "\n=== full corpus ({} molecules), single in-process pass ===",
                    mols.len()
                );
                let start = Instant::now();
                for m in &mols {
                    ecfp4_rdkit_environment_experimental(m);
                }
                let baseline_total = start.elapsed();

                let start = Instant::now();
                for m in &mols {
                    let _ = apply_aromaticity_rdkit_parity_experimental(m);
                }
                let preprocessing_total = start.elapsed();

                let start = Instant::now();
                let mut errors = 0usize;
                for m in &mols {
                    if rdkit_morgan_ecfp4_experimental(m).is_err() {
                        errors += 1;
                    }
                }
                let candidate_total = start.elapsed();

                let ratio = candidate_total.as_secs_f64() / baseline_total.as_secs_f64();
                let preprocessing_share =
                    preprocessing_total.as_secs_f64() / candidate_total.as_secs_f64();
                println!(
                    "baseline_total={baseline_total:?} candidate_total={candidate_total:?} ratio={ratio:.3}x \
                     preprocessing_total={preprocessing_total:?} ({:.0}% of candidate) errors={errors}",
                    preprocessing_share * 100.0,
                );
            }
        }
    }
}
