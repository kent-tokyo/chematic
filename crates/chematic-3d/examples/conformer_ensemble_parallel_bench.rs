//! Reproducible sequential-vs-parallel conformer-ensemble microbenchmark.
//!
//! This intentionally runs both paths in alternating order and asserts exact
//! equality before reporting a speed ratio. It is a development benchmark,
//! not a universal performance claim.

use std::hint::black_box;
use std::time::{Duration, Instant};

use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::{
    ConformerDisposition, EnsembleV2Config, EnsembleV2Result, PipelineV2Config, embed_ensemble_v2,
    embed_ensemble_v2_parallel,
};
use chematic_smiles::parse;

const SMILES: &[&str] = &[
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
];

fn config() -> EnsembleV2Config {
    let mut per_conformer = PipelineV2Config::stereo_safe(ForceFieldPolicy::Mmff94BondAngleStrict);
    per_conformer.force_field_max_iterations = 200;
    let mut config = EnsembleV2Config::new(per_conformer, 10, 42);
    config.rmsd_threshold = 0.5;
    config
}

fn run(parallel: bool) -> Duration {
    let started = Instant::now();
    for source in SMILES {
        let mol = parse(source).expect("benchmark SMILES");
        let config = config();
        let result = if parallel {
            embed_ensemble_v2_parallel(black_box(&mol), black_box(&config))
        } else {
            embed_ensemble_v2(black_box(&mol), black_box(&config))
        }
        .expect("valid ensemble config");
        black_box(result);
    }
    started.elapsed()
}

fn assert_equal(left: &EnsembleV2Result, right: &EnsembleV2Result) {
    assert_eq!(left.termination, right.termination);
    assert_eq!(left.mixed_force_field, right.mixed_force_field);
    assert_eq!(left.attempts.len(), right.attempts.len());
    assert_eq!(
        left.ensemble.conformer_count(),
        right.ensemble.conformer_count()
    );
    for (a, b) in left.attempts.iter().zip(&right.attempts) {
        assert_eq!((a.attempt_index, a.seed), (b.attempt_index, b.seed));
        match (&a.outcome, &b.outcome) {
            (Ok(a), Ok(b)) => {
                assert_eq!(a.energy, b.energy);
                assert_eq!(a.actual_force_field_used, b.actual_force_field_used);
                match (&a.disposition, &b.disposition) {
                    (
                        ConformerDisposition::Kept { conformer_index: a },
                        ConformerDisposition::Kept { conformer_index: b },
                    ) => assert_eq!(a, b),
                    (
                        ConformerDisposition::PrunedAsDuplicate {
                            representative_attempt_index: ar,
                            rmsd: ad,
                            symmetric: asym,
                        },
                        ConformerDisposition::PrunedAsDuplicate {
                            representative_attempt_index: br,
                            rmsd: bd,
                            symmetric: bsym,
                        },
                    ) => assert_eq!((ar, ad.to_bits(), asym), (br, bd.to_bits(), bsym)),
                    _ => panic!("conformer dispositions differ"),
                }
            }
            (Err(_), Err(_)) => {}
            _ => panic!("attempt outcome differs"),
        }
    }
    for i in 0..left.ensemble.conformer_count() {
        assert_eq!(
            left.ensemble.get_conformer(i).unwrap().points,
            right.ensemble.get_conformer(i).unwrap().points
        );
    }
}

fn median(values: &mut [Duration]) -> Duration {
    values.sort_unstable();
    values[values.len() / 2]
}

fn main() {
    let repeats = std::env::var("SCHEMATIC_BENCH_REPEATS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(5)
        .max(1);
    for source in SMILES {
        let mol = parse(source).expect("benchmark SMILES");
        let config = config();
        assert_equal(
            &embed_ensemble_v2(&mol, &config).expect("sequential ensemble"),
            &embed_ensemble_v2_parallel(&mol, &config).expect("parallel ensemble"),
        );
    }
    let mut sequential = Vec::with_capacity(repeats);
    let mut parallel = Vec::with_capacity(repeats);
    for block in 0..repeats {
        if block % 2 == 0 {
            sequential.push(run(false));
            parallel.push(run(true));
        } else {
            parallel.push(run(true));
            sequential.push(run(false));
        }
    }
    let sequential_samples_ms: Vec<f64> = sequential
        .iter()
        .map(|duration| duration.as_secs_f64() * 1_000.0)
        .collect();
    let parallel_samples_ms: Vec<f64> = parallel
        .iter()
        .map(|duration| duration.as_secs_f64() * 1_000.0)
        .collect();
    let sequential = median(&mut sequential);
    let parallel = median(&mut parallel);
    println!(
        "{}",
        serde_json::json!({
            "molecules": SMILES.len(),
            "conformers_per_molecule": 10,
            "repeats": repeats,
            "rayon_threads": rayon::current_num_threads(),
            "sequential_samples_ms": sequential_samples_ms,
            "parallel_samples_ms": parallel_samples_ms,
            "sequential_median_ms": sequential.as_secs_f64() * 1_000.0,
            "parallel_median_ms": parallel.as_secs_f64() * 1_000.0,
            "speedup": sequential.as_secs_f64() / parallel.as_secs_f64(),
        })
    );
}
