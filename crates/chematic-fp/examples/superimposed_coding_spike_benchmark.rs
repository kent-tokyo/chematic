//! Synthetic benchmark for the `superimposed_coding` SPIKE (see
//! `crates/chematic-fp/src/superimposed_coding.rs` for module docs, scope,
//! and the CC BY-NC provenance / COI note -- this benchmark measures our own
//! synthetic data only, never a claim of reproducing that paper's numbers).
//!
//! Generates synthetic `&[(feature_hash, count)]` pairs with a controlled
//! feature-overlap fraction (using a hand-rolled splitmix64 PRNG -- no new
//! dependency, this crate doesn't otherwise depend on `rand`), computes the
//! EXACT count-Tanimoto (ground truth, from `exact_count_tanimoto`) and each
//! folding strategy's binary Tanimoto on the same pairs, then reports mean
//! absolute error vs ground truth (primary metric) plus mean bit-vector
//! fill (`popcount / n_bits`, to explain *why* a strategy over/under-shoots)
//! per strategy, per `n_bits`, per overlap-level bucket.
//!
//! `superimposed_code_counts` is swept over `codeword_weight` ∈ {1,2,3,5}.
//! `codeword_weight=1` is a deliberate sanity-check arm: at that setting
//! `superimposed_code_counts` is provably identical to
//! `count_simulation_fold` (see the `superimposed_with_codeword_weight_one_*`
//! unit test) -- the *only* structural knob this spike's design adds over
//! the RDKit-style baseline is `codeword_weight > 1` (multiple independently
//! hashed bits per unary count-layer instead of one). The sweep exists to
//! show, numerically, what that knob alone buys or costs.
//!
//! A separate pooled section at the end of each `n_bits` block draws
//! `overlap_frac` uniformly at random per pair (instead of 5 fixed buckets)
//! so ground truth varies continuously across the sample -- bucketed Pearson
//! over 80 near-identical ground-truth values per bucket is dominated by
//! noise, not signal; Pearson only means something once truth has real
//! spread to correlate against.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example superimposed_coding_spike_benchmark
//! ```

use chematic_fp::bitvec::BitVecN;
use chematic_fp::superimposed_coding::{
    SuperimposedCodingConfig, count_simulation_fold, exact_count_tanimoto, fold_presence_counts,
    superimposed_code_counts,
};

/// A single feature-hash/count fingerprint, as consumed by every folding
/// strategy under test.
type CountFp = Vec<(u64, u32)>;

/// splitmix64 -- small, deterministic, dependency-free PRNG for synthetic
/// data generation only (not used anywhere in the library itself).
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `[lo, hi]` inclusive.
    fn next_range(&mut self, lo: u32, hi: u32) -> u32 {
        lo + (self.next_u64() % (hi - lo + 1) as u64) as u32
    }

    /// Uniform float in `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Generate a synthetic pair of count fingerprints sharing `overlap_frac` of
/// a `n_features`-sized universe. Shared features get independent random
/// counts in `A` and `B`; non-shared features are split evenly between the
/// two so both fingerprints have private, non-overlapping content too.
fn gen_pair(
    rng: &mut Rng,
    n_features: usize,
    overlap_frac: f64,
    max_count: u32,
) -> (CountFp, CountFp) {
    let features: Vec<u64> = (0..n_features).map(|_| rng.next_u64()).collect();
    let n_shared = (overlap_frac * n_features as f64).round() as usize;
    let mut a = Vec::new();
    let mut b = Vec::new();
    for (i, &f) in features.iter().enumerate() {
        if i < n_shared {
            a.push((f, rng.next_range(1, max_count)));
            b.push((f, rng.next_range(1, max_count)));
        } else if i % 2 == 0 {
            a.push((f, rng.next_range(1, max_count)));
        } else {
            b.push((f, rng.next_range(1, max_count)));
        }
    }
    (a, b)
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    if var_x == 0.0 || var_y == 0.0 {
        return f64::NAN;
    }
    cov / (var_x.sqrt() * var_y.sqrt())
}

fn mae(xs: &[f64], ys: &[f64]) -> f64 {
    xs.iter()
        .zip(ys.iter())
        .map(|(&x, &y)| (x - y).abs())
        .sum::<f64>()
        / xs.len() as f64
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

/// A folding function: fold a count-fingerprint pair down to two `BitVecN`s
/// of width `n_bits`.
type FoldFn = Box<dyn Fn(&[(u64, u32)], &[(u64, u32)], usize) -> (BitVecN, BitVecN)>;

/// A named folding strategy under test.
struct Arm {
    name: &'static str,
    fold: FoldFn,
}

fn arms(seed: u64, max_repeats: u32) -> Vec<Arm> {
    vec![
        Arm {
            name: "plain",
            fold: Box::new(move |a, b, n_bits| {
                (
                    fold_presence_counts(a, n_bits),
                    fold_presence_counts(b, n_bits),
                )
            }),
        },
        Arm {
            name: "count_sim",
            fold: Box::new(move |a, b, n_bits| {
                (
                    count_simulation_fold(a, n_bits, seed, max_repeats),
                    count_simulation_fold(b, n_bits, seed, max_repeats),
                )
            }),
        },
        Arm {
            name: "superimposed(w=1)",
            fold: Box::new(move |a, b, n_bits| {
                let cfg = SuperimposedCodingConfig {
                    n_bits,
                    seed,
                    repetitions: max_repeats,
                    codeword_weight: 1,
                };
                (
                    superimposed_code_counts(a, &cfg),
                    superimposed_code_counts(b, &cfg),
                )
            }),
        },
        Arm {
            name: "superimposed(w=2)",
            fold: Box::new(move |a, b, n_bits| {
                let cfg = SuperimposedCodingConfig {
                    n_bits,
                    seed,
                    repetitions: max_repeats,
                    codeword_weight: 2,
                };
                (
                    superimposed_code_counts(a, &cfg),
                    superimposed_code_counts(b, &cfg),
                )
            }),
        },
        Arm {
            name: "superimposed(w=3)",
            fold: Box::new(move |a, b, n_bits| {
                let cfg = SuperimposedCodingConfig {
                    n_bits,
                    seed,
                    repetitions: max_repeats,
                    codeword_weight: 3,
                };
                (
                    superimposed_code_counts(a, &cfg),
                    superimposed_code_counts(b, &cfg),
                )
            }),
        },
        Arm {
            name: "superimposed(w=5)",
            fold: Box::new(move |a, b, n_bits| {
                let cfg = SuperimposedCodingConfig {
                    n_bits,
                    seed,
                    repetitions: max_repeats,
                    codeword_weight: 5,
                };
                (
                    superimposed_code_counts(a, &cfg),
                    superimposed_code_counts(b, &cfg),
                )
            }),
        },
    ]
}

#[derive(Default)]
struct ArmSample {
    approx: Vec<f64>,
    fill: Vec<f64>,
}

fn print_header() {
    println!(
        "{:<20} {:>10} {:>10} {:>10} {:>10}",
        "arm", "n", "MAE", "pearson_r", "mean_fill"
    );
}

fn print_row(name: &str, exact: &[f64], sample: &ArmSample) {
    println!(
        "{:<20} {:>10} {:>10.4} {:>10.4} {:>9.1}%",
        name,
        exact.len(),
        mae(exact, &sample.approx),
        pearson(exact, &sample.approx),
        mean(&sample.fill) * 100.0,
    );
}

fn main() {
    let n_bits_values = [512usize, 1024, 2048];
    let overlap_levels: [f64; 5] = [0.1, 0.3, 0.5, 0.7, 0.9];
    let n_pairs_per_level = 80;
    let n_pooled_pairs = 400;
    let n_features = 400;
    let max_count = 12;
    let seed_base = 20260809u64;
    let max_repeats = 4u32;

    println!(
        "superimposed_coding SPIKE synthetic benchmark -- {n_pairs_per_level} pairs/level, {n_features} features, max_count={max_count}, max_repeats={max_repeats}"
    );
    println!(
        "Primary metric is MAE (mean absolute error vs exact count-Tanimoto). Bucketed-overlap\n\
         Pearson is included for completeness but is noise-dominated (ground truth barely varies\n\
         within one fixed overlap level); see the pooled random-overlap section for a meaningful\n\
         correlation figure.\n"
    );

    for &n_bits in &n_bits_values {
        println!("======================== n_bits = {n_bits} ========================");

        // ---- bucketed overlap levels: MAE + fill per level ----
        let arm_defs = arms(seed_base, max_repeats);

        for &overlap in &overlap_levels {
            let mut rng = Rng::new(seed_base ^ overlap.to_bits() ^ (n_bits as u64));
            let mut exact = Vec::with_capacity(n_pairs_per_level);
            let mut per_arm: Vec<ArmSample> =
                arm_defs.iter().map(|_| ArmSample::default()).collect();

            for _ in 0..n_pairs_per_level {
                let (a, b) = gen_pair(&mut rng, n_features, overlap, max_count);
                exact.push(exact_count_tanimoto(&a, &b));
                for (arm, sample) in arm_defs.iter().zip(per_arm.iter_mut()) {
                    let (fa, fb) = (arm.fold)(&a, &b, n_bits);
                    sample.approx.push(fa.tanimoto(&fb));
                    sample.fill.push(
                        (fa.popcount() as f64 + fb.popcount() as f64) / (2.0 * n_bits as f64),
                    );
                }
            }

            println!("--- overlap = {overlap:.2} ---");
            print_header();
            for (arm, sample) in arm_defs.iter().zip(per_arm.iter()) {
                print_row(arm.name, &exact, sample);
            }
            println!();
        }

        // ---- pooled continuous-overlap section: meaningful Pearson ----
        let mut rng = Rng::new(seed_base ^ 0xC0FFEE ^ (n_bits as u64));
        let mut exact_pooled = Vec::with_capacity(n_pooled_pairs);
        let mut per_arm_pooled: Vec<ArmSample> =
            arm_defs.iter().map(|_| ArmSample::default()).collect();
        for _ in 0..n_pooled_pairs {
            let overlap = 0.05 + rng.next_f64() * 0.9; // uniform in [0.05, 0.95)
            let (a, b) = gen_pair(&mut rng, n_features, overlap, max_count);
            exact_pooled.push(exact_count_tanimoto(&a, &b));
            for (arm, sample) in arm_defs.iter().zip(per_arm_pooled.iter_mut()) {
                let (fa, fb) = (arm.fold)(&a, &b, n_bits);
                sample.approx.push(fa.tanimoto(&fb));
                sample
                    .fill
                    .push((fa.popcount() as f64 + fb.popcount() as f64) / (2.0 * n_bits as f64));
            }
        }
        println!("--- pooled, overlap_frac ~ Uniform[0.05, 0.95) per pair, n={n_pooled_pairs} ---");
        print_header();
        for (arm, sample) in arm_defs.iter().zip(per_arm_pooled.iter()) {
            print_row(arm.name, &exact_pooled, sample);
        }
        println!();
    }
}
