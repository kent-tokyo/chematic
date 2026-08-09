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
//! show, numerically, what that knob alone buys or costs. See
//! `superimposed_coding`'s module docs for a short proof that
//! `codeword_weight` is exactly Tanimoto-neutral in the collision-free
//! limit, so any effect visible here is collision/saturation noise, never a
//! genuine improvement.
//!
//! Every count/count_sim/superimposed arm clamps count at `max_repeats`, so
//! their best-case output (even with zero hash collisions) is the exact
//! Tanimoto of *clamped* counts, not the unclamped ground truth. A
//! `ceiling(clamp@R)` row -- `exact_count_tanimoto` on both fingerprints'
//! counts clamped to `max_repeats`, no hashing/folding at all -- is printed
//! alongside every table so the clamp-bias contribution to each arm's MAE
//! can be told apart from hashing/collision noise.
//!
//! Two regimes are covered per `n_bits`: bucketed fixed-overlap levels (MAE
//! only -- Pearson over 80 near-identical ground-truth values per bucket is
//! noise-dominated) and a pooled continuous-overlap sample (`overlap_frac ~
//! Uniform[0.05, 0.95)` per pair, meaningful Pearson). A final section
//! re-runs the pooled comparison in a smaller, clamp-free regime
//! (`n_features=60`, `max_count=max_repeats`) closer to real ECFP4 fill
//! levels (~5-10% of 2048 bits for drug-like molecules, vs. ~15-45% for the
//! `n_features=400` headline numbers) to check whether the ordering
//! observed at high, unrealistic fill survives at realistic fill.
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
/// counts in `A` and `B` (a simplification: real similar molecules would
/// have *correlated* counts on shared features, which would push ground
/// truth higher at high overlap than this generator does); non-shared
/// features are split evenly between the two so both fingerprints have
/// private, non-overlapping content too.
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

/// Clamp every count in `fp` to `cap` -- the best-case (zero-collision)
/// output shape of any `count_sim`/`superimposed` arm with `max_repeats ==
/// cap`.
fn clamp_counts(fp: &[(u64, u32)], cap: u32) -> CountFp {
    fp.iter().map(|&(f, c)| (f, c.min(cap))).collect()
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

struct BatchResult {
    exact: Vec<f64>,
    ceiling: Vec<f64>,
    per_arm: Vec<ArmSample>,
}

/// Run `n_pairs` synthetic pairs through every arm plus the clamp-bias
/// ceiling. `overlap_fn` supplies `overlap_frac` per pair (constant for a
/// fixed-overlap bucket, random for the pooled sections).
#[allow(clippy::too_many_arguments)]
fn run_batch(
    seed: u64,
    n_features: usize,
    max_count: u32,
    max_repeats: u32,
    n_bits: usize,
    n_pairs: usize,
    arm_defs: &[Arm],
    mut overlap_fn: impl FnMut(&mut Rng) -> f64,
) -> BatchResult {
    let mut rng = Rng::new(seed);
    let mut exact = Vec::with_capacity(n_pairs);
    let mut ceiling = Vec::with_capacity(n_pairs);
    let mut per_arm: Vec<ArmSample> = arm_defs.iter().map(|_| ArmSample::default()).collect();

    for _ in 0..n_pairs {
        let overlap = overlap_fn(&mut rng);
        let (a, b) = gen_pair(&mut rng, n_features, overlap, max_count);
        exact.push(exact_count_tanimoto(&a, &b));
        ceiling.push(exact_count_tanimoto(
            &clamp_counts(&a, max_repeats),
            &clamp_counts(&b, max_repeats),
        ));
        for (arm, sample) in arm_defs.iter().zip(per_arm.iter_mut()) {
            let (fa, fb) = (arm.fold)(&a, &b, n_bits);
            sample.approx.push(fa.tanimoto(&fb));
            sample
                .fill
                .push((fa.popcount() as f64 + fb.popcount() as f64) / (2.0 * n_bits as f64));
        }
    }
    BatchResult {
        exact,
        ceiling,
        per_arm,
    }
}

fn print_header() {
    println!(
        "{:<20} {:>10} {:>10} {:>10} {:>10}",
        "arm", "n", "MAE", "pearson_r", "mean_fill"
    );
}

fn print_row(name: &str, exact: &[f64], approx: &[f64], fill: Option<&[f64]>) {
    let fill_str = match fill {
        Some(f) => format!("{:>9.1}%", mean(f) * 100.0),
        None => format!("{:>10}", "n/a"),
    };
    println!(
        "{:<20} {:>10} {:>10.4} {:>10.4} {fill_str}",
        name,
        exact.len(),
        mae(exact, approx),
        pearson(exact, approx),
    );
}

fn print_batch(label: &str, result: &BatchResult, arm_defs: &[Arm]) {
    println!("{label}");
    print_header();
    print_row("ceiling(clamp@R)", &result.exact, &result.ceiling, None);
    for (arm, sample) in arm_defs.iter().zip(result.per_arm.iter()) {
        print_row(arm.name, &result.exact, &sample.approx, Some(&sample.fill));
    }
    println!();
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
         Pearson is noise-dominated (ground truth barely varies within one fixed overlap level);\n\
         see the pooled random-overlap sections for a meaningful correlation figure. The\n\
         ceiling(clamp@R) row isolates clamp bias (count clamped to max_repeats, zero hashing)\n\
         from hashing/collision noise in the arms below it.\n"
    );

    for &n_bits in &n_bits_values {
        println!("======================== n_bits = {n_bits} ========================");
        let arm_defs = arms(seed_base, max_repeats);

        for &overlap in &overlap_levels {
            let result = run_batch(
                seed_base ^ overlap.to_bits() ^ (n_bits as u64),
                n_features,
                max_count,
                max_repeats,
                n_bits,
                n_pairs_per_level,
                &arm_defs,
                |_| overlap,
            );
            print_batch(
                &format!("--- overlap = {overlap:.2} ---"),
                &result,
                &arm_defs,
            );
        }

        let pooled = run_batch(
            seed_base ^ 0xC0FFEE ^ (n_bits as u64),
            n_features,
            max_count,
            max_repeats,
            n_bits,
            n_pooled_pairs,
            &arm_defs,
            |rng| 0.05 + rng.next_f64() * 0.9,
        );
        print_batch(
            &format!(
                "--- pooled, overlap_frac ~ Uniform[0.05, 0.95) per pair, n={n_pooled_pairs} ---"
            ),
            &pooled,
            &arm_defs,
        );
    }

    // ---- realistic-fill, clamp-free regime ----
    // n_features=60 puts single-fingerprint popcount in the same ballpark as
    // real ECFP4 (tens to ~100 active features), and max_count=max_repeats
    // means no arm clamps, so the ceiling row above should read ~0 MAE here
    // -- isolating pure hashing/collision effects from clamp bias.
    println!("======================== realistic-fill, clamp-free regime ========================");
    println!("n_features=60, max_count=max_repeats={max_repeats} (no clamping), pooled overlap\n");
    let small_n_features = 60;
    let matched_max_count = max_repeats;
    for &n_bits in &n_bits_values {
        let arm_defs = arms(seed_base, max_repeats);
        let result = run_batch(
            seed_base ^ 0xBEEF ^ (n_bits as u64),
            small_n_features,
            matched_max_count,
            max_repeats,
            n_bits,
            n_pooled_pairs,
            &arm_defs,
            |rng| 0.05 + rng.next_f64() * 0.9,
        );
        print_batch(
            &format!("--- n_bits = {n_bits}, n={n_pooled_pairs} ---"),
            &result,
            &arm_defs,
        );
    }
}
