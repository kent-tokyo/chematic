//! Superimposed count-fingerprint coding — a bounded reproducibility SPIKE.
//!
//! **Status: research spike, not a production fingerprint.** This module is
//! deliberately standalone: it is not called from [`crate::ecfp`],
//! [`crate::rdkit_morgan_ecfp4`], or any other production fingerprint path,
//! and it does not change either engine's output. See `crates/chematic-fp`'s
//! entry in the repo root `CLAUDE.md` ("Fingerprints: legacy vs.
//! RDKit-bit-exact") for why those two engines are never blended.
//!
//! **Provenance note.** The general idea explored here — approximating
//! count-aware (multiset) Tanimoto similarity on a fixed-width *binary*
//! fingerprint by superimposing several independently-hashed codewords per
//! feature count, rather than folding presence/absence alone — is discussed
//! in a paper distributed under a CC BY-NC (non-commercial) license, whose
//! authors disclosed a conflict of interest (a commercial fingerprint
//! tooling business). No code, pseudocode, or reported numbers from that
//! paper were used here. Everything below is an independent design over the
//! general concept (superimposed coding itself is a classical information-
//! retrieval technique predating that paper by decades — Mooers 1949,
//! popularized in bit-vector form by Bloom filters), built and verified
//! against synthetic data generated in this repo, not against any
//! externally reported figures.
//!
//! Three encodings of a *count fingerprint* — `&[(feature_hash, count)]`, the
//! natural output shape of e.g. [`crate::ecfp::morgan_fp_counts`] — down to a
//! fixed-width binary vector are provided for comparison:
//!
//! - [`fold_presence_counts`]: plain presence/absence folding. Counts are
//!   discarded entirely; each feature sets exactly one bit.
//! - [`count_simulation_fold`]: RDKit-style "count simulation". Each feature
//!   sets up to `max_repeats` bits, one per unit of count, at deterministic
//!   pseudo-random positions (each unit of count = a different salt).
//! - [`superimposed_code_counts`]: this spike's technique. Count is encoded
//!   as a *unary/thermometer* ladder of up to `repetitions` layers (layer
//!   `l` is active iff `count > l`), and each active layer contributes its
//!   own `codeword_weight`-bit sparse codeword (not a single bit) into the
//!   shared vector. Superimposing multi-bit codewords per layer, instead of
//!   one bit per unit of count, is the one deliberate structural difference
//!   from `count_simulation_fold` that this spike measures the effect of.
//!
//! [`exact_count_tanimoto`] computes the exact multiset (count-aware)
//! Tanimoto similarity directly on the unfolded `&[(feature_hash, count)]`
//! pairs — the ground truth the three folding strategies are compared
//! against in `examples/superimposed_coding_spike_benchmark.rs`.
//!
//! **`codeword_weight` is provably Tanimoto-neutral in the collision-free
//! limit — it cannot help, only risk collisions.** In the limit `n_bits →
//! ∞` (no two distinct `(feature_hash, layer, slot)` triples ever hash to
//! the same bit), a feature `f` with `L_f` active layers contributes
//! exactly `codeword_weight` distinct bits per layer, so
//! `|A| = codeword_weight · Σ_f L_f^A`. Two fingerprints only ever share a
//! bit at address `coded_bit(f, l, slot, seed, n_bits)` when *both* have
//! layer `l` active for the same `f` (the hash is a pure function of its
//! inputs — nothing else can land there), and when that happens all
//! `codeword_weight` slots of that layer collide identically on both
//! sides. So `|A ∩ B| = codeword_weight · Σ_f min(L_f^A, L_f^B)` and,
//! via `max = a + b - min`, `|A ∪ B| = codeword_weight · Σ_f max(L_f^A,
//! L_f^B)`. The `codeword_weight` factor cancels in the ratio:
//! `Tanimoto = Σ_f min(L_f^A, L_f^B) / Σ_f max(L_f^A, L_f^B)`, independent
//! of `codeword_weight`. On any *finite* `n_bits`, `codeword_weight > 1`
//! only adds more chances for unrelated `(feature, layer, slot)` triples to
//! collide on the same bit, which can only inflate `|A ∩ B|` and `|A ∪ B|`
//! spuriously relative to this ideal ratio — never improve on it. The
//! benchmark's `codeword_weight` sweep confirms this empirically: MAE rises
//! and Pearson correlation falls monotonically as `codeword_weight`
//! increases, tracking rising mean bit-vector fill (saturation), at every
//! `n_bits` tested.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::bitvec::BitVecN;
use crate::ecfp::fnv1a as fnv1a_hash;

/// Configuration for [`superimposed_code_counts`].
#[derive(Clone, Debug)]
pub struct SuperimposedCodingConfig {
    /// Width of the output bit vector.
    pub n_bits: usize,
    /// Seed mixed into every hash call; same seed + same input ⇒ identical
    /// output (see the `deterministic_*` tests below).
    pub seed: u64,
    /// Max number of unary "count layers" encoded per feature (layer `l` is
    /// active iff `count > l`). Counts above `repetitions` are clamped —
    /// this bounds the number of bits any single feature can contribute.
    pub repetitions: u32,
    /// Number of independently-hashed bit positions set per active layer.
    /// `1` degenerates to one-bit-per-layer — provably identical to
    /// [`count_simulation_fold`] when `repetitions` matches that function's
    /// `max_repeats` and both use the same `seed` (see the
    /// `superimposed_with_codeword_weight_one_equals_count_simulation_fold`
    /// test). Values `>1` are the actual "superimposed coding" part this
    /// spike measures, and — per the module-doc proof above — cannot ever
    /// beat `1` on Tanimoto fidelity, only match it (collision-free) or
    /// underperform it (any real, finite `n_bits`). Swept in the benchmark;
    /// kept configurable so the sweep is possible, not because `>1` is ever
    /// recommended.
    pub codeword_weight: u32,
}

impl Default for SuperimposedCodingConfig {
    fn default() -> Self {
        Self {
            n_bits: 2048,
            seed: 0,
            repetitions: 4,
            // ponytail: 1, not >1 -- the module-doc proof plus this spike's
            // own benchmark (crates/chematic-fp/examples/
            // superimposed_coding_spike_benchmark.rs) both show
            // codeword_weight>1 can only match or underperform 1, never
            // improve on it. Defaulting to the proven-best setting so
            // nobody cargo-cults `::default()` into the worst arm.
            codeword_weight: 1,
        }
    }
}

/// Hash `(feature_hash, layer, slot, seed)` to a bit position in `0..n_bits`.
fn coded_bit(feature_hash: u64, layer: u32, slot: u32, seed: u64, n_bits: usize) -> usize {
    let mut buf = [0u8; 24];
    buf[0..8].copy_from_slice(&feature_hash.to_le_bytes());
    buf[8..12].copy_from_slice(&layer.to_le_bytes());
    buf[12..16].copy_from_slice(&slot.to_le_bytes());
    buf[16..24].copy_from_slice(&seed.to_le_bytes());
    (fnv1a_hash(&buf) % n_bits as u64) as usize
}

/// Superimposed coding: fold a count fingerprint down to `config.n_bits`
/// bits by OR-ing, per feature, one `codeword_weight`-bit sparse codeword
/// per active unary count-layer (see module docs and
/// [`SuperimposedCodingConfig`]).
///
/// Deterministic: identical `counts` + `config` always produce an identical
/// [`BitVecN`].
pub fn superimposed_code_counts(
    counts: &[(u64, u32)],
    config: &SuperimposedCodingConfig,
) -> BitVecN {
    let mut bv = BitVecN::new(config.n_bits);
    for &(feature_hash, count) in counts {
        let layers = count.min(config.repetitions);
        for layer in 0..layers {
            for slot in 0..config.codeword_weight {
                let bit = coded_bit(feature_hash, layer, slot, config.seed, config.n_bits);
                bv.set(bit);
            }
        }
    }
    bv
}

/// Plain folding: presence/absence only, counts discarded. Each feature with
/// `count > 0` sets exactly one bit at `feature_hash % n_bits`. This is the
/// cheapest baseline and the one most likely to under-approximate
/// count-aware similarity for high-count features.
pub fn fold_presence_counts(counts: &[(u64, u32)], n_bits: usize) -> BitVecN {
    let mut bv = BitVecN::new(n_bits);
    for &(feature_hash, count) in counts {
        if count > 0 {
            bv.set((feature_hash % n_bits as u64) as usize);
        }
    }
    bv
}

/// RDKit-style "count simulation": each feature with count `c` sets up to
/// `min(c, max_repeats)` bits, one per unit of count, each at a
/// deterministic pseudo-random position derived from `(feature_hash, unit,
/// seed)`. `max_repeats` bounds the number of bits any single feature can
/// contribute (the literal "repeat forever" reading is unbounded and
/// impractical, so this spike caps it like any real implementation would).
pub fn count_simulation_fold(
    counts: &[(u64, u32)],
    n_bits: usize,
    seed: u64,
    max_repeats: u32,
) -> BitVecN {
    let mut bv = BitVecN::new(n_bits);
    for &(feature_hash, count) in counts {
        let reps = count.min(max_repeats);
        for unit in 0..reps {
            // slot=0: a single bit per unit of count, unlike superimposed
            // coding's multi-bit-per-layer codeword.
            let bit = coded_bit(feature_hash, unit, 0, seed, n_bits);
            bv.set(bit);
        }
    }
    bv
}

/// Exact count-aware (multiset) Tanimoto similarity on raw, unfolded
/// `&[(feature_hash, count)]` pairs:
///
/// `sum(min(a_i, b_i)) / sum(max(a_i, b_i))` over the union of features,
/// with the standard empty-union convention of `1.0`.
///
/// This is the ground truth the three folding strategies above are compared
/// against — it never folds or hashes anything.
///
/// Precondition: `feature_hash` values within `a` (and within `b`) must be
/// unique. Duplicates are silently last-wins (each slice is collected into a
/// map keyed by `feature_hash`), unlike the folding functions above, which
/// treat every `(feature_hash, count)` entry as an independent contribution.
pub fn exact_count_tanimoto(a: &[(u64, u32)], b: &[(u64, u32)]) -> f64 {
    let a_map: FxHashMap<u64, u32> = a.iter().copied().collect();
    let b_map: FxHashMap<u64, u32> = b.iter().copied().collect();
    let mut keys: FxHashSet<u64> = FxHashSet::default();
    keys.extend(a_map.keys());
    keys.extend(b_map.keys());

    let mut min_sum = 0u64;
    let mut max_sum = 0u64;
    for k in keys {
        let av = *a_map.get(&k).unwrap_or(&0) as u64;
        let bv = *b_map.get(&k).unwrap_or(&0) as u64;
        min_sum += av.min(bv);
        max_sum += av.max(bv);
    }
    if max_sum == 0 {
        1.0
    } else {
        min_sum as f64 / max_sum as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── exact_count_tanimoto ────────────────────────────────────────────

    #[test]
    fn exact_tanimoto_identical_inputs_is_one() {
        let a = [(1u64, 4u32), (2, 7), (3, 1)];
        assert_eq!(exact_count_tanimoto(&a, &a), 1.0);
    }

    #[test]
    fn exact_tanimoto_both_empty_is_one() {
        assert_eq!(exact_count_tanimoto(&[], &[]), 1.0);
    }

    #[test]
    fn exact_tanimoto_hand_computed_example() {
        // a: {1:4, 2:2}; b: {1:2, 2:2, 3:5}
        // min sum = min(4,2) + min(2,2) + min(0,5) = 2 + 2 + 0 = 4
        // max sum = max(4,2) + max(2,2) + max(0,5) = 4 + 2 + 5 = 11
        // tanimoto = 4/11
        let a = [(1u64, 4u32), (2, 2)];
        let b = [(1u64, 2u32), (2, 2), (3, 5)];
        let t = exact_count_tanimoto(&a, &b);
        assert!((t - 4.0 / 11.0).abs() < 1e-12, "got {t}");
        // symmetric
        assert_eq!(t, exact_count_tanimoto(&b, &a));
    }

    #[test]
    fn exact_tanimoto_disjoint_is_zero() {
        let a = [(1u64, 3u32)];
        let b = [(2u64, 5u32)];
        assert_eq!(exact_count_tanimoto(&a, &b), 0.0);
    }

    // ── fold_presence_counts ────────────────────────────────────────────

    #[test]
    fn plain_folding_hand_computed_positions() {
        // feature hashes chosen so `% 16` is trivial by hand.
        let counts = [(5u64, 3u32), (7u64, 1u32), (21u64, 9u32)]; // 21 % 16 == 5, collides with feature 5
        let bv = fold_presence_counts(&counts, 16);
        assert!(bv.get(5), "5 % 16 == 5");
        assert!(bv.get(7), "7 % 16 == 7");
        // popcount is 2, not 3: features 5 and 21 collide on bit 5.
        assert_eq!(bv.popcount(), 2);
    }

    #[test]
    fn plain_folding_ignores_count_magnitude() {
        let small = [(9u64, 1u32)];
        let large = [(9u64, 1000u32)];
        assert_eq!(
            fold_presence_counts(&small, 512),
            fold_presence_counts(&large, 512),
            "presence folding must be blind to count value"
        );
    }

    #[test]
    fn plain_folding_zero_count_sets_no_bit() {
        let counts = [(9u64, 0u32)];
        assert_eq!(fold_presence_counts(&counts, 512).popcount(), 0);
    }

    // ── count_simulation_fold ───────────────────────────────────────────

    #[test]
    fn count_simulation_bounds_popcount_by_min_count_max_repeats() {
        let counts = [(11u64, 5u32), (13u64, 2u32)];
        let bv = count_simulation_fold(&counts, 1024, 0, 3);
        // feature 11: min(5,3)=3 attempted bits; feature 13: min(2,3)=2.
        // Hash collisions can only reduce popcount, never exceed the sum.
        assert!(bv.popcount() <= 5);
    }

    #[test]
    fn count_simulation_single_bit_universe_collapses_to_one_bit() {
        // n_bits=1: every hash mod 1 == 0, so any nonzero count sets
        // exactly bit 0 and nothing else -- fully hand-verifiable.
        let counts = [(42u64, 7u32)];
        let bv = count_simulation_fold(&counts, 1, 0, 3);
        assert_eq!(bv.popcount(), 1);
        assert!(bv.get(0));
    }

    #[test]
    fn count_simulation_zero_count_sets_no_bit() {
        let counts = [(42u64, 0u32)];
        assert_eq!(count_simulation_fold(&counts, 512, 0, 3).popcount(), 0);
    }

    // ── superimposed_code_counts ────────────────────────────────────────

    #[test]
    fn superimposed_single_bit_universe_collapses_to_one_bit() {
        // n_bits=1: every coded_bit call maps to bit 0 regardless of layer
        // or slot, so any nonzero count sets exactly bit 0.
        let counts = [(42u64, 7u32)];
        let cfg = SuperimposedCodingConfig {
            n_bits: 1,
            seed: 0,
            repetitions: 4,
            codeword_weight: 3,
        };
        let bv = superimposed_code_counts(&counts, &cfg);
        assert_eq!(bv.popcount(), 1);
        assert!(bv.get(0));
    }

    #[test]
    fn superimposed_popcount_bounded_by_layers_times_codeword_weight() {
        let counts = [(11u64, 10u32), (13u64, 1u32)];
        let cfg = SuperimposedCodingConfig {
            n_bits: 4096,
            seed: 7,
            repetitions: 4,
            codeword_weight: 3,
        };
        let bv = superimposed_code_counts(&counts, &cfg);
        // feature 11: min(10,4)=4 layers * 3 slots = 12; feature 13: 1*3=3.
        assert!(bv.popcount() <= 15);
    }

    #[test]
    fn superimposed_zero_count_sets_no_bit() {
        let counts = [(42u64, 0u32)];
        let cfg = SuperimposedCodingConfig::default();
        assert_eq!(superimposed_code_counts(&counts, &cfg).popcount(), 0);
    }

    #[test]
    fn superimposed_more_repeated_count_never_shrinks_popcount() {
        // Monotonicity sanity check: going from count=1 to count=5 (more
        // active layers) can only add bits, never remove them, since layer
        // l < old count is a strict subset of layer l < new count.
        let cfg = SuperimposedCodingConfig {
            n_bits: 4096,
            seed: 3,
            repetitions: 8,
            codeword_weight: 2,
        };
        let low = superimposed_code_counts(&[(99u64, 1u32)], &cfg);
        let high = superimposed_code_counts(&[(99u64, 5u32)], &cfg);
        assert!(high.popcount() >= low.popcount());
    }

    // ── structural relationship to count_simulation_fold ───────────────

    #[test]
    fn superimposed_with_codeword_weight_one_equals_count_simulation_fold() {
        // codeword_weight=1 means each active layer contributes exactly one
        // bit via coded_bit(feature_hash, layer, slot=0, seed, n_bits) --
        // the same formula, same salts, as count_simulation_fold's per-unit
        // bit. So at codeword_weight=1 and repetitions == max_repeats, the
        // two functions must produce byte-identical output. This is the
        // load-bearing fact behind this spike's real finding: the
        // `codeword_weight` knob is the ONLY structural difference between
        // "superimposed coding" as implemented here and the RDKit-style
        // baseline -- see the benchmark's codeword_weight sweep.
        let counts: Vec<(u64, u32)> = (0u64..40)
            .map(|i| (i * 97 + 3, (i % 9) as u32 + 1))
            .collect();
        let cfg = SuperimposedCodingConfig {
            n_bits: 1024,
            seed: 7,
            repetitions: 4,
            codeword_weight: 1,
        };
        let sup = superimposed_code_counts(&counts, &cfg);
        let sim = count_simulation_fold(&counts, 1024, 7, 4);
        assert_eq!(
            sup, sim,
            "codeword_weight=1 must degenerate to exactly count_simulation_fold"
        );
    }

    // ── determinism ─────────────────────────────────────────────────────

    #[test]
    fn deterministic_same_seed_identical_output_all_three_strategies() {
        let counts: Vec<(u64, u32)> = (0u64..50)
            .map(|i| (i * 97 + 3, (i % 9) as u32 + 1))
            .collect();
        let cfg = SuperimposedCodingConfig {
            n_bits: 1024,
            seed: 12345,
            repetitions: 5,
            codeword_weight: 3,
        };

        let sup1 = superimposed_code_counts(&counts, &cfg);
        let sup2 = superimposed_code_counts(&counts, &cfg);
        assert_eq!(sup1, sup2, "superimposed_code_counts must be deterministic");

        let sim1 = count_simulation_fold(&counts, 1024, 12345, 5);
        let sim2 = count_simulation_fold(&counts, 1024, 12345, 5);
        assert_eq!(sim1, sim2, "count_simulation_fold must be deterministic");

        let plain1 = fold_presence_counts(&counts, 1024);
        let plain2 = fold_presence_counts(&counts, 1024);
        assert_eq!(plain1, plain2, "fold_presence_counts must be deterministic");
    }

    #[test]
    fn deterministic_different_seed_usually_differs() {
        let counts: Vec<(u64, u32)> = (0u64..50)
            .map(|i| (i * 97 + 3, (i % 9) as u32 + 1))
            .collect();
        let cfg_a = SuperimposedCodingConfig {
            n_bits: 1024,
            seed: 1,
            repetitions: 5,
            codeword_weight: 3,
        };
        let cfg_b = SuperimposedCodingConfig {
            seed: 2,
            ..cfg_a.clone()
        };
        let a = superimposed_code_counts(&counts, &cfg_a);
        let b = superimposed_code_counts(&counts, &cfg_b);
        assert_ne!(a, b, "different seeds should (almost always) diverge");
    }
}
