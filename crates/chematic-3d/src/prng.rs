/// Minimal xorshift64 PRNG — replaces the `fastrand` crate dependency.
///
/// Quality is sufficient for ETKDG torsion noise and MD velocity initialization;
/// cryptographic quality is not required here.
///
/// Each `Prng::new()` call advances a global atomic Weyl counter so that
/// successive calls (e.g. inside a conformer-ensemble loop) receive distinct
/// seeds and produce independent random streams.
pub(crate) struct Prng(u64);

/// Global Weyl-sequence counter for seeding.  The additive step (golden-ratio
/// based) ensures the period is 2^64 and every seed value is eventually visited.
static PRNG_COUNTER: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0x517c_c1b7_2722_0a95);

impl Prng {
    /// Create a new PRNG with a unique seed derived from a shared counter.
    pub fn new() -> Self {
        let seed =
            PRNG_COUNTER.fetch_add(0x9e37_79b9_7f4a_7c15, std::sync::atomic::Ordering::Relaxed) | 1; // low-bit set ensures seed is always non-zero (xorshift64 absorbs at 0)
        Self(seed)
    }

    /// Create a PRNG from an explicit, caller-controlled seed.
    ///
    /// Unlike [`Prng::new`] (process-global atomic counter, no reproducibility contract),
    /// this is fully call-local and deterministic: the same `seed` on the same
    /// target/thread-count always produces the same output stream (this codebase's
    /// existing reproducibility convention -- not a cross-platform bit-exactness claim).
    /// Added for `distance_geometry_v2`'s `EmbedParameters.random_seed`; does not change
    /// [`Prng::new`]'s existing behavior or any of its current callers.
    ///
    /// Runs `seed` through a SplitMix64 finalizer before using it as xorshift64 state.
    /// A prior version used `seed | 1` (forcing the low bit on to avoid the
    /// xorshift64 zero-absorbing state), which silently aliased every even seed `n`
    /// with odd seed `n + 1` (both map to the same state) -- verified broken via
    /// `same_seed_reproducible`-style checks at seeds (0,1), (2,3), (10,11), all
    /// "identical=true" before this fix. SplitMix64's full avalanche (flipping any
    /// one input bit flips ~half the output bits) fixes this; the same three pairs
    /// are checked non-colliding in this module's tests below.
    pub fn from_seed(seed: u64) -> Self {
        let mixed = splitmix64(seed);
        // mixed == 0 only if splitmix64's specific bit pattern lands there (not
        // reachable by any of the pairs above, but guard explicitly rather than
        // trust it can never happen): xorshift64 absorbs at zero state forever.
        Self(if mixed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            mixed
        })
    }

    /// Return a uniform f64 in [0, 1).
    pub fn f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        // Use the top 53 bits to fill the mantissa.
        (x >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Return a standard-normal sample via Box-Muller transform.
    ///
    /// Resamples the rare u1=0 case (probability 2^-53) to avoid ln(0).
    pub fn gaussian_f64(&mut self) -> f64 {
        use std::f64::consts::PI;
        let u1 = loop {
            let v = self.f64();
            if v != 0.0 {
                break v;
            }
        };
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

/// SplitMix64 finalizer (Steele, Lea & Flood 2014 / Vigna's public-domain
/// reference): a bijective, full-avalanche bit mixer used here only to turn a
/// small/adjacent-valued caller seed into well-distributed xorshift64 state --
/// not used as a standalone generator.
fn splitmix64(seed: u64) -> u64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_seed_adjacent_seeds_do_not_collide() {
        // The exact pairs the previous `seed | 1` bug aliased (verified broken
        // externally before this fix): every even/odd pair produced an identical
        // stream. Confirm the fix by checking the internal state directly rather
        // than the downstream f64 stream, so this test fails immediately if the
        // mixing regresses, without depending on embed()-level behavior.
        for (a, b) in [
            (0u64, 1u64),
            (2, 3),
            (10, 11),
            (42, 43),
            (100, 101),
            (u64::MAX - 1, u64::MAX),
        ] {
            let sa = Prng::from_seed(a).0;
            let sb = Prng::from_seed(b).0;
            assert_ne!(
                sa, sb,
                "seeds {a} and {b} collided to the same internal state"
            );
        }
    }

    #[test]
    fn from_seed_never_zero() {
        // xorshift64 absorbs at state 0 forever -- from_seed must never produce it,
        // for any input including the one seed value (0) most likely to matter to
        // callers using a default/unset seed.
        for seed in [0u64, 1, u64::MAX, 0x9E37_79B9_7F4A_7C15] {
            assert_ne!(
                Prng::from_seed(seed).0,
                0,
                "seed {seed} produced zero state"
            );
        }
    }

    #[test]
    fn from_seed_same_seed_reproducible() {
        let mut a = Prng::from_seed(777);
        let mut b = Prng::from_seed(777);
        for _ in 0..10 {
            assert_eq!(a.f64(), b.f64());
        }
    }

    #[test]
    fn from_seed_sequential_seeds_have_no_small_bit_pattern() {
        // Broader sweep than the specific pairs above: 200 sequential seeds should
        // produce 200 distinct internal states (SplitMix64 is a bijection on u64,
        // so this must hold exactly, not just "very likely").
        use std::collections::HashSet;
        let states: HashSet<u64> = (0u64..200).map(|s| Prng::from_seed(s).0).collect();
        assert_eq!(
            states.len(),
            200,
            "expected 200 distinct states from 200 sequential seeds"
        );
    }
}
