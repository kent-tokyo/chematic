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
    pub fn from_seed(seed: u64) -> Self {
        Self(seed | 1) // low-bit set ensures seed is always non-zero (xorshift64 absorbs at 0)
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
