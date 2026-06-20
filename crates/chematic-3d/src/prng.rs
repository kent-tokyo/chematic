/// Minimal xorshift64 PRNG — replaces the `fastrand` crate dependency.
///
/// Quality is sufficient for ETKDG torsion noise and MD velocity initialization;
/// cryptographic quality is not required here.
pub(crate) struct Prng(u64);

impl Prng {
    /// Create a new PRNG with a fixed non-zero seed.
    pub fn new() -> Self {
        Self(0x517c_c1b7_2722_0a95)
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
}
