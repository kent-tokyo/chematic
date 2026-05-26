//! Fixed-size 2048-bit bitvector backed by 32 × u64 words.

/// A 2048-bit bitvector.
///
/// The vector is stored as 32 little-endian 64-bit words.
/// Bit `n` resides in word `n / 64`, at position `n % 64`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitVec2048 {
    words: [u64; 32],
}

impl Default for BitVec2048 {
    fn default() -> Self {
        Self::new()
    }
}

impl BitVec2048 {
    /// Create a new all-zero bitvector.
    pub fn new() -> Self {
        Self { words: [0u64; 32] }
    }

    /// Set bit `bit` to 1.
    ///
    /// # Panics
    /// Panics if `bit >= 2048`.
    pub fn set(&mut self, bit: usize) {
        assert!(bit < 2048, "bit index {bit} out of range for BitVec2048 (max 2047)");
        self.words[bit / 64] |= 1u64 << (bit % 64);
    }

    /// Return the value of bit `bit`.
    ///
    /// # Panics
    /// Panics if `bit >= 2048`.
    pub fn get(&self, bit: usize) -> bool {
        assert!(bit < 2048, "bit index {bit} out of range for BitVec2048 (max 2047)");
        (self.words[bit / 64] >> (bit % 64)) & 1 == 1
    }

    /// Count the number of bits set to 1 (popcount / Hamming weight).
    pub fn popcount(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Bitwise AND of two bitvectors.
    pub fn and(&self, other: &Self) -> Self {
        let mut result = Self::new();
        for i in 0..32 {
            result.words[i] = self.words[i] & other.words[i];
        }
        result
    }

    /// Bitwise OR of two bitvectors.
    pub fn or(&self, other: &Self) -> Self {
        let mut result = Self::new();
        for i in 0..32 {
            result.words[i] = self.words[i] | other.words[i];
        }
        result
    }

    /// Tanimoto similarity: `|A & B| / |A | B|`.
    ///
    /// Returns `1.0` when both vectors are all-zero (the empty-set convention).
    pub fn tanimoto(&self, other: &Self) -> f64 {
        let intersection = self.and(other).popcount() as f64;
        let a = self.popcount() as f64;
        let b = other.popcount() as f64;
        let union = a + b - intersection;
        if union == 0.0 {
            1.0
        } else {
            intersection / union
        }
    }

    /// Dice similarity: `2|A & B| / (|A| + |B|)`.
    ///
    /// Returns `1.0` when both vectors are all-zero.
    pub fn dice(&self, other: &Self) -> f64 {
        let intersection = self.and(other).popcount() as f64;
        let a = self.popcount() as f64;
        let b = other.popcount() as f64;
        let denom = a + b;
        if denom == 0.0 {
            1.0
        } else {
            2.0 * intersection / denom
        }
    }

    /// XOR-fold this 2048-bit vector to `bits` bits.
    ///
    /// The upper half is XOR-folded into the lower half, and the upper half is zeroed.
    /// This is applied recursively until the target width is reached.
    ///
    /// `bits` must be 1024, 512, or 256.
    ///
    /// # Panics
    /// Panics if `bits` is not one of 1024, 512, or 256.
    pub fn fold(&self, bits: usize) -> Self {
        assert!(
            matches!(bits, 256 | 512 | 1024),
            "fold target must be 1024, 512, or 256; got {bits}"
        );

        let mut current = self.clone();
        let mut current_bits = 2048usize;

        while current_bits > bits {
            let half_words = current_bits / 64 / 2; // number of 64-bit words in half
            let mut folded = Self::new();
            for i in 0..half_words {
                folded.words[i] = current.words[i] ^ current.words[i + half_words];
            }
            current = folded;
            current_bits /= 2;
        }

        current
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_bitvec_is_all_zero() {
        let bv = BitVec2048::new();
        assert_eq!(bv.popcount(), 0, "new bitvec must have popcount 0");
    }

    #[test]
    fn set_and_get_boundary_bits() {
        let mut bv = BitVec2048::new();
        bv.set(0);
        bv.set(2047);
        assert_eq!(bv.popcount(), 2);
        assert!(bv.get(0));
        assert!(!bv.get(1));
        assert!(bv.get(2047));
    }

    #[test]
    fn tanimoto_identical_nonzero() {
        let mut bv = BitVec2048::new();
        bv.set(42);
        bv.set(100);
        assert_eq!(bv.tanimoto(&bv.clone()), 1.0);
    }

    #[test]
    fn tanimoto_disjoint_is_zero() {
        let mut a = BitVec2048::new();
        a.set(10);
        let mut b = BitVec2048::new();
        b.set(20);
        assert_eq!(a.tanimoto(&b), 0.0);
    }

    #[test]
    fn dice_identical_nonzero() {
        let mut bv = BitVec2048::new();
        bv.set(7);
        assert_eq!(bv.dice(&bv.clone()), 1.0);
    }

    #[test]
    fn fold_1024_popcount_leq_original() {
        // XOR-fold can only reduce or maintain the popcount.
        let mut bv = BitVec2048::new();
        for i in (0..2048).step_by(3) {
            bv.set(i);
        }
        let folded = bv.fold(1024);
        assert!(
            folded.popcount() <= bv.popcount(),
            "folded popcount ({}) must be <= original ({})",
            folded.popcount(),
            bv.popcount()
        );
    }

    #[test]
    fn and_or_basic_correctness() {
        let mut a = BitVec2048::new();
        a.set(5);
        a.set(10);

        let mut b = BitVec2048::new();
        b.set(10);
        b.set(15);

        let and = a.and(&b);
        assert!(and.get(10), "AND: shared bit 10 should be set");
        assert!(!and.get(5), "AND: bit 5 only in A should be clear");
        assert!(!and.get(15), "AND: bit 15 only in B should be clear");

        let or = a.or(&b);
        assert!(or.get(5), "OR: bit 5 should be set");
        assert!(or.get(10), "OR: bit 10 should be set");
        assert!(or.get(15), "OR: bit 15 should be set");
        assert!(!or.get(0), "OR: bit 0 should be clear");
    }
}
