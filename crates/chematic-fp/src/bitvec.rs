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
        assert!(
            bit < 2048,
            "bit index {bit} out of range for BitVec2048 (max 2047)"
        );
        self.words[bit / 64] |= 1u64 << (bit % 64);
    }

    /// Return the value of bit `bit`.
    ///
    /// # Panics
    /// Panics if `bit >= 2048`.
    pub fn get(&self, bit: usize) -> bool {
        assert!(
            bit < 2048,
            "bit index {bit} out of range for BitVec2048 (max 2047)"
        );
        (self.words[bit / 64] >> (bit % 64)) & 1 == 1
    }

    // SIMD note: #![forbid(unsafe_code)] + wasm32 target means hand-rolled
    // intrinsics (std::arch) are not an option. LLVM autovectorizes these loops
    // to AVX2/NEON with --release; #[inline] lets it see the fixed trip count
    // (32 words) and apply unrolling. u64::count_ones() → hardware POPCNT.

    /// Count the number of bits set to 1 (popcount / Hamming weight).
    #[inline]
    pub fn popcount(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Bitwise AND of two bitvectors.
    #[inline]
    pub fn and(&self, other: &Self) -> Self {
        let mut result = Self::new();
        for (out, (a, b)) in result
            .words
            .iter_mut()
            .zip(self.words.iter().zip(other.words.iter()))
        {
            *out = a & b;
        }
        result
    }

    /// Bitwise OR of two bitvectors.
    #[inline]
    pub fn or(&self, other: &Self) -> Self {
        let mut result = Self::new();
        for (out, (a, b)) in result
            .words
            .iter_mut()
            .zip(self.words.iter().zip(other.words.iter()))
        {
            *out = a | b;
        }
        result
    }

    /// Count the number of bits set in `A & B` without allocating a temporary.
    ///
    /// This is the fundamental primitive for bulk Tanimoto computation:
    /// hoisting the query popcount and reusing this per db entry avoids
    /// the per-pair allocation of `self.and(other).popcount()`.
    #[inline]
    pub fn intersection_popcount(&self, other: &Self) -> u32 {
        self.words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| (a & b).count_ones())
            .sum()
    }

    /// Tanimoto similarity given precomputed popcounts, returning `f32`.
    ///
    /// Use when calling one query against many db entries: hoist
    /// `self_popcount = query.popcount()` outside the loop, compute
    /// `other_popcount = db_entry.popcount()` once per entry, then call
    /// this instead of `tanimoto()` to avoid redundant allocation and
    /// recomputation.
    ///
    /// Returns `1.0` when union is zero (both all-zero convention).
    #[inline]
    pub fn tanimoto_with_counts(
        &self,
        other: &Self,
        self_popcount: u32,
        other_popcount: u32,
    ) -> f32 {
        let inter = self.intersection_popcount(other) as f32;
        let union = self_popcount as f32 + other_popcount as f32 - inter;
        if union == 0.0 { 1.0 } else { inter / union }
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

    /// Convert to a variable-length BitVecN with the same bits.
    pub fn to_bitvecn(&self) -> BitVecN {
        BitVecN {
            words: self.words.to_vec(),
            bits: 2048,
        }
    }
}

/// A variable-length bitvector with dynamic bit width.
///
/// Unlike `BitVec2048` which is fixed at 2048 bits, `BitVecN` allows
/// arbitrary bit widths (512, 1024, 2048, 4096, etc.).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitVecN {
    words: Vec<u64>,
    bits: usize,
}

impl BitVecN {
    /// Create a new all-zero bitvector with specified bit width.
    ///
    /// # Panics
    /// Panics if `bits == 0`.
    pub fn new(bits: usize) -> Self {
        assert!(bits > 0, "BitVecN must have at least 1 bit");
        let num_words = bits.div_ceil(64);
        Self {
            words: vec![0u64; num_words],
            bits,
        }
    }

    /// Return the bit width.
    pub fn bit_width(&self) -> usize {
        self.bits
    }

    /// Set bit `bit` to 1.
    ///
    /// # Panics
    /// Panics if `bit >= self.bits`.
    pub fn set(&mut self, bit: usize) {
        assert!(
            bit < self.bits,
            "bit index {bit} out of range for BitVecN (max {})",
            self.bits - 1
        );
        self.words[bit / 64] |= 1u64 << (bit % 64);
    }

    /// Return the value of bit `bit`.
    ///
    /// # Panics
    /// Panics if `bit >= self.bits`.
    pub fn get(&self, bit: usize) -> bool {
        assert!(
            bit < self.bits,
            "bit index {bit} out of range for BitVecN (max {})",
            self.bits - 1
        );
        (self.words[bit / 64] >> (bit % 64)) & 1 == 1
    }

    /// Count the number of bits set to 1.
    pub fn popcount(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Tanimoto similarity with another BitVecN.
    ///
    /// Both vectors must have the same bit width.
    ///
    /// # Panics
    /// Panics if the two vectors have different bit widths.
    pub fn tanimoto(&self, other: &Self) -> f64 {
        assert_eq!(
            self.bits, other.bits,
            "BitVecN tanimoto requires same bit width"
        );
        let mut intersection = 0u32;
        for (a, b) in self.words.iter().zip(other.words.iter()) {
            intersection += (a & b).count_ones();
        }
        let intersection = intersection as f64;
        let a = self.popcount() as f64;
        let b = other.popcount() as f64;
        let union = a + b - intersection;
        if union == 0.0 {
            1.0
        } else {
            intersection / union
        }
    }

    /// Convert from a BitVec2048 (truncates to same 2048 bits).
    pub fn from_bitvec2048(bv: &BitVec2048) -> Self {
        BitVecN {
            words: bv.words.to_vec(),
            bits: 2048,
        }
    }

    /// Convert to a BitVec2048 if the bit width is exactly 2048.
    ///
    /// Returns None if the bit width is not 2048.
    pub fn to_bitvec2048(&self) -> Option<BitVec2048> {
        if self.bits != 2048 {
            return None;
        }
        let mut arr = [0u64; 32];
        for (i, &w) in self.words.iter().enumerate() {
            arr[i] = w;
        }
        Some(BitVec2048 { words: arr })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== BitVec2048 Tests =====

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

    // ===== BitVecN Tests =====

    #[test]
    fn bitvecn_new_creates_zero_vector() {
        let bv = BitVecN::new(512);
        assert_eq!(bv.bit_width(), 512);
        assert_eq!(bv.popcount(), 0);
    }

    #[test]
    fn bitvecn_set_get_basic() {
        let mut bv = BitVecN::new(1024);
        bv.set(0);
        bv.set(512);
        bv.set(1023);
        assert!(bv.get(0));
        assert!(bv.get(512));
        assert!(bv.get(1023));
        assert!(!bv.get(1));
        assert_eq!(bv.popcount(), 3);
    }

    #[test]
    fn bitvecn_tanimoto_identical() {
        let mut bv = BitVecN::new(256);
        bv.set(10);
        bv.set(50);
        assert_eq!(bv.tanimoto(&bv.clone()), 1.0);
    }

    #[test]
    fn bitvecn_tanimoto_disjoint() {
        let mut a = BitVecN::new(256);
        a.set(5);
        let mut b = BitVecN::new(256);
        b.set(10);
        assert_eq!(a.tanimoto(&b), 0.0);
    }

    #[test]
    fn bitvecn_conversion_to_from_2048() {
        let mut bv2048 = BitVec2048::new();
        bv2048.set(42);
        bv2048.set(100);

        let bvn = BitVecN::from_bitvec2048(&bv2048);
        assert_eq!(bvn.bit_width(), 2048);
        assert!(bvn.get(42));
        assert!(bvn.get(100));
        assert_eq!(bvn.popcount(), 2);

        let bv2048_back = bvn.to_bitvec2048();
        assert_eq!(bv2048_back, Some(bv2048));
    }

    #[test]
    fn bitvecn_arbitrary_sizes() {
        for size in [128, 256, 512, 1024, 2048, 4096] {
            let mut bv = BitVecN::new(size);
            bv.set(0);
            bv.set(size - 1);
            assert_eq!(bv.popcount(), 2);
            assert_eq!(bv.bit_width(), size);
        }
    }
}
