//! OR folding with corrected target-density and minimum-width conditions.
//! Motivated by RDKit #9715; independent implementation on native bitvectors.

use crate::BitVecN;

/// Maximum accepted width (128 KiB of packed input).
pub const MAX_FOLD_DENSITY_BITS: usize = 1_048_576;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FoldDensityError {
    EmptyFingerprint,
    InvalidDensity,
    InvalidMinimumLength,
    ResourceLimit { bits: usize, limit: usize },
    PackedLength { actual: usize, expected: usize },
    NonzeroPadding,
}

impl std::fmt::Display for FoldDensityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyFingerprint => write!(f, "fold density requires at least one bit"),
            Self::InvalidDensity => write!(f, "fold density target must be finite and in [0, 1]"),
            Self::InvalidMinimumLength => write!(f, "fold density minimum length must be positive"),
            Self::ResourceLimit { bits, limit } => {
                write!(f, "fold density width {bits} exceeds limit {limit}")
            }
            Self::PackedLength { actual, expected } => {
                write!(f, "fold density expected {expected} bytes, got {actual}")
            }
            Self::NonzeroPadding => write!(f, "fold density unused high padding bits must be zero"),
        }
    }
}

impl std::error::Error for FoldDensityError {}

/// A folded copy. Reaching the minimum width need not reach the target density.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FoldDensityResult {
    pub fingerprint: BitVecN,
    pub original_bits: usize,
    pub folds: usize,
    pub target_reached: bool,
}

impl FoldDensityResult {
    pub fn density(&self) -> f64 {
        f64::from(self.fingerprint.popcount()) / self.fingerprint.bit_width() as f64
    }
}

fn validate(bits: usize, density: f64, min_length: usize) -> Result<(), FoldDensityError> {
    if bits == 0 {
        return Err(FoldDensityError::EmptyFingerprint);
    }
    if bits > MAX_FOLD_DENSITY_BITS {
        return Err(FoldDensityError::ResourceLimit {
            bits,
            limit: MAX_FOLD_DENSITY_BITS,
        });
    }
    if !density.is_finite() || !(0.0..=1.0).contains(&density) {
        return Err(FoldDensityError::InvalidDensity);
    }
    if min_length == 0 {
        return Err(FoldDensityError::InvalidMinimumLength);
    }
    Ok(())
}

/// Fold by two while below the target and the next width is at least the minimum.
///
/// Each step has width `floor(old_width / 2)`, with every on bit OR-combined at
/// `old_index % new_width`. Collisions never cancel. Odd widths are supported.
/// Equality with either the density target or minimum width is inclusive.
/// The input is never mutated; a no-op also returns an independent copy.
/// An input already shorter than `min_length` is returned without folding.
/// Finite targets in [0, 1] and positive minimum lengths are required.
/// Input width is capped before copying; at most 20 folds are possible.
pub fn fold_to_target_density(
    fp: &BitVecN,
    density: f64,
    min_length: usize,
) -> Result<FoldDensityResult, FoldDensityError> {
    validate(fp.bit_width(), density, min_length)?;
    let mut current = fp.clone();
    let mut folds = 0;
    while f64::from(current.popcount()) / (current.bit_width() as f64) < density
        && current.bit_width() / 2 >= min_length
    {
        let next_width = current.bit_width() / 2;
        let mut next = BitVecN::new(next_width);
        for bit in 0..current.bit_width() {
            if current.get(bit) {
                next.set(bit % next_width);
            }
        }
        current = next;
        folds += 1;
    }
    let target_reached = f64::from(current.popcount()) / current.bit_width() as f64 >= density;
    Ok(FoldDensityResult {
        fingerprint: current,
        original_bits: fp.bit_width(),
        folds,
        target_reached,
    })
}

/// Packed-byte entry point, using LSB-first bits as in [`BitVecN::to_le_bytes`].
///
/// Requires exactly `ceil(bit_width / 8)` bytes, with zero high padding bits.
/// Validation, including the width cap, runs before decoding or copying input.
pub fn fold_packed_to_target_density(
    packed: &[u8],
    bit_width: usize,
    density: f64,
    min_length: usize,
) -> Result<FoldDensityResult, FoldDensityError> {
    validate(bit_width, density, min_length)?;
    let expected = bit_width.div_ceil(8);
    if packed.len() != expected {
        return Err(FoldDensityError::PackedLength {
            actual: packed.len(),
            expected,
        });
    }
    let remainder = bit_width % 8;
    if remainder != 0 && packed[expected - 1] >> remainder != 0 {
        return Err(FoldDensityError::NonzeroPadding);
    }
    let mut fp = BitVecN::new(bit_width);
    for bit in 0..bit_width {
        if packed[bit / 8] & (1 << (bit % 8)) != 0 {
            fp.set(bit);
        }
    }
    fold_to_target_density(&fp, density, min_length)
}
