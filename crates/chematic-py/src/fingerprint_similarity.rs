//! Shared byte-fingerprint similarity contract for the Python bindings.
//!
//! Python exposes fingerprints as arbitrary byte strings, unlike the fixed-size
//! `BitVec2048` used by the Rust API. Keep validation and the empty-set
//! convention here so scalar, matrix, slice, and nearest-neighbour entry points
//! cannot silently drift apart.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FingerprintLengthMismatch {
    pub left: usize,
    pub right: usize,
}

impl std::fmt::Display for FingerprintLengthMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "fingerprints must be the same length ({} vs {})",
            self.left, self.right
        )
    }
}

pub(crate) fn popcount(fp: &[u8]) -> u64 {
    fp.iter().map(|byte| u64::from(byte.count_ones())).sum()
}

pub(crate) fn tanimoto_bytes(left: &[u8], right: &[u8]) -> Result<f64, FingerprintLengthMismatch> {
    tanimoto_bytes_with_counts(left, right, popcount(left), popcount(right))
}

pub(crate) fn tanimoto_bytes_with_counts(
    left: &[u8],
    right: &[u8],
    left_popcount: u64,
    right_popcount: u64,
) -> Result<f64, FingerprintLengthMismatch> {
    if left.len() != right.len() {
        return Err(FingerprintLengthMismatch {
            left: left.len(),
            right: right.len(),
        });
    }

    let intersection: u64 = left
        .iter()
        .zip(right)
        .map(|(a, b)| u64::from((a & b).count_ones()))
        .sum();
    let union = left_popcount + right_popcount - intersection;
    Ok(if union == 0 {
        // Match BitVec2048 and RDKit: two empty sets are identical.
        1.0
    } else {
        intersection as f64 / union as f64
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_fingerprints_are_identical() {
        assert_eq!(tanimoto_bytes(&[], &[]), Ok(1.0));
        assert_eq!(tanimoto_bytes(&[0, 0], &[0, 0]), Ok(1.0));
    }

    #[test]
    fn mismatched_lengths_are_rejected() {
        assert_eq!(
            tanimoto_bytes(&[1], &[1, 0]),
            Err(FingerprintLengthMismatch { left: 1, right: 2 })
        );
    }

    #[test]
    fn ordinary_similarity_is_unchanged() {
        assert_eq!(tanimoto_bytes(&[0b0011], &[0b0101]), Ok(1.0 / 3.0));
    }
}
