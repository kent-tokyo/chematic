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
    let mut chunks = fp.chunks_exact(8);
    let mut total: u64 = (&mut chunks)
        .map(|c| u64::from(u64::from_le_bytes(c.try_into().unwrap()).count_ones()))
        .sum();
    total += chunks
        .remainder()
        .iter()
        .map(|byte| u64::from(byte.count_ones()))
        .sum::<u64>();
    total
}

fn intersection_count(left: &[u8], right: &[u8]) -> u64 {
    let mut lc = left.chunks_exact(8);
    let mut rc = right.chunks_exact(8);
    let mut total: u64 = (&mut lc)
        .zip(&mut rc)
        .map(|(a, b)| {
            let a = u64::from_le_bytes(a.try_into().unwrap());
            let b = u64::from_le_bytes(b.try_into().unwrap());
            u64::from((a & b).count_ones())
        })
        .sum();
    total += lc
        .remainder()
        .iter()
        .zip(rc.remainder())
        .map(|(a, b)| u64::from((a & b).count_ones()))
        .sum::<u64>();
    total
}

/// A Python fingerprint argument: `bytes` objects are borrowed without
/// copying; any other object goes through the historical `Vec<u8>`
/// extraction (same accepted inputs and errors as before).
pub(crate) enum FpArg<'py> {
    Borrowed(pyo3::Bound<'py, pyo3::types::PyBytes>),
    Owned(Vec<u8>),
}

impl FpArg<'_> {
    pub(crate) fn as_slice(&self) -> &[u8] {
        use pyo3::types::PyBytesMethods;
        match self {
            FpArg::Borrowed(b) => b.as_bytes(),
            FpArg::Owned(v) => v,
        }
    }
}

/// Convert a list of fingerprint objects, borrowing `bytes` zero-copy.
/// `arg_name` reproduces PyO3's argument-extraction error prefix.
pub(crate) fn fp_list<'py>(
    items: Vec<pyo3::Bound<'py, pyo3::PyAny>>,
    arg_name: &str,
) -> pyo3::PyResult<Vec<FpArg<'py>>> {
    use pyo3::prelude::*;
    items
        .into_iter()
        .map(|item| match item.cast_into::<pyo3::types::PyBytes>() {
            Ok(bytes) => Ok(FpArg::Borrowed(bytes)),
            Err(err) => err
                .into_inner()
                .extract::<Vec<u8>>()
                .map(FpArg::Owned)
                .map_err(|e| {
                    pyo3::exceptions::PyTypeError::new_err(format!("argument '{arg_name}': {e}"))
                }),
        })
        .collect()
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

    let intersection = intersection_count(left, right);
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
