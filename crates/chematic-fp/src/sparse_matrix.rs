//! Checked condensed Tanimoto matrices over sparse on-bit positions.

pub const MAX_SPARSE_FINGERPRINTS: usize = 1024;
pub const MAX_SPARSE_BITS_PER_FINGERPRINT: usize = 65_536;
pub const MAX_SPARSE_INPUT_IDS: usize = 1_000_000;
pub const MAX_SPARSE_PAIR_WORK: usize = 50_000_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SparseTanimotoMetric {
    Similarity,
    Distance,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SparseMatrixError {
    RowLimit,
    BitsPerFingerprint { row: usize },
    TotalIds,
    WorkLimit,
    Allocation,
}
impl std::fmt::Display for SparseMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RowLimit => write!(
                f,
                "sparse matrix supports at most {MAX_SPARSE_FINGERPRINTS} fingerprints"
            ),
            Self::BitsPerFingerprint { row } => write!(
                f,
                "fingerprint {row} exceeds {MAX_SPARSE_BITS_PER_FINGERPRINT} input bit IDs"
            ),
            Self::TotalIds => write!(
                f,
                "sparse matrix exceeds {MAX_SPARSE_INPUT_IDS} total input bit IDs"
            ),
            Self::WorkLimit => write!(
                f,
                "sparse matrix exceeds {MAX_SPARSE_PAIR_WORK} pair-bit work units"
            ),
            Self::Allocation => f.write_str("sparse matrix allocation failed"),
        }
    }
}
impl std::error::Error for SparseMatrixError {}

fn reserved<T>(count: usize) -> Result<Vec<T>, SparseMatrixError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| SparseMatrixError::Allocation)?;
    Ok(values)
}

fn normalize(rows: &[Vec<u32>]) -> Result<Vec<Vec<u32>>, SparseMatrixError> {
    if rows.len() > MAX_SPARSE_FINGERPRINTS {
        return Err(SparseMatrixError::RowLimit);
    }
    let mut total = 0;
    for (row, ids) in rows.iter().enumerate() {
        if ids.len() > MAX_SPARSE_BITS_PER_FINGERPRINT {
            return Err(SparseMatrixError::BitsPerFingerprint { row });
        }
        total += ids.len();
        if total > MAX_SPARSE_INPUT_IDS {
            return Err(SparseMatrixError::TotalIds);
        }
    }
    let mut normalized = reserved(rows.len())?;
    for row in rows {
        let mut ids = reserved(row.len())?;
        ids.extend_from_slice(row);
        ids.sort_unstable();
        ids.dedup();
        normalized.push(ids);
    }
    Ok(normalized)
}

fn similarity(a: &[u32], b: &[u32]) -> f64 {
    let (mut i, mut j, mut common) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }
    let union = a.len() + b.len() - common;
    if union == 0 {
        1.0
    } else {
        common as f64 / union as f64
    }
}

/// Binary set Tanimoto over unsigned 32-bit on-bit IDs, without dense expansion.
/// Duplicates are idempotent and input order is ignored; source arrays are not
/// mutated. IDs share a caller-defined universe with no declared bit length.
/// This is not count-weighted similarity, MinHash or a fingerprint generator.
///
/// Returns the strict lower triangle `(1,0), (2,0), (2,1), ...` as f64 values.
/// Empty input/singletons return empty output after validation. Empty/empty
/// similarity is one; empty/nonempty is zero. Distance is one minus similarity.
///
/// At most 1024 fingerprints, 65536 raw IDs per fingerprint, one million raw
/// IDs overall, and `(n-1)*sum(unique ID counts)` <= 50 million work units.
/// Raw caps apply before sorting/deduplication; sorting is bounded separately
/// by those input caps. Numeric pair work is conservative, not a timeout.
/// Buffers are fallibly allocated and errors return no partial output.
pub fn sparse_tanimoto_matrix_checked(
    fingerprints: &[Vec<u32>],
    metric: SparseTanimotoMetric,
) -> Result<Vec<f64>, SparseMatrixError> {
    let rows = normalize(fingerprints)?;
    let n = rows.len();
    let unique: usize = rows.iter().map(Vec::len).sum();
    if n.saturating_sub(1) * unique > MAX_SPARSE_PAIR_WORK {
        return Err(SparseMatrixError::WorkLimit);
    }
    let mut output = reserved(n * n.saturating_sub(1) / 2)?;
    for i in 1..n {
        for j in 0..i {
            let value = similarity(&rows[i], &rows[j]);
            output.push(match metric {
                SparseTanimotoMetric::Similarity => value,
                SparseTanimotoMetric::Distance => 1.0 - value,
            });
        }
    }
    Ok(output)
}
