//! Checked Euclidean distances between supplied descriptor vectors.

pub const MAX_DESCRIPTOR_ROWS: usize = 1024;
pub const MAX_DESCRIPTOR_FEATURES: usize = 256;
pub const MAX_DESCRIPTOR_PAIR_FEATURES: usize = 50_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum DescriptorDistanceError {
    RowLimit,
    FeatureLimit,
    Ragged {
        row: usize,
        expected: usize,
        actual: usize,
    },
    Nonfinite {
        row: usize,
        column: usize,
    },
    WorkLimit,
    Allocation,
    NumericalRange {
        first: usize,
        second: usize,
    },
}

impl std::fmt::Display for DescriptorDistanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RowLimit => write!(
                f,
                "descriptor matrix supports at most {MAX_DESCRIPTOR_ROWS} rows"
            ),
            Self::FeatureLimit => write!(
                f,
                "descriptor vectors support at most {MAX_DESCRIPTOR_FEATURES} features"
            ),
            Self::Ragged {
                row,
                expected,
                actual,
            } => write!(
                f,
                "descriptor row {row} has {actual} features; expected {expected}"
            ),
            Self::Nonfinite { row, column } => {
                write!(f, "non-finite descriptor at row {row}, column {column}")
            }
            Self::WorkLimit => write!(
                f,
                "descriptor distances exceed {MAX_DESCRIPTOR_PAIR_FEATURES} pair-feature operations"
            ),
            Self::Allocation => f.write_str("descriptor distance output allocation failed"),
            Self::NumericalRange { first, second } => write!(
                f,
                "distance between descriptor rows {first} and {second} exceeds finite float range"
            ),
        }
    }
}
impl std::error::Error for DescriptorDistanceError {}

fn distance(a: &[f64], b: &[f64]) -> Option<f64> {
    // Scaled sum of squares: no intermediate squaring of an absolute
    // difference. Tiny and huge representable norms remain representable.
    let (mut scale, mut sum) = (0.0_f64, 1.0_f64);
    for (&x, &y) in a.iter().zip(b) {
        let value = (x - y).abs();
        if !value.is_finite() {
            return None;
        }
        if value != 0.0 {
            if value > scale {
                let ratio = scale / value;
                sum = 1.0 + sum * ratio * ratio;
                scale = value;
            } else {
                let ratio = value / scale;
                sum += ratio * ratio;
            }
        }
    }
    let result = scale * sum.sqrt();
    result.is_finite().then_some(result)
}

/// Strict lower triangle in row-major order: `(1,0), (2,0), (2,1), ...`.
/// Entry `(i,j)` for `i > j` has index `i*(i-1)/2+j`; the diagonal is omitted.
/// Values use the caller's descriptor units, without scaling or normalization.
///
/// Requires rectangular finite vectors, at most 1024 rows / 256 features and
/// 50 million pair-feature operations. The output has at most 523776 values.
/// Empty input or a singleton returns empty output after input validation;
/// zero-feature vectors have zero distances. A scaled norm avoids intermediate
/// overflow/underflow; an unrepresentable distance refuses the whole result.
/// Allocation is fallible. Finite-precision error remains; no bit parity or
/// descriptor generation, molecule mapping, clustering or feature weighting
/// is implied.
pub fn euclidean_distance_matrix_checked<T: AsRef<[f64]>>(
    rows: &[T],
) -> Result<Vec<f64>, DescriptorDistanceError> {
    let count = rows.len();
    if count > MAX_DESCRIPTOR_ROWS {
        return Err(DescriptorDistanceError::RowLimit);
    }
    let features = rows.first().map_or(0, |r| r.as_ref().len());
    if features > MAX_DESCRIPTOR_FEATURES {
        return Err(DescriptorDistanceError::FeatureLimit);
    }
    let pairs = count * count.saturating_sub(1) / 2;
    if pairs
        .checked_mul(features)
        .is_none_or(|n| n > MAX_DESCRIPTOR_PAIR_FEATURES)
    {
        return Err(DescriptorDistanceError::WorkLimit);
    }
    for (row, values) in rows.iter().enumerate() {
        let values = values.as_ref();
        if values.len() != features {
            return Err(DescriptorDistanceError::Ragged {
                row,
                expected: features,
                actual: values.len(),
            });
        }
        if let Some(column) = values.iter().position(|v| !v.is_finite()) {
            return Err(DescriptorDistanceError::Nonfinite { row, column });
        }
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(pairs)
        .map_err(|_| DescriptorDistanceError::Allocation)?;
    for first in 1..count {
        for second in 0..first {
            output.push(
                distance(rows[first].as_ref(), rows[second].as_ref())
                    .ok_or(DescriptorDistanceError::NumericalRange { first, second })?,
            );
        }
    }
    Ok(output)
}
