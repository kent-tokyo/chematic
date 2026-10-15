//! Checked Cartesian-to-Z-matrix conversion and Gaussian-style output.
//!
//! Geometry only: references are caller-selected, zero-based, distinct earlier
//! atoms. They need not be chemically bonded. No topology or electronic-state
//! perception, dummy atoms, reordering, or arbitrary torsions are introduced.
use chematic_core::Element;

mod geometry;
mod writer;
pub use geometry::zmatrix_to_cartesian;
pub use writer::{GzmatWriteOptions, write_gzmat};

/// Hard cap keeps conversion and output linear and bounded.
pub const MAX_ZMATRIX_ATOMS: usize = 100_000;
/// Minimum sine of a defined angle or reference-plane angle.
pub const ZMATRIX_MIN_SINE: f64 = 1e-10;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ZMatrixError {
    /// Zero-based atom row; None denotes a document/options error.
    pub atom: Option<usize>,
    pub detail: String,
}
impl core::fmt::Display for ZMatrixError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Z-matrix {:?}: {}", self.atom, self.detail)
    }
}
impl std::error::Error for ZMatrixError {}
fn error(atom: Option<usize>, detail: impl Into<String>) -> ZMatrixError {
    ZMatrixError {
        atom,
        detail: detail.into(),
    }
}

/// One row of internal coordinates; angstrom distance, degree angle/torsion.
/// The first min(row_index, 3) earlier atoms define these values in order.
#[derive(Debug, Clone, PartialEq)]
pub struct ZMatrixRow {
    pub element: Element,
    pub references: Vec<usize>,
    pub values: Vec<f64>,
}

fn validate_references(row: usize, references: &[usize]) -> Result<(), ZMatrixError> {
    if references.len() != row.min(3) {
        return Err(error(Some(row), "expected min(row_index, 3) references"));
    }
    for (k, &reference) in references.iter().enumerate() {
        if reference >= row || references[..k].contains(&reference) {
            return Err(error(
                Some(row),
                "references must be distinct earlier atoms",
            ));
        }
    }
    Ok(())
}

/// Validate publicly mutable rows, including ranges and reference structure.
pub fn validate_zmatrix(rows: &[ZMatrixRow]) -> Result<(), ZMatrixError> {
    validate_count(rows.len())?;
    for (i, row) in rows.iter().enumerate() {
        validate_references(i, &row.references)?;
        if row.values.len() != row.references.len() {
            return Err(error(Some(i), "reference and value counts differ"));
        }
        validate_values(i, &row.values)?;
    }
    Ok(())
}
fn validate_count(count: usize) -> Result<(), ZMatrixError> {
    if count == 0 || count > MAX_ZMATRIX_ATOMS {
        return Err(error(None, "atom count must be 1..=100000"));
    }
    Ok(())
}
fn validate_values(i: usize, values: &[f64]) -> Result<(), ZMatrixError> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err(error(Some(i), "non-finite internal coordinate"));
    }
    if values.first().is_some_and(|v| *v <= 0.0) {
        return Err(error(Some(i), "distance must be positive"));
    }
    if values
        .get(1)
        .is_some_and(|v| *v <= 0.0 || *v >= 180.0 || v.to_radians().sin() <= ZMATRIX_MIN_SINE)
    {
        return Err(error(Some(i), "angle is linear or near-linear"));
    }
    if values.get(2).is_some_and(|v| !(-180.0..=180.0).contains(v)) {
        return Err(error(Some(i), "torsion must be in [-180,180] degrees"));
    }
    Ok(())
}

/// Convert indexed Cartesian angstrom coordinates without choosing a frame.
/// Undefined/near-linear planes fail instead of receiving an arbitrary torsion.
/// At most three references per atom; all work and memory are O(atom count).
pub fn cartesian_to_zmatrix(
    elements: &[Element],
    positions_angstrom: &[[f64; 3]],
    references: &[Vec<usize>],
) -> Result<Vec<ZMatrixRow>, ZMatrixError> {
    validate_count(elements.len())?;
    if elements.len() != positions_angstrom.len() || elements.len() != references.len() {
        return Err(error(
            None,
            "element, position and reference row counts must match",
        ));
    }
    for (i, p) in positions_angstrom.iter().enumerate() {
        if p.iter().any(|v| !v.is_finite()) {
            return Err(error(Some(i), "non-finite Cartesian coordinate"));
        }
        validate_references(i, &references[i])?;
    }
    let mut rows = Vec::with_capacity(elements.len());
    for (i, &element) in elements.iter().enumerate() {
        let values = geometry::internal_values(i, positions_angstrom, &references[i])?;
        validate_values(i, &values)?;
        rows.push(ZMatrixRow {
            element,
            references: references[i].clone(),
            values,
        });
    }
    Ok(rows)
}
