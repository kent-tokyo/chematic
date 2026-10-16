//! Construction-time external partial charges for prepared MMFF fields.
use super::{MinimizerError, MmffVariant, Molecule, RdkitMmffField};

#[derive(Debug)]
#[non_exhaustive]
pub enum MmffPartialChargeError {
    Count { expected: usize, actual: usize },
    Nonfinite { index: usize },
    Coordinates,
    Threshold,
    CoefficientRange,
    Preparation(MinimizerError),
}

impl core::fmt::Display for MmffPartialChargeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Count { expected, actual } => write!(
                f,
                "partial charges must have one entry per atom: expected {expected}, got {actual}"
            ),
            Self::Nonfinite { index } => write!(f, "non-finite partial charge at atom {index}"),
            Self::Coordinates => f.write_str("coordinates must have one finite triple per atom"),
            Self::Threshold => f.write_str("nonbonded threshold must be finite and nonnegative"),
            Self::CoefficientRange => {
                f.write_str("partial charges exceed electrostatic coefficient range")
            }
            Self::Preparation(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for MmffPartialChargeError {}

impl RdkitMmffField {
    /// Prepare MMFF94/MMFF94s using caller-supplied partial charges in units
    /// of elementary charge, in input atom order (including explicit H).
    ///
    /// Overrides only electrostatics, retaining MMFF typing, bonded/vdW terms,
    /// topology exclusions, 1-4 scaling and construction-time pair selection.
    /// No charge renormalization, net-charge check or default substitution is
    /// performed. Exact zero disables that atom's electrostatics; finite small
    /// nonzero charges are retained rather than using the default near-zero
    /// cutoff. Nonrepresentable included-pair coefficients are refused.
    ///
    /// Rebuild the field to change charges. Subsequent coordinate, dielectric
    /// and minimizer operations use these prepared coefficients. This hybrid
    /// charge model is opt-in and is not a standard MMFF parametrization.
    pub fn with_terms_and_partial_charges(
        mol: &Molecule,
        coords: &[[f64; 3]],
        variant: MmffVariant,
        non_bonded_thresh: f64,
        ignore_interfrag: bool,
        terms: [bool; 7],
        charges: &[f64],
    ) -> Result<Self, MmffPartialChargeError> {
        if charges.len() != mol.atom_count() {
            return Err(MmffPartialChargeError::Count {
                expected: mol.atom_count(),
                actual: charges.len(),
            });
        }
        if let Some(index) = charges.iter().position(|q| !q.is_finite()) {
            return Err(MmffPartialChargeError::Nonfinite { index });
        }
        if coords.len() != mol.atom_count() || coords.iter().flatten().any(|v| !v.is_finite()) {
            return Err(MmffPartialChargeError::Coordinates);
        }
        if !non_bonded_thresh.is_finite() || non_bonded_thresh < 0.0 {
            return Err(MmffPartialChargeError::Threshold);
        }
        let field = Self::prepare(
            mol,
            coords,
            variant,
            non_bonded_thresh,
            ignore_interfrag,
            terms,
            Some(charges),
        )
        .map_err(MmffPartialChargeError::Preparation)?;
        if field.nonbonded.iter().any(|t| {
            t.ele
                .is_some_and(|(base, _, _)| base == 0.0 || !(332.0716 * base).is_finite())
        }) {
            return Err(MmffPartialChargeError::CoefficientRange);
        }
        Ok(field)
    }
}
