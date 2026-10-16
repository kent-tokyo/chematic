//! Mutable constant-dielectric screening for an already prepared MMFF field.
use super::RdkitMmffField;

/// A refused dielectric update leaves the prepared field unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MmffDielectricError {
    /// A dielectric constant must be finite and strictly positive.
    InvalidConstant,
    /// The requested screening would overflow a prepared electrostatic coefficient.
    CoefficientOverflow,
}

impl core::fmt::Display for MmffDielectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidConstant => "MMFF dielectric constant must be finite and positive",
            Self::CoefficientOverflow => {
                "MMFF dielectric constant overflows electrostatic coefficients"
            }
        })
    }
}
impl std::error::Error for MmffDielectricError {}

impl RdkitMmffField {
    /// Constant relative dielectric used in the buffered Coulomb term.
    /// Defaults to 1.0. This is not a distance-dependent dielectric model.
    pub fn dielectric_constant(&self) -> f64 {
        self.dielectric_constant
    }

    /// Update electrostatic screening without rebuilding types, charges or pairs.
    /// Energy, gradient and subsequent minimization immediately use this value.
    /// All other term classes and current coordinates are unaffected. Validation
    /// is atomic: a refused update leaves the previous value and terms intact.
    pub fn set_dielectric_constant(&mut self, value: f64) -> Result<(), MmffDielectricError> {
        if !value.is_finite() || value <= 0.0 {
            return Err(MmffDielectricError::InvalidConstant);
        }
        if self.nonbonded.iter().any(|t| {
            t.ele
                .is_some_and(|(charge, _, _)| !(332.0716 * (charge / value)).is_finite())
        }) {
            return Err(MmffDielectricError::CoefficientOverflow);
        }
        for term in &mut self.nonbonded {
            if let Some((base, screened, _)) = &mut term.ele {
                *screened = *base / value;
            }
        }
        self.dielectric_constant = value;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::MmffDielectricError;

    #[test]
    fn dielectric_errors_have_stable_messages() {
        assert_eq!(
            MmffDielectricError::InvalidConstant.to_string(),
            "MMFF dielectric constant must be finite and positive"
        );
        assert_eq!(
            MmffDielectricError::CoefficientOverflow.to_string(),
            "MMFF dielectric constant overflows electrostatic coefficients"
        );
    }
}
