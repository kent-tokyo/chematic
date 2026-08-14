//! Shared, crate-internal validation helpers.
//!
//! Every module that accepts raw `f64`/`[f64; 3]` values from a caller
//! (rather than deriving them from an already-validated value) routes them
//! through here, so "reject NaN/Infinity before it enters a public type" is
//! enforced in one place instead of re-implemented per module.

use crate::error::CrystalError;

/// Reject a non-finite scalar.
pub(crate) fn require_finite(value: f64, field: &'static str) -> Result<(), CrystalError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CrystalError::NonFinite { field })
    }
}

/// Reject a non-finite 3-vector (any component `NaN`/`Infinity`).
pub(crate) fn require_finite3(value: [f64; 3], field: &'static str) -> Result<(), CrystalError> {
    if value.iter().all(|c| c.is_finite()) {
        Ok(())
    } else {
        Err(CrystalError::NonFinite { field })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_finite_accepts_finite() {
        assert!(require_finite(1.5, "x").is_ok());
        assert!(require_finite(0.0, "x").is_ok());
        assert!(require_finite(-3.0, "x").is_ok());
    }

    #[test]
    fn require_finite_rejects_nan_and_infinity() {
        assert_eq!(
            require_finite(f64::NAN, "x"),
            Err(CrystalError::NonFinite { field: "x" })
        );
        assert_eq!(
            require_finite(f64::INFINITY, "x"),
            Err(CrystalError::NonFinite { field: "x" })
        );
        assert_eq!(
            require_finite(f64::NEG_INFINITY, "x"),
            Err(CrystalError::NonFinite { field: "x" })
        );
    }

    #[test]
    fn require_finite3_rejects_any_bad_component() {
        assert!(require_finite3([1.0, 2.0, 3.0], "v").is_ok());
        assert!(require_finite3([f64::NAN, 2.0, 3.0], "v").is_err());
        assert!(require_finite3([1.0, f64::INFINITY, 3.0], "v").is_err());
        assert!(require_finite3([1.0, 2.0, f64::NEG_INFINITY], "v").is_err());
    }
}
