use chematic_core::{AtomIdx, BondIdx, BondOrder};
use chematic_fp::{RdkitMorganError, fold_density::FoldDensityError};

#[test]
fn fold_density_errors_have_stable_messages() {
    let cases = [
        FoldDensityError::EmptyFingerprint,
        FoldDensityError::InvalidDensity,
        FoldDensityError::InvalidMinimumLength,
        FoldDensityError::ResourceLimit {
            bits: 65,
            limit: 64,
        },
        FoldDensityError::PackedLength {
            actual: 7,
            expected: 8,
        },
        FoldDensityError::NonzeroPadding,
    ];
    let messages: Vec<_> = cases.iter().map(ToString::to_string).collect();
    assert_eq!(
        messages,
        [
            "fold density requires at least one bit",
            "fold density target must be finite and in [0, 1]",
            "fold density minimum length must be positive",
            "fold density width 65 exceeds limit 64",
            "fold density expected 8 bytes, got 7",
            "fold density unused high padding bits must be zero",
        ]
    );
}

#[test]
fn rdkit_morgan_errors_preserve_the_failed_boundary() {
    let cases = [
        RdkitMorganError::UnsupportedBondOrder {
            bond_idx: BondIdx(3),
            order: BondOrder::QueryAny,
        },
        RdkitMorganError::UnsupportedCoordinationSanitization {
            atom_idx: AtomIdx(4),
            atomic_number: 26,
            degree: 6,
        },
        RdkitMorganError::InternalInvariantViolation {
            reason: "missing environment".into(),
        },
    ];
    for error in &cases {
        assert!(!error.to_string().is_empty());
    }
}
