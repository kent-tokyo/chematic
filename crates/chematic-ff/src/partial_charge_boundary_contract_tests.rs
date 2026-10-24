use crate::rdkit_mmff::{MmffPartialChargeError, MmffVariant, RdkitMmffField};

#[test]
fn external_charge_errors_are_typed_before_force_field_construction() {
    let diagnostics = [
        MmffPartialChargeError::Count {
            expected: 2,
            actual: 1,
        },
        MmffPartialChargeError::Nonfinite { index: 3 },
        MmffPartialChargeError::Coordinates,
        MmffPartialChargeError::Threshold,
        MmffPartialChargeError::CoefficientRange,
    ];
    for error in diagnostics {
        assert!(!error.to_string().is_empty());
    }

    let molecule = chematic_smiles::parse("CC").unwrap();
    let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
    let terms = [true; 7];
    assert!(matches!(
        RdkitMmffField::with_terms_and_partial_charges(
            &molecule,
            &coords,
            MmffVariant::Mmff94,
            100.0,
            true,
            terms,
            &[0.0],
        ),
        Err(MmffPartialChargeError::Count { .. })
    ));
    assert!(matches!(
        RdkitMmffField::with_terms_and_partial_charges(
            &molecule,
            &coords,
            MmffVariant::Mmff94,
            100.0,
            true,
            terms,
            &[f64::NAN, 0.0],
        ),
        Err(MmffPartialChargeError::Nonfinite { index: 0 })
    ));
    assert!(matches!(
        RdkitMmffField::with_terms_and_partial_charges(
            &molecule,
            &[[0.0; 3]],
            MmffVariant::Mmff94,
            100.0,
            true,
            terms,
            &[0.0, 0.0],
        ),
        Err(MmffPartialChargeError::Coordinates)
    ));
    assert!(matches!(
        RdkitMmffField::with_terms_and_partial_charges(
            &molecule,
            &coords,
            MmffVariant::Mmff94,
            -1.0,
            true,
            terms,
            &[0.0, 0.0],
        ),
        Err(MmffPartialChargeError::Threshold)
    ));
}
