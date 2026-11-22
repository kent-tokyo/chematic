use chematic_3d::{RdkitUsrError, rdkit_usr_descriptors};

#[test]
fn usr_errors_explain_input_and_arithmetic_failures() {
    let too_few = rdkit_usr_descriptors(&[[0.; 3]; 2]).unwrap_err();
    assert_eq!(too_few, RdkitUsrError::TooFewAtoms { count: 2 });
    assert!(too_few.to_string().contains("at least three atoms"));

    let nonfinite = rdkit_usr_descriptors(&[[0.; 3], [f64::NAN, 0., 0.], [0.; 3]]).unwrap_err();
    assert_eq!(
        nonfinite,
        RdkitUsrError::NonFiniteCoordinate { atom_index: 1 }
    );
    assert!(nonfinite.to_string().contains("coordinate 1"));

    let overflow = rdkit_usr_descriptors(&[[f64::MAX; 3]; 3]).unwrap_err();
    assert_eq!(overflow, RdkitUsrError::NonFiniteResult);
    assert!(overflow.to_string().contains("moments are not finite"));
}
