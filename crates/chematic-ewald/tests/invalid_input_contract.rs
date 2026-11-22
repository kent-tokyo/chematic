use chematic_ewald::{BoxVectors, EwaldError, PmeConfig, spme_energy};

#[test]
fn invalid_box_and_configuration_report_public_error() {
    let coords = [[0.0, 0.0, 0.0]];
    let charges = [1.0];
    let config = PmeConfig::default();
    let error = spme_energy(&coords, &charges, &BoxVectors::cubic(0.0), &config).unwrap_err();
    assert_eq!(error, EwaldError::SingularBoxMatrix);
    assert_eq!(
        error.to_string(),
        "invalid PME input, configuration, or singular box matrix"
    );
    let invalid = PmeConfig {
        r_cut: 0.0,
        ..config
    };
    assert_eq!(
        spme_energy(&coords, &charges, &BoxVectors::cubic(10.0), &invalid),
        Err(EwaldError::SingularBoxMatrix)
    );
}
