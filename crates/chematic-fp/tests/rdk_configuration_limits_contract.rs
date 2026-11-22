use chematic_fp::{RdkitRdkConfig, RdkitRdkConfigError, rdkit_rdk_config};
use chematic_smiles::parse;

#[test]
fn configured_rdk_rejects_invalid_limits_with_typed_errors() {
    let molecule = parse("CCO").unwrap();
    let cases = [
        (
            RdkitRdkConfig {
                min_path: 0,
                ..Default::default()
            },
            RdkitRdkConfigError::PathRange,
        ),
        (
            RdkitRdkConfig {
                nbits: 0,
                ..Default::default()
            },
            RdkitRdkConfigError::FingerprintSize,
        ),
        (
            RdkitRdkConfig {
                count_bounds: Some(vec![]),
                ..Default::default()
            },
            RdkitRdkConfigError::CountBounds,
        ),
        (
            RdkitRdkConfig {
                atom_invariants: Some(vec![1]),
                ..Default::default()
            },
            RdkitRdkConfigError::AtomInvariants,
        ),
        (
            RdkitRdkConfig {
                max_visits: 0,
                ..Default::default()
            },
            RdkitRdkConfigError::Budget,
        ),
    ];
    for (config, expected) in cases {
        assert_eq!(rdkit_rdk_config(&molecule, &config).unwrap_err(), expected);
    }
    let many_atoms = std::iter::repeat_n("C", 1025).collect::<Vec<_>>().join(".");
    assert_eq!(
        rdkit_rdk_config(&parse(&many_atoms).unwrap(), &Default::default()).unwrap_err(),
        RdkitRdkConfigError::Graph
    );
}
