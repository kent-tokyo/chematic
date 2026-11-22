use chematic_fp::{
    RdkitAtomPairConfig as Config, RdkitAtomPairError as Error, rdkit_atom_pair_config,
};

#[test]
fn invalid_atom_pair_configuration_identifies_the_failed_boundary() {
    let molecule = chematic_smiles::parse("CCO").unwrap();
    let cases = [
        (
            Config {
                min_distance: 3,
                max_distance: 2,
                ..Default::default()
            },
            Error::InvalidDistances,
            "Atom Pair requires 0 <= min_distance <= max_distance <= 30",
        ),
        (
            Config {
                nbits: 0,
                ..Default::default()
            },
            Error::InvalidBitWidth,
            "Atom Pair bit width must be between 1 and 1048576",
        ),
        (
            Config {
                count_bounds: Some(vec![]),
                ..Default::default()
            },
            Error::InvalidCountBounds,
            "Atom Pair count simulation needs nonempty bounds shorter than nbits and at most 2147483647",
        ),
    ];
    for (config, expected, message) in cases {
        let error = rdkit_atom_pair_config(&molecule, &config).unwrap_err();
        assert_eq!(error, expected);
        assert_eq!(error.to_string(), message);
    }
}
