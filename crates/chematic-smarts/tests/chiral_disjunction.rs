use chematic_smarts::{
    MatchConfig, RdkitParityConfig, find_matches_rdkit_parity, find_matches_with_config,
    parse_smarts,
};

#[test]
fn disjunction_accepts_both_enantiomers_but_not_unspecified_stereo() {
    let query = parse_smarts("F[C@,C@@](Cl)(Br)I").unwrap();
    let config = MatchConfig {
        use_chirality: true,
        ..Default::default()
    };
    for (smiles, expected) in [
        ("F[C@](Cl)(Br)I", true),
        ("F[C@@](Cl)(Br)I", true),
        ("FC(Cl)(Br)I", false),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        assert_eq!(
            !find_matches_with_config(&query, &mol, &config).is_empty(),
            expected,
            "{smiles}"
        );
        let parity = RdkitParityConfig {
            base: config.clone(),
            ..Default::default()
        };
        let (matches, exhausted) = find_matches_rdkit_parity(&query, &mol, &parity).unwrap();
        assert!(!exhausted);
        assert_eq!(!matches.is_empty(), expected, "{smiles}");
    }
}
