use chematic_smarts::{
    MatchConfig, RdkitParityConfig, find_matches_rdkit_parity, find_matches_with_config,
    parse_smarts,
};

#[test]
fn atom_chirality_negation_complements_positive_constraint_after_mapping() {
    let positive = parse_smarts("F[C;@](Cl)(Br)I").unwrap();
    let negative = parse_smarts("F[C;!@](Cl)(Br)I").unwrap();
    let config = MatchConfig {
        use_chirality: true,
        ..Default::default()
    };
    for smiles in ["F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I", "FC(Cl)(Br)I"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let yes = !find_matches_with_config(&positive, &mol, &config).is_empty();
        let no = !find_matches_with_config(&negative, &mol, &config).is_empty();
        assert_ne!(yes, no, "{smiles}");
        let parity = RdkitParityConfig {
            base: config.clone(),
            ..Default::default()
        };
        let (yes_parity, exhausted) = find_matches_rdkit_parity(&positive, &mol, &parity).unwrap();
        assert!(!exhausted);
        let (no_parity, exhausted) = find_matches_rdkit_parity(&negative, &mol, &parity).unwrap();
        assert!(!exhausted);
        assert_eq!(!yes_parity.is_empty(), yes, "{smiles}");
        assert_eq!(!no_parity.is_empty(), no, "{smiles}");
    }
}
