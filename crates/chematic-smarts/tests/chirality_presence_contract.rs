use chematic_smarts::{MatchConfig, find_matches_with_config, parse_smarts};
use chematic_smiles::parse;

#[test]
fn underconstrained_chiral_query_checks_presence_not_orientation() {
    // RDKit 2025.09.3 HasSubstructMatch(..., useChirality=True).
    let config = MatchConfig {
        use_chirality: true,
        ..Default::default()
    };
    for smarts in ["[C@]", "[C@@]"] {
        let query = parse_smarts(smarts).unwrap();
        for (smiles, expected) in [
            ("F[C@](Cl)(Br)I", true),
            ("F[C@@](Cl)(Br)I", true),
            ("FC(Cl)(Br)I", false),
        ] {
            assert_eq!(
                !find_matches_with_config(&query, &parse(smiles).unwrap(), &config).is_empty(),
                expected,
                "{smarts} against {smiles}"
            );
        }
    }
}
