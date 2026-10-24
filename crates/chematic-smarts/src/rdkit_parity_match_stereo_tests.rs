use super::*;
use crate::parser::parse_smarts;
use chematic_smiles::parse;

#[test]
fn parity_matcher_enforces_directional_double_bond_stereo() {
    let query = parse_smarts("F/C=C/F").unwrap();
    let config = RdkitParityConfig {
        base: MatchConfig {
            use_chirality: true,
            ..MatchConfig::default()
        },
        ..RdkitParityConfig::default()
    };
    for (target, expected) in [("F/C=C/F", true), ("F/C=C\\F", false), ("FC=CF", false)] {
        let (matches, exhausted) =
            find_matches_rdkit_parity(&query, &parse(target).unwrap(), &config).unwrap();
        assert!(!exhausted);
        assert_eq!(!matches.is_empty(), expected, "{target}");
    }
}
