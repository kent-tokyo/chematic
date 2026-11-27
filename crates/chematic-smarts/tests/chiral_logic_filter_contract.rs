use chematic_smarts::{MatchConfig, find_matches_with_config, parse_smarts};
use chematic_smiles::parse;

#[test]
fn chirality_enabled_filter_preserves_logical_or_and_and() {
    let config = MatchConfig {
        use_chirality: true,
        ..Default::default()
    };
    for (pattern, target, expected) in [
        ("[#6,#7]", "C", true),
        ("[#6,#7]", "N", true),
        ("[#6,#7]", "O", false),
        ("[#6;!#7]", "C", true),
        ("[#6;!#7]", "N", false),
        ("[#6@,#8]", "O", true),
        ("[#6@,#8]", "C", false),
        ("[#6@,#8]", "[C@H](F)(Cl)Br", true),
    ] {
        let query = parse_smarts(pattern).unwrap();
        let mol = parse(target).unwrap();
        assert_eq!(
            !find_matches_with_config(&query, &mol, &config).is_empty(),
            expected,
            "{pattern} on {target}"
        );
    }
}
