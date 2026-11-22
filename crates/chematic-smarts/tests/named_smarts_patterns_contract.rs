use chematic_smarts::{SmartsCache, named_pattern};

#[test]
fn named_patterns_distinguish_common_functional_groups() {
    let cases = [
        ("aromatic_ring", "c1ccccc1", "C1CCCCC1"),
        ("amine_tertiary", "CN(C)C", "CNC"),
        ("ether", "COC", "CCO"),
        ("aromatic_n", "n1ccccc1", "C1CCNCC1"),
        ("sulfonamide", "CS(=O)(=O)N", "CS(=O)(=O)O"),
        ("acceptor_strict", "n1ccccc1", "C1CCCCC1"),
    ];
    let mut cache = SmartsCache::new(cases.len());
    for (name, positive, negative) in cases {
        let pattern = named_pattern(name).unwrap();
        assert!(
            cache
                .has_match(pattern, &chematic_smiles::parse(positive).unwrap())
                .unwrap(),
            "{name} should match {positive}"
        );
        assert!(
            !cache
                .has_match(pattern, &chematic_smiles::parse(negative).unwrap())
                .unwrap(),
            "{name} should not match {negative}"
        );
    }
    assert!(named_pattern("unknown-pattern").is_none());
}
