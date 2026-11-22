use chematic_smiles::{SmilesError, SmilesParseLimits, parse_with_limits};

#[test]
fn hard_atom_limit_applies_to_chains_and_disconnected_fragments() {
    let limits = SmilesParseLimits {
        max_atoms: 100_001,
        ..Default::default()
    };
    for text in [
        "C".repeat(100_001),
        std::iter::repeat_n("C", 100_001)
            .collect::<Vec<_>>()
            .join("."),
    ] {
        let error = parse_with_limits(&text, &limits)
            .err()
            .expect("hard atom limit");
        assert!(
            matches!(error, SmilesError::InvalidBracketAtom { ref detail, .. } if detail.contains("maximum atom count 100000")),
            "{error}"
        );
    }
}
