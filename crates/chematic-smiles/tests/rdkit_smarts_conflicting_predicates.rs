use chematic_smiles::rdkit_smarts_to_smarts;

#[test]
fn contradictory_smarts_predicates_remain_explicit_in_written_query() {
    // RDKit 2025.09.3 MolToSmarts(MolFromSmarts(...)).
    for (input, expected) in [
        ("[C;H1;H2]", "[C&H1&H2]"),
        ("[N;+1;-1]", "[N&+&-]"),
        ("[C;H1;+1;H2;-1]", "[C&H1&+&H2&-]"),
        ("[CH1H2]", "[C&H1&H2]"),
        ("[N+1-1]", "[N&+&-]"),
        ("[CH1H2+1-1]", "[C&H1&H2&+&-]"),
    ] {
        assert_eq!(rdkit_smarts_to_smarts(input).unwrap(), expected, "{input}");
    }
}
