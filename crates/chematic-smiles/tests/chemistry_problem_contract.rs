use chematic_smiles::{
    parse, parse_template, rdkit_detect_chemistry_problems,
    rdkit_detect_chemistry_problems_sanitized,
};

#[test]
fn chemistry_problems_identify_valence_and_aromaticity_failures() {
    for (smiles, kind, atoms) in [
        ("c", "AtomKekulizeException", vec![0]),
        ("c1cccc1", "KekulizeException", vec![0, 1, 2, 3, 4]),
        ("C(F)(F)(F)(F)F", "AtomValenceException", vec![0]),
        ("O(C)(C)C", "AtomValenceException", vec![0]),
    ] {
        let molecule = parse(smiles).or_else(|_| parse_template(smiles)).unwrap();
        let problems = rdkit_detect_chemistry_problems(&molecule).unwrap();
        assert_eq!(problems.len(), 1, "{smiles}");
        assert_eq!(problems[0].kind, kind, "{smiles}");
        assert_eq!(problems[0].atoms, atoms, "{smiles}");
        assert!(
            rdkit_detect_chemistry_problems_sanitized(&molecule).is_err(),
            "{smiles}"
        );
    }
    let valid = parse("C1=CC=CC=C1").unwrap();
    assert!(rdkit_detect_chemistry_problems(&valid).unwrap().is_empty());
    assert!(
        rdkit_detect_chemistry_problems_sanitized(&valid)
            .unwrap()
            .is_empty()
    );
}
