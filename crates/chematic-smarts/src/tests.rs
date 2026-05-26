//! Integration tests for chematic-smarts: parser + VF2 matcher.

#[cfg(test)]
mod integration_tests {
    use chematic_smiles::parse as parse_mol;

    use crate::find_matches;
    use crate::parse_smarts;
    use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

    // -----------------------------------------------------------------------
    // Helper
    // -----------------------------------------------------------------------

    /// Count the number of matches found for `smarts` in `smiles`.
    fn match_count(smarts: &str, smiles: &str) -> usize {
        let query = parse_smarts(smarts).expect("parse_smarts failed");
        let mol = parse_mol(smiles).expect("parse_mol failed");
        find_matches(&query, &mol).len()
    }

    // -----------------------------------------------------------------------
    // Parser tests (tests 1–12)
    // -----------------------------------------------------------------------

    /// Test 1: `C` → 1 atom, Symbol("C"), aliphatic.
    #[test]
    fn test_parser_aliphatic_c() {
        let mol = parse_smarts("C").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.bonds.len(), 0);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
        );
        assert_eq!(mol.atoms[0].query, expected, "C should parse as aliphatic carbon");
    }

    /// Test 2: `c` → 1 atom, aromatic C.
    #[test]
    fn test_parser_aromatic_c() {
        let mol = parse_smarts("c").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
        );
        assert_eq!(mol.atoms[0].query, expected, "c should parse as aromatic carbon");
    }

    /// Test 3: `[#6]` → 1 atom, AtomicNum(6).
    #[test]
    fn test_parser_atomic_num() {
        let mol = parse_smarts("[#6]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(6))
        );
    }

    /// Test 4: `[!C]` → Not(Symbol("C")).
    #[test]
    fn test_parser_not() {
        let mol = parse_smarts("[!C]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Not(Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol(
                "C".to_string()
            ))))
        );
    }

    /// Test 5: `[a]` → Aromatic(true).
    #[test]
    fn test_parser_aromatic_primitive() {
        let mol = parse_smarts("[a]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Aromatic(true))
        );
    }

    /// Test 6: `[D3]` → Degree(3).
    #[test]
    fn test_parser_degree() {
        let mol = parse_smarts("[D3]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Degree(3))
        );
    }

    /// Test 7: `[r5]` → RingSize(5).
    #[test]
    fn test_parser_ring_size() {
        let mol = parse_smarts("[r5]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::RingSize(5))
        );
    }

    /// Test 8: `[H2]` → HCount(2).
    #[test]
    fn test_parser_hcount() {
        let mol = parse_smarts("[H2]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::HCount(2))
        );
    }

    /// Test 9: `CC` → 2 atoms, 1 bond (Any).
    #[test]
    fn test_parser_cc_implicit_bond() {
        let mol = parse_smarts("CC").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].query, BondQuery::Any);
    }

    /// Test 10: `C=C` → 2 atoms, 1 Double bond.
    #[test]
    fn test_parser_double_bond() {
        let mol = parse_smarts("C=C").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(
            mol.bonds[0].query,
            BondQuery::Primitive(BondPrimitive::Double)
        );
    }

    /// Test 11: `C(=O)O` → 3 atoms, 2 bonds (C=O and C~O).
    #[test]
    fn test_parser_branch() {
        let mol = parse_smarts("C(=O)O").unwrap();
        assert_eq!(mol.atoms.len(), 3, "should have 3 atoms (C, O, O)");
        assert_eq!(mol.bonds.len(), 2, "should have 2 bonds");
        // First bond: C=O (explicit double)
        assert_eq!(
            mol.bonds[0].query,
            BondQuery::Primitive(BondPrimitive::Double),
            "first bond should be Double"
        );
        // Second bond: C-O (implicit Any)
        assert_eq!(
            mol.bonds[1].query,
            BondQuery::Any,
            "second bond should be implicit Any"
        );
    }

    /// Test 12: `c1ccccc1` → 6 aromatic C atoms, 6 ring bonds.
    #[test]
    fn test_parser_benzene_ring() {
        let mol = parse_smarts("c1ccccc1").unwrap();
        assert_eq!(mol.atoms.len(), 6, "benzene ring has 6 atoms");
        assert_eq!(mol.bonds.len(), 6, "benzene ring has 6 bonds");
        for atom in &mol.atoms {
            let expected = AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
            );
            assert_eq!(atom.query, expected, "each atom should be aromatic C");
        }
    }

    // -----------------------------------------------------------------------
    // Matching tests (tests 13–22)
    // -----------------------------------------------------------------------

    /// Test 13: `C` matches all aliphatic C atoms in ethane → 2 matches.
    #[test]
    fn test_match_aliphatic_c_in_ethane() {
        assert_eq!(
            match_count("C", "CC"),
            2,
            "C should match 2 times in ethane"
        );
    }

    /// Test 14: `[#6]` matches all 6 atoms in benzene → 6 matches.
    #[test]
    fn test_match_atomic_num_in_benzene() {
        assert_eq!(
            match_count("[#6]", "c1ccccc1"),
            6,
            "[#6] should match 6 atoms in benzene"
        );
    }

    /// Test 15: `[a]` matches all 6 aromatic atoms in benzene → 6 matches.
    #[test]
    fn test_match_aromatic_in_benzene() {
        assert_eq!(
            match_count("[a]", "c1ccccc1"),
            6,
            "[a] should match 6 aromatic atoms in benzene"
        );
    }

    /// Test 16: `[A]` matches 0 atoms in benzene (all are aromatic) → 0 matches.
    #[test]
    fn test_match_aliphatic_in_benzene() {
        assert_eq!(
            match_count("[A]", "c1ccccc1"),
            0,
            "[A] (aliphatic) should match 0 atoms in benzene"
        );
    }

    /// Test 17: `C=O` matches 2 times in aspirin → 2 matches.
    #[test]
    fn test_match_co_in_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        // Two C=O groups: acetyl and carboxyl.
        let n = match_count("C=O", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n, 2, "C=O should match 2 times in aspirin, got {n}");
    }

    /// Test 18: `[OH]` matches 1 time in aspirin (the carboxylic OH) → 1 match.
    #[test]
    fn test_match_oh_in_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        // The terminal O in C(=O)O has 1 implicit H; the ester O has 0 H.
        let n = match_count("[OH]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n, 1, "[OH] should match 1 time in aspirin, got {n}");
    }

    /// Test 19: `c1ccccc1` (benzene ring pattern) matches the ring in toluene → >= 1 match.
    #[test]
    fn test_match_benzene_ring_in_toluene() {
        let n = match_count("c1ccccc1", "Cc1ccccc1");
        assert!(
            n >= 1,
            "benzene ring pattern should match in toluene, got {n} matches"
        );
    }

    /// Test 20: `[R]` (any ring atom) matches 6 atoms in benzene → 6 matches.
    #[test]
    fn test_match_ring_membership_in_benzene() {
        assert_eq!(
            match_count("[R]", "c1ccccc1"),
            6,
            "[R] should match 6 ring atoms in benzene"
        );
    }

    /// Test 21: `[D1]` (terminal atom) matches 1 atom in toluene → 1 match (the methyl C).
    #[test]
    fn test_match_terminal_in_toluene() {
        let n = match_count("[D1]", "Cc1ccccc1");
        assert_eq!(
            n, 1,
            "[D1] should match 1 terminal atom in toluene, got {n}"
        );
    }

    /// Test 22: `*` (wildcard) matches all atoms in ethane → 2 matches.
    #[test]
    fn test_match_wildcard_in_ethane() {
        assert_eq!(
            match_count("*", "CC"),
            2,
            "* should match all 2 atoms in ethane"
        );
    }
}
