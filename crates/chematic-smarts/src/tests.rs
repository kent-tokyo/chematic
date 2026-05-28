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

    // -----------------------------------------------------------------------
    // Recursive SMARTS `$(...)` tests (tests 23–30)
    // -----------------------------------------------------------------------

    /// Test 23: `[$(C(=O)O)]` — carboxyl C. Parses to 1 atom with Recursive inner.
    #[test]
    fn test_parser_recursive_smarts_structure() {
        let q = parse_smarts("[$(C(=O)O)]").unwrap();
        assert_eq!(q.atoms.len(), 1, "outer query has 1 atom");
        match &q.atoms[0].query {
            AtomQuery::Primitive(AtomPrimitive::Recursive(inner)) => {
                assert_eq!(inner.atoms.len(), 3, "inner C(=O)O has 3 atoms");
                assert_eq!(inner.bonds.len(), 2, "inner C(=O)O has 2 bonds");
            }
            other => panic!("expected Recursive, got {:?}", other),
        }
    }

    /// Test 24: `[$(C(=O)O)]` matches the carboxyl C in acetic acid but not in acetone.
    #[test]
    fn test_recursive_carboxylic_acid_carbon() {
        // Acetic acid: CH3-C(=O)-OH; the carbonyl C is bonded to both =O and -OH
        assert_eq!(
            match_count("[$(C(=O)O)]", "CC(=O)O"),
            1,
            "acetic acid should have 1 carboxylic C"
        );
        // Acetone: CH3-C(=O)-CH3; no adjacent OH oxygen
        assert_eq!(
            match_count("[$(C(=O)O)]", "CC(=O)C"),
            0,
            "acetone should have 0 carboxylic C"
        );
    }

    /// Test 25: single-atom recursive `[$([OH])]` matches the O in methanol.
    #[test]
    fn test_recursive_single_atom_oh() {
        // Methanol CO: O has 1 implicit H
        assert_eq!(
            match_count("[$([OH])]", "CO"),
            1,
            "methanol O should match [$([OH])]"
        );
        // Dimethyl ether COC: O has 0 implicit H
        assert_eq!(
            match_count("[$([OH])]", "COC"),
            0,
            "dimethyl ether O should not match [$([OH])]"
        );
    }

    /// Test 26: `[N;!$(NC=O)]` — non-amide nitrogen.
    #[test]
    fn test_recursive_non_amide_nitrogen() {
        // Methylamine CN: N is not in an amide → should match
        assert_eq!(
            match_count("[N;!$(NC=O)]", "CN"),
            1,
            "methylamine N should match [N;!$(NC=O)]"
        );
        // Acetamide CC(=O)N: N is bonded to C=O → should NOT match
        assert_eq!(
            match_count("[N;!$(NC=O)]", "CC(=O)N"),
            0,
            "acetamide N should not match [N;!$(NC=O)]"
        );
    }

    /// Test 27: `[$(CC)]` matches aliphatic C bonded to another C in propane.
    #[test]
    fn test_recursive_cc_in_propane() {
        // Propane CCC: middle C and both terminal C are bonded to C
        let n = match_count("[$(CC)]", "CCC");
        assert_eq!(n, 3, "all 3 carbons in propane match [$(CC)] — each is bonded to a C");
    }

    /// Test 28: nested recursive `[$([C;$(C(=O))])]` matches carbonyl C.
    #[test]
    fn test_recursive_nested_parse() {
        // Should parse without error
        let q = parse_smarts("[$([C;$(C(=O))])]").unwrap();
        assert_eq!(q.atoms.len(), 1, "outer has 1 atom");
        // Verify matching: carbonyl C in acetone CC(=O)C should match
        let mol = parse_mol("CC(=O)C").unwrap();
        let matches = find_matches(&q, &mol);
        assert_eq!(matches.len(), 1, "acetone has 1 carbonyl C");
    }

    /// Test 29: `[$(c1ccccc1)]` matches all 6 atoms in benzene.
    #[test]
    fn test_recursive_benzene_pattern() {
        // Every atom in benzene is part of a benzene ring
        let n = match_count("[$(c1ccccc1)]", "c1ccccc1");
        assert_eq!(n, 6, "all 6 benzene atoms match [$(c1ccccc1)]");
    }

    // -----------------------------------------------------------------------
    // New SMARTS primitives: [v], [x], [^], [+0] (tests 31–42)
    // -----------------------------------------------------------------------

    /// Test 31: `[v4]` matches atoms with total valence 4.
    #[test]
    fn test_valence_methane_carbon() {
        // Methane C: 4 implicit H → valence 4
        assert_eq!(match_count("[v4]", "C"), 1, "[v4] should match CH4 carbon");
    }

    /// Test 32: `[v2]` matches oxygen in water (2 H → valence 2).
    #[test]
    fn test_valence_water_oxygen() {
        assert_eq!(match_count("[v2]", "O"), 1, "[v2] should match water O");
    }

    /// Test 33: `[v3]` matches nitrogen in methylamine (1 C + 2 H → valence 3).
    #[test]
    fn test_valence_methylamine_nitrogen() {
        assert_eq!(match_count("[v3]", "CN"), 1, "[v3] should match methylamine N");
    }

    /// Test 34: `[x0]` matches atoms with 0 ring bonds (acyclic).
    #[test]
    fn test_ring_bond_count_zero_in_ethanol() {
        // Ethanol CCO: all atoms are acyclic → x0 for all 3
        assert_eq!(match_count("[x0]", "CCO"), 3, "[x0] should match all 3 atoms in ethanol");
    }

    /// Test 35: `[x2]` matches atoms that have exactly 2 ring bonds in benzene.
    #[test]
    fn test_ring_bond_count_two_benzene() {
        // Each C in benzene has 2 ring bonds
        assert_eq!(match_count("[x2]", "c1ccccc1"), 6, "[x2] should match all 6 benzene atoms");
    }

    /// Test 36: `[x0]` matches 0 atoms in benzene (all have 2 ring bonds).
    #[test]
    fn test_ring_bond_count_zero_not_in_benzene() {
        assert_eq!(match_count("[x0]", "c1ccccc1"), 0, "[x0] should match 0 atoms in benzene");
    }

    /// Test 37: `[^3]` matches sp3 atoms (aliphatic C with no double/triple bonds).
    #[test]
    fn test_hybridization_sp3_ethane() {
        assert_eq!(match_count("[^3]", "CC"), 2, "[^3] should match both sp3 C atoms in ethane");
    }

    /// Test 38: `[^2]` matches sp2 atoms — both C and O in ethylene/carbonyl are sp2.
    #[test]
    fn test_hybridization_sp2_ethylene() {
        // Ethylene C=C: both carbons are sp2 (double bond)
        assert_eq!(match_count("[^2]", "C=C"), 2, "[^2] should match both sp2 C atoms in ethylene");
    }

    /// Test 39: `[^2]` matches all 6 aromatic C atoms in benzene.
    #[test]
    fn test_hybridization_sp2_benzene() {
        assert_eq!(match_count("[^2]", "c1ccccc1"), 6, "[^2] should match all 6 aromatic C in benzene");
    }

    /// Test 40: `[^1]` matches sp atoms in acetylene (HC≡CH).
    #[test]
    fn test_hybridization_sp_acetylene() {
        assert_eq!(match_count("[^1]", "C#C"), 2, "[^1] should match both sp C atoms in acetylene");
    }

    /// Test 41: `[+0]` matches neutral atoms (explicit zero charge).
    #[test]
    fn test_charge_explicit_zero() {
        // Trimethylammonium: CC[N+](C)C — the N has +1 charge
        // [+0] should match the 4 carbons (neutral) but NOT the N+ (charge +1)
        assert_eq!(match_count("[+0]", "CC[N+](C)C"), 4, "[+0] should match 4 neutral C atoms");
    }

    /// Test 42: `[-0]` is also a valid explicit-zero-charge query (same as [+0]).
    #[test]
    fn test_charge_negative_zero() {
        // [-0] means charge == -(0) == 0, so same as [+0]
        assert_eq!(match_count("[-0]", "CCO"), 3, "[-0] should match all 3 neutral atoms in ethanol");
    }

    /// Test 30: `[$(C(=O)O)]` finds both ester and carboxyl C in aspirin.
    #[test]
    fn test_recursive_carbonyl_with_oxygen_in_aspirin() {
        // Aspirin CC(=O)Oc1ccccc1C(=O)O has two C(=O)-O groups:
        // the acetyl ester C and the carboxylic acid C → both match [$(C(=O)O)].
        // To match ONLY the carboxylic acid, use [$(C(=O)[OH])].
        let n_cooh = match_count("[$(C(=O)O)]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n_cooh, 2, "aspirin has 2 C atoms with C(=O)O (ester + carboxylic), got {n_cooh}");
        // [$(C(=O)[OH])] is more specific: only the COOH group matches.
        let n_acid = match_count("[$(C(=O)[OH])]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n_acid, 1, "aspirin has 1 carboxylic acid C [$(C(=O)[OH])], got {n_acid}");
    }
}
