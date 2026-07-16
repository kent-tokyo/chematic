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
    // Parser tests
    // -----------------------------------------------------------------------

    /// `C` → 1 atom, Symbol("C"), aliphatic.
    #[test]
    fn test_parser_aliphatic_c() {
        let mol = parse_smarts("C").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.bonds.len(), 0);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
        );
        assert_eq!(
            mol.atoms[0].query, expected,
            "C should parse as aliphatic carbon"
        );
    }

    /// `c` → 1 atom, aromatic C.
    #[test]
    fn test_parser_aromatic_c() {
        let mol = parse_smarts("c").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
        );
        assert_eq!(
            mol.atoms[0].query, expected,
            "c should parse as aromatic carbon"
        );
    }

    /// `[#6]` → 1 atom, AtomicNum(6).
    #[test]
    fn test_parser_atomic_num() {
        let mol = parse_smarts("[#6]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(6))
        );
    }

    /// `[!C]` → Not(Symbol("C")).
    #[test]
    fn test_parser_not() {
        // [C] inside brackets means aliphatic C (uppercase → Aromatic(false) added).
        let mol = parse_smarts("[!C]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Not(Box::new(AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
            )))
        );
    }

    /// `[a]` → Aromatic(true).
    #[test]
    fn test_parser_aromatic_primitive() {
        let mol = parse_smarts("[a]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Aromatic(true))
        );
    }

    /// `[D3]` → Degree(3).
    #[test]
    fn test_parser_degree() {
        let mol = parse_smarts("[D3]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Degree(3))
        );
    }

    /// `[r5]` → MinRingSize(5) -- smallest-ring semantics (SMARTS-R1).
    #[test]
    fn test_parser_ring_size() {
        let mol = parse_smarts("[r5]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::MinRingSize(5))
        );
    }

    /// `[H2]` → HCount(2).
    #[test]
    fn test_parser_hcount() {
        let mol = parse_smarts("[H2]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::HCount(2))
        );
    }

    /// `CC` → 2 atoms, 1 bond (Any).
    #[test]
    fn test_parser_cc_implicit_bond() {
        let mol = parse_smarts("CC").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].query, BondQuery::Any);
    }

    /// `C=C` → 2 atoms, 1 Double bond.
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

    /// `C(=O)O` → 3 atoms, 2 bonds (C=O and C~O).
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

    /// `c1ccccc1` → 6 aromatic C atoms, 6 ring bonds.
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
    // Matching tests
    // -----------------------------------------------------------------------

    /// `C` matches all aliphatic C atoms in ethane → 2 matches.
    #[test]
    fn test_match_aliphatic_c_in_ethane() {
        assert_eq!(
            match_count("C", "CC"),
            2,
            "C should match 2 times in ethane"
        );
    }

    /// `[#6]` matches all 6 atoms in benzene → 6 matches.
    #[test]
    fn test_match_atomic_num_in_benzene() {
        assert_eq!(
            match_count("[#6]", "c1ccccc1"),
            6,
            "[#6] should match 6 atoms in benzene"
        );
    }

    /// `[a]` matches all 6 aromatic atoms in benzene → 6 matches.
    #[test]
    fn test_match_aromatic_in_benzene() {
        assert_eq!(
            match_count("[a]", "c1ccccc1"),
            6,
            "[a] should match 6 aromatic atoms in benzene"
        );
    }

    /// `[A]` matches 0 atoms in benzene (all are aromatic) → 0 matches.
    #[test]
    fn test_match_aliphatic_in_benzene() {
        assert_eq!(
            match_count("[A]", "c1ccccc1"),
            0,
            "[A] (aliphatic) should match 0 atoms in benzene"
        );
    }

    /// `C=O` matches 2 times in aspirin → 2 matches.
    #[test]
    fn test_match_co_in_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        // Two C=O groups: acetyl and carboxyl.
        let n = match_count("C=O", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n, 2, "C=O should match 2 times in aspirin, got {n}");
    }

    /// `[OH]` matches 1 time in aspirin (the carboxylic OH) → 1 match.
    #[test]
    fn test_match_oh_in_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        // The terminal O in C(=O)O has 1 implicit H; the ester O has 0 H.
        let n = match_count("[OH]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(n, 1, "[OH] should match 1 time in aspirin, got {n}");
    }

    /// `c1ccccc1` (benzene ring pattern) matches the ring in toluene → >= 1 match.
    #[test]
    fn test_match_benzene_ring_in_toluene() {
        let n = match_count("c1ccccc1", "Cc1ccccc1");
        assert!(
            n >= 1,
            "benzene ring pattern should match in toluene, got {n} matches"
        );
    }

    /// `[R]` (any ring atom) matches 6 atoms in benzene → 6 matches.
    #[test]
    fn test_match_ring_membership_in_benzene() {
        assert_eq!(
            match_count("[R]", "c1ccccc1"),
            6,
            "[R] should match 6 ring atoms in benzene"
        );
    }

    /// `[D1]` (terminal atom) matches 1 atom in toluene → 1 match (the methyl C).
    #[test]
    fn test_match_terminal_in_toluene() {
        let n = match_count("[D1]", "Cc1ccccc1");
        assert_eq!(
            n, 1,
            "[D1] should match 1 terminal atom in toluene, got {n}"
        );
    }

    /// `*` (wildcard) matches all atoms in ethane → 2 matches.
    #[test]
    fn test_match_wildcard_in_ethane() {
        assert_eq!(
            match_count("*", "CC"),
            2,
            "* should match all 2 atoms in ethane"
        );
    }

    // -----------------------------------------------------------------------
    // Recursive SMARTS `$(...)` tests
    // -----------------------------------------------------------------------

    /// `[$(C(=O)O)]` — carboxyl C. Parses to 1 atom with Recursive inner.
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

    /// `[$(C(=O)O)]` matches the carboxyl C in acetic acid but not in acetone.
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

    /// Single-atom recursive `[$([OH])]` matches the O in methanol.
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

    /// `[N;!$(NC=O)]` — non-amide nitrogen.
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

    /// `[$(CC)]` matches aliphatic C bonded to another C in propane.
    #[test]
    fn test_recursive_cc_in_propane() {
        // Propane CCC: middle C and both terminal C are bonded to C
        let n = match_count("[$(CC)]", "CCC");
        assert_eq!(
            n, 3,
            "all 3 carbons in propane match [$(CC)] — each is bonded to a C"
        );
    }

    /// Nested recursive `[$([C;$(C(=O))])]` matches carbonyl C.
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

    /// `[$(c1ccccc1)]` matches all 6 atoms in benzene.
    #[test]
    fn test_recursive_benzene_pattern() {
        // Every atom in benzene is part of a benzene ring
        let n = match_count("[$(c1ccccc1)]", "c1ccccc1");
        assert_eq!(n, 6, "all 6 benzene atoms match [$(c1ccccc1)]");
    }

    // -----------------------------------------------------------------------
    // SMARTS primitives: [v], [x], [^], [+0]
    // -----------------------------------------------------------------------

    /// `[v4]` matches atoms with total valence 4.
    #[test]
    fn test_valence_methane_carbon() {
        // Methane C: 4 implicit H → valence 4
        assert_eq!(match_count("[v4]", "C"), 1, "[v4] should match CH4 carbon");
    }

    /// `[v2]` matches oxygen in water (2 H → valence 2).
    #[test]
    fn test_valence_water_oxygen() {
        assert_eq!(match_count("[v2]", "O"), 1, "[v2] should match water O");
    }

    /// `[v3]` matches nitrogen in methylamine (1 C + 2 H → valence 3).
    #[test]
    fn test_valence_methylamine_nitrogen() {
        assert_eq!(
            match_count("[v3]", "CN"),
            1,
            "[v3] should match methylamine N"
        );
    }

    /// `[x0]` matches atoms with 0 ring bonds (acyclic).
    #[test]
    fn test_ring_bond_count_zero_in_ethanol() {
        // Ethanol CCO: all atoms are acyclic → x0 for all 3
        assert_eq!(
            match_count("[x0]", "CCO"),
            3,
            "[x0] should match all 3 atoms in ethanol"
        );
    }

    /// `[x2]` matches atoms that have exactly 2 ring bonds in benzene.
    #[test]
    fn test_ring_bond_count_two_benzene() {
        // Each C in benzene has 2 ring bonds
        assert_eq!(
            match_count("[x2]", "c1ccccc1"),
            6,
            "[x2] should match all 6 benzene atoms"
        );
    }

    /// `[x0]` matches 0 atoms in benzene (all have 2 ring bonds).
    #[test]
    fn test_ring_bond_count_zero_not_in_benzene() {
        assert_eq!(
            match_count("[x0]", "c1ccccc1"),
            0,
            "[x0] should match 0 atoms in benzene"
        );
    }

    /// `[^3]` matches sp3 atoms (aliphatic C with no double/triple bonds).
    #[test]
    fn test_hybridization_sp3_ethane() {
        assert_eq!(
            match_count("[^3]", "CC"),
            2,
            "[^3] should match both sp3 C atoms in ethane"
        );
    }

    /// `[^2]` matches sp2 atoms — both C and O in ethylene/carbonyl are sp2.
    #[test]
    fn test_hybridization_sp2_ethylene() {
        // Ethylene C=C: both carbons are sp2 (double bond)
        assert_eq!(
            match_count("[^2]", "C=C"),
            2,
            "[^2] should match both sp2 C atoms in ethylene"
        );
    }

    /// `[^2]` matches all 6 aromatic C atoms in benzene.
    #[test]
    fn test_hybridization_sp2_benzene() {
        assert_eq!(
            match_count("[^2]", "c1ccccc1"),
            6,
            "[^2] should match all 6 aromatic C in benzene"
        );
    }

    /// `[^1]` matches sp atoms in acetylene (HC≡CH).
    #[test]
    fn test_hybridization_sp_acetylene() {
        assert_eq!(
            match_count("[^1]", "C#C"),
            2,
            "[^1] should match both sp C atoms in acetylene"
        );
    }

    /// `[+0]` matches neutral atoms (explicit zero charge).
    #[test]
    fn test_charge_explicit_zero() {
        // Trimethylammonium: CC[N+](C)C — the N has +1 charge
        // [+0] should match the 4 carbons (neutral) but NOT the N+ (charge +1)
        assert_eq!(
            match_count("[+0]", "CC[N+](C)C"),
            4,
            "[+0] should match 4 neutral C atoms"
        );
    }

    /// `[-0]` is also a valid explicit-zero-charge query (same as [+0]).
    #[test]
    fn test_charge_negative_zero() {
        // [-0] means charge == -(0) == 0, so same as [+0]
        assert_eq!(
            match_count("[-0]", "CCO"),
            3,
            "[-0] should match all 3 neutral atoms in ethanol"
        );
    }

    /// `[$(C(=O)O)]` finds both ester and carboxyl C in aspirin.
    #[test]
    fn test_recursive_carbonyl_with_oxygen_in_aspirin() {
        // Aspirin CC(=O)Oc1ccccc1C(=O)O has two C(=O)-O groups:
        // the acetyl ester C and the carboxylic acid C → both match [$(C(=O)O)].
        // To match ONLY the carboxylic acid, use [$(C(=O)[OH])].
        let n_cooh = match_count("[$(C(=O)O)]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(
            n_cooh, 2,
            "aspirin has 2 C atoms with C(=O)O (ester + carboxylic), got {n_cooh}"
        );
        // [$(C(=O)[OH])] is more specific: only the COOH group matches.
        let n_acid = match_count("[$(C(=O)[OH])]", "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(
            n_acid, 1,
            "aspirin has 1 carboxylic acid C [$(C(=O)[OH])], got {n_acid}"
        );
    }

    // -----------------------------------------------------------------------
    // [XN] total connectivity tests
    // -----------------------------------------------------------------------

    /// `[X4]` matches methane carbon (D=0, implH=4 → X=4).
    #[test]
    fn test_total_connectivity_methane() {
        assert_eq!(match_count("[X4]", "C"), 1, "[X4] should match methane C");
        // [D4] would NOT match because methane C has degree 0
        assert_eq!(
            match_count("[D4]", "C"),
            0,
            "[D4] should NOT match methane C"
        );
    }

    /// `[X4]` matches both carbons in ethane (D=1, implH=3 → X=4).
    #[test]
    fn test_total_connectivity_ethane() {
        assert_eq!(
            match_count("[X4]", "CC"),
            2,
            "[X4] should match both C in ethane"
        );
    }

    /// `[X2]` matches oxygen in ethanol (D=1, implH=1 → X=2), not the carbons (X=4).
    #[test]
    fn test_total_connectivity_ethanol_oxygen() {
        assert_eq!(
            match_count("[X2]", "CCO"),
            1,
            "[X2] should match O in ethanol"
        );
        assert_eq!(
            match_count("[X4]", "CCO"),
            2,
            "[X4] should match 2 C atoms in ethanol"
        );
    }

    // -----------------------------------------------------------------------
    // Compound bond expressions: OR (`=,:`), AND+NOT (`=!@`)
    // -----------------------------------------------------------------------

    /// `[#6]=,:[#6]` matches both double bonds and aromatic bonds to carbon.
    #[test]
    fn test_bond_or_double_or_aromatic_in_benzene() {
        // Benzene: 6 aromatic C-C bonds; with uniquify=true (default), returns 6 unique mappings
        // (not 12 from both A→B and B→A directions)
        assert_eq!(
            match_count("[#6]=,:[#6]", "c1ccccc1"),
            6,
            "=,: should match all 6 aromatic bonds in benzene (deduplicated)"
        );
        // Pure single-bond query should NOT match aromatic bonds
        assert_eq!(
            match_count("[#6]-[#6]", "c1ccccc1"),
            0,
            "explicit - should not match aromatic bonds"
        );
    }

    /// `[#6]=,:[#6]` also matches a double bond in ethylene.
    #[test]
    fn test_bond_or_double_or_aromatic_in_ethylene() {
        // Ethylene C=C: 1 match (with uniquify=true, one unique atom set mapping)
        let c = match_count("[#6]=,:[#6]", "C=C");
        assert_eq!(
            c, 1,
            "=,: should match the double bond in ethylene (deduplicated)"
        );
    }

    /// `C=!@C` matches a non-ring double bond (e.g. ethylene) but NOT cyclohexene ring double bond.
    #[test]
    fn test_bond_non_ring_double_ethylene() {
        assert_eq!(
            match_count("C=!@C", "C=C"),
            1,
            "=!@ should match the non-ring double bond in ethylene (deduplicated)"
        );
    }

    /// `C=!@C` does NOT match the double bond inside a ring (cyclohexene).
    #[test]
    fn test_bond_non_ring_double_cyclohexene() {
        // Cyclohexene: the C=C is in-ring → =!@ should NOT match
        assert_eq!(
            match_count("C=!@C", "C1=CCCCC1"),
            0,
            "=!@ should not match in-ring double bond"
        );
    }

    /// `[#6]-,:c` (single OR aromatic bond to aromatic C) matches anisole Ar-O bond.
    #[test]
    fn test_bond_single_or_aromatic_anisole() {
        // Anisole: the O-c bond is a single bond in the parsed graph → -,: matches
        let c = match_count("[#8]-,:[c]", "COc1ccccc1");
        assert!(c >= 1, "-,: should match O-c bond in anisole");
    }

    /// Recursive SMARTS `$(…)` beyond depth 8 is rejected with RecursionDepthExceeded.
    #[test]
    fn test_recursive_smarts_depth_limit() {
        use crate::parser::SmartsError;
        // Build 9-level nested: start with "C", wrap 9 times as [$(<inner>)]
        // Each wrapping adds one recursion level when parsed.
        let mut p = "C".to_string();
        for _ in 0..9 {
            p = format!("[$({})]", p);
        }
        let err = parse_smarts(&p).unwrap_err();
        assert_eq!(
            err,
            SmartsError::RecursionDepthExceeded,
            "depth > 8 should return RecursionDepthExceeded, got {err:?}"
        );
    }

    /// Recursive SMARTS at depth ≤ 8 is accepted.
    #[test]
    fn test_recursive_smarts_shallow_ok() {
        // 1-level recursive: [$([C])] — should parse successfully.
        let mol = parse_smarts("[$([C])]");
        assert!(mol.is_ok(), "1-level recursive SMARTS should parse ok");
    }

    // -- atom map number `:N` -------------------------------------------------

    #[test]
    fn test_atom_map_parses_without_error() {
        // [O;D1;H0:3] — semicolon-AND, degree, H count, atom map — must not error.
        let q = parse_smarts("[O;D1;H0:3]").expect("[O;D1;H0:3] should parse");
        assert_eq!(q.atom_count(), 1);
        assert_eq!(q.atoms[0].atom_map, Some(3));
    }

    #[test]
    fn test_atom_map_matches_correctly() {
        use chematic_smiles::parse;
        let q = parse_smarts("[O;D1;H0:3]").unwrap();
        // Acetic acid: the ester oxygen (D1, H0) should match.
        let acetic = parse("CC(=O)O").unwrap();
        let matches = find_matches(&q, &acetic);
        assert!(
            !matches.is_empty(),
            "carbonyl O in acetic acid should match [O;D1;H0:3]"
        );
    }

    #[test]
    fn test_atom_map_stored_on_each_atom() {
        // Multi-atom SMARTS with maps: [C:1][O:2]
        let q = parse_smarts("[C:1][O:2]").expect("[C:1][O:2] should parse");
        assert_eq!(q.atom_count(), 2);
        assert_eq!(q.atoms[0].atom_map, Some(1));
        assert_eq!(q.atoms[1].atom_map, Some(2));
    }

    #[test]
    fn test_atom_map_zero_stored() {
        // :0 is a valid (if unusual) map number.
        let q = parse_smarts("[C:0]").expect("[C:0] should parse");
        assert_eq!(q.atoms[0].atom_map, Some(0));
    }

    #[test]
    fn test_no_atom_map_is_none() {
        let q = parse_smarts("[O;D1]").unwrap();
        assert_eq!(q.atoms[0].atom_map, None);
    }
}
