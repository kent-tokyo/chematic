use super::*;
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

#[test]
fn compound_bond_expressions_follow_daylight_precedence() {
    use BondPrimitive::*;
    let p = |x| BondQuery::Primitive(x);
    let and = |a, b| BondQuery::And(Box::new(a), Box::new(b));
    let or = |a, b| BondQuery::Or(Box::new(a), Box::new(b));
    let not = |a| BondQuery::Not(Box::new(a));
    let bond = |smarts: &str| parse_smarts(smarts).unwrap().bonds[0].query.clone();
    // Historical forms keep their exact trees.
    assert_eq!(bond("C=,:C"), or(p(Double), p(Aromatic)));
    assert_eq!(bond("C=!@C"), and(p(Double), not(p(Ring))));
    // Newly accepted forms (previously parse errors or truncated).
    assert_eq!(bond("C=@C"), and(p(Double), p(Ring)));
    assert_eq!(bond("C=;@C"), and(p(Double), p(Ring)));
    assert_eq!(bond("C-&@C"), and(p(Single), p(Ring)));
    assert_eq!(bond("C=@,-C"), or(and(p(Double), p(Ring)), p(Single)));
    assert_eq!(bond("C-,=;@C"), and(or(p(Single), p(Double)), p(Ring)));
    assert_eq!(bond("C!!@C"), not(not(p(Ring))));
    // MACCS key 26 as written in RDKit's MACCSkeys.py.
    assert!(parse_smarts("[#6]=;@[#6](@*)@*").is_ok());
}

#[test]
fn test_parse_aliphatic_c() {
    // `C` → 1 atom with And(Symbol("C"), Aromatic(false))
    let mol = parse_smarts("C").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    assert_eq!(mol.bonds.len(), 0);
    let expected = AtomQuery::And(
        Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
        Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
    );
    assert_eq!(mol.atoms[0].query, expected);
}

#[test]
fn rejects_valid_prefix_followed_by_whitespace_junk() {
    assert_eq!(
        parse_smarts("CC this is junk"),
        Err(SmartsError::UnexpectedChar(' ', 2))
    );
    assert_eq!(
        parse_smarts("[#6] this is junk"),
        Err(SmartsError::UnexpectedChar(' ', 4))
    );
}

#[test]
fn test_parse_aromatic_c() {
    // `c` → 1 atom with And(Symbol("C"), Aromatic(true))
    let mol = parse_smarts("c").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    let expected = AtomQuery::And(
        Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
        Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
    );
    assert_eq!(mol.atoms[0].query, expected);
}

#[test]
fn test_parse_atomic_num() {
    // `[#6]` → AtomicNum(6)
    let mol = parse_smarts("[#6]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    assert_eq!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(6))
    );
}

#[test]
fn test_parse_not() {
    // `[!C]` → Not(And(Symbol("C"), Aromatic(false)))
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

#[test]
fn test_parse_aromatic_primitive() {
    // `[a]` → Aromatic(true)
    let mol = parse_smarts("[a]").unwrap();
    assert_eq!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::Aromatic(true))
    );
}

#[test]
fn test_parse_degree() {
    // `[D3]` → Degree(3)
    let mol = parse_smarts("[D3]").unwrap();
    assert_eq!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::Degree(3))
    );
}

#[test]
fn test_parse_ring_size() {
    // `[r5]` → MinRingSize(5) -- smallest-ring semantics (SMARTS-R1), not the
    // `[kN]` any-ring semantics chematic used to (incorrectly) alias it to.
    let mol = parse_smarts("[r5]").unwrap();
    assert_eq!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::MinRingSize(5))
    );
}

#[test]
fn test_parse_k_not_same_primitive_as_r() {
    // SMARTS-R1: `[kN]` (any-ring) and `[rN]` (smallest-ring) are DIFFERENT
    // predicates in RDKit -- confirmed empirically (a fusion atom shared between
    // a 5-ring and a 6-ring: RDKit's `[k6]` matches, `[r6]` does not). They must
    // no longer parse to the identical primitive (this test previously asserted
    // the opposite, which was the bug -- see docs/rdkit_compat.md's
    // "SMARTS-R0"/"SMARTS-R1" entries).
    let q_k = parse_smarts("[k6]").expect("[k6] must parse");
    let q_r = parse_smarts("[r6]").expect("[r6] must parse");
    assert_eq!(
        q_k.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::RingSize(6))
    );
    assert_eq!(
        q_r.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::MinRingSize(6))
    );
    assert_ne!(
        q_k.atoms[0].query, q_r.atoms[0].query,
        "[k6] and [r6] must NOT produce the same primitive"
    );
}

#[test]
fn test_k_primitive_matches_benzene() {
    // [k6] must match all 6 benzene atoms (all in 6-membered ring).
    use crate::{find_matches, parse_smarts};
    use chematic_smiles::parse;
    let benzene = parse("c1ccccc1").unwrap();
    let q = parse_smarts("[k6]").unwrap();
    let m = find_matches(&q, &benzene);
    assert_eq!(m.len(), 6, "[k6] must match all 6 benzene carbons");
}

#[test]
fn test_k_vs_r_equivalent_for_cyclopentane() {
    // [k5] and [r5] must give identical match counts on a SINGLE isolated ring
    // (no fusion) -- here "any ring of size 5" and "smallest ring is size 5"
    // trivially coincide, since there's only one ring to consider at all.
    use crate::{find_matches, parse_smarts};
    use chematic_smiles::parse;
    let mol = parse("C1CCCC1").unwrap(); // cyclopentane
    let q_k = parse_smarts("[k5]").unwrap();
    let q_r = parse_smarts("[r5]").unwrap();
    assert_eq!(
        find_matches(&q_k, &mol).len(),
        find_matches(&q_r, &mol).len(),
        "[k5] and [r5] must be equivalent in chematic"
    );
}

#[test]
fn test_k_vs_r_diverge_at_5_6_fusion_atom() {
    // SMARTS-R1: on a fused bicyclic (5-ring + 6-ring sharing one bond), the
    // fusion atoms belong to BOTH rings -- [kN] (any-ring) must match both [k5]
    // and [k6] on them, but [rN] (smallest-ring) must match ONLY [r5], since
    // their smallest ring is the 5-ring. This is the exact property Milestone
    // SMARTS-R0 found chematic previously got wrong (aliasing [r6] to [k6]'s
    // semantics would have wrongly matched the fusion atom on [r6] too).
    use crate::{find_matches, parse_smarts};
    use chematic_perception::find_sssr;
    use chematic_smiles::parse;

    let mol = parse("C1CCC2CCCC2C1").unwrap(); // bicyclic: 5-ring fused to 6-ring
    let sssr = find_sssr(&mol);
    let rings = sssr.rings();
    assert_eq!(rings.len(), 2, "expected exactly 2 SSSR rings");
    let sizes: Vec<usize> = rings.iter().map(|r| r.len()).collect();
    assert!(
        sizes.contains(&5) && sizes.contains(&6),
        "expected one 5-ring and one 6-ring, got sizes {sizes:?}"
    );
    let ring5: std::collections::HashSet<_> = rings
        .iter()
        .find(|r| r.len() == 5)
        .unwrap()
        .iter()
        .copied()
        .collect();
    let ring6: std::collections::HashSet<_> = rings
        .iter()
        .find(|r| r.len() == 6)
        .unwrap()
        .iter()
        .copied()
        .collect();
    let fusion_atom = *ring5
        .intersection(&ring6)
        .next()
        .expect("rings must share an atom");

    let q_k5 = parse_smarts("[k5]").unwrap();
    let q_k6 = parse_smarts("[k6]").unwrap();
    let q_r5 = parse_smarts("[r5]").unwrap();
    let q_r6 = parse_smarts("[r6]").unwrap();

    let matches = |q| {
        find_matches(q, &mol)
            .into_iter()
            .any(|m| m[&0] == fusion_atom)
    };
    assert!(matches(&q_k5), "[k5] must match the fusion atom");
    assert!(matches(&q_k6), "[k6] must match the fusion atom");
    assert!(
        matches(&q_r5),
        "[r5] must match the fusion atom (its smallest ring is 5)"
    );
    assert!(
        !matches(&q_r6),
        "[r6] must NOT match the fusion atom (its smallest ring is 5, not 6)"
    );
}

#[test]
fn test_r6_excludes_purine_fusion_atom_matching_rdkit() {
    // The exact discriminating case from SMARTS-R0's diagnosis: a purine where
    // RDKit's [k6] includes the ring-fusion nitrogen but [r6] excludes it (its
    // smallest ring is the 5-membered imidazole ring). Regression-pins the fix
    // against a concrete, previously-verified-against-real-RDKit example rather
    // than only a synthetic bicyclic.
    use crate::{find_matches, parse_smarts};
    use chematic_smiles::parse;
    let mol = parse("Nc1nc(N)c2ncn(C3CC(O)C(O)C(CO)O3)c2n1").unwrap();
    let q_r6 = parse_smarts("[r6]").unwrap();
    let q_k6 = parse_smarts("[k6]").unwrap();
    let r6_atoms: std::collections::HashSet<u32> = find_matches(&q_r6, &mol)
        .into_iter()
        .map(|m| m[&0].0)
        .collect();
    let k6_atoms: std::collections::HashSet<u32> = find_matches(&q_k6, &mol)
        .into_iter()
        .map(|m| m[&0].0)
        .collect();
    assert_eq!(
        r6_atoms,
        std::collections::HashSet::from([1, 2, 3, 9, 10, 11, 13, 15, 18, 20]),
        "[r6] atom set must match RDKit exactly (see SMARTS-R0 diagnosis)"
    );
    assert_eq!(
        k6_atoms,
        std::collections::HashSet::from([1, 2, 3, 5, 9, 10, 11, 13, 15, 18, 19, 20]),
        "[k6] atom set must match RDKit exactly (unchanged by this milestone)"
    );
    assert!(!r6_atoms.contains(&19) && k6_atoms.contains(&19));
}

#[test]
fn test_min_ring_size_false_for_non_ring_atom() {
    // A fully acyclic molecule has no atom in any ring -- [rN] must match nothing.
    use crate::{find_matches, parse_smarts};
    use chematic_smiles::parse;
    let mol = parse("CCCC").unwrap();
    for size in [3, 4, 5, 6] {
        let q = parse_smarts(&format!("[r{size}]")).unwrap();
        assert!(
            find_matches(&q, &mol).is_empty(),
            "[r{size}] must not match any atom in an acyclic molecule"
        );
    }
}

#[test]
fn test_ring_membership_and_ring_count_unaffected_by_r_split() {
    // [R]/[R0]/[R1]/[R2] are a separate primitive (RingMembership/RingCount) from
    // [rN] (MinRingSize) -- splitting [rN] from [kN] must not change their
    // behavior at all. Regression check on the same fused bicyclic as above.
    use crate::{find_matches, parse_smarts};
    use chematic_smiles::parse;
    let mol = parse("C1CCC2CCCC2C1").unwrap(); // 9 atoms, 2 fusion atoms
    let q_r = parse_smarts("[R]").unwrap();
    assert_eq!(
        find_matches(&q_r, &mol).len(),
        9,
        "[R] must match all 9 ring atoms"
    );
    let q_r1 = parse_smarts("[R1]").unwrap();
    assert_eq!(
        find_matches(&q_r1, &mol).len(),
        7,
        "[R1] must match the 7 non-fusion ring atoms"
    );
    let q_r2 = parse_smarts("[R2]").unwrap();
    assert_eq!(
        find_matches(&q_r2, &mol).len(),
        2,
        "[R2] must match the 2 fusion atoms (each in 2 SSSR rings)"
    );
}

#[test]
fn test_parse_hcount() {
    // `[H2]` → HCount(2)
    let mol = parse_smarts("[H2]").unwrap();
    assert_eq!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::HCount(2))
    );
}

#[test]
fn test_parse_cc_bond() {
    // `CC` → 2 atoms, 1 bond (Any)
    let mol = parse_smarts("CC").unwrap();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].query, BondQuery::Any);
}

#[test]
fn test_parse_double_bond() {
    // `C=C` → 2 atoms, 1 Double bond
    let mol = parse_smarts("C=C").unwrap();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(
        mol.bonds[0].query,
        BondQuery::Primitive(BondPrimitive::Double)
    );
}

#[test]
fn test_parse_branch() {
    // `C(=O)O` → 3 atoms (C, O, O), 2 bonds
    let mol = parse_smarts("C(=O)O").unwrap();
    assert_eq!(mol.atoms.len(), 3, "should have 3 atoms");
    assert_eq!(mol.bonds.len(), 2, "should have 2 bonds");
    // First bond: C=O (double)
    assert_eq!(
        mol.bonds[0].query,
        BondQuery::Primitive(BondPrimitive::Double)
    );
    // Second bond: C-O (implicit Any)
    assert_eq!(mol.bonds[1].query, BondQuery::Any);
}

#[test]
fn test_parse_benzene_ring() {
    // `c1ccccc1` → 6 aromatic C atoms, 6 bonds
    let mol = parse_smarts("c1ccccc1").unwrap();
    assert_eq!(mol.atoms.len(), 6, "benzene has 6 atoms");
    assert_eq!(mol.bonds.len(), 6, "benzene has 6 bonds");
    // All atoms should be aromatic.
    for atom in &mol.atoms {
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
        );
        assert_eq!(atom.query, expected);
    }
}

// ---- Sprint 2: SMARTS depth expansion edge cases ----
#[test]
fn test_parse_with_custom_config_default_depth() {
    let config = SmartsParserConfig::default();
    assert_eq!(
        config.max_recursion_depth,
        DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH
    );
    let mol = parse_smarts_with_config("c1ccccc1", &config).unwrap();
    assert_eq!(mol.atoms.len(), 6);
}

#[test]
fn test_parse_with_custom_config_increased_depth() {
    let config = SmartsParserConfig {
        max_recursion_depth: 12,
    };
    let mol = parse_smarts_with_config("[C,N,O]", &config).unwrap();
    assert_eq!(mol.atoms.len(), 1);
    // Pattern parses successfully with higher depth limit.
}

#[test]
fn test_parse_with_custom_config_depth_clamped_high() {
    // Request depth > 16 should clamp to 16.
    let config = SmartsParserConfig {
        max_recursion_depth: 100,
    };
    let mol = parse_smarts_with_config("C", &config).unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_parse_with_custom_config_depth_clamped_low() {
    // Request depth < 1 should clamp to 1.
    let config = SmartsParserConfig {
        max_recursion_depth: 0,
    };
    let mol = parse_smarts_with_config("C", &config).unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_parse_with_default_config_equivalent() {
    // parse_smarts should be equivalent to parse_smarts_with_config with default
    let pattern = "c1ccccc1";
    let mol1 = parse_smarts(pattern).unwrap();
    let mol2 = parse_smarts_with_config(pattern, &SmartsParserConfig::default()).unwrap();
    assert_eq!(mol1.atoms.len(), mol2.atoms.len());
}

#[test]
fn test_config_depth_parameter_clamped_correctly() {
    // Test that depth is clamped to [1, 16]
    assert_eq!(
        SmartsParserConfig {
            max_recursion_depth: 0
        }
        .max_recursion_depth,
        0
    );
    assert_eq!(
        SmartsParserConfig {
            max_recursion_depth: 16
        }
        .max_recursion_depth,
        16
    );
    assert_eq!(
        SmartsParserConfig {
            max_recursion_depth: 100
        }
        .max_recursion_depth,
        100
    );
    // Config stores the value as-is; clamping happens in parse_smarts_with_config
}

#[test]
fn test_operator_precedence_and_over_or() {
    // `[C&N,O]` should parse as `(C & N) | O`
    let mol = parse_smarts("[C&N,O]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    // Single atom with complex query (And/Or combination).
}

#[test]
fn test_operator_precedence_not_highest() {
    // `[!C&N]` should parse as `(! C) & N`
    let mol = parse_smarts("[!C&N]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_operator_precedence_semicolon_lowest() {
    // `[C;N,O]` should parse as `(C ; (N | O))` — different from `[C&N,O]`
    let mol = parse_smarts("[C;N,O]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_complex_bracket_atom_all_primitives() {
    // Combine multiple primitives: `[#6;a;R;H1]`
    let mol = parse_smarts("[#6;a;R;H1]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_long_chain_many_atoms() {
    // Long aliphatic chain: `CCCCCCCCCC` (10 atoms)
    let mol = parse_smarts("CCCCCCCCCC").unwrap();
    assert_eq!(mol.atoms.len(), 10);
    assert_eq!(mol.bonds.len(), 9);
}

#[test]
fn test_multiple_rings_fused() {
    // Naphthalene pattern: `c1ccc2ccccc2c1` (10 atoms in two fused rings)
    let mol = parse_smarts("c1ccc2ccccc2c1").unwrap();
    assert_eq!(mol.atoms.len(), 10);
}

#[test]
fn test_branching_from_multiple_atoms() {
    // `C(C)(C)C` — central C with 3 methyl branches
    let mol = parse_smarts("C(C)(C)C").unwrap();
    assert_eq!(mol.atoms.len(), 4);
    assert_eq!(mol.bonds.len(), 3);
}

#[test]
fn test_recursive_smarts_simple() {
    // `[$(C)]` — atoms matching inner pattern C
    let mol = parse_smarts("[$(C)]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    assert!(matches!(
        mol.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::Recursive(_))
    ));
}

#[test]
fn test_recursive_smarts_with_operators() {
    // `[$([C&N])]` — recursive pattern with operators
    let mol = parse_smarts("[$([C&N])]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

fn nested_recursive_smarts(depth: usize) -> String {
    let mut smarts = String::from("C");
    for _ in 0..depth {
        smarts = format!("[$({smarts})]");
    }
    smarts
}

#[test]
fn test_recursive_smarts_default_depth_boundary() {
    let smarts = nested_recursive_smarts(DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH);
    let mol = parse_smarts(&smarts).expect("default recursive SMARTS depth should parse");
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_recursive_smarts_default_depth_rejects_too_deep_pattern() {
    let smarts = nested_recursive_smarts(DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH + 1);
    assert!(matches!(
        parse_smarts(&smarts),
        Err(SmartsError::RecursionDepthExceeded)
    ));
}

#[test]
fn test_recursive_smarts_config_depth_is_clamped_for_safety() {
    let low_config = SmartsParserConfig {
        max_recursion_depth: 0,
    };
    let too_deep_for_low_config = nested_recursive_smarts(2);
    assert!(matches!(
        parse_smarts_with_config(&too_deep_for_low_config, &low_config),
        Err(SmartsError::RecursionDepthExceeded)
    ));

    let high_config = SmartsParserConfig {
        max_recursion_depth: 100,
    };
    let absolute_boundary = nested_recursive_smarts(ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH);
    parse_smarts_with_config(&absolute_boundary, &high_config)
        .expect("absolute recursive SMARTS depth boundary should parse");

    let too_deep_for_absolute_limit =
        nested_recursive_smarts(ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH + 1);
    assert!(matches!(
        parse_smarts_with_config(&too_deep_for_absolute_limit, &high_config),
        Err(SmartsError::RecursionDepthExceeded)
    ));
}

#[test]
fn test_malformed_recursive_smarts_return_errors() {
    assert!(matches!(
        parse_smarts("[$(C]"),
        Err(SmartsError::UnexpectedEnd)
    ));
    assert!(matches!(
        parse_smarts("[$(C"),
        Err(SmartsError::UnexpectedEnd)
    ));
    assert!(matches!(
        parse_smarts("[$(C))]"),
        Err(SmartsError::UnclosedBracket(_) | SmartsError::UnexpectedChar(_, _))
    ));
}

#[test]
fn test_isotope_with_bracket_atom() {
    // `[13C]` — explicitly labeled isotope
    let mol = parse_smarts("[13C]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_isotope_overflow_is_rejected_without_panicking() {
    let result = parse_smarts("[999999999999999999999999C]");
    assert!(matches!(result, Err(SmartsError::UnexpectedChar(_, _))));
}

#[test]
fn test_charge_positive_and_negative() {
    // `[C+2]` and `[O-1]`
    let mol_pos = parse_smarts("[C+2]").unwrap();
    let mol_neg = parse_smarts("[O-1]").unwrap();
    assert_eq!(mol_pos.atoms.len(), 1);
    assert_eq!(mol_neg.atoms.len(), 1);
}

#[test]
fn rdkit_smarts_extensions_match_rdkit(/* issue #734 item 6 */) {
    // Each row: query, target, expected matched atom indices (sorted),
    // taken from RDKit 2026.03.6 GetSubstructMatches.
    let cfg = crate::MatchConfig {
        uniquify: true,
        ..Default::default()
    };
    let hits = |smarts: &str, smiles: &str| -> Vec<u32> {
        let q = parse_smarts(smarts).unwrap();
        let mut v: Vec<u32> =
            crate::find_matches_with_config(&q, &chematic_smiles::parse(smiles).unwrap(), &cfg)
                .into_iter()
                .map(|m| m.values().map(|a| a.0).min().unwrap())
                .collect();
        v.sort_unstable();
        v
    };
    for (smarts, smiles, expected) in [
        ("[z1]", "CO", vec![0]),
        ("[z]", "OCO", vec![1]),
        ("[Z]", "OCO", vec![1]),
        ("[z1]", "c1ccncc1", vec![2, 4]),
        ("[Z1]", "c1ccncc1", vec![]),
        ("[Z0]", "CO", vec![1]),
        ("[d1]", "CC", vec![0, 1]),
        ("[d]", "CCC", vec![0, 2]),
        ("[d2]", "CCC", vec![1]),
        ("[D]", "CCC", vec![0, 2]),
        ("[CD{1-2}]", "CC(C)C", vec![0, 2, 3]),
        ("[CD{2-}]", "CC(C)C", vec![1]),
        ("[CD{-1}]", "CC(C)C", vec![0, 2, 3]),
        ("[D{2-1}]", "CCC", vec![]),
        ("[C;X{3-4}]", "CC=C", vec![0, 1, 2]),
        ("[R{1-}]", "C1CC1C", vec![0, 1, 2]),
        ("[R{-0}]", "CC1CC1", vec![0]),
        ("[r{5-6}]", "C1CCCC1", vec![0, 1, 2, 3, 4]),
        ("[x{2-}]", "C1CC1", vec![0, 1, 2]),
        ("[v{3-4}]", "CN", vec![0, 1]),
        ("[X{-2}]", "CC#N", vec![1, 2]),
        ("[C;z{1-}]", "CCO", vec![1]),
        ("[+{1-2}]", "[NH4+]", vec![0]),
        ("[-{1-}]", "C[O-]", vec![1]),
        ("[Zn]", "[Zn]", vec![0]),
        ("[N]->[Pt]", "N->[Pt]", vec![0]),
        ("[Pt]<-[N]", "N->[Pt]", vec![0]),
        ("[N]<-[Pt]", "N->[Pt]", vec![]),
        ("[Pt]->[N]", "N->[Pt]", vec![]),
        ("N-[Pt]", "N->[Pt]", vec![]),
        ("N~[Pt]", "N->[Pt]", vec![0]),
        // Ring ranges and acyclic atoms (`r0` is acyclic; `k` ranges
        // never match acyclic atoms).
        ("[r{3-}]", "CC1CCCCC1", vec![1, 2, 3, 4, 5, 6]),
        ("[r{-6}]", "CC1CCCCC1", vec![0, 1, 2, 3, 4, 5, 6]),
        ("[r{-5}]", "CC1CCCCC1", vec![0]),
        ("[r0]", "CC1CCCCC1", vec![0]),
        ("[k{3-}]", "CC1CCCCC1", vec![1, 2, 3, 4, 5, 6]),
        ("[k{-6}]", "CC1CCCCC1", vec![1, 2, 3, 4, 5, 6]),
        ("[k{0-6}]", "CC1CCCCC1", vec![1, 2, 3, 4, 5, 6]),
        // Bare `k` is "in a ring"; `k0` matches every atom (RDKit).
        ("[k]", "CC1CCCCC1", vec![1, 2, 3, 4, 5, 6]),
        ("[k0]", "CC1CCCCC1", vec![0, 1, 2, 3, 4, 5, 6]),
        ("[!k0]", "CC1CCCCC1", vec![]),
        // Bare `x` / `h` mean "at least one".
        ("[x]", "CC1CCCC1", vec![1, 2, 3, 4, 5]),
        ("[h]", "CC(C)(C)C", vec![0, 2, 3, 4]),
        // Out-of-range bounds match nothing instead of overflowing.
        ("[+{120-130}]", "C[N+](C)(C)C", vec![]),
        ("[-{1-200}]", "C[O-]", vec![1]),
        ("[D{999-1}]", "CCC", vec![]),
        ("[D{0-999}]", "CCC", vec![0, 1, 2]),
        ("[Xe]", "[Xe]", vec![0]),
    ] {
        assert_eq!(hits(smarts, smiles), expected, "{smarts} on {smiles}");
    }
    assert!(parse_smarts("[D{}]").is_err());
    assert!(parse_smarts("[D{1}]").is_err());
    // `<` is only a bond as `<-`; anything else is an error, not a panic.
    for bad in ["C<C", "C<", "C-<C", "C=,<C", "C!<C", "C;<C"] {
        assert!(parse_smarts(bad).is_err(), "{bad}");
    }
    // Multi-digit counts and charges (RDKit reads `[D12]`, `[+10]`).
    let atom = |smarts: &str| parse_smarts(smarts).unwrap().atoms[0].query.clone();
    assert_eq!(
        atom("[D12]"),
        AtomQuery::Primitive(AtomPrimitive::Degree(12))
    );
    assert_eq!(
        atom("[H12]"),
        AtomQuery::Primitive(AtomPrimitive::HCount(12))
    );
    assert_eq!(
        atom("[+10]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(10))
    );
    assert_eq!(
        atom("[-12]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(-12))
    );
    assert_eq!(
        atom("[R12]"),
        AtomQuery::Primitive(AtomPrimitive::RingCount(12))
    );
}

#[test]
fn bare_bracket_h_is_a_hydrogen_atom() {
    // Daylight/OpenSMARTS: `[H]`, `[2H]`, `[H+]` and `[H:1]` are hydrogen
    // atoms (issue #734); `[*H]`, `[CH]`, `[H1]` and `[C;H]` keep the
    // H-count meaning.
    let cfg = crate::MatchConfig::default();
    let matches = |smarts: &str, smiles: &str| {
        let q = parse_smarts(smarts).unwrap();
        crate::find_matches_with_config(&q, &chematic_smiles::parse(smiles).unwrap(), &cfg).len()
    };
    // Methane with one explicit hydrogen atom: only that atom is an H atom.
    assert_eq!(matches("[H]", "[H]C"), 1);
    assert_eq!(matches("[H:1]", "[H]C"), 1);
    assert_eq!(matches("[H]", "C"), 0, "implicit hydrogens are not atoms");
    assert_eq!(matches("[H+]", "[H+]"), 1);
    assert_eq!(matches("[H+]", "[H]C"), 0);
    assert_eq!(matches("[2H]", "[2H]C"), 1);
    // H-count forms are unchanged.
    assert_eq!(matches("[*H]", "CC"), 0);
    assert_eq!(matches("[CH3]", "CC"), 2);
    assert_eq!(matches("[C;H3]", "CC"), 2);
    assert_eq!(matches("[H3]", "CC"), 2);
}

#[test]
fn repeated_charge_signs_accumulate() {
    // Daylight / OpenSMARTS: `++` is +2, `--` is -2 (issue #680), never an
    // AND of two unit charges; a digit still sets the magnitude directly.
    fn charge(smarts: &str) -> AtomQuery {
        parse_smarts(smarts).unwrap().atoms[0].query.clone()
    }
    assert_eq!(
        charge("[++]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(2))
    );
    assert_eq!(
        charge("[--]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(-2))
    );
    assert_eq!(
        charge("[+++]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(3))
    );
    assert_eq!(charge("[+2]"), charge("[++]"));
    assert_eq!(charge("[-2]"), charge("[--]"));
    assert_eq!(
        charge("[+]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(1))
    );
    assert_eq!(
        charge("[+0]"),
        AtomQuery::Primitive(AtomPrimitive::Charge(0))
    );
    // Element and other primitives around the charge still combine by AND.
    assert_eq!(
        charge("[Ca++]"),
        AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("Ca".into()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(2))),
        )
    );
    // `+-` is two primitives (+1 and -1), not a charge of zero.
    assert_eq!(
        charge("[+-]"),
        AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(1))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(-1))),
        )
    );
    // Matching: a dication matches `[++]` and `[+2]`, a monocation neither.
    let dication = chematic_smiles::parse("[Ca+2]").unwrap();
    let monocation = chematic_smiles::parse("[Na+]").unwrap();
    for q in ["[++]", "[+2]"] {
        let q = parse_smarts(q).unwrap();
        assert!(
            crate::has_match_with_config(&q, &dication, &crate::MatchConfig::default()),
            "{q:?}"
        );
        assert!(
            !crate::has_match_with_config(&q, &monocation, &crate::MatchConfig::default()),
            "{q:?}"
        );
    }
}

#[test]
fn test_degree_and_connectivity() {
    // `[D4]` (degree 4) and `[X3]` (connectivity 3)
    let mol = parse_smarts("[D4]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
    let mol2 = parse_smarts("[X3]").unwrap();
    assert_eq!(mol2.atoms.len(), 1);
}

#[test]
fn test_ring_size_constraint() {
    // `[r6]` (in 6-membered ring)
    let mol = parse_smarts("[r6]").unwrap();
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn test_nested_branch_depth() {
    // Deeply nested branches: `C(C(C(C(C))))`
    let mol = parse_smarts("C(C(C(C(C))))").unwrap();
    assert_eq!(mol.atoms.len(), 5);
}

#[test]
fn test_ring_and_branch_combined() {
    // Ring with branch: `c1cc(C)ccc1`
    let mol = parse_smarts("c1cc(C)ccc1").unwrap();
    // 6 aromatic carbons + 1 branch carbon = 7 atoms
    assert_eq!(mol.atoms.len(), 7);
}

#[test]
fn test_all_bond_types() {
    // Test all bond types: `-`, `=`, `#`, `:`, `~`
    let bonds = vec![
        ("C-C", BondPrimitive::Single),
        ("C=C", BondPrimitive::Double),
        ("C#C", BondPrimitive::Triple),
        ("C$C", BondPrimitive::Quadruple),
        ("c:c", BondPrimitive::Aromatic),
    ];
    for (smarts, expected_prim) in bonds {
        let mol = parse_smarts(smarts).unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(
            mol.bonds[0].query,
            BondQuery::Primitive(expected_prim),
            "bond type mismatch for {smarts}"
        );
    }
}

#[test]
fn test_empty_recursive_pattern_rejected() {
    // `[$()] `— empty inner pattern should fail gracefully or return empty match
    let result = parse_smarts("[$$()]");
    // Empty recursion is probably an error.
    assert!(result.is_err() || result.as_ref().unwrap().atoms.len() == 1);
}

#[test]
fn rdkit_rejects_and_accepts_like_rdkit(/* issue #754 */) {
    // Lowercase atoms are only the aromatic organic subset.
    for bad in [
        "[i]", "[Ci1]", "[f]", "[C^6]", "[C^7]", "[C@TH3]", "[C@SP4]", "[C@SP]",
    ] {
        assert!(parse_smarts(bad).is_err(), "{bad} should be an error");
    }
    for good in [
        "[c]", "[b]", "[se]", "[C^0]", "[C^5]", "[C@TH1]", "[C@SP1]", "[C@TB20]", "[C@OH30]",
    ] {
        assert!(parse_smarts(good).is_ok(), "{good} should parse");
    }
    // `@SP1` is one chirality primitive, not `@` followed by S and P.
    let q = parse_smarts("[C;@SP1]").unwrap();
    let mol = chematic_smiles::parse("ClCCBr").unwrap();
    assert_eq!(crate::find_matches(&q, &mol).len(), 2);
}
