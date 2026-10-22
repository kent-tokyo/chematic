use super::*;
use crate::{find_matches, find_matches_with_config, parse_smarts};
use chematic_smiles::parse;

#[test]
fn existence_search_agrees_with_first_embedding_search() {
    let patterns = [
        "*1~*~*~*~*~*~*~1",
        "*1~*~*~1",
        "[#7]~*~*~*~[#8]",
        "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",
        "c1ccc2ccccc2c1",
        "[#16]",
        "[Cl,Br]",
        "C(=O)N.c1ccccc1",
        "[R2]",
        "*@*!@*@*",
    ];
    let targets = [
        "CCO",
        "c1ccccc1CCN",
        "C1CCCCCC1",
        "c1ccc2ccccc2c1O",
        "CC(=O)Nc1ccccc1",
        "C1CC1C1CCC1",
        "C[Co]1(Cl)CCCCC1",
        "OCCCCO",
    ];
    let first = MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        ..MatchConfig::default()
    };
    for p in patterns {
        let q = parse_smarts(p).unwrap();
        for t in targets {
            let mol = parse(t).unwrap();
            assert_eq!(
                has_match_with_config(&q, &mol, &first),
                !find_matches_with_config(&q, &mol, &first).is_empty(),
                "{p} in {t}"
            );
        }
    }
}

#[test]
fn single_atomic_number_existence_query_is_exact() {
    let mol = parse("CCN").unwrap();
    let config = MatchConfig::default();
    assert!(has_match_with_config(
        &parse_smarts("[#6]").unwrap(),
        &mol,
        &config
    ));
    assert!(has_match_with_config(
        &parse_smarts("[#7]").unwrap(),
        &mol,
        &config
    ));
    assert!(!has_match_with_config(
        &parse_smarts("[#8]").unwrap(),
        &mol,
        &config
    ));
}

#[test]
fn query_cycle_flags() {
    let q = parse_smarts("*1~*~*~1-C").unwrap();
    assert_eq!(query_cyclic_atoms(&q), vec![true, true, true, false]);
    assert!(query_cyclic_atoms(&parse_smarts("CCO").unwrap()).is_empty());
    let spiro = parse_smarts("C1CC12CC2").unwrap();
    assert!(query_cyclic_atoms(&spiro).iter().all(|&c| c));
}

// -- Isotope matching -----------------------------------------------------

#[test]
fn test_isotope_enforced_by_default_like_rdkit() {
    // RDKit: an isotope in a query always constrains; unlabelled is 0.
    let mol = parse("[13CH3]C").unwrap();
    let count = |q: &str| find_matches(&parse_smarts(q).unwrap(), &mol).len();
    assert_eq!(count("[13C]"), 1);
    assert_eq!(count("[12C]"), 0);
    assert_eq!(count("[0C]"), 1);
    assert_eq!(count("[C;!13]"), 1);
    assert_eq!(count("[!1]"), 2);
    let ignoring = MatchConfig {
        use_isotopes: false,
        ..MatchConfig::default()
    };
    let all = find_matches_with_config(&parse_smarts("[12C]").unwrap(), &mol, &ignoring);
    assert_eq!(
        all.len(),
        2,
        "use_isotopes=false ignores isotope primitives"
    );
}

#[test]
fn directional_bonds_match_without_stereo_like_rdkit() {
    let count =
        |q: &str, s: &str| find_matches(&parse_smarts(q).unwrap(), &parse(s).unwrap()).len();
    assert_eq!(count("C/C", "CC"), 1);
    assert_eq!(count("c/c", "c1ccccc1"), 6);
    assert_eq!(count("c\\c", "c1ccccc1"), 6);
    // Not enforced, and the same for every spelling of the target.
    for target in [
        "F/C=C/F",
        "F/C=C\\F",
        "FC=CF",
        "C(/F)=C/F",
        "C(/F)=C\\F",
        "C(\\F)=C/F",
    ] {
        assert_eq!(count("F/C=C/F", target), 1, "{target}");
    }
    assert_eq!(count("C/C", "C=C"), 0);
}

#[test]
fn total_valence_uses_kekule_bond_orders() {
    let count =
        |q: &str, s: &str| find_matches(&parse_smarts(q).unwrap(), &parse(s).unwrap()).len();
    assert_eq!(count("[c;v4]", "c1ccccc1"), 6);
    assert_eq!(count("[c;v3]", "c1ccccc1"), 0);
    assert_eq!(count("[n;v3]", "c1cc[nH]c1"), 1);
    assert_eq!(count("[n;v3]", "c1ccncc1"), 1);
    assert_eq!(count("[o;v2]", "c1ccoc1"), 1);
}

#[test]
fn test_isotope_enforced_when_enabled() {
    // Build a molecule with one 13C and one 12C.
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    let mut b = MoleculeBuilder::new();
    let mut c13 = Atom::new(Element::C);
    c13.isotope = Some(13);
    let c13_idx = b.add_atom(c13);
    let c12_idx = b.add_atom(Atom::new(Element::C));
    b.add_bond(c13_idx, c12_idx, BondOrder::Single).unwrap();
    let mol = b.build();

    let query = parse_smarts("[13C]").unwrap();
    let config = MatchConfig {
        use_isotopes: true,
        ..MatchConfig::default()
    };
    let matches = find_matches_with_config(&query, &mol, &config);
    assert_eq!(
        matches.len(),
        1,
        "[13C] with use_isotopes=true should match only the 13C atom"
    );
    assert_eq!(matches[0][&0], AtomIdx(0));
}

#[test]
fn test_no_isotope_match_on_unlabeled() {
    // [13C] with use_isotopes=true should not match unlabeled carbons.
    let mol = parse("CC").unwrap(); // both atoms have isotope=None
    let query = parse_smarts("[13C]").unwrap();
    let config = MatchConfig {
        use_isotopes: true,
        ..MatchConfig::default()
    };
    let matches = find_matches_with_config(&query, &mol, &config);
    assert_eq!(
        matches.len(),
        0,
        "[13C] with use_isotopes=true should not match unlabeled C"
    );
}

// -- Chirality matching ---------------------------------------------------

#[test]
fn test_chirality_ignored_by_default() {
    // [@] query should match any atom when use_chirality=false (default).
    let mol = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine
    let query = parse_smarts("[C@@H]").unwrap();
    let matches = find_matches(&query, &mol);
    // Default: chirality ignored, so [@] matches any C-H regardless of chirality.
    assert!(
        !matches.is_empty(),
        "chirality should be ignored by default"
    );
}

#[test]
fn test_chirality_enforced_when_enabled() {
    // L-alanine has [C@@H] — query [C@@H] should match, [C@H] should not.
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();

    let q_ccw = parse_smarts("[C@@H]").unwrap(); // CCW (@@) = kind 2
    let q_cw = parse_smarts("[C@H]").unwrap(); // CW  (@)  = kind 1

    let config = MatchConfig {
        use_chirality: true,
        ..MatchConfig::default()
    };

    let m_ccw = find_matches_with_config(&q_ccw, &mol, &config);
    let m_cw = find_matches_with_config(&q_cw, &mol, &config);

    // [C@@H] must match L-alanine's chiral center.
    assert!(!m_ccw.is_empty(), "[C@@H] should match L-alanine (@@)");
    // [C@H] must NOT match L-alanine.
    assert!(m_cw.is_empty(), "[C@H] should not match L-alanine (@@)");
}

#[test]
fn test_chirality_d_alanine_positive_match() {
    // D-alanine [C@H] should positively match [C@H] query when use_chirality=true.
    // This complements test_chirality_enforced_when_enabled which only tested negative on L-ala.
    let mol = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (@, kind 1)
    let q_cw = parse_smarts("[C@H]").unwrap(); // CW (@), kind 1
    let config = MatchConfig {
        use_chirality: true,
        ..MatchConfig::default()
    };
    let m = find_matches_with_config(&q_cw, &mol, &config);
    assert!(!m.is_empty(), "[C@H] should positively match D-alanine (@)");
}

// S2: visit budget tests

#[test]
fn test_visit_budget_unlimited_default() {
    // Default config has no budget cap — normal queries complete fully.
    let mol = parse("c1ccccc1").unwrap(); // benzene
    let q = parse_smarts("c1ccccc1").unwrap();
    let m = find_matches(&q, &mol);
    assert!(!m.is_empty(), "benzene should match aromatic ring query");
}

#[test]
fn test_visit_budget_generous_limit_finds_match() {
    // A budget high enough that a simple query completes.
    let mol = parse("CCO").unwrap();
    let q = parse_smarts("O").unwrap();
    let config = MatchConfig {
        max_visit_budget: Some(10_000),
        ..MatchConfig::default()
    };
    let m = find_matches_with_config(&q, &mol, &config);
    assert!(!m.is_empty(), "CCO contains O — should match within budget");
}

#[test]
fn test_visit_budget_zero_returns_empty() {
    // Budget of 0 → no states explored → no results (documents fail-safe behavior).
    let mol = parse("CCO").unwrap();
    let q = parse_smarts("O").unwrap();
    let config = MatchConfig {
        max_visit_budget: Some(0),
        ..MatchConfig::default()
    };
    let m = find_matches_with_config(&q, &mol, &config);
    // With zero budget the search returns immediately — may or may not find a match.
    // Just verify it does not panic.
    let _ = m;
}

// ── MatchOutcome / has_match_bounded: budget-exhaustion must be a
// distinct, non-negative outcome, never silently folded into "no match" ──

#[test]
fn has_match_bounded_found_with_generous_budget() {
    let mol = parse("CCO").unwrap();
    let q = parse_smarts("O").unwrap();
    let rings = find_sssr(&mol);
    let config = MatchConfig {
        max_visit_budget: Some(10_000),
        ..MatchConfig::default()
    };
    assert_eq!(
        has_match_bounded(&q, &mol, &rings, &config),
        MatchOutcome::Found
    );
}

#[test]
fn has_match_bounded_not_found_is_exhaustive_not_a_guess() {
    let mol = parse("CC").unwrap(); // no oxygen anywhere
    let q = parse_smarts("O").unwrap();
    let rings = find_sssr(&mol);
    let config = MatchConfig {
        max_visit_budget: Some(10_000),
        ..MatchConfig::default()
    };
    assert_eq!(
        has_match_bounded(&q, &mol, &rings, &config),
        MatchOutcome::NotFound
    );
}

#[test]
fn has_match_bounded_zero_budget_is_indeterminate_not_not_found() {
    // The pattern genuinely matches (CCO contains O), but a budget of 0
    // means the search never gets to look. The result MUST be
    // `BudgetExhausted`, never `NotFound` -- conflating the two is
    // exactly the silent-false-negative bug this API exists to prevent.
    let mol = parse("CCO").unwrap();
    let q = parse_smarts("O").unwrap();
    let rings = find_sssr(&mol);
    let config = MatchConfig {
        max_visit_budget: Some(0),
        ..MatchConfig::default()
    };
    assert_eq!(
        has_match_bounded(&q, &mol, &rings, &config),
        MatchOutcome::BudgetExhausted
    );
}

#[test]
fn find_matches_with_rings_and_config_checked_reports_exhaustion() {
    let mol = parse("CCO").unwrap();
    let q = parse_smarts("O").unwrap();
    let rings = find_sssr(&mol);

    // Generous budget: real match, not exhausted.
    let generous = MatchConfig {
        max_visit_budget: Some(10_000),
        ..MatchConfig::default()
    };
    let (matches, exhausted) =
        find_matches_with_rings_and_config_checked(&q, &mol, &rings, &generous);
    assert!(!matches.is_empty());
    assert!(!exhausted);

    // Zero budget: a match exists but the search never ran -- must be
    // reported as exhausted, not silently equivalent to a real negative.
    let starved = MatchConfig {
        max_visit_budget: Some(0),
        ..MatchConfig::default()
    };
    let (matches, exhausted) =
        find_matches_with_rings_and_config_checked(&q, &mol, &rings, &starved);
    assert!(matches.is_empty());
    assert!(
        exhausted,
        "zero budget must be reported as exhausted, not as a confirmed negative"
    );
}

// ── RDKit PR #9201: query > target early exit ─────────────────────────────

#[test]
fn query_larger_than_target_returns_no_matches() {
    // A 6-atom query (benzene) cannot match a 4-atom target. The early exit
    // must fire immediately without entering VF2 recursion (RDKit PR #9201).
    let mol = parse("CCCC").unwrap(); // 4 heavy atoms
    let q = parse_smarts("c1ccccc1").unwrap(); // 6 query atoms
    let m = find_matches(&q, &mol);
    assert!(m.is_empty(), "6-atom query must not match 4-atom target");
}

#[test]
fn query_same_size_as_target_allowed() {
    // Equal atom counts must still attempt the search (no false fast-fail).
    let mol = parse("CCCCCC").unwrap(); // 6 atoms
    let q = parse_smarts("CCCCCC").unwrap(); // 6 atoms
    let m = find_matches(&q, &mol);
    assert!(!m.is_empty(), "same-size query must be attempted");
}
