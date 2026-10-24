use super::*;
use crate::parser::parse_smarts;
use chematic_smiles::parse;

#[test]
fn query_boolean_and_bond_primitive_boundaries_match_without_panics() {
    let targets = [
        parse("[13CH3][C@H](F)Cl").unwrap(),
        parse("c1ncccc1").unwrap(),
        parse("C#N").unwrap(),
        parse("[NH4+]").unwrap(),
        parse("C1CC1").unwrap(),
    ];
    let patterns = [
        "[13#6]",
        "[C,N;!H0]",
        "[a,n]",
        "[A;D1]",
        "[v4]",
        "[R]",
        "[R1]",
        "[r3]",
        "[x2]",
        "[z0]",
        "[Z0]",
        "[d1]",
        "[+1]",
        "[^3]",
        "[C@H]",
        "[C,C@]",
        "[N,C@]",
        "[!C@]",
        "[$(C)]",
        "[$(CC)]",
        "[$(N)]",
        "C-C",
        "C=C",
        "C#N",
        "c:c",
        "C~N",
        "C@C",
        "C!@C",
        "C/C",
        "C\\C",
    ];
    for pattern in patterns {
        let query = parse_smarts(pattern).unwrap();
        for target in &targets {
            let _ = has_match_rdkit_parity_bounded(&query, target, &RdkitParityConfig::default())
                .unwrap();
        }
    }

    for (pattern, target) in [("[N]->[Cu]", "N->[Cu]"), ("[Cu]<-[N]", "[Cu]<-N")] {
        let query = parse_smarts(pattern).unwrap();
        let molecule = parse(target).unwrap();
        assert_eq!(
            has_match_rdkit_parity_bounded(&query, &molecule, &RdkitParityConfig::default())
                .unwrap(),
            MatchOutcome::Found
        );
    }

    let isotope_query = parse_smarts("[13#6]").unwrap();
    let isotope_target = parse("C").unwrap();
    let mut isotope_agnostic = RdkitParityConfig::default();
    isotope_agnostic.base.use_isotopes = false;
    assert_eq!(
        has_match_rdkit_parity_bounded(&isotope_query, &isotope_target, &isotope_agnostic).unwrap(),
        MatchOutcome::Found
    );
}

#[test]
fn candidate_tristate_and_keyed_cache_boundaries_are_explicit() {
    let molecule = parse("C").unwrap();
    let rings = chematic_perception::find_sssr(&molecule);
    let cache = TargetCache::default();
    let config = MatchConfig {
        use_chirality: true,
        ..MatchConfig::default()
    };
    let ctx = EvalCtx {
        mol: &molecule,
        rings: &rings,
        symmetrized_rings: None,
        size_model: None,
        ring_model: None,
        config: &config,
        visit_budget: std::cell::Cell::new(u64::MAX),
        budget_exhausted: std::cell::Cell::new(false),
        cache: &cache,
        state_pool: std::cell::RefCell::new(Vec::new()),
    };
    let carbon = AtomQuery::Primitive(AtomPrimitive::AtomicNum(6));
    let nitrogen = AtomQuery::Primitive(AtomPrimitive::AtomicNum(7));
    let chirality = AtomQuery::Primitive(AtomPrimitive::Chirality(1));
    let atom = AtomIdx(0);

    assert_eq!(eval_atom_query_candidate(&chirality, atom, &ctx), None);
    assert_eq!(eval_atom_query_candidate(&carbon, atom, &ctx), Some(true));
    assert_eq!(
        eval_atom_query_candidate(&nitrogen, atom, &ctx),
        Some(false)
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::And(Box::new(carbon.clone()), Box::new(carbon.clone())),
            atom,
            &ctx,
        ),
        Some(true)
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::And(Box::new(carbon.clone()), Box::new(nitrogen.clone())),
            atom,
            &ctx,
        ),
        Some(false)
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::And(Box::new(carbon.clone()), Box::new(chirality.clone())),
            atom,
            &ctx,
        ),
        None
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::Or(Box::new(carbon.clone()), Box::new(nitrogen.clone())),
            atom,
            &ctx,
        ),
        Some(true)
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::Or(Box::new(nitrogen.clone()), Box::new(nitrogen.clone())),
            atom,
            &ctx,
        ),
        Some(false)
    );
    assert_eq!(
        eval_atom_query_candidate(
            &AtomQuery::Or(Box::new(nitrogen.clone()), Box::new(chirality.clone())),
            atom,
            &ctx,
        ),
        None
    );
    assert_eq!(
        eval_atom_query_candidate(&AtomQuery::Not(Box::new(chirality)), atom, &ctx),
        None
    );

    let keyed_cell = std::sync::OnceLock::new();
    let mut local = None;
    assert_eq!(*keyed(&keyed_cell, 1_u8, &mut local, || 7_u8), 7);
    let mut other_key = None;
    assert_eq!(*keyed(&keyed_cell, 2_u8, &mut other_key, || 9_u8), 9);
    assert_eq!(other_key, Some(9));

    assert!(!eval_chirality(atom, &ctx, 1));
    let no_chirality = MatchConfig::default();
    let no_chirality_ctx = EvalCtx {
        config: &no_chirality,
        ..ctx
    };
    assert!(eval_chirality(atom, &no_chirality_ctx, 1));
}
