use crate::{
    MatchConfig, MatchOutcome, RdkitParityConfig, find_match_atom_sets_rdkit_parity,
    has_match_rdkit_parity_bounded, parse_smarts,
};
#[test]
fn parity_matcher_recursive_and_atom_primitives_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-smarts-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let smarts = row["smarts"].as_str().unwrap();
        let mol = if source.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            chematic_smiles::parse(source).unwrap()
        };
        let query = parse_smarts(smarts).unwrap();
        let expected: Vec<Vec<u32>> = serde_json::from_value(row["sets"].clone()).unwrap();
        for shared in [false, true] {
            let config = RdkitParityConfig {
                base: MatchConfig {
                    use_chirality: row["use_chirality"].as_bool().unwrap(),
                    ..Default::default()
                },
                use_shared_symmetrized_sssr: shared,
                ..Default::default()
            };
            let (actual, exhausted) =
                find_match_atom_sets_rdkit_parity(&query, &mol, &config).unwrap();
            assert!(!exhausted);
            let exists = has_match_rdkit_parity_bounded(&query, &mol, &config).unwrap();
            assert_eq!(
                exists,
                if expected.is_empty() {
                    MatchOutcome::NotFound
                } else {
                    MatchOutcome::Found
                },
                "{source}: {smarts}"
            );
            if actual != expected {
                failures.push(format!(
                    "{source}: {smarts} shared={shared}: {actual:?} != {expected:?}"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn parity_matcher_budget_exhaustion_is_distinct_from_no_match() {
    let mol = chematic_smiles::parse("c1ccccc1").unwrap();
    let query = parse_smarts("ccc").unwrap();
    let config = RdkitParityConfig {
        base: MatchConfig {
            max_visit_budget: Some(0),
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        has_match_rdkit_parity_bounded(&query, &mol, &config).unwrap(),
        MatchOutcome::BudgetExhausted
    );
    let query = parse_smarts("[13C]").unwrap();
    let mol = chematic_smiles::parse("CC").unwrap();
    assert_eq!(
        has_match_rdkit_parity_bounded(&query, &mol, &Default::default()).unwrap(),
        MatchOutcome::NotFound
    );
    let config = RdkitParityConfig {
        base: MatchConfig {
            use_isotopes: false,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        has_match_rdkit_parity_bounded(&query, &mol, &config).unwrap(),
        MatchOutcome::Found
    );
}

#[test]
fn cached_ring_matching_survives_eviction_and_applies_per_call_limits() {
    let mol = chematic_smiles::parse("c1ccncc1").unwrap();
    let rings = chematic_perception::find_sssr(&mol);
    let mut cache = crate::SmartsCache::new(2);
    for round in 0..40 {
        for (query, expected) in [("c", 5), ("n", 1), ("[R0]", 0), ("[r6]", 6)] {
            assert_eq!(
                cache
                    .find_matches_with_rings(query, &mol, &rings)
                    .unwrap()
                    .len(),
                expected
            );
            let config = MatchConfig {
                max_matches: Some(1),
                ..Default::default()
            };
            assert_eq!(
                cache
                    .find_matches_with_config(query, &mol, &config)
                    .unwrap()
                    .len(),
                expected.min(1),
                "round {round}: {query}"
            );
            assert!(cache.len() <= 2);
        }
        assert!(
            cache
                .find_matches_with_config("[", &mol, &Default::default())
                .is_err()
        );
    }
}
