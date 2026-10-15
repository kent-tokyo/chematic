use chematic_smarts::{
    MatchConfig, QueryHydrogenError, find_matches_with_config, merge_query_hydrogens, parse_smarts,
};
use chematic_smiles::parse;

fn matches(pattern: &str, target: &str, chirality: bool) -> bool {
    let query = diagnostic(merge_query_hydrogens(&parse_smarts(pattern).unwrap())).unwrap();
    !find_matches_with_config(
        &query,
        &parse(target).unwrap(),
        &MatchConfig {
            use_chirality: chirality,
            ..Default::default()
        },
    )
    .is_empty()
}

#[test]
fn issue9675_matches_intended_enantiomer_after_h_folding() {
    let q = "[C@@]([H])(C)(N)C(=O)O";
    assert!(matches(q, "C[C@H](N)C(=O)O", true));
    assert!(!matches(q, "C[C@@H](N)C(=O)O", true));
    assert!(matches(q, "C[C@@H](N)C(=O)O", false));
}

#[test]
fn virtual_h_parity_tracks_every_neighbor_position() {
    let equivalent = [
        "[H][C@@](C)(N)C(=O)O",
        "[C@@]([H])(C)(N)C(=O)O",
        "[C@](C)([H])(N)C(=O)O",
        "[C@@](C)(N)([H])C(=O)O",
        "[C@](C)(N)(C(=O)O)[H]",
        "C[C@]([H])(N)C(=O)O",
        "C[C@@](N)([H])C(=O)O",
    ];
    for q in equivalent {
        assert!(matches(q, "C[C@H](N)C(=O)O", true), "{q}");
        assert!(!matches(q, "C[C@@H](N)C(=O)O", true), "{q}");
    }
}

#[test]
fn explicit_h_targets_keep_same_stereo_result() {
    for target in [
        "[H][C@@](C)(N)C(=O)O",
        "C[C@]([H])(N)C(=O)O",
        "[2H][C@@](C)(N)C(=O)O",
    ] {
        assert!(matches("[C@@]([H])(C)(N)C(=O)O", target, true), "{target}");
    }
    assert!(!matches(
        "[C@@]([H])(C)(N)C(=O)O",
        "[H][C@](C)(N)C(=O)O",
        true
    ));
}

#[test]
fn h_count_is_minimum_not_exact() {
    for target in ["C", "CC", "CCC", "[CH2]", "[CH3]"] {
        assert!(matches("C([H])[H]", target, false), "{target}");
    }
    for target in ["[C]", "[CH]", "C(F)(F)(F)F"] {
        assert!(!matches("C([H])[H]", target, false), "{target}");
    }
    assert!(!matches("[CH0][H]", "C", false));
}

#[test]
fn isotopic_mapped_charged_and_complex_h_are_retained() {
    for q in [
        "[2H]C", "[H:7]C", "[H+]C", "[#1;D1]C", "[H]~C", "[H]=C", "[H][H]", "[H][H+]", "[H]",
    ] {
        let original = parse_smarts(q).unwrap();
        assert_eq!(
            diagnostic(merge_query_hydrogens(&original)).unwrap(),
            original,
            "{q}"
        );
    }
    assert_eq!(
        diagnostic(merge_query_hydrogens(&parse_smarts("[#1]C").unwrap()))
            .unwrap()
            .atom_count(),
        1
    );
}

#[test]
fn output_is_idempotent_and_input_is_immutable() {
    let input = parse_smarts("[H][C@](C)(N)C(=O)O").unwrap();
    let copy = input.clone();
    let output = diagnostic(merge_query_hydrogens(&input)).unwrap();
    assert_eq!(input, copy);
    assert_eq!(output.atom_count(), input.atom_count() - 1);
    assert_eq!(diagnostic(merge_query_hydrogens(&output)).unwrap(), output);
}

#[test]
fn unsafe_partial_frames_recursive_queries_and_bad_graphs_fail() {
    for q in ["[C@]([H])(N)C", "[C@]([H])([H])(N)C"] {
        assert!(matches!(
            diagnostic(merge_query_hydrogens(&parse_smarts(q).unwrap())),
            Err(QueryHydrogenError::UnsupportedStereo { .. })
        ));
    }
    assert_eq!(
        diagnostic(merge_query_hydrogens(&parse_smarts("[$(C[H])]").unwrap())),
        Err(QueryHydrogenError::RecursiveQuery)
    );
    let mut q = parse_smarts("C[H]").unwrap();
    q.adj[0][0].1 = 99;
    assert_eq!(
        diagnostic(merge_query_hydrogens(&q)),
        Err(QueryHydrogenError::InvalidGraph)
    );
}

#[test]
fn ring_openings_reserve_order_for_both_closure_orders() {
    for (q, same, opposite) in [
        (
            "[C@]12([H])CC(C)CC1NC2",
            "[C@H]12CC(C)CC1NC2",
            "[C@@H]12CC(C)CC1NC2",
        ),
        (
            "[C@]12([H])CC(C)CC2NC1",
            "[C@H]12CC(C)CC2NC1",
            "[C@@H]12CC(C)CC2NC1",
        ),
        (
            "[C@]1([H])2CC(C)CC1NC2",
            "[C@@H]12CC(C)CC1NC2",
            "[C@H]12CC(C)CC1NC2",
        ),
    ] {
        assert!(matches(q, same, true), "{q}");
        assert!(!matches(q, opposite, true), "{q}");
    }
}

#[test]
fn limits_are_checked_before_cloning_invalid_public_queries() {
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondQuery, QueryMolecule};
    let mut graph = QueryMolecule::new();
    for _ in 0..1025 {
        graph.add_atom(AtomQuery::Primitive(AtomPrimitive::Wildcard));
    }
    assert_eq!(
        diagnostic(merge_query_hydrogens(&graph)),
        Err(QueryHydrogenError::ResourceLimit)
    );
    let mut expr = AtomQuery::Primitive(AtomPrimitive::Wildcard);
    for _ in 0..33 {
        expr = AtomQuery::Not(Box::new(expr));
    }
    let mut graph = QueryMolecule::new();
    graph.add_atom(expr);
    assert_eq!(
        diagnostic(merge_query_hydrogens(&graph)),
        Err(QueryHydrogenError::ResourceLimit)
    );
    let mut graph = parse_smarts("CC").unwrap();
    let mut expr = BondQuery::Any;
    for _ in 0..33 {
        expr = BondQuery::Not(Box::new(expr));
    }
    graph.bonds[0].query = expr;
    assert_eq!(
        diagnostic(merge_query_hydrogens(&graph)),
        Err(QueryHydrogenError::ResourceLimit)
    );
}

#[test]
fn surviving_atom_maps_and_non_h_bond_queries_are_preserved() {
    let input = parse_smarts("[C:9]([H])~[N:2]").unwrap();
    let output = diagnostic(merge_query_hydrogens(&input)).unwrap();
    assert_eq!(
        output.atoms.iter().map(|a| a.atom_map).collect::<Vec<_>>(),
        [Some(9), Some(2)]
    );
    assert_eq!(output.bonds[0].query, input.bonds[1].query);
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
