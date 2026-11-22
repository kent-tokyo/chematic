use chematic_smarts::{
    find_matches, mcs_query_to_smarts, parse_smarts,
    query::{AtomPrimitive as A, AtomQuery as Q, BondPrimitive as B, BondQuery, QueryMolecule},
};
use chematic_smiles::parse;

fn atom(number: u8) -> Q {
    Q::Primitive(A::AtomicNum(number))
}

#[test]
fn mcs_smarts_preserves_isotope_and_charge_constraints() {
    let mut graph = QueryMolecule::new();
    let carbon = graph.add_atom(Q::And(
        Box::new(Q::And(
            Box::new(atom(6)),
            Box::new(Q::Primitive(A::Isotope(13))),
        )),
        Box::new(Q::Primitive(A::Charge(-1))),
    ));
    let fluorine = graph.add_atom(atom(9));
    graph.add_bond(carbon, fluorine, BondQuery::Primitive(B::Single));
    let text = mcs_query_to_smarts(&graph).unwrap();
    let query = parse_smarts(&text).unwrap();
    for (smiles, expected) in [
        ("[13C-]F", true),
        ("[12C-]F", false),
        ("[13C+]F", false),
        ("CF", false),
    ] {
        assert_eq!(
            !find_matches(&query, &parse(smiles).unwrap()).is_empty(),
            expected,
            "{text} against {smiles}"
        );
    }
}

#[test]
fn mcs_smarts_preserves_quadruple_and_any_bond_queries() {
    for (bond, yes, no) in [
        (B::Quadruple, "[Re]$[Re]", "[Re][Re]"),
        (B::Any, "[Re]$[Re]", "[Re].[Re]"),
    ] {
        let mut graph = QueryMolecule::new();
        let a = graph.add_atom(atom(75));
        let b = graph.add_atom(atom(75));
        graph.add_bond(a, b, BondQuery::Primitive(bond));
        let text = mcs_query_to_smarts(&graph).unwrap();
        let query = parse_smarts(&text).unwrap();
        assert!(
            !find_matches(&query, &parse(yes).unwrap()).is_empty(),
            "{text}"
        );
        assert!(
            find_matches(&query, &parse(no).unwrap()).is_empty(),
            "{text}"
        );
    }
}

#[test]
fn implicit_mcs_bond_keeps_single_or_aromatic_semantics() {
    let mut graph = QueryMolecule::new();
    let a = graph.add_atom(atom(6));
    let b = graph.add_atom(atom(6));
    graph.add_bond(a, b, BondQuery::Any);
    let query = parse_smarts(&mcs_query_to_smarts(&graph).unwrap()).unwrap();
    for (smiles, expected) in [("CC", true), ("c1ccccc1", true), ("C=C", false)] {
        assert_eq!(
            !find_matches(&query, &parse(smiles).unwrap()).is_empty(),
            expected,
            "{smiles}"
        );
    }
}

fn triangle_chain(count: usize) -> QueryMolecule {
    let mut graph = QueryMolecule::new();
    let mut previous = None;
    for _ in 0..count {
        let a = graph.add_atom(atom(6));
        let b = graph.add_atom(atom(6));
        let c = graph.add_atom(atom(6));
        for (left, right) in [(a, b), (b, c), (c, a)] {
            graph.add_bond(left, right, BondQuery::Primitive(B::Single));
        }
        if let Some(last) = previous {
            graph.add_bond(last, a, BondQuery::Primitive(B::Single));
        }
        previous = Some(c);
    }
    graph
}

#[test]
fn mcs_smarts_roundtrips_more_than_nine_ring_closures() {
    let graph = triangle_chain(12);
    let text = mcs_query_to_smarts(&graph).unwrap();
    assert!(text.contains("%10"), "{text}");
    let reread = parse_smarts(&text).unwrap();
    assert_eq!((reread.atoms.len(), reread.bonds.len()), (36, 47));
    let mut before: Vec<_> = graph.adj.iter().map(Vec::len).collect();
    let mut after: Vec<_> = reread.adj.iter().map(Vec::len).collect();
    before.sort_unstable();
    after.sort_unstable();
    assert_eq!(before, after);
}

#[test]
fn mcs_smarts_refuses_more_than_99_ring_closures() {
    let boundary = mcs_query_to_smarts(&triangle_chain(99)).unwrap();
    assert!(boundary.contains("%99"));
    assert!(parse_smarts(&boundary).is_ok());
    assert_eq!(
        mcs_query_to_smarts(&triangle_chain(100)).unwrap_err(),
        chematic_smarts::McsSmartsError::TooManyRingClosures
    );
}
