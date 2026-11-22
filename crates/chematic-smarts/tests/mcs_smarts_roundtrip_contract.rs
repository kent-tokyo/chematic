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
