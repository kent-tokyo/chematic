use chematic_smarts::{
    McsSmartsError, mcs_query_to_smarts,
    query::{AtomPrimitive as A, AtomQuery as Q, BondPrimitive as B, BondQuery, QueryMolecule},
};

#[test]
fn mcs_smarts_refuses_unrepresentable_query_and_stereo() {
    let mut query = QueryMolecule::new();
    query.add_atom(Q::Primitive(A::Charge(1)));
    let error = mcs_query_to_smarts(&query).unwrap_err();
    assert_eq!(error, McsSmartsError::UnsupportedAtomQuery);
    assert_eq!(error.to_string(), "UnsupportedAtomQuery");

    let mut query = QueryMolecule::new();
    let center = query.add_atom(Q::And(
        Box::new(Q::Primitive(A::AtomicNum(6))),
        Box::new(Q::Primitive(A::Chirality(1))),
    ));
    for _ in 0..5 {
        let neighbor = query.add_atom(Q::Primitive(A::AtomicNum(9)));
        query.add_bond(center, neighbor, BondQuery::Primitive(B::Single));
    }
    let error = mcs_query_to_smarts(&query).unwrap_err();
    assert_eq!(error, McsSmartsError::UnsupportedChirality);
    assert!(
        error
            .to_string()
            .contains("representable absolute tetrahedral stereo")
    );
}
