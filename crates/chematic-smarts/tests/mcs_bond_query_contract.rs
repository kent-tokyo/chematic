use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_smarts::{
    find_mcs,
    query::{BondPrimitive as P, BondQuery as Q},
};

fn p(primitive: P) -> Q {
    Q::Primitive(primitive)
}
fn either(a: P, b: P) -> Q {
    Q::Or(Box::new(p(a)), Box::new(p(b)))
}

#[test]
fn single_molecule_mcs_preserves_query_bond_alternatives() {
    let cases = [
        (BondOrder::Single, p(P::Single)),
        (BondOrder::Double, p(P::Double)),
        (BondOrder::Triple, p(P::Triple)),
        (BondOrder::Aromatic, p(P::Aromatic)),
        (BondOrder::QueryAny, p(P::Any)),
        (BondOrder::QuerySingleOrDouble, either(P::Single, P::Double)),
        (
            BondOrder::QuerySingleOrAromatic,
            either(P::Single, P::Aromatic),
        ),
        (
            BondOrder::QueryDoubleOrAromatic,
            either(P::Double, P::Aromatic),
        ),
    ];
    for (order, expected) in cases {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, b, order).unwrap();
        let mol = builder.build();
        let query = find_mcs(&[&mol]);
        assert_eq!(query.bonds.len(), 1);
        assert_eq!(query.bonds[0].query, expected, "{order:?}");
    }
}
