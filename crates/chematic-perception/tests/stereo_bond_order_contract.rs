use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::stereo_validation::{substituent_signatures_distinct, validate_stereo};

fn branches(first: BondOrder, second: BondOrder) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(Atom::new(Element::C));
    for order in [first, second] {
        let carbon = builder.add_atom(Atom::new(Element::C));
        let oxygen = builder.add_atom(Atom::new(Element::O));
        builder.add_bond(center, carbon, BondOrder::Single).unwrap();
        builder.add_bond(carbon, oxygen, order).unwrap();
    }
    builder.build()
}

#[test]
fn substituent_signatures_respect_nonstandard_bond_orders() {
    let orders = [
        BondOrder::Triple,
        BondOrder::Quadruple,
        BondOrder::Aromatic,
        BondOrder::Zero,
        BondOrder::Dative,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ];
    for order in orders {
        let identical = branches(order, order);
        assert!(
            !substituent_signatures_distinct(&identical, AtomIdx(0), &[AtomIdx(1), AtomIdx(3)]),
            "{order:?}"
        );
        assert!(validate_stereo(&identical).is_empty(), "{order:?}");
        let different = branches(BondOrder::Single, order);
        assert!(
            substituent_signatures_distinct(&different, AtomIdx(0), &[AtomIdx(1), AtomIdx(3)]),
            "{order:?}"
        );
        assert!(validate_stereo(&different).is_empty(), "{order:?}");
    }
}
