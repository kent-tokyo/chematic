use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use chematic_smiles::{parse, write, write_with_atom_order};

#[test]
fn emitted_order_tracks_reordered_branches_without_changing_spelling() {
    let mut builder = MoleculeBuilder::new();
    for element in [Element::C, Element::C, Element::O, Element::N] {
        builder.add_atom(Atom::new(element));
    }
    for (a, b) in [(0, 1), (0, 2), (1, 3)] {
        builder
            .add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Single)
            .unwrap();
    }
    let mol = builder.build();
    let (text, order) = write_with_atom_order(&mol);
    assert_eq!(text, write(&mol));
    assert_eq!(order, [AtomIdx(0), AtomIdx(1), AtomIdx(3), AtomIdx(2)]);
    let parsed = parse(&text).unwrap();
    for (idx, source) in order.iter().enumerate() {
        assert_eq!(
            parsed.atom(AtomIdx(idx as u32)).element,
            mol.atom(*source).element
        );
    }
}

#[test]
fn empty_rings_and_fragments_have_one_mapping_per_emitted_atom() {
    for smiles in ["", "C1CC1.O", "[13CH3:17]C(Cl)(Br)F", "[H][H]"] {
        let mol = if smiles.is_empty() {
            MoleculeBuilder::new().build()
        } else {
            parse(smiles).unwrap()
        };
        let (text, mut order) = write_with_atom_order(&mol);
        assert_eq!(text, write(&mol));
        assert_eq!(order.len(), mol.atom_count());
        order.sort();
        assert_eq!(
            order,
            (0..mol.atom_count())
                .map(|i| AtomIdx(i as u32))
                .collect::<Vec<_>>()
        );
    }
}
