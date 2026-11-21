use chematic_core::kekulization::{atom_must_be_matched, kekulize};
use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

#[test]
fn aromatic_silicon_anion_can_take_a_double_bond() {
    let mut b = MoleculeBuilder::new();
    let mut si = Atom::aromatic(Element::SI);
    si.charge = -1;
    let silicon = b.add_atom(si);
    let mut ring = vec![silicon];
    ring.extend((0..5).map(|_| b.add_atom(Atom::aromatic(Element::C))));
    for i in 0..6 {
        b.add_bond(ring[i], ring[(i + 1) % 6], BondOrder::Aromatic)
            .unwrap();
    }
    let mol = b.build();
    assert!(atom_must_be_matched(&mol, silicon));
    let kek = kekulize(&mol).expect("six-membered anionic silicon ring must kekulize");
    let doubles = mol
        .neighbors(silicon)
        .filter(|(_, edge)| kek.get(edge) == Some(&BondOrder::Double))
        .count();
    assert_eq!(doubles, 1);
}
