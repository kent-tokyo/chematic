use chematic_core::{Atom, Element, MoleculeBuilder};

#[test]
fn isotope_formula_counts_labeled_and_unlabeled_atoms_independently() {
    let mut builder = MoleculeBuilder::new();
    for element in [Element::H, Element::C, Element::C, Element::H, Element::O] {
        builder.add_atom(Atom::new(element));
    }
    for _ in 0..2 {
        let mut isotope = Atom::new(Element::C);
        isotope.isotope = Some(13);
        builder.add_atom(isotope);
    }
    let mol = builder.build();
    assert_eq!(mol.formula_with_isotopes(), "C2H213C2O");
    assert_eq!(mol.atom_count(), 7);
}
