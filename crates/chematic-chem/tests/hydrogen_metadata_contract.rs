use chematic_chem::{add_hydrogens, remove_hydrogens};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder, RGroupLabel};

#[test]
fn hydrogen_expansion_and_removal_preserve_heavy_atom_metadata() {
    let mut builder = MoleculeBuilder::new();
    // Place the explicit H first so removing it changes the carbon's index.
    let hydrogen = builder.add_atom(Atom::new(Element::H));
    let carbon = builder.add_atom(Atom::new(Element::C));
    builder
        .add_bond(hydrogen, carbon, BondOrder::Single)
        .unwrap();
    let placeholder = builder.add_atom(Atom::new(Element::C));
    builder
        .add_bond(carbon, placeholder, BondOrder::Single)
        .unwrap();
    let label = RGroupLabel::numbered(7).unwrap();
    builder.set_r_group(placeholder, label);
    builder.set_tag(carbon, Some(431));
    let molecule = builder.build();

    let expanded = add_hydrogens(&molecule);
    assert_eq!(expanded.atom_tag(carbon).map(|tag| tag.get()), Some(431));
    assert_eq!(expanded.r_group_label(placeholder), Some(label));
    assert_eq!(expanded.atom_count(), 5);

    let compact = remove_hydrogens(&expanded);
    assert_eq!(compact.atom_count(), 2);
    assert_eq!(compact.atom(AtomIdx(0)).element, Element::C);
    assert_eq!(compact.atom_tag(AtomIdx(0)).map(|tag| tag.get()), Some(431));
    assert_eq!(compact.r_group_label(AtomIdx(1)), Some(label));
}
