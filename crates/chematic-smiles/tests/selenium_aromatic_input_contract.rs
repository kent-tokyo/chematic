use chematic_core::{AtomIdx, Element};
use chematic_smiles::parse;

#[test]
fn unbracketed_aromatic_selenium_is_one_atom() {
    let molecule = parse("se1cccc1").unwrap();
    assert_eq!(molecule.atom_count(), 5);
    assert_eq!(molecule.atom(AtomIdx(0)).element, Element::SE);
    assert!(molecule.atom(AtomIdx(0)).aromatic);
    assert_eq!(molecule.bond_count(), 5);
}
