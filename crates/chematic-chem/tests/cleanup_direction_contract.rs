use chematic_core::{AtomIdx, BondIdx, BondOrder};

#[test]
fn metal_disconnection_preserves_direction_and_anchor_on_retained_bonds() {
    let mut molecule = chematic_smiles::parse("[Na]O.F/C=C/F").unwrap();
    molecule.set_bond_direction(BondIdx(1), BondOrder::Up);
    molecule.set_bond_direction_anchor(BondIdx(1), AtomIdx(2));
    let cleaned = chematic_chem::cleanup_unchecked(&molecule);
    assert_eq!(cleaned.bond_count(), 3);
    assert_eq!(cleaned.bond_direction(BondIdx(0)), Some(BondOrder::Up));
    assert_eq!(cleaned.bond_direction_anchor(BondIdx(0)), Some(AtomIdx(2)));
}
