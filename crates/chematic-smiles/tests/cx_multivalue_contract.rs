use chematic_core::{BondIdx, BondOrder};
use chematic_smiles::{parse_cxsmiles, write_cxsmiles};

#[test]
fn multiple_zero_bonds_survive_a_cxsmiles_roundtrip() {
    let cx = parse_cxsmiles("C~O~N |Z:0,1|").unwrap();
    assert_eq!(cx.mol.bond(BondIdx(0)).order, BondOrder::Zero);
    assert_eq!(cx.mol.bond(BondIdx(1)).order, BondOrder::Zero);
    let restored = parse_cxsmiles(&write_cxsmiles(&cx)).unwrap();
    assert_eq!(restored.mol.bond(BondIdx(0)).order, BondOrder::Zero);
    assert_eq!(restored.mol.bond(BondIdx(1)).order, BondOrder::Zero);
}

#[test]
fn multiple_wavy_bond_indices_survive_a_cxsmiles_roundtrip() {
    let cx = parse_cxsmiles("C=CC=C |w:0,2|").unwrap();
    assert_eq!(cx.wavy_bonds, vec![BondIdx(0), BondIdx(2)]);
    let restored = parse_cxsmiles(&write_cxsmiles(&cx)).unwrap();
    assert_eq!(restored.wavy_bonds, cx.wavy_bonds);
}
