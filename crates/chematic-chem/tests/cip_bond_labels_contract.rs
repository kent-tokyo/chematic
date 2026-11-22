use chematic_chem::assign_accurate_cip_with_bond_labels;
use chematic_core::{BondIdx, CipCode};

#[test]
fn accurate_cip_bond_api_returns_absolute_ez_labels() {
    for (smiles, expected) in [(r"F/C=C/F", CipCode::E), (r"F/C=C\F", CipCode::Z)] {
        let molecule = chematic_smiles::parse(smiles).unwrap();
        let (_, labels) = assign_accurate_cip_with_bond_labels(&molecule).unwrap();
        assert_eq!(labels, vec![(BondIdx(1), expected)], "{smiles}");
    }
    let molecule = chematic_smiles::parse("FC=CF").unwrap();
    let (atoms, bonds) = assign_accurate_cip_with_bond_labels(&molecule).unwrap();
    assert!(atoms.assignments.is_empty());
    assert!(atoms.unresolved.is_empty());
    assert!(bonds.is_empty());
}
