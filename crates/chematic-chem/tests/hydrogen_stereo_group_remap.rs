use chematic_chem::hydrogen::{add_hydrogens, remove_hydrogens};
use chematic_core::{AtomIdx, BondIdx, StereoGroup, StereoGroupKind};

#[test]
fn removing_hydrogens_remaps_stereo_groups_and_drops_empty_groups() {
    // The first atom and first bond disappear, shifting every surviving index.
    let mut mol = chematic_smiles::parse("[H]C(C)O").unwrap();
    mol.add_stereo_group_with_bonds(
        StereoGroup::new(
            StereoGroupKind::And(1),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(3)],
        ),
        vec![BondIdx(0), BondIdx(2)],
    );
    mol.add_stereo_group(StereoGroup::new(StereoGroupKind::Or(2), vec![AtomIdx(0)]));
    let result = remove_hydrogens(&mol);
    assert_eq!(result.atom_count(), 3);
    assert_eq!(result.stereo_groups().len(), 1);
    assert_eq!(result.stereo_groups()[0].kind, StereoGroupKind::And(1));
    assert_eq!(
        result.stereo_groups()[0].atom_indices,
        [AtomIdx(0), AtomIdx(2)]
    );
    assert_eq!(result.stereo_group_bond_indices(0), Some(&[BondIdx(1)][..]));
}

#[test]
fn adding_and_removing_hydrogens_keeps_bond_only_stereo_group() {
    let mut mol = chematic_smiles::parse("C=C").unwrap();
    mol.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::Absolute, vec![]),
        vec![BondIdx(0)],
    );
    let expanded = add_hydrogens(&mol);
    let restored = remove_hydrogens(&expanded);
    assert_eq!(restored.stereo_groups(), mol.stereo_groups());
    assert_eq!(
        restored.stereo_group_bond_indices(0),
        Some(&[BondIdx(0)][..])
    );
}
