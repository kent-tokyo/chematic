use chematic_core::AtomIdx;
use chematic_perception::{
    find_smallest_rings_bfs_with_blocked_bonds, find_smallest_rings_bfs_with_trimmed_bonds,
};

#[test]
fn leaf_trimming_preserves_rooted_cycles_and_respects_blocked_edges() {
    let mol = chematic_smiles::parse("CCC1CCCCC1CC").unwrap();
    let original_bond_count = mol.bond_count();
    let root = AtomIdx(2);
    let mut blocked = Default::default();
    let before = find_smallest_rings_bfs_with_blocked_bonds(&mol, root, &blocked);
    let trimmed = find_smallest_rings_bfs_with_trimmed_bonds(&mol, root, &blocked);
    assert_eq!(trimmed, before);
    assert_eq!(trimmed.len(), 1);
    assert_eq!(trimmed[0].len(), 6);
    let edge = mol.bond_between(AtomIdx(2), AtomIdx(3)).unwrap().0;
    blocked.insert(edge);
    assert!(find_smallest_rings_bfs_with_trimmed_bonds(&mol, root, &blocked).is_empty());
    assert_eq!(mol.bond_count(), original_bond_count);
}
