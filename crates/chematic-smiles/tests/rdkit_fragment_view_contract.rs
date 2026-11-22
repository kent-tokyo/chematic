use chematic_smiles::{parse, rdkit_mol_view};

#[test]
fn component_views_reindex_bonds_and_recompute_ring_membership() {
    let molecule = parse("C1CCCCC1.CCO").unwrap();
    let view = rdkit_mol_view(&molecule).unwrap();
    let components = view.fragments();
    assert_eq!(components, vec![vec![0, 1, 2, 3, 4, 5], vec![6, 7, 8]]);
    let ring = view.fragment(&components[0]).unwrap();
    assert_eq!((ring.num_atoms(), ring.num_bonds()), (6, 6));
    assert_eq!(ring.atom_rings.len(), 1);
    assert_eq!(ring.bond_rings.len(), 1);
    assert_eq!(ring.atom_rings[0].len(), 6);
    assert_eq!(ring.bond_rings[0].len(), 6);
    let chain = view.fragment(&components[1]).unwrap();
    assert_eq!((chain.num_atoms(), chain.num_bonds()), (3, 2));
    assert!(chain.atom_rings.is_empty());
    assert!(chain.bond_between(0, 1).is_some());
    assert!(chain.bond_between(1, 2).is_some());
    assert_eq!((view.num_atoms(), view.num_bonds()), (9, 8));
}
