use super::*;
use crate::derived_cache::{DerivedScalar, DerivedSlot};

fn rich_graph() -> Molecule {
    let mut b = MoleculeBuilder::new();
    let a = b.add_atom(Atom::new(Element::C));
    let c = b.add_atom(Atom::new(Element::O));
    b.add_bond(a, c, BondOrder::Single).unwrap();
    b.set_tag(a, Some(29));
    b.set_r_group(a, RGroupLabel::numbered(17).unwrap());
    b.set_stereo_neighbor_order(a, vec![c.0, STEREO_H_SENTINEL]);
    b.set_bond_direction(BondIdx(0), BondOrder::Up);
    b.set_bond_direction_anchor(BondIdx(0), a);
    b.set_smiles_ring_closure(BondIdx(0), 7, 1);
    let mut mol = b.build();
    mol.set_stereo_groups_with_bonds(vec![(
        StereoGroup::new(crate::StereoGroupKind::And(3), vec![a]),
        vec![BondIdx(0)],
    )]);
    mol
}

#[test]
fn clear_removes_all_graph_and_source_metadata() {
    let mut mol = rich_graph();
    mol.clear();
    assert!(mol.atoms.is_empty());
    assert!(mol.atom_tags.is_empty());
    assert!(mol.bonds.is_empty());
    assert!(mol.adjacency.is_empty());
    assert!(mol.r_groups.is_empty());
    assert!(mol.stereo_groups.is_empty());
    assert!(mol.stereo_group_bonds.is_empty());
    assert!(mol.stereo_neighbor_order.is_empty());
    assert!(mol.bond_directions.is_empty());
    assert!(mol.bond_direction_anchors.is_empty());
    assert!(mol.smiles_ring_closure_keys.is_empty());
    assert_eq!(mol.formula(), "");
}

#[test]
fn clear_drops_all_cache_kinds_but_retained_derived_views_remain_valid() {
    let mut mol = rich_graph();
    // Diagnostic scalar values test cache lifecycle, not descriptor chemistry.
    for slot in [
        DerivedScalar::Tpsa,
        DerivedScalar::RdkitTpsa,
        DerivedScalar::LabuteAsa,
        DerivedScalar::Qed,
    ] {
        assert_eq!(mol.derived_scalar(slot, || 7.0), 7.0);
    }
    assert_eq!(mol.derived_rdkit_num_rings(|| 9), 9);
    let view = mol.derived(DerivedSlot::RdkitScalarModel, || Some(mol.clone()));
    assert!(view.is_some());
    mol.clear();
    assert!(
        mol.derived_if_computed::<Option<Molecule>>(DerivedSlot::RdkitScalarModel)
            .is_none()
    );
    for slot in [
        DerivedScalar::Tpsa,
        DerivedScalar::RdkitTpsa,
        DerivedScalar::LabuteAsa,
        DerivedScalar::Qed,
    ] {
        assert_eq!(mol.derived_scalar(slot, || 0.0), 0.0);
    }
    assert_eq!(mol.derived_rdkit_num_rings(|| 0), 0);
    assert_eq!(view.as_ref().as_ref().unwrap().atom_count(), 2);
}

#[test]
fn clear_keeps_clones_independent_and_allows_reuse_at_zero() {
    let mut mol = rich_graph();
    let snapshot = mol.clone();
    for _ in 0..10 {
        mol.clear();
        mol.clear();
        assert_eq!(mol.add_atom(Atom::new(Element::N)), AtomIdx(0));
        assert_eq!(mol.atom_tag(AtomIdx(0)), None);
        assert_eq!(mol.r_group_label(AtomIdx(0)), None);
        assert_eq!(mol.stereo_neighbor_order(AtomIdx(0)), None);
    }
    assert_eq!(snapshot.atom_count(), 2);
    assert_eq!(snapshot.atom_tag(AtomIdx(0)).unwrap().get(), 29);
    assert_eq!(
        snapshot.r_group_label(AtomIdx(0)),
        RGroupLabel::numbered(17)
    );
    assert_eq!(snapshot.stereo_groups.len(), 1);
}
