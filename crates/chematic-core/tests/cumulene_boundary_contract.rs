use chematic_core::{
    Atom, AtomIdx, BondOrder, Element, MoleculeBuilder, cumulene::planar_cumulene_at,
};

fn branched_axis(
    reverse_substituent_order: bool,
) -> (chematic_core::Molecule, chematic_core::BondIdx) {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..8 {
        builder.add_atom(Atom::new(Element::C));
    }
    for (a, b) in [(0, 1), (1, 2), (2, 3), (5, 7)] {
        builder
            .add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Double)
            .unwrap();
    }
    let mut reference = None;
    let edges = if reverse_substituent_order {
        [(0, 5), (0, 4)]
    } else {
        [(0, 4), (0, 5)]
    };
    for (a, b) in edges {
        let bond = builder
            .add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Up)
            .unwrap();
        if b == 4 {
            reference = Some(bond);
        }
    }
    builder
        .add_bond(AtomIdx(3), AtomIdx(6), BondOrder::Up)
        .unwrap();
    (builder.build(), reference.unwrap())
}

#[test]
fn cumulene_terminal_prefers_nonconjugated_direction_reference_in_either_bond_order() {
    for reverse in [false, true] {
        let (mol, expected_reference) = branched_axis(reverse);
        let axis = planar_cumulene_at(&mol, AtomIdx(0)).unwrap();
        assert_eq!(axis.atoms, [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)]);
        assert_eq!(axis.ends[0].reference, expected_reference);
    }
}

#[test]
fn tetrahedral_remap_reports_duplicate_in_canonical_order() {
    use chematic_core::stereo_geometry::{StereoGeometryError, remap_tetrahedral_parity};
    assert_eq!(
        remap_tetrahedral_parity([1, 2, 3, 4], [1, 2, 2, 4]),
        Err(StereoGeometryError::DuplicateSlotId(2))
    );
}

#[test]
fn cumulene_end_requires_consistent_markers_and_at_most_two_substituents() {
    fn molecule(second_marker: BondOrder, extra_substituent: bool) -> chematic_core::Molecule {
        let mut builder = MoleculeBuilder::new();
        for _ in 0..8 {
            builder.add_atom(Atom::new(Element::C));
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            builder
                .add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Double)
                .unwrap();
        }
        builder
            .add_bond(AtomIdx(0), AtomIdx(4), BondOrder::Up)
            .unwrap();
        builder
            .add_bond(AtomIdx(0), AtomIdx(5), second_marker)
            .unwrap();
        builder
            .add_bond(AtomIdx(3), AtomIdx(6), BondOrder::Up)
            .unwrap();
        if extra_substituent {
            builder
                .add_bond(AtomIdx(0), AtomIdx(7), BondOrder::Single)
                .unwrap();
        }
        builder.build()
    }
    assert!(planar_cumulene_at(&molecule(BondOrder::Down, false), AtomIdx(0)).is_some());
    assert!(planar_cumulene_at(&molecule(BondOrder::Up, false), AtomIdx(0)).is_none());
    assert!(planar_cumulene_at(&molecule(BondOrder::Down, true), AtomIdx(0)).is_none());
}
