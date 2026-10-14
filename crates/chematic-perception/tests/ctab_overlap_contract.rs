use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder};
use chematic_perception::{
    StereoRejectionReason, apply_ctab_local_parity_from_wedges_with_diagnostics,
    local_parity_from_wedges,
};

fn drawing(angle: f64, scale: f64, rotation: f64) -> (Molecule, Vec<(f64, f64)>) {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(Atom::new(Element::C));
    for (index, element) in [Element::F, Element::CL, Element::BR, Element::I]
        .into_iter()
        .enumerate()
    {
        let neighbor = builder.add_atom(Atom::new(element));
        builder
            .add_bond(
                center,
                neighbor,
                if index == 0 {
                    BondOrder::Up
                } else {
                    BondOrder::Single
                },
            )
            .unwrap();
    }
    let angle = angle.to_radians();
    let (s, c) = rotation.sin_cos();
    let coords = [
        (0.0, 0.0),
        (0.0, -1.0),
        (1.0, 0.0),
        (2.0 * angle.cos(), 2.0 * angle.sin()),
        (-1.0, 1.0),
    ]
    .into_iter()
    .map(|(x, y)| (scale * (c * x - s * y), scale * (s * x + c * y)))
    .collect();
    (builder.build(), coords)
}

#[test]
fn almost_parallel_flat_bonds_are_reported_without_assigning_ctab_chirality() {
    for scale in [0.5, 1.0, 10.0] {
        for rotation in [0.0, 0.7, 2.5] {
            let (mut mol, coords) = drawing(1.0, scale, rotation);
            let diagnostic =
                apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
            assert_eq!(diagnostic.len(), 1);
            assert_eq!(diagnostic[0].atom, AtomIdx(0));
            assert_eq!(
                diagnostic[0].reason,
                StereoRejectionReason::DegenerateGeometry
            );
            assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
            assert!(mol.stereo_neighbor_order(AtomIdx(0)).is_none());
        }
    }
}

#[test]
fn separated_flat_bonds_and_generic_geometry_keep_their_previous_parity() {
    let (mut mol, coords) = drawing(3.0, 1.0, 0.0);
    let expected = local_parity_from_wedges(&mol, &coords, AtomIdx(0)).unwrap();
    assert!(apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords).is_empty());
    assert_eq!(mol.atom(AtomIdx(0)).chirality, expected.0);
    assert_eq!(
        mol.stereo_neighbor_order(AtomIdx(0)),
        Some(expected.1.as_slice())
    );
    let (mol, coords) = drawing(1.0, 1.0, 0.0);
    assert!(local_parity_from_wedges(&mol, &coords, AtomIdx(0)).is_some());
}

#[test]
fn wedged_and_flat_neighbors_may_overlap_in_the_drawing() {
    let (mut mol, _) = drawing(3.0, 1.0, 0.0);
    let coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0)];
    assert!(apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords).is_empty());
    assert_ne!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
}
