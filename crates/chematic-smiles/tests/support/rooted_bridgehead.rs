use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
};

pub fn bridgehead(explicit_h: bool, flips: u8, rotation: usize, reverse: bool) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for (i, element) in [
        Element::C,
        Element::C,
        Element::C,
        Element::O,
        Element::C,
        Element::C,
        Element::C,
    ]
    .into_iter()
    .enumerate()
    {
        let mut atom = Atom::new(element);
        atom.atom_map = Some((i + 1) as u16);
        if i == 2 || i == 4 {
            let ccw = (i == 2) ^ (flips & (if i == 2 { 1 } else { 2 }) != 0);
            atom.chirality = if ccw {
                Chirality::CounterClockwise
            } else {
                Chirality::Clockwise
            };
            if explicit_h {
                atom.hydrogen_count = Some(1);
            }
        }
        builder.add_atom(atom);
    }
    let mut edges = vec![
        (2, 3, 1),
        (3, 4, 1),
        (2, 1, 1),
        (4, 0, 1),
        (1, 0, 2),
        (5, 6, 1),
        (2, 6, 1),
        (4, 5, 1),
    ];
    edges.rotate_left(rotation % 8);
    if reverse {
        edges.reverse();
    }
    for (a, b, order) in edges {
        builder
            .add_bond(
                AtomIdx(a),
                AtomIdx(b),
                if order == 2 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            )
            .unwrap();
    }
    let mut mol = builder.build();
    for i in [2, 4] {
        // A declared programmatic frame: inserted neighbor slots, then virtual H.
        // This is independent of the traversal that will write the SMILES.
        let mut order: Vec<_> = mol.neighbors(AtomIdx(i)).map(|(nb, _)| nb.0).collect();
        order.push(STEREO_H_SENTINEL);
        mol.set_stereo_neighbor_order(AtomIdx(i), order);
    }
    mol
}
