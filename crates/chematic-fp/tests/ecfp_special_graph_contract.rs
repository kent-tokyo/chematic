use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder, SquarePlanarPermutation,
};
use chematic_fp::{EcfpConfig, ecfp};

#[test]
fn ecfp_distinguishes_special_bond_orders() {
    let orders = [
        BondOrder::Single,
        BondOrder::Quadruple,
        BondOrder::Zero,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ];
    let config = EcfpConfig {
        radius: 1,
        ..Default::default()
    };
    let fingerprints: Vec<_> = orders
        .iter()
        .map(|&order| {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(Atom::new(Element::C));
            builder.add_atom(Atom::new(Element::N));
            builder.add_bond(AtomIdx(0), AtomIdx(1), order).unwrap();
            ecfp(&builder.build(), &config)
        })
        .collect();
    for i in 0..orders.len() {
        for j in i + 1..orders.len() {
            assert_ne!(
                fingerprints[i], fingerprints[j],
                "{:?} vs {:?}",
                orders[i], orders[j]
            );
        }
    }
}

#[test]
fn ecfp_stereo_flag_controls_non_tetrahedral_tags() {
    let tags = [
        Chirality::None,
        Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        Chirality::SquarePlanarUnnumbered,
        Chirality::TrigonalBipyramidal(1),
        Chirality::Octahedral(1),
    ];
    let mut fingerprints = Vec::new();
    let mut ignored = Vec::new();
    for tag in tags {
        let mut atom = Atom::new(Element::P);
        atom.chirality = tag;
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(atom);
        let molecule = builder.build();
        fingerprints.push(ecfp(
            &molecule,
            &EcfpConfig {
                radius: 0,
                use_chirality: true,
                ..Default::default()
            },
        ));
        ignored.push(ecfp(
            &molecule,
            &EcfpConfig {
                radius: 0,
                use_chirality: false,
                ..Default::default()
            },
        ));
    }
    for fp in &ignored {
        assert_eq!(fp, &ignored[0]);
    }
    for i in 0..fingerprints.len() {
        for j in i + 1..fingerprints.len() {
            assert_ne!(fingerprints[i], fingerprints[j]);
        }
    }
}
