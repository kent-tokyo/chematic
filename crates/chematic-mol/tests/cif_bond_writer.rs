//! Open Babel #2726: explicit bond geometry with exact atom-site references.
use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_mol::{
    CifBondWriteError, UnitCell, parse_cif, write_cif_checked, write_cif_with_bonds_checked,
};

fn graph(elements: &[Element], edges: &[(u32, u32, BondOrder)]) -> Molecule {
    let mut b = MoleculeBuilder::new();
    for &element in elements {
        b.add_atom(Atom::new(element));
    }
    for &(a, c, order) in edges {
        b.add_bond(AtomIdx(a), AtomIdx(c), order).unwrap();
    }
    b.build()
}

fn rows(text: &str) -> Vec<Vec<&str>> {
    text.split("_ccdc_geom_bond_type\n")
        .nth(1)
        .unwrap_or("")
        .lines()
        .map(|line| line.split_whitespace().collect())
        .collect()
}

#[test]
fn water_labels_endpoints_and_cartesian_distances() {
    let mol = graph(
        &[Element::O, Element::H, Element::H],
        &[(2, 0, BondOrder::Single), (0, 1, BondOrder::Single)],
    );
    let coords = [(0., 0., 0.), (0.9, 0., 0.), (0., 1.1, 0.)];
    let text = write_cif_with_bonds_checked(&mol, &coords, None).unwrap();
    assert_eq!(
        rows(&text),
        [
            vec!["H3", "O1", "1.1", ".", ".", "S"],
            vec!["O1", "H2", "0.9", ".", ".", "S"]
        ]
    );
    assert!(text.starts_with(&write_cif_checked(&mol, &coords, None).unwrap()));
    let parsed = parse_cif(&text).unwrap();
    assert_eq!(parsed.coords, coords);
    assert_eq!(
        parsed.mol.bond_count(),
        0,
        "plain CIF reader remains positions-only"
    );
}

#[test]
fn single_double_triple_aromatic_are_distinct() {
    let mol = graph(
        &[Element::C; 5],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Triple),
            (3, 4, BondOrder::Aromatic),
        ],
    );
    let coords: Vec<_> = (0..5).map(|i| (f64::from(i), 0., 0.)).collect();
    let text = write_cif_with_bonds_checked(&mol, &coords, None).unwrap();
    assert_eq!(
        rows(&text).iter().map(|r| r[5]).collect::<Vec<_>>(),
        ["S", "D", "T", "A"]
    );
}

#[test]
fn cells_change_coordinate_columns_but_not_bond_distance_or_images() {
    let mol = graph(&[Element::C; 2], &[(0, 1, BondOrder::Single)]);
    let coords = [(0.1, 0., 0.), (9.9, 0., 0.)];
    let plain = write_cif_with_bonds_checked(&mol, &coords, None).unwrap();
    for cell in [
        UnitCell {
            a: 10.,
            b: 10.,
            c: 10.,
            ..Default::default()
        },
        UnitCell {
            a: 12.,
            b: 14.,
            c: 15.,
            alpha: 75.,
            beta: 83.,
            gamma: 64.,
        },
    ] {
        let text = write_cif_with_bonds_checked(&mol, &coords, Some(&cell)).unwrap();
        assert!(text.contains("_atom_site_fract_x"));
        assert_eq!(rows(&text), rows(&plain));
        assert_eq!(rows(&text)[0][2].parse::<f64>().unwrap(), 9.8);
    }
}

#[test]
fn stable_norm_does_not_overflow_or_underflow_when_distance_is_representable() {
    let mol = graph(&[Element::C; 2], &[(0, 1, BondOrder::Single)]);
    for length in [1e308, 1e-300, f64::from_bits(1)] {
        let text =
            write_cif_with_bonds_checked(&mol, &[(0., 0., 0.), (length, 0., 0.)], None).unwrap();
        assert_eq!(rows(&text)[0][2].parse::<f64>().unwrap(), length);
    }
    let text = write_cif_with_bonds_checked(&mol, &[(0., 0., 0.), (3., 4., 12.)], None).unwrap();
    assert_eq!(rows(&text)[0][2], "13");
}

#[test]
fn nonfinite_and_zero_bond_distances_refuse() {
    let mol = graph(&[Element::C; 2], &[(0, 1, BondOrder::Single)]);
    for coords in [
        [(0., 0., 0.), (0., 0., 0.)],
        [(1e308, 0., 0.), (-1e308, 0., 0.)],
    ] {
        assert!(matches!(
            write_cif_with_bonds_checked(&mol, &coords, None),
            Err(CifBondWriteError::InvalidDistance { bond: 0, .. })
        ));
    }
    let self_bond = graph(&[Element::C], &[(0, 0, BondOrder::Single)]);
    assert!(matches!(
        write_cif_with_bonds_checked(&self_bond, &[(0., 0., 0.)], None),
        Err(CifBondWriteError::InvalidDistance { .. })
    ));
}

#[test]
fn unsupported_orders_never_become_single_bonds() {
    for order in [
        BondOrder::Zero,
        BondOrder::Dative,
        BondOrder::Quadruple,
        BondOrder::Up,
        BondOrder::Down,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mol = graph(&[Element::C; 2], &[(0, 1, order)]);
        assert!(
            matches!(write_cif_with_bonds_checked(&mol, &[(0.,0.,0.), (1.,0.,0.)], None), Err(CifBondWriteError::UnsupportedBond { bond:0, order: actual }) if actual == order)
        );
    }
}

#[test]
fn atom_and_cell_validation_is_shared_with_existing_checked_writer() {
    let mol = graph(&[Element::C], &[]);
    for coords in [vec![], vec![(f64::NAN, 0., 0.)], vec![(0., 0., 0.); 2]] {
        assert!(matches!(
            write_cif_with_bonds_checked(&mol, &coords, None),
            Err(CifBondWriteError::InvalidInput(_))
        ));
    }
    let bad = UnitCell {
        b: 0.,
        ..Default::default()
    };
    assert!(matches!(
        write_cif_with_bonds_checked(&mol, &[(0., 0., 0.)], Some(&bad)),
        Err(CifBondWriteError::InvalidInput(_))
    ));
}

#[test]
fn empty_and_disconnected_graphs_do_not_gain_bonds() {
    for elements in [vec![], vec![Element::NA, Element::CL]] {
        let mol = graph(&elements, &[]);
        let coords = vec![(0., 0., 0.); elements.len()];
        let text = write_cif_with_bonds_checked(&mol, &coords, None).unwrap();
        assert_eq!(text, write_cif_checked(&mol, &coords, None).unwrap());
        assert!(!text.contains("_geom_bond"));
    }
}
