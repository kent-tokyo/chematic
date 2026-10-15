use super::*;
use chematic_core::{Element, MoleculeBuilder};

#[test]
fn numeric_edges_and_missing_coordinates_are_byte_identical() {
    let mut b = MoleculeBuilder::new();
    let mut atom = Atom::new(Element::C);
    atom.charge = -2;
    atom.isotope = Some(13);
    atom.hydrogen_count = Some(0);
    atom.atom_map = Some(u16::MAX);
    b.add_atom(atom);
    let mol = b.build();
    let metadata = MolMetadata::default();
    let values = [
        0.0,
        -0.0,
        0.00005,
        -0.00005,
        1.23445,
        -999.99995,
        1e11,
        f64::MAX,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    for x in values {
        let coords = [(x, -x)];
        assert_eq!(
            write_mol_v3000(&mol, &metadata, &coords),
            reference::write_mol_v3000(&mol, &metadata, &coords)
        );
        let conformer = Coords3D {
            points: vec![Point3::new(x, -x, x)],
        };
        assert_eq!(
            write_mol_v3000_with_conformer(&mol, &metadata, &conformer),
            reference::write_mol_v3000_with_conformer(&mol, &metadata, &conformer)
        );
    }
    assert_eq!(
        write_mol_v3000(&mol, &metadata, &[]),
        reference::write_mol_v3000(&mol, &metadata, &[])
    );
    let empty = Coords3D::new_zeroed(0);
    assert_eq!(
        write_mol_v3000_with_conformer(&mol, &metadata, &empty),
        reference::write_mol_v3000_with_conformer(&mol, &metadata, &empty)
    );
}

#[test]
fn growing_records_and_opaque_properties_are_byte_identical() {
    for count in [0, 1, 2, 8, 64, 256, 1024] {
        let mut b = MoleculeBuilder::new();
        for i in 0..count {
            let mut atom = Atom::new(Element::C);
            if i % 7 == 0 {
                atom.atom_map = Some(i as u16);
            }
            b.add_atom(atom);
            if i > 0 {
                b.add_bond(
                    AtomIdx((i - 1) as u32),
                    AtomIdx(i as u32),
                    BondOrder::Single,
                )
                .unwrap();
            }
        }
        let mol = b.build();
        let metadata = MolMetadata {
            name: "テスト name".into(),
            comment: "opaque bytes".into(),
            v3000_atom_properties: vec![
                (1, "FOO=(2 1 2)  LABEL=abc".into()),
                (1, "ignored duplicate".into()),
            ],
            v3000_bond_properties: vec![
                (1, "  CFG=2  TOPO=1  ".into()),
                (1, "ignored duplicate".into()),
            ],
            ..Default::default()
        };
        let coords: Vec<_> = (0..count).map(|i| (i as f64 * 1.2345, 0.12345)).collect();
        let conformer = Coords3D {
            points: coords.iter().map(|&(x, y)| Point3::new(x, y, -y)).collect(),
        };
        assert_eq!(
            write_mol_v3000(&mol, &metadata, &coords),
            reference::write_mol_v3000(&mol, &metadata, &coords)
        );
        assert_eq!(
            write_mol_v3000_with_conformer(&mol, &metadata, &conformer),
            reference::write_mol_v3000_with_conformer(&mol, &metadata, &conformer)
        );
        assert_eq!(
            write_mol_v3000_with_conformer_checked(&mol, &metadata, &conformer).unwrap(),
            reference::write_mol_v3000_with_conformer_checked(&mol, &metadata, &conformer).unwrap()
        );
    }
}

#[test]
fn cfg_filter_retains_first_match_and_historical_whitespace() {
    for properties in [
        "",
        "   ",
        "CFG=2",
        "CFG=2 CFG=1",
        " CFG=2  TOPO=1\tENDPTS=(2 1 2) ",
        "cfg=2 CFGX=3 X=5",
        "LABEL=\"a  b\" CFG=2",
    ] {
        let metadata = MolMetadata {
            v3000_bond_properties: vec![
                (9, "other".into()),
                (1, properties.into()),
                (1, "ignored duplicate".into()),
            ],
            ..Default::default()
        };
        for written in [false, true] {
            for id in [1, 9, 99] {
                let mut actual = String::new();
                append_bond_properties(&mut actual, &metadata, id, written);
                assert_eq!(
                    actual,
                    reference::properties_for_test(&metadata, id, written)
                );
            }
        }
    }
}
