use super::*;
use chematic_core::{BondIdx, Element, MoleculeBuilder, RGroupLabel, StereoGroup, StereoGroupKind};

#[test]
fn numeric_edges_and_missing_coordinates_are_byte_identical() {
    let mut b = MoleculeBuilder::new();
    let mut atom = Atom::new(Element::C);
    atom.charge = -2;
    atom.isotope = Some(13);
    // Keep this coordinate-format comparison outside the intentional
    // zero-valence CTfile change exercised by v3000_bracket_valence.
    atom.hydrogen_count = Some(1);
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

#[test]
fn reference_writer_covers_query_bonds_rgroups_and_enhanced_stereo() {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(Atom::new(Element::C));
    let numbered = builder.add_atom(Atom::new(Element::C));
    builder.set_r_group(numbered, RGroupLabel::numbered(17).unwrap());
    let unnumbered = builder.add_atom(Atom::new(Element::C));
    builder.set_r_group(unnumbered, RGroupLabel::unnumbered());
    let mut wildcard = Atom::new(Element::C);
    wildcard.wildcard = true;
    let wildcard = builder.add_atom(wildcard);

    let orders = [
        BondOrder::Zero,
        BondOrder::Single,
        BondOrder::Up,
        BondOrder::Down,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Aromatic,
        BondOrder::Quadruple,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
        BondOrder::QueryAny,
        BondOrder::Dative,
    ];
    let mut bonds = Vec::new();
    for order in orders {
        let leaf = builder.add_atom(Atom::new(Element::C));
        bonds.push(builder.add_bond(center, leaf, order).unwrap());
    }
    builder.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::Absolute, vec![numbered]),
        vec![bonds[0]],
    );
    builder.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::Or(2), vec![unnumbered]),
        vec![bonds[1]],
    );
    builder.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::And(3), vec![wildcard]),
        vec![bonds[2]],
    );
    let mol = builder.build();
    let metadata = MolMetadata {
        v3000_sgroups: vec!["1 DAT 0 ATOMS=(1 1) FIELDNAME=source".into()],
        v3000_atom_properties: vec![(1, "INVRET=1".into())],
        v3000_bond_properties: vec![(1, "CFG=2 TOPO=1".into())],
        ..Default::default()
    };
    let conformer = Coords3D {
        points: (0..mol.atom_count())
            .map(|index| Point3::new(index as f64, -(index as f64), 0.5))
            .collect(),
    };

    let expected_2d = write_mol_v3000(&mol, &metadata, &[]);
    let reference_2d = reference::write_mol_v3000(&mol, &metadata, &[]);
    assert_eq!(reference_2d, expected_2d);
    assert!(reference_2d.contains("R#"));
    assert!(reference_2d.contains(" RGROUPS=(1 17)"));
    assert!(reference_2d.contains(" MDLV30/STEABS ATOMS="));
    assert!(reference_2d.contains(" MDLV30/STEBREL2 BONDS="));
    assert!(reference_2d.contains(" MDLV30/STERAC3 ATOMS="));

    assert_eq!(
        reference::write_mol_v3000_with_conformer(&mol, &metadata, &conformer),
        write_mol_v3000_with_conformer(&mol, &metadata, &conformer)
    );
}

#[test]
fn reference_writer_exercises_depiction_and_checked_failure_paths() {
    let chiral = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
    let metadata = MolMetadata::default();
    let (reference_block, reference_loss) =
        reference::write_mol_v3000_with_stereo_report(&chiral, &metadata, &[]);
    let (optimized_block, optimized_loss) =
        write_mol_v3000_with_stereo_report(&chiral, &metadata, &[]);
    assert_eq!(reference_loss, optimized_loss);
    assert!(reference_block.contains(" CFG="));
    let reference_roundtrip = read_mol_v3000_with_diagnostics(&reference_block).unwrap();
    let optimized_roundtrip = read_mol_v3000_with_diagnostics(&optimized_block).unwrap();
    assert_eq!(reference_roundtrip.mol.atom_count(), chiral.atom_count());
    assert_eq!(reference_roundtrip.mol.bond_count(), chiral.bond_count());
    assert_eq!(
        reference_roundtrip.mol.atom(AtomIdx(1)).chirality,
        optimized_roundtrip.mol.atom(AtomIdx(1)).chirality
    );

    let square_planar = chematic_smiles::parse("[Pt@SP1](Cl)(Cl)(N)N").unwrap();
    let empty = Coords3D::new_zeroed(0);
    assert_eq!(
        reference::write_mol_v3000_with_conformer_checked(&square_planar, &metadata, &empty,)
            .unwrap_err(),
        write_mol_v3000_with_conformer_checked(&square_planar, &metadata, &empty).unwrap_err()
    );

    // The first bond is also the first enhanced-stereo bond member above.
    assert_eq!(BondIdx(0).0, 0);
}
