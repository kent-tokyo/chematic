use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder as B, Chirality, Element, Molecule, MoleculeBuilder,
    StereoGroup, StereoGroupKind,
};
use chematic_mol::{CjsonWriteError as E, parse_cjson_checked, write_cjson, write_cjson_checked};

fn pair(order: B) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut carbon = Atom::new(Element::C);
    carbon.isotope = Some(13);
    carbon.charge = -1;
    let a = builder.add_atom(carbon);
    let b = builder.add_atom(Atom::new(Element::N));
    builder.add_bond(a, b, order).unwrap();
    builder.build()
}

#[test]
fn accepted_graph_and_coordinate_bits_roundtrip_without_perception() {
    for order in [B::Single, B::Double, B::Triple, B::Aromatic] {
        let mol = pair(order);
        let coords = [
            (-0.0, 1e-200, f64::MAX),
            (1.2345678901234567, -1e200, f64::MIN_POSITIVE),
        ];
        for points in [&[][..], &coords[..]] {
            let output = write_cjson_checked(&mol, points).unwrap();
            let (parsed, actual) = parse_cjson_checked(&output).unwrap();
            assert_eq!(mol.atom_count(), parsed.atom_count());
            for ((index, atom), (other, decoded)) in mol.atoms().zip(parsed.atoms()) {
                assert_eq!(index, other);
                assert_eq!(atom, decoded);
            }
            for ((index, bond), (other, decoded)) in mol.bonds().zip(parsed.bonds()) {
                assert_eq!(index, other);
                assert_eq!(bond, decoded);
            }
            assert_eq!(actual.len(), points.len());
            for (expected, actual) in points.iter().zip(actual) {
                assert_eq!(
                    [
                        expected.0.to_bits(),
                        expected.1.to_bits(),
                        expected.2.to_bits()
                    ],
                    [actual.0.to_bits(), actual.1.to_bits(), actual.2.to_bits()]
                );
            }
        }
    }
    assert!(write_cjson_checked(&MoleculeBuilder::new().build(), &[]).is_ok());
}

#[test]
fn refuses_every_coerced_or_ambiguous_bond_order() {
    for order in [
        B::Zero,
        B::Dative,
        B::Quadruple,
        B::Up,
        B::Down,
        B::QueryAny,
        B::QuerySingleOrDouble,
        B::QuerySingleOrAromatic,
        B::QueryDoubleOrAromatic,
    ] {
        let mol = pair(order);
        assert_eq!(
            write_cjson_checked(&mol, &[]),
            Err(E::BondOrder {
                bond_index: 0,
                order
            })
        );
        let (legacy, _) = parse_cjson_checked(&write_cjson(&mol, &[])).unwrap();
        assert_ne!(legacy.bond(BondIdx(0)).order, order);
    }
}

#[test]
fn refuses_unsupported_atom_properties_individually() {
    let mut cases = Vec::new();
    for (property, change) in [
        (
            "explicit hydrogen count",
            (|a: &mut Atom| a.hydrogen_count = Some(0)) as fn(&mut Atom),
        ),
        ("atom chirality", |a: &mut Atom| {
            a.chirality = Chirality::Clockwise
        }),
        ("wildcard identity", |a: &mut Atom| a.wildcard = true),
        ("atom aromaticity", |a: &mut Atom| a.aromatic = true),
        ("atom map", |a: &mut Atom| a.atom_map = Some(0)),
        ("CIP descriptor", |a: &mut Atom| {
            a.cip_code = Some(chematic_core::CipCode::R)
        }),
        ("isotope mass zero", |a: &mut Atom| a.isotope = Some(0)),
    ] {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(Atom::new(Element::C));
        let mut atom = Atom::new(Element::N);
        change(&mut atom);
        builder.add_atom(atom);
        let mol = builder.build();
        cases.push((mol, property));
    }
    let mut tagged = pair(B::Single);
    tagged.set_tag(AtomIdx(1), Some(7));
    cases.push((tagged, "atom tag"));
    for (mol, property) in cases {
        assert_eq!(
            write_cjson_checked(&mol, &[]),
            Err(E::AtomProperty {
                atom_index: 1,
                property
            })
        );
    }
}

#[test]
fn refuses_self_bonds_allowed_by_the_core_builder() {
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(Atom::new(Element::C));
    builder.add_bond(atom, atom, B::Single).unwrap();
    assert_eq!(
        write_cjson_checked(&builder.build(), &[]),
        Err(E::SelfBond {
            bond_index: 0,
            atom_index: 0
        })
    );
}

#[test]
fn refuses_stereo_side_channels_even_without_atom_chirality() {
    let mut mol = pair(B::Aromatic);
    mol.set_bond_direction(BondIdx(0), B::Up);
    assert_eq!(
        write_cjson_checked(&mol, &[]),
        Err(E::BondDirection { bond_index: 0 })
    );
    mol.clear_bond_direction(BondIdx(0));
    mol.set_stereo_groups(vec![StereoGroup::new(
        StereoGroupKind::Or(1),
        vec![AtomIdx(0)],
    )]);
    assert_eq!(write_cjson_checked(&mol, &[]), Err(E::StereoGroups));
}

#[test]
fn refuses_incomplete_or_nonfinite_geometry_without_partial_output() {
    let mol = pair(B::Single);
    assert_eq!(
        write_cjson_checked(&mol, &[(0.0, 0.0, 0.0)]),
        Err(E::CoordinateCount {
            expected: 2,
            actual: 1
        })
    );
    assert!(matches!(
        write_cjson_checked(&mol, &[(0.0, 0.0, 0.0); 3]),
        Err(E::CoordinateCount {
            expected: 2,
            actual: 3
        })
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for axis in 0..3 {
            let mut row = [0.0; 3];
            row[axis] = value;
            assert_eq!(
                write_cjson_checked(&mol, &[(0.0, 0.0, 0.0), (row[0], row[1], row[2])]),
                Err(E::NonfiniteCoordinate {
                    atom_index: 1,
                    axis
                })
            );
        }
    }
}
