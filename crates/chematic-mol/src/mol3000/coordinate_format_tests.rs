//! Differential checks against the frozen original complete 3D atom writer.
use super::*;

fn original_atoms(out: &mut String, mol: &Molecule, metadata: &MolMetadata, conformer: &Coords3D) {
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms() {
        let point = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        write!(
            out,
            "M  V30 {} {} {:.4} {:.4} {:.4} {}",
            idx.0 + 1,
            v3000_atom_symbol(mol, idx, atom),
            point.x,
            point.y,
            point.z,
            atom.atom_map.unwrap_or(0)
        )
        .expect("writing to String cannot fail");
        append_atom_properties(out, mol, metadata, idx, atom);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
}

fn compare(mol: &Molecule, metadata: &MolMetadata, conformer: &Coords3D) {
    let mut actual = String::from("retained prefix\n");
    let mut expected = actual.clone();
    append_atoms_3d(&mut actual, mol, metadata, conformer);
    original_atoms(&mut expected, mol, metadata, conformer);
    assert_eq!(actual, expected);
}

#[test]
fn complete_3d_rows_preserve_f64_edges_and_map_boundaries() {
    let mut b = MoleculeBuilder::new();
    for map in [
        None,
        Some(0),
        Some(1),
        Some(9),
        Some(99),
        Some(999),
        Some(1000),
        Some(u16::MAX),
    ] {
        let mut atom = Atom::new(Element::C);
        atom.atom_map = map;
        atom.charge = -1;
        atom.isotope = Some(13);
        atom.hydrogen_count = Some(0);
        b.add_atom(atom);
    }
    let mol = b.build();
    let metadata = MolMetadata {
        v3000_atom_properties: vec![
            (1, "LABEL=\"unicode 日本語\"  X=7".into()),
            (1, "ignored duplicate".into()),
        ],
        ..Default::default()
    };
    let edges = [
        0.0,
        -0.0,
        0.00005,
        -0.00005,
        1.23445,
        -1.23445,
        999.99995,
        -999.99995,
        1e11,
        -1e11,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    for x in edges {
        compare(
            &mol,
            &metadata,
            &Coords3D {
                points: vec![Point3::new(x, -x, x); mol.atom_count()],
            },
        );
    }
    for length in 0..=mol.atom_count() + 1 {
        compare(
            &mol,
            &metadata,
            &Coords3D {
                points: vec![Point3::new(-0.0, 1.0, -2.0); length],
            },
        );
    }
}

#[test]
fn deterministic_float_bits_and_rounding_ties_match_original() {
    let mut b = MoleculeBuilder::new();
    b.add_atom(Atom::new(Element::C));
    let mol = b.build();
    let metadata = MolMetadata::default();
    let mut bits = 0x63a5f19b8c2de407u64;
    for _ in 0..20000 {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        let x = f64::from_bits(bits);
        compare(
            &mol,
            &metadata,
            &Coords3D {
                points: vec![Point3::new(x, -x, x)],
            },
        );
    }
    for i in -10000..10000 {
        let x = (f64::from(i) + 0.5) / 10000.0;
        for value in [x, x.next_up(), x.next_down()] {
            compare(
                &mol,
                &metadata,
                &Coords3D {
                    points: vec![Point3::new(value, -value, value)],
                },
            );
        }
    }
}

#[test]
fn growing_atom_ids_and_pseudoatoms_match_original() {
    let mut b = MoleculeBuilder::new();
    for i in 0..1025 {
        let mut atom = Atom::new(if i % 2 == 0 { Element::C } else { Element::CL });
        atom.wildcard = i % 11 == 0;
        b.add_atom(atom);
    }
    let mut mol = b.build();
    mol.set_r_group(AtomIdx(0), RGroupLabel::numbered(7).unwrap());
    let metadata = MolMetadata::default();
    compare(
        &mol,
        &metadata,
        &Coords3D {
            points: (0..1025)
                .map(|i| Point3::new(f64::from(i) * 1.2345, -1.0, 0.1234))
                .collect(),
        },
    );
    compare(
        &MoleculeBuilder::new().build(),
        &metadata,
        &Coords3D::new_zeroed(0),
    );
}
