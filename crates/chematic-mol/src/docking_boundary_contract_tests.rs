use crate::mol2_tripos::{Mol2ParseLimits, parse_mol2, parse_mol2_with_limits, write_mol2};
use crate::pdbqt::{
    PdbqtParseLimits, autodock_atom_type, parse_pdbqt, parse_pdbqt_with_limits, write_pdbqt,
};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};

#[test]
fn tripos_atom_and_bond_failures_report_the_record_location() {
    let make = |atom: &str, bond: &str| {
        format!(
            "@<TRIPOS>MOLECULE\nfixture\n2 1 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n{atom}\n2 O 1 0 0 O.2\n@<TRIPOS>BOND\n{bond}\n"
        )
    };
    for atom in [
        "1 C 0 0 0",
        "x C 0 0 0 C.3",
        "1 C x 0 0 C.3",
        "1 C 0 x 0 C.3",
        "1 C 0 0 x C.3",
        "1 Xx 0 0 0 Xx",
    ] {
        let error = parse_mol2(&make(atom, "1 1 2 1"))
            .err()
            .expect("invalid input must fail");
        assert!(error.to_string().contains("MOL2"));
        assert!(error.to_string().contains("7"), "{error}");
    }
    for bond in ["1 1 2", "1 x 2 1", "1 1 x 1", "1 99 2 1", "1 1 99 1"] {
        let error = parse_mol2(&make("1 C 0 0 0 C.3", bond))
            .err()
            .expect("invalid input must fail");
        assert!(error.to_string().contains("10"), "{error}");
    }
    let text = make("1 C 0 0 0 C.3", "1 1 2 1");
    for limits in [
        Mol2ParseLimits {
            max_lines: 1,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_sections: 1,
            ..Default::default()
        },
    ] {
        assert!(
            parse_mol2_with_limits(&text, &limits)
                .err()
                .expect("invalid input must fail")
                .to_string()
                .contains("limit")
        );
    }
    for order in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
        BondOrder::Aromatic,
        BondOrder::Zero,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let o = b.add_atom(Atom::new(Element::O));
        b.add_bond(c, o, order).unwrap();
        let mol = b.build();
        let text = write_mol2(&mol, &[(0., 0., 0.), (1., 2., 3.)]);
        let (read, coords) = parse_mol2(&text).unwrap();
        assert_eq!(coords, vec![(0., 0., 0.), (1., 2., 3.)]);
        let expected = match order {
            BondOrder::QuerySingleOrDouble
            | BondOrder::QuerySingleOrAromatic
            | BondOrder::QueryDoubleOrAromatic => BondOrder::QueryAny,
            _ => order,
        };
        assert_eq!(read.bonds().next().unwrap().1.order, expected);
    }
}

#[test]
fn docking_types_preserve_elements_coordinates_and_partial_charges() {
    for (source, expected) in [
        ("N", "NA"),
        ("[NH4+]", "N"),
        ("S", "SA"),
        ("[S-]", "S"),
        ("P", "P"),
        ("F", "F"),
        ("Cl", "Cl"),
        ("Br", "Br"),
        ("I", "I"),
        ("[Mg+2]", "Mg"),
        ("[Mn+2]", "Mn"),
        ("[Zn+2]", "Zn"),
        ("[Ca+2]", "Ca"),
        ("[Fe+2]", "Fe"),
        ("[H]N", "HD"),
        ("[H]C", "H"),
    ] {
        let mol = chematic_smiles::parse(source).unwrap();
        assert_eq!(autodock_atom_type(&mol, AtomIdx(0)), expected, "{source}");
        let coords = vec![(1.25, -2.5, 3.75); mol.atom_count()];
        let charges = vec![0.125; mol.atom_count()];
        let text = write_pdbqt(&mol, &coords, &charges, "boundary");
        let (read, positions, q) =
            parse_pdbqt(&text).unwrap_or_else(|e| panic!("{source}: {e}: {text}"));
        assert_eq!(read.atom(AtomIdx(0)).element, mol.atom(AtomIdx(0)).element);
        assert_eq!(positions, coords);
        assert_eq!(q, charges);
    }
}

#[test]
fn docking_malformed_numeric_fields_and_limits_are_errors() {
    let mol = chematic_smiles::parse("C").unwrap();
    let text = write_pdbqt(&mol, &[(1., 2., 3.)], &[0.25], "fixture");
    let line = text.lines().find(|l| l.starts_with("ATOM")).unwrap();
    for (offset, width, value) in [
        (30, 8, "bad"),
        (30, 8, "NaN"),
        (38, 8, "inf"),
        (46, 8, "NaN"),
        (70, 6, "bad"),
        (70, 6, "NaN"),
        (77, 2, "Xx"),
    ] {
        let mut changed = line.as_bytes().to_vec();
        let value = format!("{value:>width$}");
        changed[offset..offset + width].copy_from_slice(value.as_bytes());
        let changed = String::from_utf8(changed).unwrap();
        let error = parse_pdbqt(&changed)
            .err()
            .expect("invalid input must fail");
        assert!(!error.to_string().is_empty());
    }
    assert!(
        parse_pdbqt("ATOM")
            .err()
            .expect("invalid input must fail")
            .to_string()
            .contains("short")
    );
    assert!(
        parse_pdbqt_with_limits(
            &text,
            &PdbqtParseLimits {
                max_lines: 1,
                ..Default::default()
            }
        )
        .err()
        .expect("invalid input must fail")
        .to_string()
        .contains("limit")
    );
}
