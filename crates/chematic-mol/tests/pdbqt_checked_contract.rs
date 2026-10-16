use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_mol::{PdbqtWriteError as Error, parse_pdbqt, write_pdbqt, write_pdbqt_checked};
use chematic_smiles::parse;

#[test]
fn ordinary_output_keeps_legacy_bytes_and_fixed_columns() {
    let mol = parse("[H]OC(Cl)Br").unwrap();
    let coords: Vec<_> = (0..mol.atom_count())
        .map(|i| (i as f64, -2.345, 6.789))
        .collect();
    let charges = vec![-0.123; mol.atom_count()];
    let checked = write_pdbqt_checked(&mol, &coords, &charges, "explicit H").unwrap();
    assert_eq!(checked, write_pdbqt(&mol, &coords, &charges, "explicit H"));
    for (i, line) in checked
        .lines()
        .filter(|s| s.starts_with("ATOM"))
        .enumerate()
    {
        assert_eq!(line.len(), 79);
        assert_eq!(line[6..11].trim().parse::<usize>().unwrap(), i + 1);
        assert_eq!(line[30..38].trim().parse::<f64>().unwrap(), coords[i].0);
        assert_eq!(&line[38..46], "  -2.345");
        assert_eq!(&line[46..54], "   6.789");
        assert_eq!(&line[70..76], "-0.123");
        assert_eq!(line[77..79].trim(), ["HD", "OA", "C", "Cl", "Br"][i]);
    }
    let (back, xyz, q) = parse_pdbqt(&checked).unwrap();
    assert_eq!(
        back.atoms().map(|(_, a)| a.element).collect::<Vec<_>>(),
        mol.atoms().map(|(_, a)| a.element).collect::<Vec<_>>()
    );
    assert_eq!(xyz, coords);
    assert_eq!(q, charges);
}

#[test]
fn overflow_is_checked_after_rounding_in_every_numeric_field() {
    let mol = parse("CC").unwrap();
    for (field, width, bad) in [
        ("x", 8, 10000.0),
        ("y", 8, -1000.0),
        ("z", 8, 2261608.729),
        ("x", 8, 9999.9999),
        ("y", 8, -999.9999),
        ("charge", 6, 100.0),
        ("charge", 6, -10.0),
        ("charge", 6, 99.9999),
        ("charge", 6, -9.9999),
    ] {
        let mut xyz = vec![(0.0, 0.0, 0.0); 2];
        let mut q = vec![0.0; 2];
        match field {
            "x" => xyz[1].0 = bad,
            "y" => xyz[1].1 = bad,
            "z" => xyz[1].2 = bad,
            _ => q[1] = bad,
        }
        assert_eq!(
            write_pdbqt_checked(&mol, &xyz, &q, ""),
            Err(Error::FieldOverflow {
                atom: 1,
                field,
                width
            })
        );
    }
}

#[test]
fn representable_edges_and_negative_zero_keep_their_columns() {
    let mol = parse("CCC").unwrap();
    let xyz = vec![
        (9999.999, -999.999, -0.0001),
        (9999.9994, -999.9994, 0.0001),
        (0.0, 0.0, 0.0),
    ];
    let text = write_pdbqt_checked(&mol, &xyz, &[99.999, -9.999, -0.0001], "").unwrap();
    let lines: Vec<_> = text.lines().filter(|s| s.starts_with("ATOM")).collect();
    assert_eq!(&lines[0][30..54], "9999.999-999.999  -0.000");
    assert_eq!(&lines[0][70..76], "99.999");
    assert_eq!(&lines[1][70..76], "-9.999");
    assert_eq!(&lines[2][70..76], "-0.000");
    assert!(lines.iter().all(|s| s.len() == 79));
}

#[test]
fn nonfinite_numbers_never_reach_the_formatter() {
    let mol = parse("C").unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for field in ["x", "y", "z", "charge"] {
            let mut p = (0.0, 0.0, 0.0);
            let mut q = 0.0;
            match field {
                "x" => p.0 = bad,
                "y" => p.1 = bad,
                "z" => p.2 = bad,
                _ => q = bad,
            }
            assert_eq!(
                write_pdbqt_checked(&mol, &[p], &[q], "x"),
                Err(Error::NonFinite { atom: 0, field })
            );
        }
    }
}

#[test]
fn counts_include_explicit_hydrogens_and_do_not_default_or_drop_values() {
    let empty = MoleculeBuilder::new().build();
    assert_eq!(
        write_pdbqt_checked(&empty, &[], &[], ""),
        Err(Error::EmptyMolecule)
    );
    let mol = parse("[H]O[H]").unwrap();
    for values in [0, 1, 2, 4] {
        assert_eq!(
            write_pdbqt_checked(&mol, &vec![(0.0, 0.0, 0.0); values], &[0.0; 3], ""),
            Err(Error::CountMismatch {
                field: "coordinates",
                atoms: 3,
                values
            })
        );
        assert_eq!(
            write_pdbqt_checked(&mol, &[(0.0, 0.0, 0.0); 3], &vec![0.0; values], ""),
            Err(Error::CountMismatch {
                field: "charges",
                atoms: 3,
                values
            })
        );
    }
    assert!(write_pdbqt(&mol, &[], &[], "").contains("ATOM"));
}

#[test]
fn names_cannot_inject_records_and_unknown_types_are_rejected() {
    let mol = parse("C").unwrap();
    for name in ["\nATOM", "a\rb", "\0", "\t", "\u{7f}", "分子"] {
        assert_eq!(
            write_pdbqt_checked(&mol, &[(0.0, 0.0, 0.0)], &[0.0], name),
            Err(Error::InvalidName)
        );
    }
    assert!(write_pdbqt_checked(&mol, &[(0.0, 0.0, 0.0)], &[0.0], &"x".repeat(1024)).is_ok());
    assert_eq!(
        write_pdbqt_checked(&mol, &[(0.0, 0.0, 0.0)], &[0.0], &"x".repeat(1025)),
        Err(Error::InvalidName)
    );
    assert_eq!(
        write_pdbqt_checked(&parse("*").unwrap(), &[(0.0, 0.0, 0.0)], &[0.0], ""),
        Err(Error::WildcardAtom { atom: 0 })
    );
    for (smiles, atomic_number) in [("[B]", 5), ("[Si]", 14)] {
        assert_eq!(
            write_pdbqt_checked(&parse(smiles).unwrap(), &[(0.0, 0.0, 0.0)], &[0.0], ""),
            Err(Error::UnsupportedAtom {
                atom: 0,
                atomic_number
            })
        );
    }
}

#[test]
fn serial_boundary_is_bounded_without_expanding_the_field() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..99_999 {
        b.add_atom(Atom::new(Element::C));
    }
    let mol = b.build();
    let text =
        write_pdbqt_checked(&mol, &vec![(0.0, 0.0, 0.0); 99_999], &vec![0.0; 99_999], "").unwrap();
    let last = text.lines().rev().nth(1).unwrap();
    assert_eq!(&last[6..11], "99999");
    assert_eq!(last.len(), 79);
    let mut b = MoleculeBuilder::new();
    for _ in 0..100_000 {
        b.add_atom(Atom::new(Element::C));
    }
    assert_eq!(
        write_pdbqt_checked(&b.build(), &[], &[], ""),
        Err(Error::AtomCount {
            count: 100_000,
            limit: 99_999
        })
    );
}

#[test]
fn typed_errors_have_field_and_atom_diagnostics() {
    let errors = [
        Error::EmptyMolecule,
        Error::WildcardAtom { atom: 0 },
        Error::AtomCount {
            count: 100_000,
            limit: 99_999,
        },
        Error::CountMismatch {
            field: "coordinates",
            atoms: 2,
            values: 1,
        },
        Error::InvalidName,
        Error::NonFinite {
            atom: 1,
            field: "x",
        },
        Error::FieldOverflow {
            atom: 1,
            field: "charge",
            width: 6,
        },
        Error::UnsupportedAtom {
            atom: 1,
            atomic_number: 5,
        },
    ];
    for e in errors {
        assert!(e.to_string().starts_with("PDBQT"));
    }
}
