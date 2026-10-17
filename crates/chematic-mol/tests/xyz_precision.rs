use chematic_core::Element;
use chematic_mol::{XyzAtom, XyzNotation, XyzWriteOptions, parse_xyz, write_xyz_checked};
fn atom(x: f64, y: f64, z: f64) -> XyzAtom {
    XyzAtom {
        element: Element::O,
        x,
        y,
        z,
    }
}

#[test]
fn requested_fixed_format_has_exact_fraction_and_minimum_width() {
    let a = atom(1.2345678901234, -2.3456789012345, -0.0);
    let text = write_xyz_checked(&[a], "original geometry", &XyzWriteOptions::default()).unwrap();
    assert_eq!(
        text,
        "1\noriginal geometry\nO        1.234567890123      -2.345678901234      -0.000000000000\n"
    );
    let frame = parse_xyz(&text).unwrap();
    assert!((frame.atoms[0].x - a.x).abs() <= 5e-13);
    assert!((frame.atoms[0].y - a.y).abs() <= 5.01e-13);
    assert_eq!(frame.atoms[0].z.to_bits(), (-0.0f64).to_bits());
}
#[test]
fn seventeen_significant_scientific_digits_roundtrip_finite_binary_values() {
    let options = XyzWriteOptions {
        notation: XyzNotation::Scientific,
        decimal_places: 16,
        min_width: 0,
        ..Default::default()
    };
    for value in [
        0.0,
        -0.0,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
        1.0 / 3.0,
        1e-200,
        1e200,
    ] {
        let text = write_xyz_checked(&[atom(value, value, value)], "", &options).unwrap();
        let parsed = parse_xyz(&text).unwrap();
        assert_eq!(parsed.atoms[0].x.to_bits(), value.to_bits());
        assert_eq!(parsed.atoms[0].y.to_bits(), value.to_bits());
        assert_eq!(parsed.atoms[0].z.to_bits(), value.to_bits());
    }
}
#[test]
fn large_fixed_values_expand_and_tokens_stay_separate() {
    let options = XyzWriteOptions {
        min_width: 1,
        decimal_places: 17,
        ..Default::default()
    };
    let text = write_xyz_checked(&[atom(f64::MAX, -f64::MAX, 1.0)], "", &options).unwrap();
    assert_eq!(text.lines().nth(2).unwrap().split_whitespace().count(), 4);
    assert_eq!(parse_xyz(&text).unwrap().atoms[0].x, f64::MAX);
}
#[test]
fn rounded_scientific_overflow_is_refused() {
    for sign in [1.0, -1.0] {
        let options = XyzWriteOptions {
            notation: XyzNotation::Scientific,
            decimal_places: 0,
            ..Default::default()
        };
        let e = write_xyz_checked(&[atom(sign * f64::MAX, 0.0, 0.0)], "", &options).unwrap_err();
        assert_eq!(e.atom, Some(0));
        assert!(e.detail.contains("rounded"));
    }
}
#[test]
fn nonfinite_coordinates_are_rejected_in_every_axis() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for xyz in [[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]] {
            assert_eq!(
                write_xyz_checked(
                    &[atom(0.0, 0.0, 0.0), atom(xyz[0], xyz[1], xyz[2])],
                    "",
                    &Default::default()
                )
                .unwrap_err()
                .atom,
                Some(1)
            );
        }
    }
}
#[test]
fn comments_and_options_cannot_inject_frames_or_extxyz() {
    for comment in [
        "bad\nline",
        "bad\rline",
        "\t",
        "Properties=species:S:1:pos:R:3",
        "Lattice=cell",
        "café",
        "\0",
    ] {
        assert!(write_xyz_checked(&[], comment, &Default::default()).is_err());
    }
    assert!(write_xyz_checked(&[], &"x".repeat(4097), &Default::default()).is_err());
    for options in [
        XyzWriteOptions {
            decimal_places: 18,
            ..Default::default()
        },
        XyzWriteOptions {
            min_width: 65,
            ..Default::default()
        },
        XyzWriteOptions {
            max_output_bytes: 0,
            ..Default::default()
        },
        XyzWriteOptions {
            max_output_bytes: (64 << 20) + 1,
            ..Default::default()
        },
    ] {
        assert!(write_xyz_checked(&[], "", &options).is_err());
    }
}
#[test]
fn aggregate_output_cap_includes_padding_headers_and_exact_endpoint() {
    let atoms = [atom(1.0, 2.0, 3.0); 2];
    let expected = write_xyz_checked(&atoms, "header", &Default::default()).unwrap();
    let mut options = XyzWriteOptions {
        max_output_bytes: expected.len(),
        ..Default::default()
    };
    assert_eq!(
        write_xyz_checked(&atoms, "header", &options).unwrap(),
        expected
    );
    options.max_output_bytes -= 1;
    assert!(
        write_xyz_checked(&atoms, "header", &options)
            .unwrap_err()
            .detail
            .contains("byte limit")
    );
    assert_eq!(
        write_xyz_checked(&[], "", &Default::default()).unwrap(),
        "0\n\n"
    );
}
#[test]
fn all_elements_and_explicit_atom_order_survive() {
    let atoms: Vec<_> = (1..=118)
        .map(|n| XyzAtom {
            element: Element::from_atomic_number(n).unwrap(),
            ..atom(f64::from(n), 0.0, 0.0)
        })
        .collect();
    let text = write_xyz_checked(&atoms, "all elements", &Default::default()).unwrap();
    let frame = parse_xyz(&text).unwrap();
    assert_eq!(frame.atoms.len(), 118);
    for (actual, expected) in frame.atoms.iter().zip(atoms) {
        assert_eq!(*actual, expected);
    }
}
#[test]
fn atom_cap_is_checked() {
    let mut atoms = vec![atom(0.0, 0.0, 0.0); 100_000];
    assert!(write_xyz_checked(&atoms, "", &Default::default()).is_ok());
    atoms.push(atom(0.0, 0.0, 0.0));
    assert!(
        write_xyz_checked(&atoms, "", &Default::default())
            .unwrap_err()
            .detail
            .contains("atom count")
    );
}

#[test]
fn deterministic_finite_bit_patterns_roundtrip_in_native_reader() {
    let mut state = 0x123456789abcdef0u64;
    let mut atoms = Vec::new();
    while atoms.len() < 4096 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let value = f64::from_bits(state);
        if value.is_finite() {
            atoms.push(atom(value, -value, value));
        }
    }
    let options = XyzWriteOptions {
        notation: XyzNotation::Scientific,
        decimal_places: 16,
        min_width: 0,
        ..Default::default()
    };
    let text = write_xyz_checked(&atoms, "finite bit patterns", &options).unwrap();
    let frame = parse_xyz(&text).unwrap();
    assert_eq!(frame.atoms.len(), atoms.len());
    for (actual, expected) in frame.atoms.iter().zip(atoms) {
        assert_eq!(actual.x.to_bits(), expected.x.to_bits());
        assert_eq!(actual.y.to_bits(), expected.y.to_bits());
        assert_eq!(actual.z.to_bits(), expected.z.to_bits());
    }
}
