use chematic_core::Element;
use chematic_mol::{GjfCoordinateUnits, GjfWriteOptions, parse_gjf_geometry, write_gjf_cartesian};

#[test]
fn geometry_rewrite_retains_numbers_freeze_codes_and_job_state() {
    let input = "# HF/STO-3G\n\noriginal\n\n-1 2\nC -1 1.2345678901234567 -2.5e-20 3.75\nO 0 -4.5 5.25 -6.75\nH 7.25 -8.5 9.75\n\n";
    let g = parse_gjf_geometry(input).unwrap();
    let elements: Vec<_> = g.mol.atoms().map(|(_, a)| a.element).collect();
    let coords: Vec<_> = g.coords.iter().map(|&(x, y, z)| [x, y, z]).collect();
    for (units, keyword) in [
        (GjfCoordinateUnits::Angstrom, "Ang"),
        (GjfCoordinateUnits::Bohr, "AU"),
    ] {
        let mut options = GjfWriteOptions::new(g.charge, g.multiplicity, units);
        options.route = "#p HF/STO-3G".into();
        options.title = "明示座標 geometry".into();
        let output = write_gjf_cartesian(&elements, &coords, &g.freeze_codes, &options).unwrap();
        assert!(output.starts_with(&format!(
            "#p HF/STO-3G Units={keyword}\n\n明示座標 geometry\n\n-1 2\n"
        )));
        let back = parse_gjf_geometry(&output).unwrap();
        assert_eq!(back.coords, g.coords);
        assert_eq!(back.freeze_codes, g.freeze_codes);
        assert_eq!((back.charge, back.multiplicity), (g.charge, g.multiplicity));
    }
}

#[test]
fn finite_f64_values_survive_bitwise_including_signed_zero_and_subnormals() {
    let mut bits = 0x9876_5432_1234_5678u64;
    let mut values = vec![
        0.0,
        -0.0,
        f64::MAX,
        -f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -f64::from_bits(1),
        1e-200,
        1.2345678901234567,
    ];
    for _ in 0..4096 {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        let value = f64::from_bits(bits);
        if value.is_finite() {
            values.push(value);
        }
    }
    let coords: Vec<_> = values.iter().map(|&v| [v, -v, 0.0]).collect();
    let output = write_gjf_cartesian(
        &vec![Element::C; coords.len()],
        &coords,
        &vec![None; coords.len()],
        &GjfWriteOptions::new(0, 1, GjfCoordinateUnits::Angstrom),
    )
    .unwrap();
    let back = parse_gjf_geometry(&output).unwrap();
    for (&(x, y, z), position) in back.coords.iter().zip(coords) {
        assert_eq!(
            [x.to_bits(), y.to_bits(), z.to_bits()],
            position.map(f64::to_bits)
        );
    }
}

#[test]
fn mismatched_arrays_and_nonfinite_coordinates_never_produce_output() {
    let options = GjfWriteOptions::new(0, 1, GjfCoordinateUnits::Angstrom);
    for (coords, codes) in [
        (vec![], vec![None]),
        (vec![[0.; 3]; 2], vec![None]),
        (vec![[0.; 3]], vec![]),
        (vec![[0.; 3]], vec![None; 2]),
    ] {
        assert!(write_gjf_cartesian(&[Element::C], &coords, &codes, &options).is_err());
    }
    assert!(write_gjf_cartesian(&[], &[], &[], &options).is_err());
    for nonfinite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for axis in 0..3 {
            let mut coords = [[0.; 3]];
            coords[0][axis] = nonfinite;
            assert!(write_gjf_cartesian(&[Element::C], &coords, &[Some(-1)], &options).is_err());
        }
    }
}

#[test]
fn headers_cannot_insert_records_comments_or_conflicting_units() {
    for field in [
        "",
        " ",
        "title\n0 1",
        "title\r0 1",
        "title\0",
        "a\u{2028}b",
        "a\u{2029}b",
    ] {
        let mut options = GjfWriteOptions::new(0, 1, GjfCoordinateUnits::Angstrom);
        options.title = field.into();
        assert!(write_gjf_cartesian(&[Element::C], &[[0.; 3]], &[None], &options).is_err());
    }
    for route in [
        "HF/STO-3G",
        "# HF\n\n0 1",
        "# ! swallowed",
        "## HF",
        "# Units=AU",
        "# uNiTs(Bohr)",
        "# UN=AU",
        "# U=AU",
        "# HF 漢",
    ] {
        let mut options = GjfWriteOptions::new(0, 1, GjfCoordinateUnits::Angstrom);
        options.route = route.into();
        assert!(
            write_gjf_cartesian(&[Element::C], &[[0.; 3]], &[None], &options).is_err(),
            "{route}"
        );
    }
    let mut options = GjfWriteOptions::new(0, 0, GjfCoordinateUnits::Angstrom);
    assert!(write_gjf_cartesian(&[Element::C], &[[0.; 3]], &[None], &options).is_err());
    options.multiplicity = 1;
    options.title = "x".repeat(1025);
    assert!(write_gjf_cartesian(&[Element::C], &[[0.; 3]], &[None], &options).is_err());
}

#[test]
fn element_profile_and_atom_cap_are_enforced() {
    let options = GjfWriteOptions::new(0, 1, GjfCoordinateUnits::Angstrom);
    let elements: Vec<_> = (1..=118)
        .map(|n| Element::from_atomic_number(n).unwrap())
        .collect();
    let output =
        write_gjf_cartesian(&elements, &vec![[0.; 3]; 118], &vec![None; 118], &options).unwrap();
    assert_eq!(parse_gjf_geometry(&output).unwrap().mol.atom_count(), 118);
    if let Some(dummy) = Element::from_atomic_number(0) {
        assert!(write_gjf_cartesian(&[dummy], &[[0.; 3]], &[None], &options).is_err());
    }
    let elements = vec![Element::C; GjfWriteOptions::MAX_ATOMS];
    let coords = vec![[f64::MAX, -f64::MAX, f64::from_bits(1)]; elements.len()];
    let codes = vec![Some(i32::MIN); elements.len()];
    assert!(
        write_gjf_cartesian(&elements, &coords, &codes, &options)
            .unwrap()
            .len()
            < 16 * 1024 * 1024
    );
    assert!(
        write_gjf_cartesian(
            &vec![Element::C; GjfWriteOptions::MAX_ATOMS + 1],
            &[],
            &[],
            &options
        )
        .is_err()
    );
}
