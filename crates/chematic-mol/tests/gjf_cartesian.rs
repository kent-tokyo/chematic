use chematic_core::AtomIdx;
use chematic_mol::{
    GaussianError, GaussianParseLimits, parse_gjf, parse_gjf_geometry,
    parse_gjf_geometry_with_limits,
};

fn input(rows: &str) -> String {
    format!("%chk=test.chk\n# HF/STO-3G\n\nCartesian probe\n\n-1 2\n{rows}\n\n")
}

#[test]
fn freeze_codes_do_not_shift_coordinates_or_disappear() {
    let text = input("C1 -1 1.25 2.5 3.75\n8 0 -4.5 5.25 -6.75\nCl12 7.25 -8.5 9.75");
    let g = parse_gjf_geometry(&text).unwrap();
    assert_eq!(
        g.coords,
        [(1.25, 2.5, 3.75), (-4.5, 5.25, -6.75), (7.25, -8.5, 9.75)]
    );
    assert_eq!(g.freeze_codes, [Some(-1), Some(0), None]);
    assert_eq!((g.charge, g.multiplicity), (-1, 2));
    assert_eq!(g.mol.atom_count(), 3);
    assert_eq!(g.mol.bond_count(), 0);
    for (i, symbol) in ["C", "O", "Cl"].into_iter().enumerate() {
        assert_eq!(g.mol.atom(AtomIdx(i as u32)).element.symbol(), symbol);
        assert_eq!(g.mol.atom(AtomIdx(i as u32)).charge, 0);
    }
    let (_, coords, charge, multiplicity) = parse_gjf(&text).unwrap();
    assert_eq!(coords, g.coords);
    assert_eq!((charge, multiplicity), (-1, 2));
}

#[test]
fn commas_comments_and_numeric_values_are_retained() {
    let text = input(
        "C, -1, 1.2345678901234567, -2.5e-20, 3.75 ! frozen\n! ignored comment\nH 0 2.5 3.5 4.5",
    )
    .replace("-1 2\n", "-1,2 ! job\n");
    let g = parse_gjf_geometry(&text).unwrap();
    assert_eq!(
        g.coords,
        [(1.2345678901234567, -2.5e-20, 3.75), (2.5, 3.5, 4.5)]
    );
    assert_eq!(g.freeze_codes, [Some(-1), Some(0)]);
    // This reader returns input numeric units. It must not silently scale Bohr input.
    let bohr = parse_gjf_geometry(&text.replace("# HF/STO-3G", "# HF/STO-3G Units=Bohr")).unwrap();
    assert_eq!(g.coords, bohr.coords);
}

#[test]
fn incomplete_or_extra_rows_cannot_return_partial_geometry() {
    for bad in [
        "H",
        "H 1",
        "H 1 2",
        "H -1 1 2 3 4",
        "H 0.5 1 2 3",
        "H x 1 2 3",
        "H 2147483648 1 2 3",
        "Variables:",
        "Constants:",
        "C 1 1.2 2 90.0 3 180.0",
        "C 0 1 2 3 H",
    ] {
        for rows in [
            format!("{bad}\nC 0 0 0"),
            format!("C 0 0 0\n{bad}"),
            format!("C 0 0 0\n{bad}\nO 1 2 3"),
        ] {
            assert!(
                parse_gjf_geometry(&input(&rows)).is_err(),
                "accepted {rows}"
            );
            assert!(parse_gjf(&input(&rows)).is_err(), "legacy accepted {rows}");
        }
    }
}

#[test]
fn all_three_coordinates_are_checked_even_after_a_freeze_code() {
    for row in [
        "C -1 NaN 2 3",
        "C -1 1 inf 3",
        "C -1 1 2 -inf",
        "C 1 2 NaN",
        "C 0 1 2 1e999",
    ] {
        assert!(
            matches!(
                parse_gjf_geometry(&input(row)),
                Err(GaussianError::NonFiniteCoordinate(_))
            ),
            "{row}"
        );
    }
}

#[test]
fn multiple_jobs_and_fragment_spin_are_refused() {
    let text = input("C 1 2 3");
    assert!(parse_gjf_geometry(&text.replace("-1 2\n", "-1 2 0 1\n")).is_err());
    assert!(parse_gjf_geometry(&text.replace("-1 2\n", "-1 0\n")).is_err());
    for marker in ["--Link1--", "  --link1--  "] {
        assert!(parse_gjf_geometry(&format!("{text}{marker}\n{text}")).is_err());
    }
}

#[test]
fn resource_limits_apply_to_the_new_geometry_reader() {
    let text = input("C 0 1 2 3\nH 4 5 6");
    for limits in [
        GaussianParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        GaussianParseLimits {
            max_input_bytes: text.len() - 1,
            ..Default::default()
        },
        GaussianParseLimits {
            max_sections: 2,
            ..Default::default()
        },
        GaussianParseLimits {
            max_lines: 2,
            ..Default::default()
        },
        GaussianParseLimits {
            max_line_bytes: 2,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            parse_gjf_geometry_with_limits(&text, &limits),
            Err(GaussianError::ResourceLimit { .. })
        ));
    }
    let limits = GaussianParseLimits {
        max_atoms: 2,
        max_input_bytes: text.len(),
        ..Default::default()
    };
    assert_eq!(
        parse_gjf_geometry_with_limits(&text, &limits)
            .unwrap()
            .coords
            .len(),
        2
    );
}
