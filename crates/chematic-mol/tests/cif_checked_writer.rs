//! Open Babel #2818: invalid cells must fail before emitting NaN coordinates.
use chematic_mol::{CifError, UnitCell, parse_cif, parse_xyz, write_cif, write_cif_checked};

fn silicon() -> chematic_core::Molecule {
    parse_xyz("1\nsilicon\nSi 0 0 0\n").unwrap().to_molecule()
}

#[test]
fn singular_poscar_equivalent_cell_is_rejected() {
    // Issue #2818: scaled vectors (4,0,0), (0,0,0), (0,0,4).
    let cell = UnitCell {
        a: 4.,
        b: 0.,
        c: 4.,
        ..UnitCell::default()
    };
    assert!(write_cif(&silicon(), &[(0., 0., 0.)], Some(&cell)).contains("NaN"));
    assert!(matches!(
        write_cif_checked(&silicon(), &[(0., 0., 0.)], Some(&cell)),
        Err(CifError::InvalidCellParameters(_))
    ));
}

#[test]
fn invalid_cell_lengths_and_angles_are_rejected() {
    for key in ["a", "b", "c", "alpha", "beta", "gamma"] {
        let values: &[f64] = if key.len() == 1 {
            &[0., -1., f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
        } else {
            &[
                0.,
                -90.,
                180.,
                270.,
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ]
        };
        for &value in values {
            let mut cell = UnitCell::default();
            match key {
                "a" => cell.a = value,
                "b" => cell.b = value,
                "c" => cell.c = value,
                "alpha" => cell.alpha = value,
                "beta" => cell.beta = value,
                _ => cell.gamma = value,
            }
            assert!(
                matches!(cell.validate(), Err(CifError::InvalidCellParameters(_))),
                "{cell:?}"
            );
            assert!(write_cif_checked(&silicon(), &[(0., 0., 0.)], Some(&cell)).is_err());
        }
    }
}

#[test]
fn coplanar_impossible_and_numerically_singular_angles_are_rejected() {
    for (alpha, beta, gamma) in [(60., 60., 120.), (10., 10., 150.), (90., 90., 0.000001)] {
        let cell = UnitCell {
            alpha,
            beta,
            gamma,
            ..UnitCell::default()
        };
        assert!(
            matches!(cell.validate(), Err(CifError::InvalidCellParameters(_))),
            "{cell:?}"
        );
    }
    for length in [1e200, 1e-200] {
        let cell = UnitCell {
            a: length,
            b: length,
            c: length,
            ..UnitCell::default()
        };
        assert!(cell.validate().is_err(), "{cell:?}");
    }
}

#[test]
fn checked_writer_round_trips_orthogonal_triclinic_and_small_cells() {
    let mol = silicon();
    for cell in [
        UnitCell {
            a: 4.,
            b: 6.,
            c: 8.,
            ..UnitCell::default()
        },
        UnitCell {
            a: 8.123456789,
            b: 9.,
            c: 10.,
            alpha: 75.,
            beta: 83.,
            gamma: 64.,
        },
        UnitCell {
            a: 0.00000123456789,
            b: 0.000002,
            c: 0.000003,
            ..UnitCell::default()
        },
    ] {
        let expected_frac = (0.125, -0.75, 1.25);
        let cart = cell.frac_to_cart(expected_frac.0, expected_frac.1, expected_frac.2);
        let out = write_cif_checked(&mol, &[cart], Some(&cell)).unwrap();
        let parsed = parse_cif(&out).unwrap();
        assert_eq!(parsed.cell, Some(cell));
        let actual = parsed.coords[0];
        assert!((actual.0 - cart.0).abs() < 1e-12);
        assert!((actual.1 - cart.1).abs() < 1e-12);
        assert!((actual.2 - cart.2).abs() < 1e-12);
    }
}

#[test]
fn no_cell_keeps_full_precision_cartesian_coordinates() {
    let coords = [(2.123456789123456, -3.000000123456789, 4.75)];
    let text = write_cif_checked(&silicon(), &coords, None).unwrap();
    assert!(text.contains("_atom_site_Cartn_x"));
    let parsed = parse_cif(&text).unwrap();
    assert!(parsed.cell.is_none());
    assert_eq!(parsed.coords, coords);
}

#[test]
fn missing_and_excess_coordinates_are_rejected() {
    for coords in [vec![], vec![(0., 0., 0.); 2]] {
        assert!(matches!(
            write_cif_checked(&silicon(), &coords, None),
            Err(CifError::CoordinateCountMismatch { expected: 1, actual }) if actual == coords.len()
        ));
    }
}

#[test]
fn nonfinite_input_and_conversion_overflow_are_rejected() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for coords in [(value, 0., 0.), (0., value, 0.), (0., 0., value)] {
            for cell in [None, Some(UnitCell::default())] {
                assert!(matches!(
                    write_cif_checked(&silicon(), &[coords], cell.as_ref()),
                    Err(CifError::InvalidCoordinate(_))
                ));
            }
        }
    }
    let cell = UnitCell {
        a: 1e-100,
        b: 1e-100,
        c: 1e-100,
        ..UnitCell::default()
    };
    cell.validate().unwrap();
    assert!(matches!(
        write_cif_checked(&silicon(), &[(1e308, 0., 0.)], Some(&cell)),
        Err(CifError::InvalidCoordinate(_))
    ));
}

#[test]
fn fractional_reader_uses_same_cell_validation() {
    for (a, b, alpha, beta, gamma) in [
        (-4., -4., 90., 90., 90.),
        (4., 4., 270., 90., 90.),
        (4., 4., 60., 60., 120.),
    ] {
        let input = format!(
            "data_x\n_cell_length_a {a}\n_cell_length_b {b}\n_cell_length_c 4\n_cell_angle_alpha {alpha}\n_cell_angle_beta {beta}\n_cell_angle_gamma {gamma}\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nSi 0 0 0\n"
        );
        assert!(matches!(
            parse_cif(&input),
            Err(CifError::InvalidCellParameters(_))
        ));
    }
}
