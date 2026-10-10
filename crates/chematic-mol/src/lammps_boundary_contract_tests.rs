use crate::*;

fn data(section: &str, rows: &str) -> String {
    format!(
        "contract\n\n2 atoms\n1 atom types\n1 bonds\n1 bond types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n2 1 1 0 0\n\n{section}\n\n{rows}\n"
    )
}

#[test]
fn malformed_lammps_sections_have_useful_diagnostics() {
    for (section, rows, diagnostic) in [
        ("Masses", "1", "row"),
        ("Masses", "1 nope", "row"),
        ("Masses", "1 NaN", "finite"),
        ("Masses", "1.5 12", "row"),
        ("Velocities", "1 0 0", "row"),
        ("Velocities", "1 0 x 0", "row"),
        ("Velocities", "1 0 inf 0", "finite"),
        ("Bonds", "1 1 1", "row"),
        ("Bonds", "1 1 1 x", "row"),
        ("Bonds", "1 1 3 2", "atom"),
        ("Bonds", "1 1 1 3", "atom"),
        ("Atom Type Labels", "1 carbon", "label"),
    ] {
        let error = parse_lammps_data(&data(section, rows), LammpsAtomStyle::Atomic)
            .err()
            .unwrap();
        assert!(
            error.to_string().to_lowercase().contains(diagnostic),
            "{section} {rows}: {error}"
        );
    }
    let valid = data("Bonds", "1 1 1 2");
    for (from, to) in [
        ("2 atoms", "3 atoms"),
        ("2 atoms", "-1 atoms"),
        ("2 atoms\n", ""),
        ("0 10 xlo xhi", "NaN 10 xlo xhi"),
        ("0 10 xlo xhi", "10 0 xlo xhi"),
        ("1 1 0 0 0", "1 1 0 0"),
        ("2 1 1 0 0", "1 1 1 0 0"),
        ("2 1 1 0 0", "2 1 NaN 0 0"),
    ] {
        let error = parse_lammps_data(&valid.replace(from, to), LammpsAtomStyle::Atomic)
            .err()
            .unwrap();
        assert!(!error.to_string().is_empty());
    }
    assert!(matches!(
        parse_lammps_data(&valid, LammpsAtomStyle::Other("hybrid".into())),
        Err(LammpsDataError::UnsupportedAtomStyle { .. })
    ));
}

#[test]
fn all_lammps_resource_limits_apply_to_their_sections() {
    let text = format!(
        "{}\nMasses\n\n1 12.0\n\nVelocities\n\n1 0 0 0\n\nPair Coeffs\n\n1 0.5 1.2\n",
        data("Bonds", "1 1 1 2")
    );
    for limits in [
        LammpsDataParseLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_line_bytes: 1,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_header_counts: 1,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_masses: 0,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_velocities: 0,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_opaque_section_bytes: 1,
            ..Default::default()
        },
        LammpsDataParseLimits {
            max_sections: 1,
            ..Default::default()
        },
    ] {
        let error = parse_lammps_data_with_limits(&text, LammpsAtomStyle::Atomic, &limits)
            .err()
            .unwrap();
        assert!(
            matches!(error, LammpsDataError::ResourceLimit { .. }),
            "{error}"
        );
        assert!(error.to_string().contains("limit"));
    }
}

#[test]
fn lammps_boxes_reject_nonfinite_bounds_and_tilts() {
    for axis in 0..3 {
        for (lo, hi, diagnostic) in [
            (f64::NAN, 10., "finite"),
            (0., f64::INFINITY, "finite"),
            (10., 0., "greater"),
        ] {
            let mut b = LammpsBox {
                lo: [0.; 3],
                hi: [10.; 3],
                tilt: None,
            };
            b.lo[axis] = lo;
            b.hi[axis] = hi;
            assert!(b.validate().unwrap_err().contains(diagnostic));
        }
        let mut tilt = [0.; 3];
        tilt[axis] = f64::NAN;
        assert!(
            LammpsBox {
                lo: [0.; 3],
                hi: [10.; 3],
                tilt: Some(tilt)
            }
            .validate()
            .unwrap_err()
            .contains("tilt")
        );
    }
}
