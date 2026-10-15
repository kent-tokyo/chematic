//! Synthetic MIT OR Apache-2.0 regression inputs; manifest: validation/openbabel_issues/turbomole/lattice.json.
use chematic_core::Element;
use chematic_mol::{
    TURBOMOLE_BOHR_ANGSTROM as BOHR, TurbomoleParseLimits, parse_turbomole,
    parse_turbomole_with_limits, write_turbomole,
};

const PERIODIC: &str = "$periodic 3\n$lattice angs\n 4 1 0\n -1 5 2\n 1 0 -6\n$coord angs\n 1 -2 3 si f\n -4 5 -6 cl\n$end\n";
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0), "{a} != {b}");
}
fn triple(a: [f64; 3], b: [f64; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        close(a, b);
    }
}

#[test]
fn preserves_oriented_left_handed_lattice_and_atom_order() {
    let g = diagnostic(parse_turbomole(PERIODIC)).unwrap();
    assert_eq!(
        g.lattice_angstrom.unwrap(),
        [[4., 1., 0.], [-1., 5., 2.], [1., 0., -6.]]
    );
    assert_eq!(g.periodicity(), 3);
    assert_eq!(g.atoms[0].element, Element::SI);
    assert!(g.atoms[0].frozen);
    assert!(!g.atoms[1].frozen);
    assert_eq!(g.atoms[1].element, Element::CL);
    assert_eq!(g.atoms[1].position_angstrom, [-4., 5., -6.]);
    assert_eq!(g.source_text, PERIODIC);
}
#[test]
fn write_read_keeps_all_coordinates_flags_and_lattice_components() {
    let g = diagnostic(parse_turbomole(PERIODIC)).unwrap();
    let text = diagnostic(write_turbomole(&g)).unwrap();
    assert!(text.starts_with("$periodic 3\n$lattice\n"));
    let back = diagnostic(parse_turbomole(&text)).unwrap();
    for (a, b) in g.atoms.iter().zip(&back.atoms) {
        assert_eq!(a.element, b.element);
        assert_eq!(a.frozen, b.frozen);
        triple(a.position_angstrom, b.position_angstrom);
    }
    for (a, b) in g
        .lattice_angstrom
        .unwrap()
        .into_iter()
        .zip(back.lattice_angstrom.unwrap())
    {
        triple(a, b);
    }
}
#[test]
fn fractional_uses_rows_without_wrapping_or_basis_rotation() {
    let text = PERIODIC.replace("$coord angs", "$coord frac");
    let g = diagnostic(parse_turbomole(&text)).unwrap();
    triple(g.atoms[0].position_angstrom, [9., -9., -22.]);
    triple(g.atoms[1].position_angstrom, [-27., 21., 46.]);
}
#[test]
fn default_bohr_and_fortran_exponents() {
    let g = diagnostic(parse_turbomole("$coord\n 1D0 -2d0 3e0 si\n$end\n")).unwrap();
    triple(g.atoms[0].position_angstrom, [BOHR, -2. * BOHR, 3. * BOHR]);
    assert_eq!(g.periodicity(), 0);
}
#[test]
fn triclinic_cell_reuses_validated_unit_cell_conversion() {
    let text = "$coord frac\n 0.5 0.25 -0.5 c\n$periodic 3\n$cell angs\n 4 5 6 90 60 90\n$end\n";
    let g = diagnostic(parse_turbomole(text)).unwrap();
    triple(
        g.atoms[0].position_angstrom,
        [0.5, 1.25, -3. * (3.0f64).sqrt() / 2.],
    );
}
#[test]
fn header_order_and_comment_line_endings_are_independent() {
    let text = "# synthetic\r\n$coord angs # units\r\n 0 0 0 SI f # fixed\r\n$lattice angs\r\n 1 0 0\r\n 0 2 0\r\n 0 0 3\r\n$periodic 3\r\n$end\r\n# tail\r\n";
    let g = diagnostic(parse_turbomole(text)).unwrap();
    assert_eq!(g.source_text, text);
    assert!(g.atoms[0].frozen);
    assert_eq!(g.periodicity(), 3);
}
#[test]
fn rejects_ambiguous_or_unsupported_documents() {
    let cases = [
        "$coord\n0 0 0 c\n",
        "$coord\n$end",
        "$end",
        "$coord frac\n0 0 0 c\n$end",
        "$coord file=coord\n0 0 0 c\n$end",
        "$coord\n0 0 0 c\n$end\n$coord\n0 0 0 c\n$end",
        "$coord\n0 0 0 c\n$coord\n0 0 0 c\n$end",
        "$coord\n0 0 0 c\n$periodic 2\n$end",
        "$coord\n0 0 0 c\n$periodic 3\n$end",
        "$coord\n0 0 0 c\n$eht charge=1\n$end",
        "$coord\n0 0 0 c x\n$end",
        "$coord\nNaN 0 0 c\n$end",
        "$coord\n0 0 0 13c\n$end",
        "$coord\n0 0 0 xx\n$end",
        "$coord\n0 0 0 c\n$end stuff",
        "$coord\n0 0 0 c\n$end\ntrailing",
        "$coord\n0 0 0 c\n$periodic\n3\n$end",
        "$coord ang\n0 0 0 c\n$end",
        "$coord\n0 0 0 c extra f\n$end",
        "$coord\n0 0 0 🧪\n$end",
    ];
    for text in cases {
        assert!(
            diagnostic(parse_turbomole(text)).is_err(),
            "accepted {text}"
        );
    }
}
#[test]
fn refuses_singular_near_singular_and_nonfinite_cells() {
    for rows in [
        "1 0 0\n2 0 0\n0 0 1",
        "1 0 0\n1 1e-8 0\n0 0 1",
        "inf 0 0\n0 1 0\n0 0 1",
        "0 0 0\n0 1 0\n0 0 1",
        "1e300 0 0\n0 1e300 0\n0 0 1e300",
    ] {
        let text = format!("$periodic 3\n$lattice angs\n{rows}\n$coord\n0 0 0 c\n$end");
        assert!(diagnostic(parse_turbomole(&text)).is_err(), "{rows}");
    }
}
#[test]
fn refuses_missing_periodicity_and_conflicting_cell_forms() {
    assert!(diagnostic(parse_turbomole(&PERIODIC.replace("$periodic 3\n", ""))).is_err());
    assert!(
        diagnostic(parse_turbomole(
            &PERIODIC.replace("$periodic 3", "$periodic 0")
        ))
        .is_err()
    );
    assert!(
        diagnostic(parse_turbomole(
            &PERIODIC.replace("$end", "$cell angs\n4 5 6 90 90 90\n$end")
        ))
        .is_err()
    );
}
#[test]
fn cell_width_angles_and_lengths_are_checked() {
    for row in [
        "4 5 6 90 90",
        "4 5 6 0 90 90",
        "4 5 6 10 10 170",
        "-4 5 6 90 90 90",
    ] {
        let text = format!("$coord\n0 0 0 c\n$periodic 3\n$cell\n{row}\n$end");
        assert!(diagnostic(parse_turbomole(&text)).is_err());
    }
}
#[test]
fn input_line_and_atom_limits_are_applied() {
    let base = TurbomoleParseLimits::default();
    for limits in [
        TurbomoleParseLimits {
            max_input_bytes: 1,
            ..base
        },
        TurbomoleParseLimits {
            max_line_len: 3,
            ..base
        },
        TurbomoleParseLimits {
            max_atoms: 1,
            ..base
        },
    ] {
        assert!(diagnostic(parse_turbomole_with_limits(PERIODIC, &limits)).is_err());
    }
    assert!(
        diagnostic(parse_turbomole_with_limits(
            PERIODIC,
            &TurbomoleParseLimits {
                max_atoms: 2,
                ..base
            }
        ))
        .is_ok()
    );
}
#[test]
fn writer_validates_mutated_geometry_instead_of_replaying_source() {
    let mut g = diagnostic(parse_turbomole(PERIODIC)).unwrap();
    g.atoms[0].position_angstrom[0] = 8.0;
    close(
        diagnostic(parse_turbomole(&diagnostic(write_turbomole(&g)).unwrap()))
            .unwrap()
            .atoms[0]
            .position_angstrom[0],
        8.,
    );
    g.atoms[0].position_angstrom[0] = f64::NAN;
    assert!(diagnostic(write_turbomole(&g)).is_err());
    g.atoms.clear();
    assert!(diagnostic(write_turbomole(&g)).is_err());
}
#[test]
fn writer_refuses_unrepresentable_bohr_coordinate() {
    let mut g = diagnostic(parse_turbomole(PERIODIC)).unwrap();
    g.atoms[0].position_angstrom[0] = f64::MAX;
    assert!(diagnostic(write_turbomole(&g)).is_err());
}
#[test]
fn row_field_counts_and_external_cell_options_are_checked() {
    for (from, to) in [
        ("4 1 0", "4 1"),
        ("4 1 0", "4 1 0 9"),
        ("$lattice angs", "$lattice file=cell"),
        ("1 0 -6\n", ""),
        ("$periodic 3", "$periodic 3 extra"),
    ] {
        assert!(diagnostic(parse_turbomole(&PERIODIC.replace(from, to))).is_err());
    }
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
