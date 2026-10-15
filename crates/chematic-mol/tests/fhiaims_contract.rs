//! Original MIT OR Apache-2.0 fixtures; provenance/hashes are in the oracle manifest.
use chematic_mol::{
    FhiaimsCoordinates, FhiaimsParseLimits, parse_fhiaims, parse_fhiaims_with_limits, write_fhiaims,
};
const TEXT: &str = "# original\nlattice_vector 4 1 0\nconstrain_relaxation x\nlattice_vector -1 5 2\nlattice_vector 1 0 -6\natom_frac 0.5 0.25 -0.5 Si\nconstrain_relaxation y\ninitial_charge 0.25\ninitial_moment -1.5\natom 1 2 3 Cl\n";

#[test]
fn fractional_basis_and_metadata_are_preserved() {
    let g = diagnostic(parse_fhiaims(TEXT)).unwrap();
    assert_eq!(g.source_text, TEXT);
    assert_eq!(
        g.positions_angstrom().unwrap(),
        [[1.25, 1.75, 3.5], [1., 2., 3.]]
    );
    assert_eq!(
        g.atoms[0].coordinates,
        FhiaimsCoordinates::Fractional([0.5, 0.25, -0.5])
    );
    assert_eq!(g.atoms[0].initial_charge, Some(0.25));
    assert_eq!(g.atoms[0].initial_moment, Some(-1.5));
    assert_eq!(g.atoms[0].constraints, [false, true, false]);
    assert_eq!(g.lattice_vectors[0].constraints, [true, false, false]);
    let back = diagnostic(parse_fhiaims(&diagnostic(write_fhiaims(&g)).unwrap())).unwrap();
    assert_eq!(g.atoms, back.atoms);
    assert_eq!(g.lattice_vectors, back.lattice_vectors);
}

#[test]
fn zero_to_three_periodic_directions_keep_exact_basis() {
    for n in 0..=3 {
        let mut text = String::new();
        for row in ["2 1 0", "0 3 1", "0 0 -4"].iter().take(n) {
            text.push_str(&format!("lattice_vector {row}\n"));
        }
        text.push_str("atom 1 2 3 C\n");
        let g = diagnostic(parse_fhiaims(&text)).unwrap();
        assert_eq!(g.periodicity(), n);
        assert_eq!(
            g.lattice_vectors,
            diagnostic(parse_fhiaims(&diagnostic(write_fhiaims(&g)).unwrap()))
                .unwrap()
                .lattice_vectors
        );
    }
}

#[test]
fn all_constraint_masks_and_zero_vs_absent_initial_values() {
    for mask in 0..8 {
        let mut g = diagnostic(parse_fhiaims("atom 0 0 0 C\ninitial_charge 0\n")).unwrap();
        g.atoms[0].constraints = [mask & 1 != 0, mask & 2 != 0, mask & 4 != 0];
        let back = diagnostic(parse_fhiaims(&diagnostic(write_fhiaims(&g)).unwrap())).unwrap();
        assert_eq!(g.atoms, back.atoms);
        assert_eq!(back.atoms[0].initial_charge, Some(0.));
        assert_eq!(back.atoms[0].initial_moment, None);
    }
}

#[test]
fn numbers_comments_crlf_and_out_of_order_lattice() {
    let text = "atom_frac 0.5 0.25 -0.5 Si\r\ninitial_charge 2.5D-1 # inline\r\nlattice_vector 4 1 0\r\nlattice_vector -1 5 2\r\nlattice_vector 1 0 -6\r\n";
    let g = diagnostic(parse_fhiaims(text)).unwrap();
    assert_eq!(g.positions_angstrom().unwrap(), [[1.25, 1.75, 3.5]]);
    assert_eq!(g.atoms[0].initial_charge, Some(0.25));
    let g = diagnostic(parse_fhiaims(
        "atom -0.0 4.9406564584124654e-324 1.7976931348623157e308 H",
    ))
    .unwrap();
    let back = diagnostic(parse_fhiaims(&diagnostic(write_fhiaims(&g)).unwrap())).unwrap();
    assert_eq!(g.atoms, back.atoms);
    assert_eq!(
        back.positions_angstrom().unwrap()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn malformed_unsupported_and_detached_metadata_fail() {
    for text in [
        "",
        "atom 0 0 C",
        "atom 0 0 0 fe",
        "atom 0 0 0 Fe_custom",
        "atom NaN 0 0 C",
        "empty 0 0 0 C",
        "atom 0 0 0 C\nvelocity 1 2 3",
        "initial_charge 0.25\natom 0 0 0 C",
        "lattice_vector 1 0 0\ninitial_moment 1\natom 0 0 0 C",
        "constrain_relaxation x\natom 0 0 0 C",
        "atom 0 0 0 C\nconstrain_relaxation x y",
        "atom 0 0 0 C\nconstrain_relaxation x\nconstrain_relaxation .false.",
        "atom 0 0 0 C\ninitial_charge 0\ninitial_charge 1",
        "atom_frac 0 0 0 C",
        "atom 0 0 0 C\ninitial_moment inf",
    ] {
        assert!(diagnostic(parse_fhiaims(text)).is_err(), "{text}");
    }
}

#[test]
fn invalid_cell_and_fractional_overflow_fail() {
    for text in [
        "lattice_vector 0 0 0\natom 0 0 0 C",
        "lattice_vector 1 0 0\nlattice_vector 2 0 0\natom 0 0 0 C",
        "lattice_vector 1 0 0\nlattice_vector 0 1 0\nlattice_vector 1 1 1e-12\natom 0 0 0 C",
        "lattice_vector 1 0 0\nlattice_vector 0 1 0\natom_frac 0 0 0 C",
        "lattice_vector 2 0 0\nlattice_vector 0 2 0\nlattice_vector 0 0 2\natom_frac 1e308 0 0 C",
        "lattice_vector 1 0 0\nlattice_vector 0 1 0\nlattice_vector 0 0 1\nlattice_vector 1 1 1\natom 0 0 0 C",
    ] {
        assert!(diagnostic(parse_fhiaims(text)).is_err(), "{text}");
    }
}

#[test]
fn limits_and_source_lines_are_explicit() {
    for limits in [
        FhiaimsParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        FhiaimsParseLimits {
            max_input_bytes: 2,
            ..Default::default()
        },
        FhiaimsParseLimits {
            max_line_len: 5,
            ..Default::default()
        },
    ] {
        assert!(diagnostic(parse_fhiaims_with_limits(TEXT, &limits)).is_err());
    }
    assert_eq!(
        diagnostic(parse_fhiaims("atom 0 0 0 C\ninitial_charge NaN"))
            .unwrap_err()
            .line,
        2
    );
}

#[test]
fn writing_current_data_validates_and_never_replays_source() {
    let mut g = diagnostic(parse_fhiaims(TEXT)).unwrap();
    g.atoms[0].coordinates = FhiaimsCoordinates::Cartesian([0.125, 2., 3.]);
    g.atoms[0].initial_charge = None;
    let back = diagnostic(parse_fhiaims(&diagnostic(write_fhiaims(&g)).unwrap())).unwrap();
    assert_eq!(g.atoms, back.atoms);
    g.atoms[0].initial_moment = Some(f64::NAN);
    assert!(diagnostic(write_fhiaims(&g)).is_err());
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
