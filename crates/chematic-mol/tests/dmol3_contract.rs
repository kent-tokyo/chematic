//! Original synthetic MIT OR Apache-2.0 fixtures; see the DMol3 manifest.
use chematic_core::Element;
use chematic_mol::{
    DMOL3_BOHR_ANGSTROM, Dmol3Atom, Dmol3ParseLimits, parse_dmol3, parse_dmol3_with_limits,
    write_dmol3,
};
const TEXT: &str = "$cell vectors\n4 1 0\n-1 5 2\n1 0 -6\n$coordinates\nSi .125 -1.25 2.5\nCl -1.5 .125 .5\n$end\n";

#[test]
fn raw_bohr_basis_orientation_and_projection() {
    let g = diagnostic(parse_dmol3(TEXT)).unwrap();
    assert_eq!(g.source_text, TEXT);
    assert_eq!(g.periodicity(), 3);
    assert_eq!(g.atoms[0].position_bohr, [0.125, -1.25, 2.5]);
    assert_eq!(
        g.positions_angstrom().unwrap()[0],
        [0.125, -1.25, 2.5].map(|v| v * DMOL3_BOHR_ANGSTROM)
    );
    assert_eq!(g.cell_bohr.unwrap()[2], [1., 0., -6.]);
    assert_eq!(
        g.cell_angstrom().unwrap().unwrap()[2],
        [1., 0., -6.].map(|v| v * DMOL3_BOHR_ANGSTROM)
    );
    let back = diagnostic(parse_dmol3(&diagnostic(write_dmol3(&g)).unwrap())).unwrap();
    assert_eq!(back.atoms, g.atoms);
    assert_eq!(back.cell_bohr, g.cell_bohr);
}

#[test]
fn all_canonical_elements_and_nonperiodic_documents() {
    let mut g = diagnostic(parse_dmol3("$coordinates\nH 0 0 0\n$end")).unwrap();
    g.atoms = (1..=118)
        .map(|n| Dmol3Atom {
            element: Element::from_atomic_number(n).unwrap(),
            position_bohr: [f64::from(n), 0., 0.],
        })
        .collect();
    assert_eq!(g.periodicity(), 0);
    assert_eq!(g.cell_angstrom().unwrap(), None);
    let back = diagnostic(parse_dmol3(&diagnostic(write_dmol3(&g)).unwrap())).unwrap();
    assert_eq!(back.atoms, g.atoms);
}

#[test]
fn precision_crlf_blank_lines_and_fortran_exponents() {
    let text = "\r\n$coordinates\r\n\r\nH -0.0 4.9406564584124654D-324 1.7976931348623157e308\r\n$end\r\n\r\n";
    let g = diagnostic(parse_dmol3(text)).unwrap();
    assert_eq!(g.source_text, text);
    let back = diagnostic(parse_dmol3(&diagnostic(write_dmol3(&g)).unwrap())).unwrap();
    for (a, b) in g.atoms[0]
        .position_bohr
        .iter()
        .zip(back.atoms[0].position_bohr)
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert!(
        g.positions_angstrom().unwrap()[0]
            .iter()
            .all(|v| v.is_finite())
    );
}

#[test]
fn malformed_unknown_and_multiple_sections_fail() {
    for text in [
        "",
        "$coordinates\n$end",
        "$coordinates\nC 0 0 0",
        "$coordinates\nC 0 0 0\n$end\n$end",
        "preamble\n$coordinates\nC 0 0 0\n$end",
        "$coordinates angstrom\nC 0 0 0\n$end",
        "$coordinates\nC 0 0 0 0.5\n$end",
        "$coordinates\nc 0 0 0\n$end",
        "$coordinates\nC1 0 0 0\n$end",
        "$coordinates\n6 0 0 0\n$end",
        "$coordinates\nC NaN 0 0\n$end",
        "$coordinates\nC 1e309 0 0\n$end",
        "$cell vectors\n1 0 0\n0 1 0\n$coordinates\nC 0 0 0\n$end",
        "$coordinates\nC 0 0 0\n$cell vectors\n$end",
        "$coordinates\n# comment\nC 0 0 0\n$end",
    ] {
        assert!(diagnostic(parse_dmol3(text)).is_err(), "{text}");
    }
}

#[test]
fn invalid_and_nearly_dependent_cells_are_refused() {
    for cell in [
        "0 0 0\n0 1 0\n0 0 1",
        "1 0 0\n2 0 0\n0 0 1",
        "1 0 0\n0 1 0\n1 1 1e-12",
        "1e308 0 0\n0 1e308 0\n0 0 1e308",
        "NaN 0 0\n0 1 0\n0 0 1",
    ] {
        assert!(
            diagnostic(parse_dmol3(&format!(
                "$cell vectors\n{cell}\n$coordinates\nC 0 0 0\n$end"
            )))
            .is_err()
        );
    }
}

#[test]
fn caller_limits_and_physical_lines() {
    for limits in [
        Dmol3ParseLimits {
            max_input_bytes: 2,
            ..Default::default()
        },
        Dmol3ParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        Dmol3ParseLimits {
            max_line_len: 5,
            ..Default::default()
        },
    ] {
        assert!(diagnostic(parse_dmol3_with_limits(TEXT, &limits)).is_err());
    }
    assert_eq!(
        diagnostic(parse_dmol3("$coordinates\n\nC bad 0 0\n$end"))
            .unwrap_err()
            .line,
        3
    );
    let text = format!("{}$coordinates\nC 0 0 0\n$end", " ".repeat(1025));
    assert!(
        diagnostic(parse_dmol3_with_limits(
            &text,
            &Dmol3ParseLimits {
                max_line_len: usize::MAX,
                ..Default::default()
            }
        ))
        .is_err()
    );
}

#[test]
fn writer_and_projections_validate_mutated_data() {
    let mut g = diagnostic(parse_dmol3(TEXT)).unwrap();
    g.atoms[0].position_bohr[0] = 0.12345678901234566;
    g.cell_bohr = None;
    let text = diagnostic(write_dmol3(&g)).unwrap();
    assert_ne!(text, TEXT);
    let back = diagnostic(parse_dmol3(&text)).unwrap();
    assert_eq!(back.atoms, g.atoms);
    assert!(back.cell_bohr.is_none());
    g.atoms[0].position_bohr[0] = f64::NAN;
    assert!(diagnostic(write_dmol3(&g)).is_err());
    assert!(g.positions_angstrom().is_err());
    assert!(g.cell_angstrom().is_err());
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
