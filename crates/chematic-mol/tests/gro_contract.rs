use chematic_core::Element;
use chematic_mol::{GroError, GroParseLimits, parse_gro, parse_gro_with_limits};
const FIXTURE: &str =
    include_str!("../../../validation/openbabel_issues/gro/two-letter-elements.gro");
fn row(name: &str, precision: usize, velocity: bool) -> String {
    let width = precision + 5;
    let mut row = format!("{:5}{:<5}{:>5}{:5}", 99999, "ALA", name, 0);
    for value in [-12.25, 0.5, 4.75] {
        row.push_str(&format!("{value:>width$.precision$}"));
    }
    if velocity {
        let precision = precision + 1;
        for value in [-0.125, 0.25, 0.5] {
            row.push_str(&format!("{value:>width$.precision$}"));
        }
    }
    row
}
fn frame(row: &str, box_line: &str) -> String {
    format!("synthetic\n1\n{row}\n{box_line}\n")
}
#[test]
fn names_are_preserved_without_first_letter_element_inference() {
    let d = diagnostic(parse_gro(FIXTURE)).unwrap();
    assert_eq!(
        d.atoms
            .iter()
            .map(|a| a.atom_name.as_str())
            .collect::<Vec<_>>(),
        vec!["BE", "CE", "FE"]
    );
    assert!(d.atoms.iter().all(|a| a.element.is_none()));
    assert_eq!(d.atoms[1].position_nm, [2.0, 0.0, 0.0]);
    assert_eq!(d.source_text, FIXTURE);
}
#[test]
fn supplied_elements_resolve_two_letter_species_on_independent_copy() {
    let d = diagnostic(parse_gro(FIXTURE)).unwrap();
    let before = d.clone();
    let resolved = d
        .with_elements(&[Element::BE, Element::CE, Element::FE])
        .unwrap();
    assert_eq!(
        resolved
            .atoms
            .iter()
            .map(|a| a.element.unwrap().symbol())
            .collect::<Vec<_>>(),
        vec!["Be", "Ce", "Fe"]
    );
    assert_eq!(d, before);
    assert_eq!(resolved.source_text, d.source_text);
    assert_eq!(resolved.box_vectors_nm, d.box_vectors_nm);
    assert!(matches!(
        d.with_elements(&[Element::B]),
        Err(GroError::ElementCount {
            atoms: 3,
            elements: 1
        })
    ));
}
#[test]
fn ca_is_not_silently_assumed_carbon_or_calcium() {
    let d = diagnostic(parse_gro(&frame(&row("CA", 3, false), "1 1 1"))).unwrap();
    assert_eq!(d.atoms[0].element, None);
    assert_eq!(
        d.with_elements(&[Element::C]).unwrap().atoms[0].element,
        Some(Element::C)
    );
    assert_eq!(
        d.with_elements(&[Element::CA]).unwrap().atoms[0].element,
        Some(Element::CA)
    );
}
#[test]
fn positive_decimal_precision_is_inferred_from_adjacent_decimal_points() {
    for precision in 2..=12 {
        let d = diagnostic(parse_gro(&frame(&row("CA", precision, false), "1 2 3"))).unwrap();
        assert_eq!(d.atoms[0].position_nm, [-12.25, 0.5, 4.75]);
        assert_eq!(d.atoms[0].velocity_nm_per_ps, None);
        assert_eq!(d.atoms[0].residue_number, 99999);
        assert_eq!(d.atoms[0].atom_number, 0);
    }
}
#[test]
fn optional_velocities_keep_nm_per_ps_without_unit_rescaling() {
    for precision in 2..=8 {
        let d = diagnostic(parse_gro(&frame(&row("OW", precision, true), "0 0 0"))).unwrap();
        let velocity = d.atoms[0].velocity_nm_per_ps.unwrap();
        assert_eq!(velocity, [-0.125, 0.25, 0.5]);
    }
}
#[test]
fn nine_box_values_map_to_row_vectors_in_documented_order() {
    let d = diagnostic(parse_gro(&frame(
        &row("CA", 3, false),
        "1 2 3 0 0 0.4 0 0.5 0.6",
    )))
    .unwrap();
    assert_eq!(
        d.box_vectors_nm,
        [[1.0, 0.0, 0.0], [0.4, 2.0, 0.0], [0.5, 0.6, 3.0]]
    );
    assert_eq!(d.box_values_nm.len(), 9);
    let d = diagnostic(parse_gro(&frame(&row("CA", 3, false), "1 2 3"))).unwrap();
    assert_eq!(
        d.box_vectors_nm,
        [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    );
    assert_eq!(d.box_values_nm.len(), 3);
}
#[test]
fn zero_box_does_not_invent_a_periodic_cell_and_empty_frame_is_retained() {
    let text = "\n0\n0 0 0\n";
    let d = diagnostic(parse_gro(text)).unwrap();
    assert!(d.atoms.is_empty());
    assert_eq!(d.title, "");
    assert_eq!(d.box_vectors_nm, [[0.0; 3]; 3]);
    assert_eq!(d.with_elements(&[]).unwrap(), d);
}
#[test]
fn crlf_unicode_title_and_trailing_blank_lines_keep_exact_source() {
    let text = frame(&row("CA", 3, false), "1 1 1")
        .replace("synthetic", "構造 t= 2.5")
        .replace('\n', "\r\n")
        + "\r\n";
    let d = diagnostic(parse_gro(&text)).unwrap();
    assert_eq!(d.source_text, text);
    assert_eq!(d.title, "構造 t= 2.5");
}
#[test]
fn wrapped_identifiers_are_literal_not_unique_row_keys() {
    let text = format!(
        "ids\n2\n{}\n{}\n1 1 1\n",
        row("CA", 3, false),
        row("CA", 3, false)
    );
    let d = diagnostic(parse_gro(&text)).unwrap();
    assert_eq!(d.atoms.len(), 2);
    assert_eq!(d.atoms[0], d.atoms[1]);
}
#[test]
fn additional_frames_and_nonblank_trailing_data_are_rejected() {
    for tail in [FIXTURE, "extra\n"] {
        assert!(matches!(
            diagnostic(parse_gro(&(FIXTURE.to_owned() + tail))),
            Err(GroError::TrailingData { .. })
        ));
    }
}
#[test]
fn malformed_short_nonascii_or_name_fields_fail_without_panics() {
    let good = row("CA", 3, false);
    let variants = [
        "short".into(),
        good.replacen("ALA", "原子", 1),
        good.replacen("   CA", "     ", 1),
        good.replacen("99999", "-0001", 1),
        good.replacen("    0", "   -1", 1),
    ];
    for r in variants {
        assert!(diagnostic(parse_gro(&frame(&r, "1 1 1"))).is_err(), "{r}");
    }
}
#[test]
fn malformed_decimal_precision_partial_velocities_and_bad_box_fail() {
    let good = row("CA", 3, true);
    for r in [
        good[..good.len() - 8].into(),
        good.replacen("-12.250", "-12e250", 1),
        good.replacen("0.500", "0..00", 1),
    ] {
        assert!(diagnostic(parse_gro(&frame(&r, "1 1 1"))).is_err(), "{r}");
    }
    for b in ["1 2", "1 2 3 4", "NaN 1 1", "1 inf 1", "1e999 1 1"] {
        assert!(diagnostic(parse_gro(&frame(&row("CA", 3, false), b))).is_err());
    }
}
#[test]
fn declared_counts_and_missing_lines_fail_before_large_reservations() {
    for text in [
        "",
        "title\n",
        "title\n-1\n",
        "title\nabc\n",
        "title\n2000000\n0 0 0\n",
        "title\n1\n",
    ] {
        assert!(diagnostic(parse_gro(text)).is_err());
    }
}
#[test]
fn input_atom_and_line_limits_are_enforced() {
    let limits = GroParseLimits {
        max_input_bytes: 8,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(parse_gro_with_limits(FIXTURE, &limits)),
        Err(GroError::InputTooLarge { .. })
    ));
    let limits = GroParseLimits {
        max_atoms: 2,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(parse_gro_with_limits(FIXTURE, &limits)),
        Err(GroError::TooManyAtoms { .. })
    ));
    let limits = GroParseLimits {
        max_line_len: 8,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(parse_gro_with_limits(FIXTURE, &limits)),
        Err(GroError::LineTooLong { line: 1, .. })
    ));
}

#[test]
fn one_decimal_coordinate_fields_are_supported() {
    let r = format!("{:5}{:<5}{:>5}{:5}  -1.2   0.5   4.8", 1, "ALA", "CA", 1);
    let d = diagnostic(parse_gro(&frame(&r, "1 1 1"))).unwrap();
    assert_eq!(d.atoms[0].position_nm, [-1.2, 0.5, 4.8]);
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
