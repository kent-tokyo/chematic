use super::*;

fn atomic_fixture() -> &'static str {
    "comment line, ignored\n\
\n\
3 atoms\n\
1 atom types\n\
2 bonds\n\
1 bond types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Masses\n\
\n\
1 12.011\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.5 0.0 0.0\n\
3 1 3.0 0.0 0.0\n\
\n\
Bonds\n\
\n\
1 1 1 2\n\
2 1 2 3\n\
"
}

#[test]
fn parses_atomic_fixture() {
    let d = parse_lammps_data(atomic_fixture(), LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d.atoms.len(), 3);
    assert_eq!(d.bonds.len(), 2);
    assert_eq!(d.masses.len(), 1);
    assert_eq!(d.atoms[0].molecule_id, None);
    assert_eq!(d.atoms[0].charge, None);
    assert_eq!(d.count("atoms"), Some(3));
    assert_eq!(d.count("bond types"), Some(1));
}

#[test]
fn round_trip_atomic_style() {
    let d = parse_lammps_data(atomic_fixture(), LammpsAtomStyle::Atomic).unwrap();
    let written = write_lammps_data(&d);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d, d2);
}

fn charge_fixture() -> String {
    "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
-5.0 5.0 xlo xhi\n\
-5.0 5.0 ylo yhi\n\
-5.0 5.0 zlo zhi\n\
\n\
Atoms # charge\n\
\n\
1 1 -0.5 0.0 0.0 0.0\n\
2 1 0.5 1.0 0.0 0.0\n\
"
    .to_string()
}

#[test]
fn round_trip_charge_style() {
    let d = parse_lammps_data(&charge_fixture(), LammpsAtomStyle::Charge).unwrap();
    assert_eq!(d.atoms[0].charge, Some(-0.5));
    assert_eq!(d.atoms[0].molecule_id, None);
    let written = write_lammps_data(&d);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Charge).unwrap();
    assert_eq!(d, d2);
}

fn molecular_fixture() -> String {
    "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # molecular\n\
\n\
1 1 1 0.0 0.0 0.0\n\
2 1 1 1.0 0.0 0.0\n\
"
    .to_string()
}

#[test]
fn round_trip_molecular_style() {
    let d = parse_lammps_data(&molecular_fixture(), LammpsAtomStyle::Molecular).unwrap();
    assert_eq!(d.atoms[0].molecule_id, Some(1));
    assert_eq!(d.atoms[0].charge, None);
    let written = write_lammps_data(&d);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Molecular).unwrap();
    assert_eq!(d, d2);
}

fn full_fixture_with_images() -> String {
    "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # full\n\
\n\
1 1 1 -0.834 0.0 0.0 0.0 1 0 0\n\
2 1 1 0.417 1.0 0.0 0.0 0 -1 2\n\
"
    .to_string()
}

#[test]
fn round_trip_full_style_with_image_flags() {
    let d = parse_lammps_data(&full_fixture_with_images(), LammpsAtomStyle::Full).unwrap();
    assert_eq!(d.atoms[0].molecule_id, Some(1));
    assert_eq!(d.atoms[0].charge, Some(-0.834));
    assert_eq!(d.atoms[0].image, Some([1, 0, 0]));
    assert_eq!(d.atoms[1].image, Some([0, -1, 2]));
    let written = write_lammps_data(&d);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Full).unwrap();
    assert_eq!(d, d2);
}

#[test]
fn round_trip_triclinic_box() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
1.0 2.0 0.5 xy xz yz\n\
\n\
Atoms # atomic\n\
\n\
1 1 5.0 5.0 5.0\n\
";
    let d = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d.simulation_box.tilt, Some([1.0, 2.0, 0.5]));
    let written = write_lammps_data(&d);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d, d2);
}

#[test]
fn opaque_sections_round_trip_as_exact_fixed_point() {
    let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
1 angle types\n\
1 dihedral types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
\n\
Angles\n\
\n\
1 1 1 2 1 # a comment on an opaque row\n\
\n\
Angle Coeffs\n\
\n\
1 100.0 109.5\n\
";
    let d1 = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(
        d1.unparsed_sections,
        vec![
            (
                "Angles".to_string(),
                "1 1 1 2 1 # a comment on an opaque row".to_string()
            ),
            ("Angle Coeffs".to_string(), "1 100.0 109.5".to_string()),
        ]
    );
    let written = write_lammps_data(&d1);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d2.unparsed_sections, d1.unparsed_sections);
    // Fixed point: parsing the rewritten file again changes nothing further.
    let written2 = write_lammps_data(&d2);
    let d3 = parse_lammps_data(&written2, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d3, d2);
}

#[test]
fn zero_row_opaque_section_is_a_fixed_point() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
\n\
Pair Coeffs\n\
\n\
";
    let d1 = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(
        d1.unparsed_sections,
        vec![("Pair Coeffs".to_string(), String::new())]
    );
    let written = write_lammps_data(&d1);
    let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d2.unparsed_sections, d1.unparsed_sections);
}

#[test]
fn atom_count_mismatch_is_typed_error() {
    let text = "comment\n\
\n\
5 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert_eq!(
        err,
        LammpsDataError::AtomCountMismatch {
            declared: 5,
            actual: 2
        }
    );
}

#[test]
fn malformed_header_is_typed_error() {
    // "not" doesn't parse as numeric, so this reads as "no header
    // lines at all, first section is named 'not a valid header line
    // at all @@@'" -- no box bounds were ever declared, so this
    // surfaces as InvalidBox (missing bound lines), not
    // TruncatedInput.
    let text = "comment\nnot a valid header line at all @@@\n";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert!(matches!(err, LammpsDataError::InvalidBox { .. }));

    // A genuinely malformed header line: a leading number that
    // matches none of the known header shapes (not a `<count>
    // <label>` count line, not a 2-float box-bound pair, not a
    // 3-float tilt line).
    let text2 = "comment\n5 3 3 3 3 3 3\n";
    let err2 = parse_lammps_data(text2, LammpsAtomStyle::Atomic).unwrap_err();
    assert!(matches!(err2, LammpsDataError::MalformedHeader { .. }));
}

#[test]
fn non_finite_value_is_typed_error() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 NaN 0.0 0.0\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert!(matches!(
        err,
        LammpsDataError::NonFiniteValue { section, .. } if section == "Atoms"
    ));
}

#[test]
fn unsupported_atom_style_is_typed_error() {
    let err = parse_lammps_data(
        atomic_fixture(),
        LammpsAtomStyle::Other("sphere".to_string()),
    )
    .unwrap_err();
    assert_eq!(
        err,
        LammpsDataError::UnsupportedAtomStyle {
            style: "sphere".to_string()
        }
    );
}

#[test]
fn duplicate_atom_id_is_typed_error() {
    let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
1 1 1.0 0.0 0.0\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    // Pin the actual source line numbers (line 12: first `1 1 0.0 0.0
    // 0.0` row; line 13: the duplicate `1 1 1.0 0.0 0.0` row) so a
    // regression back to reporting the atom's position instead of its
    // real source line is caught.
    assert_eq!(
        err,
        LammpsDataError::DuplicateAtomId {
            id: 1,
            first_line: 12,
            line: 13,
        }
    );
}

#[test]
fn bond_referencing_unknown_atom_is_typed_error() {
    let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
1 bond types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
\n\
Bonds\n\
\n\
1 1 1 99\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert_eq!(
        err,
        LammpsDataError::BondReferencesUnknownAtom {
            bond_id: 1,
            atom_id: 99
        }
    );
}

#[test]
fn type_labels_framework_is_rejected_not_misparsed() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atom Type Labels\n\
\n\
1 C\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert_eq!(
        err,
        LammpsDataError::TypeLabelsUnsupported {
            section: "Atom Type Labels".to_string()
        }
    );
}

#[test]
fn trailing_comments_are_stripped_on_typed_rows() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Masses\n\
\n\
1 12.011 # carbon\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0 # salt ion\n\
";
    let d = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(d.masses[0].mass, 12.011);
    assert_eq!(d.atoms[0].x, 0.0);
}

#[test]
fn invalid_box_is_typed_error() {
    let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
10.0 0.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
";
    let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
    assert!(matches!(err, LammpsDataError::InvalidBox { .. }));
}

#[test]
fn bounded_parser_rejects_input_and_line_limits() {
    let text = atomic_fixture();
    assert!(matches!(
        parse_lammps_data_with_limits(
            text,
            LammpsAtomStyle::Atomic,
            &LammpsDataParseLimits {
                max_input_bytes: 8,
                ..Default::default()
            }
        ),
        Err(LammpsDataError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));

    let long_line = format!("{}\n{}", "x".repeat(32), text);
    assert!(matches!(
        parse_lammps_data_with_limits(
            &long_line,
            LammpsAtomStyle::Atomic,
            &LammpsDataParseLimits {
                max_line_bytes: 16,
                ..Default::default()
            }
        ),
        Err(LammpsDataError::ResourceLimit {
            resource: "line bytes",
            ..
        })
    ));
}

#[test]
fn bounded_parser_rejects_typed_section_rows() {
    let err = parse_lammps_data_with_limits(
        atomic_fixture(),
        LammpsAtomStyle::Atomic,
        &LammpsDataParseLimits {
            max_atoms: 0,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        err,
        LammpsDataError::ResourceLimit {
            resource: "section rows",
            ..
        }
    ));
}
