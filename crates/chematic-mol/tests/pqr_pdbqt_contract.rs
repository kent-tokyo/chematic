use chematic_mol::{PqrPdbqtError, parse_pdbqt, parse_pqr, pqr_to_pdbqt};

const FIXTURE: &str = include_str!("../../../validation/openbabel_issues/pqr/charged-records.pqr");
const TYPES: &[&str] = &["N", "C", "HD", "Zn", "Cl"];

#[test]
fn charges_coords_and_identity_survive_in_source_order() {
    let source = parse_pqr(FIXTURE).unwrap();
    let before = source.clone();
    let report = diagnostic(pqr_to_pdbqt(&source, TYPES, true)).unwrap();
    assert_eq!(source, before);
    let (mol, coords, charges) = parse_pdbqt(&report.pdbqt).unwrap();
    assert_eq!(mol.atom_count(), 5);
    for (i, a) in source.atoms.iter().enumerate() {
        assert!((charges[i] - a.charge).abs() <= 0.0005);
        let c = coords[i];
        for (actual, expected) in [(c.0, a.x), (c.1, a.y), (c.2, a.z)] {
            assert!((actual - expected).abs() <= 0.0005);
        }
    }
    let rows: Vec<_> = report
        .pdbqt
        .lines()
        .filter(|s| s.starts_with("ATOM") || s.starts_with("HETATM"))
        .collect();
    for (i, (row, a)) in rows.iter().zip(&source.atoms).enumerate() {
        assert_eq!(row.len(), 79);
        assert_eq!(row[0..6].trim(), a.group_pdb);
        assert_eq!(row[6..11].trim().parse::<i64>().unwrap(), a.serial);
        assert_eq!(row[12..16].trim(), a.atom_name);
        assert_eq!(row[17..20].trim(), a.res_name);
        assert_eq!(row[21..22].trim(), a.chain_id.as_deref().unwrap_or(""));
        assert_eq!(row[22..26].trim().parse::<i64>().unwrap(), a.res_seq);
        assert_eq!(
            row[26..27].trim(),
            a.icode.map(|c| c.to_string()).unwrap_or_default()
        );
        assert_eq!(row[77..79].trim(), TYPES[i]);
    }
}

#[test]
fn sidecar_reports_original_values_and_measured_rounding() {
    let source = parse_pqr(FIXTURE).unwrap();
    let report = diagnostic(pqr_to_pdbqt(&source, TYPES, true)).unwrap();
    assert!(report.radii_omitted);
    assert_eq!(
        report.source_charges,
        source.atoms.iter().map(|a| a.charge).collect::<Vec<_>>()
    );
    assert_eq!(
        report.source_radii,
        source.atoms.iter().map(|a| a.radius).collect::<Vec<_>>()
    );
    assert_eq!(report.source_coords[0], (-0.966123, 1.523, 1.412));
    assert!((report.max_charge_rounding_error - 0.00045).abs() < 1e-12);
    assert!((report.max_coordinate_rounding_error - 0.000123).abs() < 1e-12);
    assert!(!report.pdbqt.contains("ROOT"));
    assert!(!report.pdbqt.contains("TORSDOF"));
}

#[test]
fn radius_loss_requires_explicit_consent() {
    assert_eq!(
        diagnostic(pqr_to_pdbqt(&parse_pqr(FIXTURE).unwrap(), TYPES, false)).unwrap_err(),
        PqrPdbqtError::RadiusLossRequiresConsent
    );
}

#[test]
fn empty_source_and_type_count_are_errors() {
    let mut source = parse_pqr(FIXTURE).unwrap();
    assert!(matches!(
        diagnostic(pqr_to_pdbqt(&source, &["N"], true)),
        Err(PqrPdbqtError::AtomTypeCount { .. })
    ));
    source.atoms.clear();
    assert_eq!(
        diagnostic(pqr_to_pdbqt(&source, &[], true)).unwrap_err(),
        PqrPdbqtError::EmptyInput
    );
}

#[test]
fn types_must_be_known_and_match_source_elements() {
    let source = parse_pqr(FIXTURE).unwrap();
    for invalid in ["X", "C\n", "HD ", "NAA", "H D", "Na"] {
        let types = [invalid, "C", "HD", "Zn", "Cl"];
        assert!(matches!(
            diagnostic(pqr_to_pdbqt(&source, &types, true)),
            Err(PqrPdbqtError::InvalidAtomType { atom: 0, .. })
        ));
    }
    assert!(matches!(
        diagnostic(pqr_to_pdbqt(&source, &["C", "C", "HD", "Zn", "Cl"], true)),
        Err(PqrPdbqtError::ElementMismatch { atom: 0, .. })
    ));
}

#[test]
fn rejects_silent_metadata_truncation_and_record_injection() {
    let source = parse_pqr(FIXTURE).unwrap();
    for (field, value) in [
        ("atom name", "ABCDE"),
        ("atom name", "C\n"),
        ("residue name", "A1AJ1"),
        ("chain", "AB"),
        ("chain", ""),
        ("record", "MODEL"),
    ] {
        let mut s = source.clone();
        match field {
            "atom name" => s.atoms[0].atom_name = value.into(),
            "residue name" => s.atoms[0].res_name = value.into(),
            "chain" => s.atoms[0].chain_id = Some(value.into()),
            _ => s.atoms[0].group_pdb = value.into(),
        }
        assert_eq!(
            diagnostic(pqr_to_pdbqt(&s, TYPES, true)).unwrap_err(),
            PqrPdbqtError::InvalidField { atom: 0, field }
        );
    }
}

#[test]
fn serial_and_residue_boundaries_are_checked() {
    let source = parse_pqr(FIXTURE).unwrap();
    for serial in [0, -1, 100_000] {
        let mut s = source.clone();
        s.atoms[0].serial = serial;
        assert!(matches!(
            diagnostic(pqr_to_pdbqt(&s, TYPES, true)),
            Err(PqrPdbqtError::InvalidField {
                field: "serial",
                ..
            })
        ));
    }
    let mut duplicate = source.clone();
    duplicate.atoms[1].serial = 7;
    assert!(matches!(
        diagnostic(pqr_to_pdbqt(&duplicate, TYPES, true)),
        Err(PqrPdbqtError::InvalidField {
            atom: 1,
            field: "serial"
        })
    ));
    for seq in [-1000, 10000] {
        let mut s = source.clone();
        s.atoms[0].res_seq = seq;
        assert!(matches!(
            diagnostic(pqr_to_pdbqt(&s, TYPES, true)),
            Err(PqrPdbqtError::InvalidField {
                field: "residue number",
                ..
            })
        ));
    }
    let mut s = source;
    s.atoms[0].serial = 99_999;
    s.atoms[0].res_seq = -999;
    s.atoms[0].atom_name = "ABCD".into();
    s.atoms[0].icode = None;
    assert!(diagnostic(pqr_to_pdbqt(&s, TYPES, true)).is_ok());
}

#[test]
fn rejects_numeric_overflow_including_rounding_across_field_edge() {
    let source = parse_pqr(FIXTURE).unwrap();
    for (field, value) in [
        ("x", 10000.0),
        ("x", 9999.9996),
        ("x", -1000.0),
        ("charge", 100.0),
        ("charge", -10.0),
        ("charge", 99.9996),
    ] {
        let mut s = source.clone();
        if field == "x" {
            s.atoms[0].x = value;
        } else {
            s.atoms[0].charge = value;
        }
        assert_eq!(
            diagnostic(pqr_to_pdbqt(&s, TYPES, true)).unwrap_err(),
            PqrPdbqtError::InvalidField { atom: 0, field }
        );
    }
}

#[test]
fn rejects_directly_constructed_nonfinite_values_and_negative_radii() {
    let source = parse_pqr(FIXTURE).unwrap();
    for field in ["x", "y", "z", "charge", "radius"] {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut s = source.clone();
            match field {
                "x" => s.atoms[0].x = value,
                "y" => s.atoms[0].y = value,
                "z" => s.atoms[0].z = value,
                "charge" => s.atoms[0].charge = value,
                _ => s.atoms[0].radius = value,
            }
            assert_eq!(
                diagnostic(pqr_to_pdbqt(&s, TYPES, true)).unwrap_err(),
                PqrPdbqtError::InvalidField { atom: 0, field }
            );
        }
    }
    let mut s = source;
    s.atoms[0].radius = -1.0;
    assert!(matches!(
        diagnostic(pqr_to_pdbqt(&s, TYPES, true)),
        Err(PqrPdbqtError::InvalidField {
            field: "radius",
            ..
        })
    ));
}

#[test]
fn all_supported_types_roundtrip_with_matching_element() {
    let groups: &[(&str, &[&str])] = &[
        ("C", &["C", "A"]),
        ("N", &["N", "NA"]),
        ("O", &["O", "OA", "OS"]),
        ("S", &["S", "SA"]),
        ("H", &["H", "HD"]),
        ("P", &["P"]),
        ("F", &["F"]),
        ("Cl", &["Cl"]),
        ("Br", &["Br"]),
        ("I", &["I"]),
        ("Mg", &["Mg"]),
        ("Mn", &["Mn"]),
        ("Zn", &["Zn"]),
        ("Ca", &["Ca"]),
        ("Fe", &["Fe"]),
    ];
    let mut source = parse_pqr(FIXTURE).unwrap();
    source.atoms.truncate(1);
    for (element, types) in groups {
        source.atoms[0].element = chematic_core::Element::from_symbol(element).unwrap();
        for t in *types {
            let report = diagnostic(pqr_to_pdbqt(&source, &[t], true)).unwrap();
            let (mol, _, _) = parse_pdbqt(&report.pdbqt).unwrap();
            assert_eq!(
                mol.atoms().next().unwrap().1.element,
                source.atoms[0].element
            );
        }
    }
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
