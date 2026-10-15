use chematic_mol::{MmcifPdbError, PdbNamespace, parse_mmcif_residues, write_mmcif_pdb_checked};
const GLYCAN: &str =
    include_str!("../../../validation/openbabel_issues/mmcif/branched-residues.cif");
const LONG: &str =
    include_str!("../../../validation/openbabel_issues/mmcif/long-component-name.cif");

#[test]
fn refuses_five_character_names_even_with_metadata_loss_consent() {
    let d = diagnostic(parse_mmcif_residues(LONG)).unwrap();
    for allow in [false, true] {
        assert!(matches!(
            diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, allow)),
            Err(MmcifPdbError::InvalidField {
                field: "residue name",
                ..
            })
        ));
    }
}
#[test]
fn preserves_three_glycan_residues_serials_elements_and_positions() {
    let d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    let before = d.clone();
    let r = diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).unwrap();
    assert_eq!(d, before);
    let rows: Vec<_> = r.pdb.lines().filter(|s| s.starts_with("HETATM")).collect();
    assert_eq!(rows.len(), 6);
    for (i, (line, atom)) in rows.iter().zip(&d.parsed.atoms).enumerate() {
        assert_eq!(line.len(), 80);
        assert_eq!(line[6..11].trim().parse::<i64>().unwrap(), atom.serial);
        assert_eq!(line[12..16].trim(), atom.atom_name);
        assert_eq!(line[17..20].trim(), atom.res_name);
        assert_eq!(line[21..22].trim(), "A");
        assert_eq!(line[22..26].trim().parse::<usize>().unwrap(), i / 2 + 1);
        assert_eq!(line[30..38].trim().parse::<f64>().unwrap(), atom.x);
        assert_eq!(line[76..78].trim(), atom.element.symbol());
    }
    assert!(r.losses.contains(&"identity_namespaces"));
    assert!(r.losses.contains(&"entity_identifiers"));
}
#[test]
fn loss_requires_explicit_consent() {
    let d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    assert_eq!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, false)).err(),
        Some(MmcifPdbError::MetadataLossRequiresConsent)
    );
}
#[test]
fn label_namespace_is_explicit_and_missing_labels_are_not_fabricated() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    assert!(matches!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Label, true)),
        Err(MmcifPdbError::InvalidField {
            field: "residue number",
            ..
        })
    ));
    for (i, r) in d.residues.iter_mut().enumerate() {
        r.identity.label_seq_id = Some(10 + i as i64);
    }
    let r = diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Label, true)).unwrap();
    let line = r.pdb.lines().next().unwrap();
    assert_eq!(&line[21..22], "E");
    assert_eq!(line[22..26].trim(), "10");
}
#[test]
fn alias_collisions_cannot_merge_residues_or_atoms() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.residues[1].identity.auth_seq_id = Some("1".into());
    assert!(matches!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)),
        Err(MmcifPdbError::ResidueCollision {
            first: 0,
            second: 1
        })
    ));
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.atom_ids[1].auth_atom_id = Some("C1".into());
    assert!(matches!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)),
        Err(MmcifPdbError::AtomCollision {
            first: 0,
            second: 1
        })
    ));
}
#[test]
fn guards_widths_ascii_and_identifier_normalization() {
    for (field, value) in [
        ("residue name", "LONG"),
        ("chain identifier", "AB"),
        ("chain identifier", "$"),
        ("residue number", "10000"),
        ("residue number", "01"),
        ("residue number", "+1"),
        ("insertion code", "AA"),
        ("insertion code", "*"),
        ("atom name", "ABCDE"),
        ("atom name", "C\n"),
        ("atom name", "原子"),
        ("alternate location", "AB"),
    ] {
        let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
        match field {
            "residue name" => d.residues[0].identity.auth_comp_id = Some(value.into()),
            "chain identifier" => d.residues[0].identity.auth_asym_id = Some(value.into()),
            "residue number" => d.residues[0].identity.auth_seq_id = Some(value.into()),
            "insertion code" => d.residues[0].identity.insertion_code = Some(value.into()),
            "atom name" => d.atom_ids[0].auth_atom_id = Some(value.into()),
            _ => d.atom_ids[0].alt_id = Some(value.into()),
        }
        assert_eq!(
            diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).err(),
            Some(MmcifPdbError::InvalidField { atom: 0, field })
        );
    }
}
#[test]
fn fixed_columns_preserve_altloc_insertion_charge_and_two_letter_element() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.atom_ids[0].alt_id = Some("B".into());
    d.residues[0].identity.insertion_code = Some("A".into());
    d.parsed.atoms[0].element = chematic_core::Element::ZN;
    d.parsed.atoms[0].formal_charge = Some(2);
    let r = diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).unwrap();
    let line = r.pdb.lines().next().unwrap();
    assert_eq!(&line[16..17], "B");
    assert_eq!(&line[26..27], "A");
    assert_eq!(&line[76..80], "ZN2+");
    assert_eq!(&line[12..16], "C1  ");
}
#[test]
fn reports_coordinate_occupancy_and_temperature_rounding() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.parsed.atoms[0].x = 20.123456;
    d.parsed.atoms[0].occupancy = 0.6789;
    d.parsed.atoms[0].b_iso = 12.3456;
    let r = diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).unwrap();
    assert!((r.max_coordinate_rounding_error - 0.000456).abs() < 1e-12);
    assert!((r.max_occupancy_rounding_error - 0.0011).abs() < 1e-12);
    assert!((r.max_b_iso_rounding_error - 0.0044).abs() < 1e-12);
    assert!(r.losses.contains(&"numeric_rounding"));
}
#[test]
fn nonfinite_and_rounding_overflow_are_errors() {
    for (field, value) in [
        ("x", f64::NAN),
        ("x", 9999.9996),
        ("x", -1000.0),
        ("occupancy", 1.1),
        ("occupancy", f64::INFINITY),
        ("B factor", 999.999),
    ] {
        let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
        match field {
            "x" => d.parsed.atoms[0].x = value,
            "occupancy" => d.parsed.atoms[0].occupancy = value,
            _ => d.parsed.atoms[0].b_iso = value,
        }
        assert_eq!(
            diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).err(),
            Some(MmcifPdbError::InvalidField { atom: 0, field })
        );
    }
}
#[test]
fn serials_and_formal_charges_are_not_clamped_or_renumbered() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.parsed.atoms[1].serial = 1;
    assert!(matches!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)),
        Err(MmcifPdbError::DuplicateSerial { atom: 1, serial: 1 })
    ));
    for v in [0, 100000] {
        let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
        d.parsed.atoms[0].serial = v;
        assert!(matches!(
            diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)),
            Err(MmcifPdbError::InvalidField {
                field: "serial",
                ..
            })
        ));
    }
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.parsed.atoms[0].formal_charge = Some(10);
    assert!(matches!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)),
        Err(MmcifPdbError::InvalidField {
            field: "formal charge",
            ..
        })
    ));
}
#[test]
fn source_metadata_omissions_are_reported() {
    let mut d = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    d.parsed.space_group = Some("P 1".into());
    d.parsed.unhandled_columns.push("_atom_site.extra".into());
    d.parsed.atoms[0].formal_charge = Some(0);
    let r = diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).unwrap();
    for code in [
        "cell_and_symmetry",
        "unhandled_atom_site_columns",
        "optional_formal_charge_presence",
    ] {
        assert!(r.losses.contains(&code));
    }
}
#[test]
fn inconsistent_public_documents_and_unsupported_models_fail_safely() {
    let original = diagnostic(parse_mmcif_residues(GLYCAN)).unwrap();
    let mut variants = Vec::new();
    let mut d = original.clone();
    d.atom_to_residue.clear();
    variants.push(d);
    let mut d = original.clone();
    d.atom_ids.clear();
    variants.push(d);
    let mut d = original.clone();
    d.residues[0].atom_indices.push(100);
    variants.push(d);
    let mut d = original.clone();
    d.residues[0].atom_indices.push(0);
    variants.push(d);
    let mut d = original.clone();
    d.atom_to_residue[0] = 100;
    variants.push(d);
    for d in variants {
        assert_eq!(
            diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).err(),
            Some(MmcifPdbError::InvalidDocument)
        );
    }
    let mut d = original;
    d.residues[0].identity.model_num = 2;
    assert_eq!(
        diagnostic(write_mmcif_pdb_checked(&d, PdbNamespace::Author, true)).err(),
        Some(MmcifPdbError::UnsupportedModels)
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
