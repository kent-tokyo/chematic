//! Hand-derived MOL2 preservation regressions. Source IDs deliberately differ
//! from graph indices; charges cross rounding boundaries without implying ions.

use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{
    Mol2Error, Mol2ParseLimits,
    mol2_source::{
        parse_mol2_source as parse_mol2, parse_mol2_source_record as parse_mol2_record,
        parse_mol2_source_record_with_limits as parse_mol2_record_with_limits,
        write_mol2_source as write_mol2, write_mol2_source_record as write_mol2_record,
    },
};

const AMIDE: &str = include_str!("../../../validation/mol2_preservation/amide.mol2");

#[test]
fn record_preserves_charges_types_ids_and_opaque_source() {
    let record = parse_mol2_record(AMIDE).unwrap();
    assert_eq!(record.name(), Some("acetamide"));
    assert_eq!(record.molecule_type(), Some("SMALL"));
    assert_eq!(record.charge_type(), Some("USER_CHARGES"));
    assert_eq!(record.coordinates()[0], (0.123456789, 0.0, 0.0));
    assert_eq!(record.atoms()[2].source_id, 90);
    assert_eq!(record.atoms()[2].name, "oxygen");
    assert_eq!(record.atoms()[2].atom_type, "O.2");
    assert_eq!(record.atoms()[2].substructure_id, Some(1));
    assert_eq!(record.atoms()[2].substructure_name.as_deref(), Some("ACM"));
    assert_eq!(record.atoms()[2].partial_charge, Some(-0.70));
    assert_eq!(record.atoms()[0].partial_charge, Some(0.70));
    assert!(record.molecule().atoms().all(|(_, atom)| atom.charge == 0));
    assert_eq!(record.bonds()[2].source_id, 12);
    assert_eq!(record.bonds()[2].atom_ids, [7, 3]);
    assert_eq!(record.bonds()[2].bond_type, "am");
    assert_eq!(
        record.molecule().bonds().nth(2).unwrap().1.order,
        BondOrder::Single
    );
    assert_eq!(record.opaque_sections(), ["SUBSTRUCTURE", "VENDOR_DATA"]);
    assert_eq!(record.source(), AMIDE);
    assert_eq!(write_mol2_record(&record), AMIDE);
    let reparsed = parse_mol2_record(&write_mol2_record(&record)).unwrap();
    assert_eq!(reparsed.atoms(), record.atoms());
    assert_eq!(reparsed.bonds(), record.bonds());
    assert_eq!(reparsed.coordinates(), record.coordinates());
}

#[test]
fn source_writer_retains_crlf_and_empty_title() {
    let source = AMIDE.replace("acetamide\n", "\n").replace('\n', "\r\n");
    let record = parse_mol2_record(&source).unwrap();
    assert_eq!(record.name(), Some(""));
    assert_eq!(record.charge_type(), Some("USER_CHARGES"));
    assert_eq!(write_mol2_record(&record), source);
}

#[test]
fn inspection_schema_separates_partial_formal_and_opaque_information() {
    let record = parse_mol2_record(AMIDE).unwrap();
    let json = record.to_json();
    assert_eq!(json["schema"], "chematic.mol2-record.v1");
    assert_eq!(json["atoms"][2]["partial_charge"], -0.70);
    assert_eq!(json["formal_charges"], serde_json::json!([0, 0, 0, 0]));
    assert_eq!(
        json["graph_projection"]["opaque_sections"],
        "not_interpreted"
    );
    assert_eq!(json["source"], AMIDE);
}

#[test]
fn absent_partial_charge_remains_absent() {
    let input = "@<TRIPOS>MOLECULE\natom\n1 0\nSMALL\nNO_CHARGES\n\
@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n";
    let record = parse_mol2_record(input).unwrap();
    assert_eq!(record.atoms()[0].partial_charge, None);
    assert_eq!(record.atoms()[0].substructure_id, None);
    assert_eq!(record.atoms()[0].substructure_name, None);
    assert_eq!(record.charge_type(), Some("NO_CHARGES"));
    assert_eq!(write_mol2_record(&record), input);
}

#[test]
fn n4_formal_charge_comes_from_type_not_partial_charge() {
    let input = "@<TRIPOS>MOLECULE\nion\n1 0\nSMALL\nUSER_CHARGES\n\
@<TRIPOS>ATOM\n1 N1 0 0 0 N.4 1 LIG -0.70\n";
    let record = parse_mol2_record(input).unwrap();
    assert_eq!(record.molecule().atom(AtomIdx(0)).charge, 1);
    assert_eq!(record.atoms()[0].partial_charge, Some(-0.70));
}

#[test]
fn malformed_charges_and_nonfinite_geometry_are_rejected() {
    for value in ["NaN", "inf", "-inf", "1e309", "bad"] {
        for input in [
            AMIDE.replace("ACM -0.70", &format!("ACM {value}")),
            AMIDE.replace("0.123456789", value),
        ] {
            assert!(
                matches!(parse_mol2(&input), Err(Mol2Error::InvalidAtomLine { .. })),
                "{value}"
            );
            assert!(parse_mol2_record(&input).is_err());
        }
    }
}

#[test]
fn duplicate_ids_edges_self_bonds_and_dangling_references_are_rejected() {
    for input in [
        AMIDE.replace("90 oxygen", "50 oxygen"),
        AMIDE.replace("90 oxygen", "0 oxygen"),
        AMIDE.replace("12 7 3", "8 7 3"),
        AMIDE.replace("12 7 3", "0 7 3"),
        AMIDE.replace("12 7 3", "12 50 7"),
        AMIDE.replace("12 7 3", "12 7 7"),
        AMIDE.replace("12 7 3", "12 7 999"),
    ] {
        assert!(parse_mol2(&input).is_err(), "{input}");
        assert!(parse_mol2_record(&input).is_err());
    }
}

#[test]
fn truncation_and_multiple_records_cannot_succeed_as_one_molecule() {
    for input in [
        AMIDE.replace("4 3 1 0 0", "5 3 1 0 0"),
        AMIDE.replace("4 3 1 0 0", "4 4 1 0 0"),
        AMIDE.replace("4 3 1 0 0", "4 -1 1 0 0"),
        AMIDE.replace("4 3 1 0 0", "broken"),
        format!("{AMIDE}{AMIDE}"),
        format!("{AMIDE}@<TRIPOS>ATOM\n1 C1 0 0 0 C\n"),
    ] {
        assert!(parse_mol2(&input).is_err(), "{input}");
        assert!(parse_mol2_record(&input).is_err());
    }
}

#[test]
fn bond_tokens_keep_distinct_source_meaning() {
    for (token, order) in [
        ("1", BondOrder::Single),
        ("2", BondOrder::Double),
        ("3", BondOrder::Triple),
        ("ar", BondOrder::Aromatic),
        ("am", BondOrder::Single),
        ("un", BondOrder::QueryAny),
        ("du", BondOrder::Zero),
        ("nc", BondOrder::Zero),
    ] {
        let input = AMIDE.replace("12 7 3 am", &format!("12 7 3 {token}"));
        let record = parse_mol2_record(&input).unwrap();
        assert_eq!(record.bonds()[2].bond_type, token);
        assert_eq!(record.molecule().bonds().nth(2).unwrap().1.order, order);
        assert_eq!(write_mol2_record(&record), input);
    }
    for token in ["1.5", "vendor", "4"] {
        assert!(parse_mol2(&AMIDE.replace("12 7 3 am", &format!("12 7 3 {token}"))).is_err());
    }
}

#[test]
fn lowercase_section_markers_follow_the_same_limits_and_projection() {
    let input = AMIDE.replace("@<TRIPOS>", "@<tripos>");
    let record = parse_mol2_record(&input).unwrap();
    assert_eq!(record.atoms().len(), 4);
    assert_eq!(record.bonds().len(), 3);
    assert_eq!(write_mol2_record(&record), input);
    assert!(matches!(
        parse_mol2_record_with_limits(
            &input,
            &Mol2ParseLimits {
                max_sections: 4,
                ..Default::default()
            }
        ),
        Err(Mol2Error::ResourceLimit {
            resource: "sections",
            ..
        })
    ));
}

#[test]
fn source_retention_obeys_input_limits() {
    assert!(matches!(
        parse_mol2_record_with_limits(
            AMIDE,
            &Mol2ParseLimits {
                max_input_bytes: AMIDE.len() - 1,
                ..Default::default()
            }
        ),
        Err(Mol2Error::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));
}

#[test]
fn regenerated_graph_output_declares_absent_partial_charges() {
    let (molecule, coords) = parse_mol2(AMIDE).unwrap();
    let output = write_mol2(&molecule, &coords);
    let reparsed = parse_mol2_record(&output).unwrap();
    assert_eq!(reparsed.charge_type(), Some("NO_CHARGES"));
    assert!(
        reparsed
            .atoms()
            .iter()
            .all(|atom| atom.partial_charge == Some(0.0))
    );
}
