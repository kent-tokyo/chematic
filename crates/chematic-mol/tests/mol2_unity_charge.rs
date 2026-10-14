//! RDKit #6716 and related MOL2 formal-charge projection boundaries.
use chematic_core::AtomIdx;
use chematic_mol::mol2_source::{parse_mol2_source, write_mol2_source};
use chematic_mol::{
    ConversionPolicy, Mol2Error, convert_with_report, parse_cjson_record, parse_mol2_record,
    parse_mol2_source_record, write_mol2_record, write_mol2_source_record,
};

const AMMONIUM: &str =
    include_str!("../../../validation/rdkit_issues/formats/ammonium_unity_charge.mol2");
const IONS: &str = "@<TRIPOS>MOLECULE\nions\n2 0\nSMALL\nUSER_CHARGES\n\n\
@<TRIPOS>ATOM\n42 O 0 0 0 O.3 1 LIG 0.75\n7 N 2 0 0 N.4 1 LIG -0.75\n\
@<TRIPOS>BOND\n@<TRIPOS>UNITY_ATOM_ATTR\n7 1\ncharge 1\n42 1\ncharge -1\n";

#[test]
fn rdkit_6716_ammonium_charge_survives_both_readers_and_writers() {
    let projected = parse_mol2_source_record(AMMONIUM).unwrap();
    assert_eq!(projected.molecule().atom_count(), 5);
    assert_eq!(projected.molecule().bond_count(), 4);
    assert_eq!(projected.molecule().atom(AtomIdx(0)).charge, 1);
    assert_eq!(write_mol2_source_record(&projected), AMMONIUM);
    let editable = parse_mol2_record(AMMONIUM).unwrap();
    let reparsed = parse_mol2_record(&write_mol2_record(&editable)).unwrap();
    assert_eq!(reparsed.molecule.atom(AtomIdx(0)).charge, 1);
}

#[test]
fn unity_uses_source_ids_and_keeps_partial_charges_separate() {
    let source = parse_mol2_source_record(IONS).unwrap();
    assert_eq!(source.molecule().atom(AtomIdx(0)).charge, -1);
    assert_eq!(source.molecule().atom(AtomIdx(1)).charge, 1);
    assert_eq!(source.atoms()[0].partial_charge, Some(0.75));
    assert_eq!(source.atoms()[1].partial_charge, Some(-0.75));
    assert!(source.opaque_sections().is_empty());
    assert_eq!(
        source.to_json()["graph_projection"]["formal_charge_policy"],
        "UNITY_charge_then_N.4"
    );
    let editable = parse_mol2_record(IONS).unwrap();
    assert_eq!(editable.molecule.atom(AtomIdx(0)).charge, -1);
    assert_eq!(editable.molecule.atom(AtomIdx(1)).charge, 1);
}

#[test]
fn graph_writer_preserves_both_signs_without_inventing_partial_charges() {
    let (mol, coords) = parse_mol2_source(IONS).unwrap();
    let text = write_mol2_source(&mol, &coords);
    let projected = parse_mol2_source_record(&text).unwrap();
    let editable = parse_mol2_record(&text).unwrap();
    assert_eq!(projected.charge_type(), Some("NO_CHARGES"));
    for (index, charge) in [(0, -1), (1, 1)] {
        assert_eq!(projected.molecule().atom(AtomIdx(index)).charge, charge);
        assert_eq!(editable.molecule.atom(AtomIdx(index)).charge, charge);
        assert_eq!(projected.atoms()[index as usize].partial_charge, Some(0.0));
    }
}

#[test]
fn explicit_zero_overrides_type_and_absence_uses_n4_fallback() {
    for (source, charge) in [
        (AMMONIUM.replace("charge 1", "charge 0"), 0),
        (
            AMMONIUM
                .split("@<TRIPOS>UNITY_ATOM_ATTR")
                .next()
                .unwrap()
                .to_owned(),
            1,
        ),
    ] {
        assert_eq!(
            parse_mol2_source_record(&source)
                .unwrap()
                .molecule()
                .atom(AtomIdx(0))
                .charge,
            charge
        );
        assert_eq!(
            parse_mol2_record(&source)
                .unwrap()
                .molecule
                .atom(AtomIdx(0))
                .charge,
            charge
        );
    }
}

#[test]
fn malformed_or_ambiguous_unity_attributes_are_rejected() {
    for attributes in [
        "999 1\ncharge 1\n",
        "0 1\ncharge 1\n",
        "7 -1\n",
        "7 2\ncharge 1\n",
        "7 1\ncharge\n",
        "7 1\ncharge NaN\n",
        "7 1\ncharge 128\n",
        "7 1\ncharge 1 extra\n",
        "7 2\ncharge 1\ncharge -1\n",
        "7 1\ncharge 1\n7 1\ncharge -1\n",
        "7 1\ncharge 1\n@<TRIPOS>UNITY_ATOM_ATTR\n42 1\ncharge -1\n",
    ] {
        let input = format!(
            "{}@<TRIPOS>UNITY_ATOM_ATTR\n{attributes}",
            IONS.split("@<TRIPOS>UNITY_ATOM_ATTR").next().unwrap()
        );
        assert!(
            matches!(
                parse_mol2_source_record(&input),
                Err(Mol2Error::InvalidAtomLine { .. })
            ),
            "source accepted {attributes}"
        );
        assert!(
            matches!(
                parse_mol2_record(&input),
                Err(Mol2Error::InvalidAtomLine { .. })
            ),
            "editable accepted {attributes}"
        );
    }
}

#[test]
fn unknown_attributes_remain_opaque_and_source_text_is_exact() {
    let input = IONS
        .replace("7 1\ncharge 1", "7 2\ncharge 1\ncolor red")
        .replace("charge -1\n", "charge -1\n999 1\nlabel detached metadata\n")
        .replace('\n', "\r\n");
    let record = parse_mol2_source_record(&input).unwrap();
    assert_eq!(record.molecule().atom(AtomIdx(1)).charge, 1);
    assert_eq!(record.opaque_sections(), ["UNITY_ATOM_ATTR"]);
    assert_eq!(write_mol2_source_record(&record), input);
    let editable = parse_mol2_record(&input).unwrap();
    assert_eq!(
        editable.unity_atom_attributes.last().unwrap().attributes,
        [("label".to_owned(), "detached metadata".to_owned())]
    );
}

#[test]
fn mol2_to_cjson_keeps_both_kinds_of_charge_through_verified_mapping() {
    let report = convert_with_report(IONS, "mol2", "cjson", ConversionPolicy::AllowLoss).unwrap();
    assert!(report.candidate_graph_verified);
    let output = parse_cjson_record(report.output.as_ref().unwrap()).unwrap();
    assert_eq!(report.candidate_atom_mapping.len(), 2);
    for mapping in &report.candidate_atom_mapping {
        let (formal, partial) = match mapping.source_id {
            42 => (-1, 0.75),
            7 => (1, -0.75),
            _ => panic!("unknown source id"),
        };
        assert_eq!(
            output
                .molecule()
                .atom(AtomIdx(mapping.candidate_index as u32))
                .charge,
            formal
        );
        assert_eq!(
            output.partial_charge_datasets()["USER_CHARGES"][mapping.candidate_index],
            partial
        );
    }
}
