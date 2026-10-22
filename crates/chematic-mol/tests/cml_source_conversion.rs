//! Original source-retention fixtures; see validation/cml_source_preservation/cases.json.
use chematic_mol::{
    ConversionError, ConversionLimits, ConversionPolicy as Policy, PreservationStatus as Status,
    convert_with_report, convert_with_report_json, convert_with_report_with_limits,
};

const SIMPLE: &str = "<molecule xmlns=\"http://www.xml-cml.org/schema\"><atom id=\"a\" elementType=\"He\"/></molecule>";

fn convert(text: &str) -> Result<chematic_mol::ConversionReport, ConversionError> {
    convert_with_report(text, "cml", "cml", Policy::Strict)
}

#[test]
fn original_documents_retain_all_bytes_under_each_policy() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cml_source_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        for policy in [Policy::Strict, Policy::AllowLoss, Policy::PreserveSource] {
            let report = convert_with_report(source, " .CML ", ".cml", policy).unwrap();
            assert_eq!(report.output.as_deref(), Some(source));
            assert_eq!(report.status, "completed");
            assert_eq!(report.source_record_index, 0);
            assert_eq!(report.verification_scope, "exact_source_bytes");
            assert!(!report.candidate_graph_verified);
            assert!(report.candidate_atom_mapping.is_empty());
            assert!(report.candidate_bond_mapping.is_empty());
            assert_eq!(report.observations.len(), 2);
            assert_eq!(report.observations[0].status, Status::Preserved);
            assert_eq!(report.observations[1].status, Status::OpaqueRetained);
        }
    }
}

#[test]
fn malformed_xml_and_unsupported_document_framing_fail_without_output() {
    for input in [
        "",
        "<molecule/>",
        "<other xmlns='http://www.xml-cml.org/schema'/>",
        "<molecule xmlns='wrong'/>",
        "<molecule xmlns='http://www.xml-cml.org/schema'><atom></molecule>",
        "<molecule xmlns='http://www.xml-cml.org/schema' id='x' id='y'/>",
        "<molecule xmlns='http://www.xml-cml.org/schema'>&unknown;</molecule>",
        "<molecule xmlns='http://www.xml-cml.org/schema'>&#0;</molecule>",
        "<molecule xmlns='http://www.xml-cml.org/schema'><v:x/></molecule>",
        "<cml xmlns='http://www.xml-cml.org/schema'><molecule/></cml>",
        "<molecule xmlns='http://www.xml-cml.org/schema'><molecule/></molecule>",
    ] {
        assert!(convert(input).is_err(), "accepted {input}");
    }
    assert!(convert(&format!("{SIMPLE}{SIMPLE}")).is_err());
    assert!(convert_with_report(SIMPLE, "cml", "smiles", Policy::AllowLoss).is_err());
}

#[test]
fn dtds_and_unsupported_declarations_are_rejected_but_literal_markup_is_retained() {
    for prolog in [
        "<!DOCTYPE molecule>",
        "<!DOCTYPE molecule SYSTEM 'file:///nonexistent-cml-source-test'>",
        "<!DOCTYPE molecule [<!ENTITY x 'value'>]>",
        "<?xml version='1.1'?>",
        "<?xml version='9.9'?>",
        "<?xml version='1&#46;0'?>",
        "<?xml version='1.0' encoding='UTF&#45;8'?>",
        "<?xml version='1.0' encoding='ISO-8859-1'?>",
    ] {
        assert!(
            convert(&format!("{prolog}{SIMPLE}")).is_err(),
            "accepted {prolog}"
        );
    }
    let literal = "\u{feff}<?xml version='1.0' encoding='UTF-8'?>\r\n<!-- <!DOCTYPE literal> --><?fixture retain?><molecule xmlns='http://www.xml-cml.org/schema'><opaque><![CDATA[<!DOCTYPE literal>]]></opaque></molecule>";
    assert_eq!(convert(literal).unwrap().output.as_deref(), Some(literal));
    let unicode = format!("<?xml version='1.0' encoding='US-ASCII'?>{SIMPLE}<!-- α -->");
    assert!(convert(&unicode).is_err());
    assert!(
        convert(&format!(
            "<?xml version='1.0' encoding='US-ASCII'?>{SIMPLE}"
        ))
        .is_ok()
    );
}

#[test]
fn retained_opaque_payload_is_never_reported_as_verified_chemistry() {
    let opaque = SIMPLE.replace(
        "elementType=\"He\"",
        "elementType=\"Xx\" formalCharge=\"oops\" atomRefs4=\"missing\"",
    );
    let report = convert(&opaque).unwrap();
    assert_eq!(report.output.as_deref(), Some(opaque.as_str()));
    assert!(!report.candidate_graph_verified);
    assert_eq!(report.observations[1].status, Status::OpaqueRetained);
    assert!(convert("<molecule xmlns='http://www.xml-cml.org/schema'/>").is_ok());
}

#[test]
fn limits_apply_to_input_xml_and_output_envelope() {
    for mode in 0..5 {
        let mut limits = ConversionLimits::default();
        match mode {
            0 => limits.mol2.max_input_bytes = SIMPLE.len() - 1,
            1 => limits.mol2.max_line_bytes = SIMPLE.len() - 1,
            2 => limits.mol2.max_lines = 0,
            3 => limits.mol2.max_atoms = 0,
            _ => limits.max_output_bytes = SIMPLE.len() - 1,
        }
        assert!(matches!(
            convert_with_report_with_limits(SIMPLE, "cml", "cml", Policy::Strict, &limits),
            Err(ConversionError::ResourceLimit(_))
        ));
    }
    let mut limits = ConversionLimits::default();
    limits.mol2.max_bonds = 0;
    let bonded = SIMPLE.replace("</molecule>", "<bond atomRefs2='a a'/></molecule>");
    assert!(matches!(
        convert_with_report_with_limits(&bonded, "cml", "cml", Policy::Strict, &limits),
        Err(ConversionError::ResourceLimit(_))
    ));
    let nested = |count| {
        SIMPLE.replace(
            "<atom id=\"a\" elementType=\"He\"/>",
            &("<x>".repeat(count) + &"</x>".repeat(count)),
        )
    };
    let shallow = nested(63);
    assert!(convert(&shallow).is_ok());
    let deep = nested(64);
    assert!(matches!(
        convert(&deep),
        Err(ConversionError::ResourceLimit(_))
    ));
    let many = SIMPLE.replace(
        "<atom id=\"a\" elementType=\"He\"/>",
        &"<x/>".repeat(200_001),
    );
    assert!(matches!(
        convert(&many),
        Err(ConversionError::ResourceLimit(_))
    ));
    limits = ConversionLimits::default();
    limits.max_output_bytes = SIMPLE.len();
    assert!(convert_with_report_with_limits(SIMPLE, "cml", "cml", Policy::Strict, &limits).is_ok());
    assert!(matches!(
        convert_with_report_json(SIMPLE, "cml", "cml", "strict", &limits),
        Err(ConversionError::ResourceLimit(_))
    ));
}
