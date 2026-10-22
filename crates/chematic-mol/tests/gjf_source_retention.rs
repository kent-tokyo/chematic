use chematic_mol::{
    GaussianError, GaussianParseLimits, parse_gjf_source_record,
    parse_gjf_source_record_with_limits, write_gjf_source_record,
};

const SOURCE: &str =
    include_str!("../../../validation/gjf-cartesian-2026-10-22/configured-source.gjf");

#[test]
fn complete_source_is_retained_separately_from_the_geometry_projection() {
    for input in [
        SOURCE.to_string(),
        SOURCE.replace('\n', "\r\n"),
        format!("\u{feff}{SOURCE}"),
        SOURCE.trim_end_matches('\n').to_string(),
    ] {
        let record = parse_gjf_source_record(&input).unwrap();
        assert_eq!(record.source_text().as_bytes(), input.as_bytes());
        assert_eq!(
            write_gjf_source_record(&record).as_bytes(),
            input.as_bytes()
        );
        let geometry = record.geometry();
        assert_eq!((geometry.charge, geometry.multiplicity), (-1, 2));
        assert_eq!(geometry.mol.atom_count(), 2);
        assert_eq!(geometry.freeze_codes, [Some(-1), Some(0)]);
        assert_eq!(geometry.coords[0].1.to_bits(), (-0.0f64).to_bits());
        assert!(record.source_text().contains("%chk=original.chk"));
        assert!(record.source_text().contains("STO-3G"));
    }
}

#[test]
fn original_source_packet_is_either_retained_exactly_or_refused() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gjf-cartesian-2026-10-22/source-cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        let result = parse_gjf_source_record(source);
        if case["expected"] == "retain" {
            assert_eq!(
                write_gjf_source_record(&result.unwrap()).as_bytes(),
                source.as_bytes(),
                "{}",
                case["id"]
            );
        } else {
            assert!(result.is_err(), "accepted {}", case["id"]);
        }
    }
}

#[test]
fn source_limits_include_original_bom_and_opaque_sections() {
    let source = format!("\u{feff}{SOURCE}");
    let exact = GaussianParseLimits {
        max_input_bytes: source.len(),
        max_atoms: 2,
        ..Default::default()
    };
    assert!(parse_gjf_source_record_with_limits(&source, &exact).is_ok());
    for limits in [
        GaussianParseLimits {
            max_input_bytes: source.len() - 1,
            ..exact
        },
        GaussianParseLimits {
            max_atoms: 1,
            ..exact
        },
        GaussianParseLimits {
            max_sections: 3,
            ..exact
        },
        GaussianParseLimits {
            max_line_bytes: 8,
            ..exact
        },
        GaussianParseLimits {
            max_lines: 3,
            ..exact
        },
    ] {
        assert!(matches!(
            parse_gjf_source_record_with_limits(&source, &limits),
            Err(GaussianError::ResourceLimit { .. })
        ));
    }
}

#[test]
fn opaque_sections_are_not_claimed_as_parsed_or_external_resources() {
    let source = format!("{SOURCE}@file-that-is-not-opened.gbs\nunknown opaque payload\n\n");
    let record = parse_gjf_source_record(&source).unwrap();
    assert_eq!(write_gjf_source_record(&record), source);
    assert_eq!(record.geometry().mol.atom_count(), 2);
}
