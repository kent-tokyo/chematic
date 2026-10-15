//! Original MIT OR Apache-2.0 RDFile fixtures, with manifest in validation/openbabel_issues/rdfile.
use chematic_mol::{
    RdfileData, RdfileParseLimits, parse_rdfile, parse_rdfile_with_limits, write_rdfile,
};
const FIXTURE: &str = include_str!("../../../validation/openbabel_issues/rdfile/mixed-records.rdf");
#[test]
fn mixed_records_keep_order_registry_spelling_and_duplicate_fields() {
    let d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    assert_eq!(d.records.len(), 5);
    assert_eq!(d.records[0].internal_registry.as_deref(), Some("0007"));
    assert_eq!(d.records[0].data.len(), 2);
    assert_eq!(d.records[0].data[0].name, d.records[0].data[1].name);
    assert_eq!(
        d.records[0].data[1].value,
        "second\nwith  two spaces\n\nlast"
    );
    assert_eq!(d.source_text, FIXTURE);
}
#[test]
fn graph_projection_preserves_reference_positions_and_reaction_roles() {
    let d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let s = d
        .records
        .iter()
        .map(|r| r.smiles().unwrap())
        .collect::<Vec<_>>();
    assert!(s[0].is_some());
    assert!(s[1].as_ref().unwrap().contains("13C"));
    assert!(s[2].as_ref().unwrap().contains(">>"));
    assert!(s[3].is_none() && s[4].is_none());
}
#[test]
fn canonical_write_reparse_retains_records_and_raw_structures() {
    let d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let back = diagnostic(parse_rdfile(&diagnostic(write_rdfile(&d)).unwrap())).unwrap();
    assert_eq!(back.records, d.records);
    assert_eq!(back.date_comment, d.date_comment);
}
#[test]
fn crlf_source_and_structure_bytes_remain_exact() {
    let text = FIXTURE.replace('\n', "\r\n");
    let d = diagnostic(parse_rdfile(&text)).unwrap();
    assert_eq!(d.source_text, text);
    assert!(d.records[0].structure.as_ref().unwrap().contains("\r\n"));
    assert_eq!(
        diagnostic(parse_rdfile(&diagnostic(write_rdfile(&d)).unwrap()))
            .unwrap()
            .records,
        d.records
    );
}
#[test]
fn writing_mutated_fields_does_not_replay_stale_source() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    d.records[0].data[0].value = "updated".into();
    let text = diagnostic(write_rdfile(&d)).unwrap();
    assert_eq!(
        diagnostic(parse_rdfile(&text)).unwrap().records[0].data[0].value,
        "updated"
    );
    assert_eq!(d.source_text, FIXTURE);
}
#[test]
fn column_81_plus_joins_wrapped_data_without_losing_newlines() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let value = format!("{}\n  continuation\n", "abc ".repeat(80));
    d.records[0].data = vec![RdfileData {
        name: "long".into(),
        value: value.clone(),
    }];
    let text = diagnostic(write_rdfile(&d)).unwrap();
    assert!(text.lines().any(|s| s.len() == 81 && s.ends_with('+')));
    assert_eq!(
        diagnostic(parse_rdfile(&text)).unwrap().records[0].data[0].value,
        value
    );
}
#[test]
fn keyword_shaped_mol_title_is_not_a_record_boundary() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let block = d.records[0].structure.as_mut().unwrap();
    *block = format!("$RFMT\n{}", block.split_once('\n').unwrap().1);
    assert_eq!(
        diagnostic(parse_rdfile(&diagnostic(write_rdfile(&d)).unwrap()))
            .unwrap()
            .records
            .len(),
        5
    );
}
#[test]
fn registry_only_documents_are_not_resolved_or_discarded() {
    let text = "$RDFILE 1\n$DATM\n$MIREG 0001\n$MEREG external molecule\n$RIREG 002\n$REREG external reaction\n";
    let d = diagnostic(parse_rdfile(text)).unwrap();
    assert_eq!(d.records.len(), 4);
    assert!(d.records.iter().all(|r| r.smiles().unwrap().is_none()));
    assert_eq!(
        diagnostic(parse_rdfile(&diagnostic(write_rdfile(&d)).unwrap()))
            .unwrap()
            .records,
        d.records
    );
}
#[test]
fn empty_library_and_empty_datum_are_valid() {
    assert!(
        diagnostic(parse_rdfile("$RDFILE 1\n$DATM\n"))
            .unwrap()
            .records
            .is_empty()
    );
    let d = diagnostic(parse_rdfile(
        "$RDFILE 1\n$DATM\n$MIREG 1\n$DTYPE empty\n$DATUM\n",
    ))
    .unwrap();
    assert_eq!(d.records[0].data[0].value, "");
}
#[test]
fn malformed_headers_identifiers_and_missing_data_fail() {
    for text in [
        "",
        "$RDFILE 2\n$DATM\n",
        "$RDFILE 1\n",
        "$RDFILE 1\nwrong\n",
        "$RDFILE 1\n$DATM\n\n",
        "$RDFILE 1\n$DATM\n$MIREG\n",
        "$RDFILE 1\n$DATM\n$MIREG -1\n",
        "$RDFILE 1\n$DATM\n$MFMT $MIREG\n",
        "$RDFILE 1\n$DATM\n$MIREG 1\n$DTYPE field\n",
        "$RDFILE 1\n$DATM\n$MIREG 1\n$DTYPE\n$DATUM a\n",
        "$RDFILE 1\n$DATM\n$MIREG 1\n$DATUM orphan\n",
    ] {
        assert!(diagnostic(parse_rdfile(text)).is_err(), "{text}");
    }
}
#[test]
fn malformed_embedded_payloads_and_counts_fail() {
    for text in [
        FIXTURE.replacen("M  END", "M  BAD", 1),
        FIXTURE.replacen("$RXN", "$RXN V3000", 1),
        FIXTURE.replace("$MOL", "$UNKNOWN"),
        FIXTURE.replace("$MFMT $MIREG 0007", "$MFMT $RIREG 0007"),
    ] {
        assert!(diagnostic(parse_rdfile(&text)).is_err());
    }
}
#[test]
fn unsupported_embedded_data_and_keyword_continuations_fail() {
    for value in ["$MFMT", "x\n$UNKNOWN anything", "x\n$DATUM orphan"] {
        let text = format!("$RDFILE 1\n$DATM\n$MIREG 1\n$DTYPE x\n$DATUM {value}\n");
        assert!(diagnostic(parse_rdfile(&text)).is_err());
    }
    let text = format!(
        "$RDFILE 1\n$DATM\n$MIREG 1\n$DTYPE x\n{}+\n",
        "a".repeat(80)
    );
    assert!(diagnostic(parse_rdfile(&text)).is_err());
}
#[test]
fn writer_refuses_delimiter_injection_and_ambiguous_fields() {
    let base = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    for value in [
        "x\n$MFMT",
        "  lost leading spaces",
        "unicode 🧪",
        "carriage\rreturn",
    ] {
        let mut d = base.clone();
        d.records[0].data[0].value = value.into();
        assert!(diagnostic(write_rdfile(&d)).is_err(), "{value}");
    }
    let mut d = base.clone();
    d.records[0].internal_registry = Some("bad".into());
    assert!(diagnostic(write_rdfile(&d)).is_err());
    let mut d = base;
    d.records[0].external_registry = Some("second".into());
    assert!(diagnostic(write_rdfile(&d)).is_err());
}
#[test]
fn writer_refuses_structure_boundary_injection() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    d.records[0]
        .structure
        .as_mut()
        .unwrap()
        .push_str("$MIREG 88\n");
    assert!(diagnostic(write_rdfile(&d)).is_err());
}
#[test]
fn all_reader_resource_limits_apply() {
    let b = RdfileParseLimits::default();
    for limits in [
        RdfileParseLimits {
            max_input_bytes: 8,
            ..b
        },
        RdfileParseLimits {
            max_records: 4,
            ..b
        },
        RdfileParseLimits {
            max_fields_per_record: 1,
            ..b
        },
        RdfileParseLimits {
            max_line_len: 4,
            ..b
        },
        RdfileParseLimits {
            max_reaction_components: 1,
            ..b
        },
    ] {
        assert!(diagnostic(parse_rdfile_with_limits(FIXTURE, &limits)).is_err());
    }
}
#[test]
fn missing_final_lf_is_accepted_and_writer_adds_it() {
    let d = diagnostic(parse_rdfile(FIXTURE.trim_end_matches('\n'))).unwrap();
    assert!(diagnostic(write_rdfile(&d)).unwrap().ends_with('\n'));
}

#[test]
fn explicit_unwrapped_profile_preserves_long_plain_data() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let value = "a".repeat(320);
    d.records[0].data = vec![RdfileData {
        name: "long".into(),
        value: value.clone(),
    }];
    let out = chematic_mol::write_rdfile_with_options(
        &d,
        chematic_mol::RdfileWriteOptions { wrap_data: false },
    )
    .unwrap();
    assert!(
        out.lines()
            .any(|s| s.len() > 81 && s.starts_with("$DATUM "))
    );
    assert_eq!(
        diagnostic(parse_rdfile(&out)).unwrap().records[0].data[0].value,
        value
    );
}
#[test]
fn writer_record_limit_and_ambiguous_literal_continuation_are_checked() {
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    let r = d.records[3].clone();
    d.records = vec![r; 10001];
    assert!(diagnostic(write_rdfile(&d)).is_err());
    let mut d = diagnostic(parse_rdfile(FIXTURE)).unwrap();
    d.records[0].data = vec![RdfileData {
        name: "x".into(),
        value: format!("{}+", "a".repeat(73)),
    }];
    assert!(
        chematic_mol::write_rdfile_with_options(
            &d,
            chematic_mol::RdfileWriteOptions { wrap_data: false }
        )
        .is_err()
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
