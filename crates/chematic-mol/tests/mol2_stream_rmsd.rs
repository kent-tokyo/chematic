//! Original hand-derived translation controls for Open Babel #2671.
use chematic_mol::{
    Mol2Error, Mol2ParseLimits, Mol2RmsdError, Mol2SourceReader, Mol2StreamError, Mol2StreamLimits,
    mol2_reference_rmsds, parse_mol2, parse_mol2_source_record,
};

const REFERENCES: &str =
    include_str!("../../../validation/openbabel_issues/mol2-rmsd/references.mol2");
const PROBE: &str = include_str!("../../../validation/openbabel_issues/mol2-rmsd/probe.mol2");

fn rmsds(reference: &str, align: bool) -> Vec<f64> {
    mol2_reference_rmsds(reference, PROBE, align, &Default::default())
        .unwrap()
        .into_iter()
        .map(|row| row.rmsd)
        .collect()
}

#[test]
fn every_record_retains_source_metadata_and_coordinates() {
    let records = Mol2SourceReader::new(REFERENCES)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 5);
    assert_eq!(
        records.iter().map(|r| r.source()).collect::<String>(),
        REFERENCES
    );
    for (i, record) in records.iter().enumerate() {
        assert_eq!(record.name(), Some(format!("reference_{i}").as_str()));
        assert_eq!(record.atoms()[0].source_id, 10);
        assert_eq!(record.atoms()[2].partial_charge, Some(-0.3));
        assert_eq!(record.opaque_sections(), ["SUBSTRUCTURE", "VENDOR_DATA"]);
        assert_eq!(record.molecule().atom_count(), 4);
    }
    assert_eq!(records[4].coordinates()[0], (3.0, 4.0, 0.0));
}

#[test]
fn preamble_crlf_case_and_no_final_newline_are_retained() {
    let input = format!(
        "# original preamble\r\n{}",
        REFERENCES
            .trim_end()
            .replace('\n', "\r\n")
            .replace("@<TRIPOS>MOLECULE", "  @<tripos>molecule  ")
    );
    let records = Mol2SourceReader::new(&input)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 5);
    assert_eq!(
        records.iter().map(|r| r.source()).collect::<String>(),
        input
    );
}

#[test]
fn embedded_marker_text_does_not_split_a_record() {
    let input = PROBE.replace(
        "original synthetic",
        "# @<TRIPOS>MOLECULE original synthetic",
    );
    let records = Mol2SourceReader::new(&input)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].source(), input);
}

#[test]
fn single_record_apis_still_reject_concatenation() {
    assert!(matches!(
        parse_mol2(REFERENCES),
        Err(Mol2Error::MultipleMolecules { count: 5 })
    ));
    assert!(parse_mol2_source_record(REFERENCES).is_err());
}

#[test]
fn malformed_middle_record_is_reported_and_next_boundary_is_available() {
    let bad = PROBE.replace("0.000000", "NaN");
    let input = format!("{PROBE}{bad}{PROBE}");
    let mut reader = Mol2SourceReader::new(&input);
    assert!(reader.next().unwrap().is_ok());
    assert!(matches!(
        reader.next().unwrap(),
        Err(Mol2StreamError::Parse {
            record_index: 1,
            ..
        })
    ));
    assert!(reader.next().unwrap().is_ok());
    assert!(reader.next().is_none());
    assert!(reader.next().is_none());
}

#[test]
fn empty_input_and_nonempty_missing_header_are_distinct() {
    for empty in ["", " \r\n\t"] {
        assert!(Mol2SourceReader::new(empty).next().is_none());
    }
    assert!(matches!(
        Mol2SourceReader::new("garbage").next().unwrap(),
        Err(Mol2StreamError::Parse {
            record_index: 0,
            ..
        })
    ));
}

#[test]
fn whole_input_byte_limit_is_exact_and_fatal() {
    let mut limits = Mol2StreamLimits {
        max_input_bytes: REFERENCES.len(),
        ..Default::default()
    };
    assert_eq!(Mol2SourceReader::with_limits(REFERENCES, limits).count(), 5);
    limits.max_input_bytes -= 1;
    let mut reader = Mol2SourceReader::with_limits(REFERENCES, limits);
    assert!(matches!(
        reader.next().unwrap(),
        Err(Mol2StreamError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));
    assert!(reader.next().is_none());
}

#[test]
fn record_count_limit_never_silently_truncates() {
    let limits = Mol2StreamLimits {
        max_records: 4,
        ..Default::default()
    };
    let mut reader = Mol2SourceReader::with_limits(REFERENCES, limits);
    for _ in 0..4 {
        assert!(reader.next().unwrap().is_ok());
    }
    assert_eq!(
        reader.next().unwrap().err().unwrap(),
        Mol2StreamError::ResourceLimit {
            resource: "records",
            actual: 5,
            limit: 4
        }
    );
    assert!(reader.next().is_none());
    let limits = Mol2StreamLimits {
        max_records: 0,
        ..Default::default()
    };
    assert!(
        Mol2SourceReader::with_limits(PROBE, limits)
            .next()
            .unwrap()
            .is_err()
    );
}

#[test]
fn per_record_limits_are_enforced_with_record_index() {
    let limits = [
        Mol2ParseLimits {
            max_input_bytes: PROBE.len() - 1,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_atoms: 3,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_bonds: 2,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_lines: 3,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_sections: 2,
            ..Default::default()
        },
        Mol2ParseLimits {
            max_line_bytes: 5,
            ..Default::default()
        },
    ];
    for record in limits {
        let mut reader = Mol2SourceReader::with_limits(
            PROBE,
            Mol2StreamLimits {
                record,
                ..Default::default()
            },
        );
        assert!(matches!(
            reader.next().unwrap(),
            Err(Mol2StreamError::Parse {
                record_index: 0,
                error: Mol2Error::ResourceLimit { .. }
            })
        ));
    }
}

#[test]
fn all_five_raw_rmsds_match_hand_derived_translations() {
    assert_eq!(rmsds(REFERENCES, false), [0.0, 1.0, 2.0, 3.0, 5.0]);
    let rows = mol2_reference_rmsds(REFERENCES, PROBE, false, &Default::default()).unwrap();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.reference_index, i);
        assert_eq!(row.reference_name, format!("reference_{i}"));
        assert_eq!(row.probe_name, "probe");
        assert_eq!(row.heavy_atom_count, 3);
    }
}

#[test]
fn optional_alignment_removes_translations_and_proper_rotation() {
    assert!(rmsds(REFERENCES, true).iter().all(|v| v.abs() < 1e-12));
    let rotated = PROBE
        .replace(
            "1.000000 0.000000 0.000000 C.3",
            "0.000000 1.000000 0.000000 C.3",
        )
        .replace(
            "1.000000 1.000000 0.000000 O.3",
            "-1.000000 1.000000 0.000000 O.3",
        );
    assert!(rmsds(&rotated, true)[0] < 1e-12);
    assert!(rmsds(&rotated, false)[0] > 1.0);
}

#[test]
fn hydrogen_coordinates_do_not_contribute_to_the_metric() {
    let input = PROBE.replace(
        "0.000000 0.000000 1.000000 H",
        "100.000000 200.000000 300.000000 H",
    );
    assert_eq!(rmsds(&input, false), [0.0]);
}

#[test]
fn alignment_does_not_reflect_a_nonplanar_heavy_atom_set() {
    let probe = PROBE
        .replace("40 H4", "40 N4")
        .replace("1.000000 H 1", "1.000000 N.3 1");
    let mirror = probe
        .replace(
            "1.000000 0.000000 0.000000 C.3",
            "-1.000000 0.000000 0.000000 C.3",
        )
        .replace(
            "1.000000 1.000000 0.000000 O.3",
            "-1.000000 1.000000 0.000000 O.3",
        );
    let rows = mol2_reference_rmsds(&mirror, &probe, true, &Default::default()).unwrap();
    assert_eq!(rows[0].heavy_atom_count, 4);
    assert!(rows[0].rmsd > 0.1);
}

#[test]
fn identity_and_indexed_bonds_are_checked_instead_of_element_sequence_only() {
    for changed in [
        PROBE.replace("10 C1", "10 other"),
        PROBE.replace("10 C1", "11 C1").replace(" 10 ", " 11 "),
        PROBE.replace("1 10 20 1", "1 10 30 1"),
        PROBE.replace("2 20 30 1", "2 20 30 2"),
        PROBE.replace("C.3", "C.2"),
        format!("{PROBE}@<TRIPOS>UNITY_ATOM_ATTR\n10 1\ncharge 1\n"),
    ] {
        assert!(parse_mol2_source_record(&changed).is_ok());
        assert_eq!(
            mol2_reference_rmsds(&changed, PROBE, false, &Default::default()).unwrap_err(),
            Mol2RmsdError::IdentityMismatch { reference_index: 0 }
        );
    }
}

#[test]
fn reordering_identical_elements_requires_an_explicit_mapping() {
    let a = "10 C1 0.000000 0.000000 0.000000 C.3 1 LIG 0.1";
    let b = "20 C2 1.000000 0.000000 0.000000 C.3 1 LIG 0.1";
    let changed = PROBE.replace(&format!("{a}\n{b}"), &format!("{b}\n{a}"));
    assert!(matches!(
        mol2_reference_rmsds(&changed, PROBE, false, &Default::default()),
        Err(Mol2RmsdError::IdentityMismatch { .. })
    ));
}

#[test]
fn bond_row_order_and_endpoint_spelling_may_differ() {
    let changed = PROBE.replace(
        "1 10 20 1\n2 20 30 1\n3 10 40 1",
        "3 40 10 1\n2 30 20 1\n1 20 10 1",
    );
    assert_eq!(rmsds(&changed, false), [0.0]);
}

#[test]
fn malformed_or_mismatched_later_reference_aborts_the_whole_result() {
    let bad = PROBE.replace("0.000000", "NaN");
    assert!(matches!(
        mol2_reference_rmsds(&format!("{PROBE}{bad}"), PROBE, false, &Default::default()),
        Err(Mol2RmsdError::Reference(Mol2StreamError::Parse {
            record_index: 1,
            ..
        }))
    ));
    let mismatch = PROBE.replace("10 C1", "10 other");
    assert!(matches!(
        mol2_reference_rmsds(
            &format!("{PROBE}{mismatch}"),
            PROBE,
            false,
            &Default::default()
        ),
        Err(Mol2RmsdError::IdentityMismatch { reference_index: 1 })
    ));
}

#[test]
fn one_probe_empty_reference_and_hydrogen_only_are_explicit_contracts() {
    assert!(matches!(
        mol2_reference_rmsds(PROBE, REFERENCES, false, &Default::default()),
        Err(Mol2RmsdError::Probe(_))
    ));
    assert_eq!(
        mol2_reference_rmsds("", PROBE, false, &Default::default()).unwrap_err(),
        Mol2RmsdError::EmptyReference
    );
    let h = "@<TRIPOS>MOLECULE\nh\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 H1 0 0 0 H\n";
    assert_eq!(
        mol2_reference_rmsds(h, h, false, &Default::default()).unwrap_err(),
        Mol2RmsdError::NoHeavyAtoms
    );
}

#[test]
fn finite_input_overflow_never_returns_infinite_rmsd() {
    let changed = PROBE.replace("0.000000", "1e308");
    assert!(parse_mol2_source_record(&changed).is_ok());
    for align in [false, true] {
        assert!(matches!(
            mol2_reference_rmsds(&changed, PROBE, align, &Default::default()),
            Err(Mol2RmsdError::NonFiniteComputation { reference_index: 0 })
        ));
    }
}

#[test]
fn rmsd_enforces_limits_on_both_inputs() {
    let limits = Mol2StreamLimits {
        max_records: 4,
        ..Default::default()
    };
    assert!(matches!(
        mol2_reference_rmsds(REFERENCES, PROBE, false, &limits),
        Err(Mol2RmsdError::Reference(_))
    ));
    let limits = Mol2StreamLimits {
        record: Mol2ParseLimits {
            max_atoms: 3,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        mol2_reference_rmsds(REFERENCES, PROBE, false, &limits),
        Err(Mol2RmsdError::Probe(_))
    ));
}
