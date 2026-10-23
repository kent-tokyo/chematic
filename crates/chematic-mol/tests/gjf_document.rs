use chematic_mol::{
    GaussianError, GjfDocumentLimits, parse_gjf_source_document,
    parse_gjf_source_document_with_limits, parse_gjf_source_record, write_gjf_source_document,
};

const SOURCE: &str = include_str!("../../../validation/gjf-document-2026-10-22/two-jobs.gjf");

#[test]
fn all_jobs_and_source_ranges_match_independent_fixture_expectations() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gjf-document-2026-10-22/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        let result = parse_gjf_source_document(source);
        if case["expected"] == "refuse" {
            assert!(result.is_err(), "{}", case["id"]);
            continue;
        }
        let document = result.unwrap();
        assert_eq!(
            write_gjf_source_document(&document).as_bytes(),
            source.as_bytes()
        );
        let actual = document.to_json();
        for (index, expected) in case["jobs"].as_array().unwrap().iter().enumerate() {
            let job = &actual["jobs"][index];
            for key in [
                "source_start",
                "source_end",
                "elements",
                "freeze_codes",
                "charge",
                "multiplicity",
            ] {
                assert_eq!(job[key], expected[key], "{} job {index} {key}", case["id"]);
            }
            for axis in 0..3 {
                assert_eq!(
                    [
                        document.jobs()[index].geometry().coords[0].0,
                        document.jobs()[index].geometry().coords[0].1,
                        document.jobs()[index].geometry().coords[0].2
                    ][axis]
                        .to_bits(),
                    expected["coords"][0][axis].as_f64().unwrap().to_bits()
                );
            }
            let range = document.jobs()[index].source_range();
            assert!(!document.source_text()[range].contains("--Link1--"));
        }
        assert_eq!(
            document.jobs().len(),
            case["jobs"].as_array().unwrap().len()
        );
    }
}

#[test]
fn cumulative_budgets_cannot_reset_at_link1() {
    let exact = GjfDocumentLimits {
        gaussian: chematic_mol::GaussianParseLimits {
            max_input_bytes: SOURCE.len(),
            max_atoms: 2,
            max_sections: 6,
            ..Default::default()
        },
        max_jobs: 2,
    };
    assert!(parse_gjf_source_document_with_limits(SOURCE, &exact).is_ok());
    let mut policies = vec![
        GjfDocumentLimits {
            max_jobs: 1,
            ..exact
        },
        GjfDocumentLimits {
            max_jobs: 0,
            ..exact
        },
    ];
    for gaussian in [
        chematic_mol::GaussianParseLimits {
            max_input_bytes: SOURCE.len() - 1,
            ..exact.gaussian
        },
        chematic_mol::GaussianParseLimits {
            max_atoms: 1,
            ..exact.gaussian
        },
        chematic_mol::GaussianParseLimits {
            max_sections: 5,
            ..exact.gaussian
        },
        chematic_mol::GaussianParseLimits {
            max_lines: SOURCE.lines().count() - 1,
            ..exact.gaussian
        },
        chematic_mol::GaussianParseLimits {
            max_line_bytes: 8,
            ..exact.gaussian
        },
    ] {
        policies.push(GjfDocumentLimits { gaussian, ..exact });
    }
    for limits in policies {
        assert!(matches!(
            parse_gjf_source_document_with_limits(SOURCE, &limits),
            Err(GaussianError::ResourceLimit { .. })
        ));
    }
}

#[test]
fn single_job_apis_still_refuse_implicit_selection() {
    assert!(parse_gjf_source_record(SOURCE).is_err());
    assert!(chematic_mol::parse_gjf_geometry(SOURCE).is_err());
    let doc = parse_gjf_source_document(SOURCE).unwrap();
    assert_eq!(doc.jobs().len(), 2);
    for job in doc.jobs() {
        assert_eq!(job.geometry().mol.bond_count(), 0);
    }
}

#[test]
fn job_limit_boundary_and_no_partial_success() {
    let one = "# HF/STO-3G\n\nHe\n\n0 1\nHe 0 0 0\n\n";
    let source = vec![one; 128].join("--Link1--\n");
    assert_eq!(
        parse_gjf_source_document(&source).unwrap().jobs().len(),
        128
    );
    assert!(parse_gjf_source_document(&format!("{source}--Link1--\n{one}")).is_err());
    assert!(parse_gjf_source_document(&format!("{source}--Link1--\ninvalid")).is_err());
}
