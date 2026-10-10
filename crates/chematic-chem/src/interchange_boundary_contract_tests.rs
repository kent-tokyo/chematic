use crate::nmr::*;

fn spectrum() -> NmrSpectrum {
    NmrSpectrum {
        schema_version: 1,
        nucleus: "1H".into(),
        solvent: None,
        frequency_mhz: Some(400.0),
        temperature_k: Some(298.0),
        source: None,
        peaks: vec![NmrPeak {
            shift: 7.2,
            intensity: 1.0,
            assignment: Some("H1".into()),
            multiplicity: Some("d".into()),
            coupling_hz: vec![7.0],
        }],
        normalization: NmrNormalization::Raw,
        vendor_metadata: serde_json::json!({"nested":{"opaque":[1,"original"]}}),
    }
}
#[test]
fn nmr_validation_accumulates_typed_paths_for_independent_invalid_fields() {
    let mut bad = spectrum();
    bad.schema_version = 2;
    bad.nucleus = "  ".into();
    bad.frequency_mhz = Some(f64::NAN);
    bad.temperature_k = Some(-1.0);
    bad.peaks[0].shift = f64::INFINITY;
    bad.peaks[0].intensity = -0.1;
    bad.peaks[0].coupling_hz = vec![-1.0, f64::NAN];
    let limits = NmrLimits {
        max_peaks: 0,
        ..Default::default()
    };
    let report = bad.validate(&limits);
    assert!(!report.valid);
    let actual: Vec<_> = report
        .diagnostics
        .iter()
        .map(|d| (d.code.as_str(), d.path.as_str()))
        .collect();
    assert_eq!(
        actual,
        vec![
            ("schema_version", "schema_version"),
            ("missing_nucleus", "nucleus"),
            ("peak_limit", "peaks"),
            ("non_finite_shift", "peaks[0].shift"),
            ("invalid_intensity", "peaks[0].intensity"),
            ("invalid_coupling", "peaks[0].coupling_hz[0]"),
            ("invalid_coupling", "peaks[0].coupling_hz[1]"),
            ("invalid_experiment_value", "frequency_mhz"),
            ("invalid_experiment_value", "temperature_k")
        ]
    );
    assert!(report.diagnostics.iter().all(|d| !d.message.is_empty()));
}
#[test]
fn nmr_json_limits_are_exact_and_normalization_metadata_round_trips() {
    for norm in [
        NmrNormalization::Raw,
        NmrNormalization::MaxOne,
        NmrNormalization::SumOne,
    ] {
        let mut original = spectrum();
        original.normalization = norm;
        let text = serde_json::to_string(&original).unwrap();
        let limits = NmrLimits {
            max_peaks: 1,
            max_serialized_bytes: text.len(),
        };
        assert!(validate_nmr_json(&text, &limits).valid);
        assert_eq!(
            serde_json::from_str::<NmrSpectrum>(&text).unwrap(),
            original
        );
        let report = validate_nmr_json(
            &text,
            &NmrLimits {
                max_serialized_bytes: text.len() - 1,
                ..limits
            },
        );
        assert_eq!(report.diagnostics[0].code, "serialized_size_limit");
        assert_eq!(report.diagnostics[0].path, "$");
    }
    for text in ["", "{", "null", "{}"] {
        let r = validate_nmr_json(text, &Default::default());
        assert!(!r.valid);
        assert_eq!(r.diagnostics[0].code, "invalid_json");
        assert_eq!(r.diagnostics[0].path, "$");
    }
}
