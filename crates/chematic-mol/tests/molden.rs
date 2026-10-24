use chematic_mol::{
    MoldenLengthUnit, MoldenLimits, molden_record_from_json, parse_molden,
    parse_molden_with_limits, write_molden,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/molden_vibrations/cases.json"
    ))
    .unwrap()
}

#[test]
fn shared_profile_records_preserve_units_vectors_modes_and_optional_intensities() {
    for c in packet()["cases"].as_array().unwrap() {
        let record = parse_molden(c["input"].as_str().unwrap()).unwrap();
        assert_eq!(
            serde_json::to_value(&record).unwrap(),
            c["expected"],
            "{}",
            c["id"]
        );
        assert_eq!(
            parse_molden(&write_molden(&record).unwrap()).unwrap(),
            record,
            "{}",
            c["id"]
        );
        assert_eq!(
            molden_record_from_json(&serde_json::to_string(&record).unwrap()).unwrap(),
            record
        );
    }
}
#[test]
fn unsupported_sections_and_malformed_input_are_not_silently_dropped() {
    for c in packet()["refusals"].as_array().unwrap() {
        assert!(
            parse_molden(c["input"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
}
#[test]
fn writer_emits_edited_values_and_keeps_absent_raman_distinct_from_zero() {
    let mut record = parse_molden(include_str!(
        "../../../validation/molden_vibrations/dimer.molden"
    ))
    .unwrap();
    record.frequencies.as_mut().unwrap()[0] = -123.45678901234567;
    record.normal_modes.as_mut().unwrap()[0].displacements_bohr[0] =
        [1e-130, 0.12345678901234567, -3.0];
    record.intensities.as_mut().unwrap()[0].ir = 0.9876543210987654;
    let original = write_molden(&record).unwrap();
    assert_eq!(parse_molden(&original).unwrap(), record);
    record.intensities.as_mut().unwrap()[0].raman = Some(0.0);
    let edited = write_molden(&record).unwrap();
    assert_ne!(original, edited);
    assert_eq!(parse_molden(&edited).unwrap(), record);
}
#[test]
fn reordered_sections_and_crlf_have_the_same_data() {
    let source = include_str!("../../../validation/molden_vibrations/dimer.molden");
    let (head, intensities) = source.split_once("[INT]\n").unwrap();
    let reordered = format!(
        "[Molden Format]\n[INT]\n{intensities}{}",
        head.trim_start_matches("[Molden Format]\n")
    );
    assert_eq!(
        parse_molden(source).unwrap(),
        parse_molden(&reordered.replace('\n', "\r\n")).unwrap()
    );
}
#[test]
fn invalid_constructed_data_and_unknown_json_keys_cannot_bypass_checks() {
    let source = parse_molden(include_str!(
        "../../../validation/molden_vibrations/dimer.molden"
    ))
    .unwrap();
    let mut r = source.clone();
    r.normal_modes.as_mut().unwrap()[0].index = 0;
    assert!(write_molden(&r).is_err());
    let mut r = source.clone();
    r.normal_modes.as_mut().unwrap()[0].displacements_bohr.pop();
    assert!(write_molden(&r).is_err());
    let mut r = source.clone();
    r.intensities.as_mut().unwrap().clear();
    assert!(write_molden(&r).is_err());
    let mut r = source.clone();
    r.fr_coords = None;
    assert!(write_molden(&r).is_err());
    for n in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut r = source.clone();
        r.frequencies.as_mut().unwrap()[0] = n;
        assert!(write_molden(&r).is_err());
        let mut r = source.clone();
        r.intensities.as_mut().unwrap()[0].raman = Some(n);
        assert!(write_molden(&r).is_err());
    }
    let mut json = serde_json::to_value(source).unwrap();
    json["orbital_coefficients"] = serde_json::json!([1.0]);
    assert!(molden_record_from_json(&json.to_string()).is_err());
    json.as_object_mut().unwrap().remove("orbital_coefficients");
    json["normal_modes"][0]["normalization"] = "mass-weighted".into();
    assert!(molden_record_from_json(&json.to_string()).is_err());
}
#[test]
fn limits_bound_total_displacements_and_not_only_each_mode() {
    let source = include_str!("../../../validation/molden_vibrations/dimer.molden");
    let defaults = MoldenLimits::default();
    for limits in [
        MoldenLimits {
            max_input_bytes: source.len() - 1,
            ..defaults
        },
        MoldenLimits {
            max_atoms: 1,
            ..defaults
        },
        MoldenLimits {
            max_modes: 0,
            ..defaults
        },
        MoldenLimits {
            max_displacements: 1,
            ..defaults
        },
        MoldenLimits {
            max_line_len: 5,
            ..defaults
        },
    ] {
        assert!(parse_molden_with_limits(source, &limits).is_err());
    }
    let c = packet()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "mode_order_and_optional_raman")
        .unwrap()
        .clone();
    assert!(
        parse_molden_with_limits(
            c["input"].as_str().unwrap(),
            &MoldenLimits {
                max_displacements: 3,
                ..defaults
            }
        )
        .is_err()
    );
    assert!(
        parse_molden_with_limits(
            c["input"].as_str().unwrap(),
            &MoldenLimits {
                max_displacements: 4,
                ..defaults
            }
        )
        .is_ok()
    );
}
#[test]
fn text_preserves_extreme_finite_values_without_converting_atom_units() {
    let mut r = parse_molden("[Molden Format]\n[Atoms] AU\nH 1 1 0 0 0\n").unwrap();
    for n in [
        f64::MAX,
        -f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -0.0,
    ] {
        r.atoms.as_mut().unwrap().atoms[0].position = [n; 3];
        let back = parse_molden(&write_molden(&r).unwrap()).unwrap();
        assert_eq!(back.atoms.as_ref().unwrap().unit, MoldenLengthUnit::Bohr);
        assert_eq!(
            back.atoms.as_ref().unwrap().atoms[0].position[0].to_bits(),
            n.to_bits()
        );
    }
}
