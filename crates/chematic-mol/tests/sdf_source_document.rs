use chematic_mol::{
    MolParseError, SdfSourceLimits, roundtrip_sdf_source, roundtrip_sdf_source_with_limits,
    single_sdf_mol_block, validate_sdf_source_with_limits,
};
use serde_json::Value;
const PACKET: &str = include_str!("../../../validation/sdf_source_document/cases.json");
fn cases() -> Vec<Value> {
    serde_json::from_str::<Value>(PACKET).unwrap()["cases"]
        .as_array()
        .unwrap()
        .clone()
}
#[test]
fn every_record_and_source_byte_survive_in_all_layouts() {
    for case in cases() {
        let source = case["input"].as_str().unwrap();
        for input in [
            source.to_owned(),
            source.replace('\n', "\r\n"),
            source.trim_end_matches('\n').to_owned(),
        ] {
            assert_eq!(roundtrip_sdf_source(&input).unwrap(), input);
            let summary =
                validate_sdf_source_with_limits(&input, &SdfSourceLimits::default()).unwrap();
            assert_eq!(summary.records as u64, case["records"].as_u64().unwrap());
            assert_eq!(
                summary.total_atoms as u64,
                case["total_atoms"].as_u64().unwrap()
            );
            if summary.records == 1 {
                let block = single_sdf_mol_block(&input).unwrap();
                assert!(block.ends_with("M  END\n") || block.ends_with("M  END\r\n"));
                assert!(!block.contains("<NOTE>"));
            } else {
                assert!(single_sdf_mol_block(&input).is_err());
            }
        }
    }
}
#[test]
fn later_failures_and_empty_records_never_return_a_valid_prefix() {
    let packet = cases();
    let source = packet[0]["input"].as_str().unwrap();
    let one = source.split("$$$$\n").next().unwrap();
    let second_bad = source.replace("    2.1250 O", "       NaN O");
    let missing_last_end = source
        .rsplit_once("M  END")
        .map(|(a, b)| format!("{a}{b}"))
        .unwrap();
    for input in [
        String::new(),
        " \n".into(),
        "$$$$\n".into(),
        format!("$$$$\n{source}"),
        format!("{one}$$$$\n$$$$\n"),
        format!("{source}junk\n"),
        missing_last_end,
        second_bad,
        source.replace("M  END", "M  END suffix"),
    ] {
        assert!(roundtrip_sdf_source(&input).is_err(), "{input}");
    }
}
#[test]
fn all_resource_limits_are_inclusive() {
    let packet = cases();
    let input = packet[0]["input"].as_str().unwrap();
    let record_bytes = input.split("$$$$\n").map(str::len).max().unwrap();
    let mut exact = SdfSourceLimits::default();
    exact.sdf.max_input_bytes = input.len();
    exact.sdf.max_record_bytes = record_bytes;
    exact.sdf.max_line_bytes = input.split_inclusive('\n').map(str::len).max().unwrap();
    exact.sdf.max_records = 2;
    exact.max_atoms_per_record = 1;
    exact.max_total_atoms = 2;
    assert_eq!(
        roundtrip_sdf_source_with_limits(input, &exact).unwrap(),
        input
    );
    let mut bad = Vec::new();
    let mut limit = exact;
    limit.sdf.max_input_bytes -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.sdf.max_record_bytes -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.sdf.max_line_bytes -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.sdf.max_records -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_atoms_per_record = 0;
    bad.push(limit);
    let mut limit = exact;
    limit.max_total_atoms -= 1;
    bad.push(limit);
    for limit in bad {
        assert!(matches!(
            roundtrip_sdf_source_with_limits(input, &limit),
            Err(MolParseError::ResourceLimit { .. })
        ));
    }
}

#[test]
fn zero_atom_records_are_not_empty_records() {
    let input = "zero\n\n\n  0  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n>  <NOTE>\nmetadata only\n\n$$$$\n";
    let summary = validate_sdf_source_with_limits(input, &SdfSourceLimits::default()).unwrap();
    assert_eq!(summary.records, 1);
    assert_eq!(summary.total_atoms, 0);
    assert_eq!(roundtrip_sdf_source(input).unwrap(), input);
    assert!(single_sdf_mol_block(&format!("{input}{input}")).is_err());
    let mut limits = SdfSourceLimits::default();
    limits.sdf.max_records = 1;
    assert!(roundtrip_sdf_source_with_limits(&format!("{input}{input}"), &limits).is_err());
}
