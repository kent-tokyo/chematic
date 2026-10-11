//! MOL/SDF parser replay of the checked-in malicious-input corpus.

use chematic_mol::sdf::parse_sdf;
use chematic_mol::{parse_mol, parse_mol_v3000};
use serde_json::Value;

const CORPUS: &str = include_str!("../../../validation/parser_security_corpus_v1.json");

#[test]
fn mol_and_sdf_security_cases_match_the_recorded_disposition_without_panicking() {
    let document: Value = serde_json::from_str(CORPUS).unwrap();
    let mut count = 0;
    for case in document["cases"].as_array().unwrap() {
        let payload = case["payload"].as_str().unwrap();
        let result = match case["format"].as_str().unwrap() {
            "mol_v2000" => parse_mol(payload).map(|_| ()),
            "mol_v3000" => parse_mol_v3000(payload).map(|_| ()),
            "sdf" => parse_sdf(payload).map(|_| ()),
            _ => continue,
        };
        count += 1;
        assert_eq!(
            result.is_ok(),
            case["expected_status"] == "accepted",
            "{}",
            case["id"]
        );
    }
    assert_eq!(count, 63);
}
