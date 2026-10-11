//! Public SMARTS parser replay of the checked-in malicious-input corpus.

use chematic_smarts::parse_smarts;
use serde_json::Value;

const CORPUS: &str = include_str!("../../../validation/parser_security_corpus_v1.json");

#[test]
fn smarts_security_cases_match_the_recorded_disposition_without_panicking() {
    let document: Value = serde_json::from_str(CORPUS).unwrap();
    let mut count = 0;
    for case in document["cases"].as_array().unwrap() {
        if case["format"] != "smarts" {
            continue;
        }
        count += 1;
        let result = parse_smarts(case["payload"].as_str().unwrap());
        assert_eq!(
            result.is_ok(),
            case["expected_status"] == "accepted",
            "{}",
            case["id"]
        );
    }
    assert_eq!(count, 20);
}
