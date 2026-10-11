//! Public SMILES parser replay of the checked-in malicious-input corpus.

use chematic_smiles::{SmilesParseLimits, parse, parse_with_limits};
use serde_json::Value;

const CORPUS: &str = include_str!("../../../validation/parser_security_corpus_v1.json");

#[test]
fn smiles_security_cases_are_bounded_and_match_the_recorded_disposition() {
    let document: Value = serde_json::from_str(CORPUS).unwrap();
    let mut count = 0;
    for case in document["cases"].as_array().unwrap() {
        if case["format"] != "smiles" {
            continue;
        }
        count += 1;
        let payload = case["payload"].as_str().unwrap();
        let result = parse(payload);
        assert_eq!(
            result.is_ok(),
            case["expected_status"] == "accepted",
            "{}",
            case["id"]
        );
        let limits = SmilesParseLimits {
            max_input_bytes: payload.len().saturating_sub(1),
            ..SmilesParseLimits::default()
        };
        let _ = parse_with_limits(payload, &limits);
    }
    assert_eq!(count, 20);
}
