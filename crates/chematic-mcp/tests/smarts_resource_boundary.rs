use serde_json::json;

fn match_call(carbons: usize, nitrogens: usize) -> serde_json::Value {
    let smiles = [vec!["C"; carbons], vec!["N"; nitrogens]]
        .concat()
        .join(".");
    let request = json!({
        "jsonrpc": "2.0", "id": 1, "method": "tools/call",
        "params": {"name": "smarts_match", "arguments": {"smiles": smiles, "smarts": "C.N"}}
    });
    chematic_mcp::handle_line(&request.to_string()).unwrap()
}

#[test]
fn smarts_tool_reports_match_limit_instead_of_returning_oversized_results() {
    let small = match_call(2, 3);
    let text = small["result"]["content"][0]["text"].as_str().unwrap();
    let result: serde_json::Value = serde_json::from_str(text).unwrap();
    assert_eq!(result["match_count"], 6);

    let large = match_call(101, 101);
    assert_eq!(large["id"], 1);
    assert_eq!(large["error"]["code"], -32000);
    assert_eq!(
        large["error"]["message"],
        "SMARTS matching produced more than 10,000 matches"
    );
}
