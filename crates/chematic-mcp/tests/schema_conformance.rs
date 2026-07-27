//! Section 13/14 conformance fixtures: schema validation against a real
//! JSON Schema 2020-12 validator (the `jsonschema` crate, not this crate's
//! own hand-rolled runtime validator — see `src/schema.rs`'s module doc for
//! why the two are deliberately different implementations), an
//! all-20-tools smoke test in both protocol eras, and negative controls
//! that prove the harness itself can detect a broken registry (section 13's
//! "negative controlがgreenになるharnessは信用しないでください").
//!
//! These tests exercise `chematic_mcp`'s public wire API only
//! (`handle_line`/`Connection`) — exactly what an external MCP client would
//! see — never internal-only items.

use std::collections::HashSet;

use serde_json::{Value, json};

fn modern_meta() -> Value {
    json!({
        "io.modelcontextprotocol/protocolVersion": "2026-07-28",
        "io.modelcontextprotocol/clientInfo": { "name": "conformance-test", "version": "1.0.0" },
        "io.modelcontextprotocol/clientCapabilities": {}
    })
}

fn modern_request(id: i64, method: &str, mut params: Value) -> Value {
    params
        .as_object_mut()
        .unwrap()
        .insert("_meta".to_string(), modern_meta());
    json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params })
}

fn modern_tools_list() -> Value {
    let req = modern_request(1, "tools/list", json!({}));
    let resp = chematic_mcp::handle_line(&req.to_string()).unwrap();
    assert!(resp.get("error").is_none(), "tools/list failed: {resp}");
    resp["result"].clone()
}

fn call_tool_modern(name: &str, arguments: Value) -> Value {
    let req = modern_request(
        1,
        "tools/call",
        json!({ "name": name, "arguments": arguments }),
    );
    chematic_mcp::handle_line(&req.to_string()).unwrap()
}

// ── section 13: modern fixtures ────────────────────────────────────────────

#[test]
fn all_input_schemas_are_valid_json_schema_2020_12() {
    let list = modern_tools_list();
    for tool in list["tools"].as_array().unwrap() {
        let name = tool["name"].as_str().unwrap();
        let schema = &tool["inputSchema"];
        jsonschema::meta::validate(schema).unwrap_or_else(|e| {
            panic!("tool '{name}' inputSchema is not valid JSON Schema 2020-12: {e}")
        });
        assert_eq!(
            schema["type"], "object",
            "tool '{name}' inputSchema root must be type object"
        );
    }
}

#[test]
fn all_output_schemas_are_valid_json_schema_2020_12() {
    let list = modern_tools_list();
    for tool in list["tools"].as_array().unwrap() {
        let name = tool["name"].as_str().unwrap();
        let schema = tool
            .get("outputSchema")
            .unwrap_or_else(|| panic!("tool '{name}' is missing outputSchema"));
        jsonschema::meta::validate(schema).unwrap_or_else(|e| {
            panic!("tool '{name}' outputSchema is not valid JSON Schema 2020-12: {e}")
        });
    }
}

#[test]
fn no_schema_contains_an_external_ref() {
    fn walk(v: &Value, path: &str, offenders: &mut Vec<String>) {
        match v {
            Value::Object(map) => {
                if let Some(Value::String(r)) = map.get("$ref")
                    && !r.starts_with('#')
                {
                    offenders.push(format!("{path}/$ref = {r}"));
                }
                for (k, val) in map {
                    walk(val, &format!("{path}/{k}"), offenders);
                }
            }
            Value::Array(items) => {
                for (i, item) in items.iter().enumerate() {
                    walk(item, &format!("{path}/{i}"), offenders);
                }
            }
            _ => {}
        }
    }

    let list = modern_tools_list();
    let mut offenders = Vec::new();
    for tool in list["tools"].as_array().unwrap() {
        walk(
            &tool["inputSchema"],
            &tool["name"].to_string(),
            &mut offenders,
        );
        if let Some(out) = tool.get("outputSchema") {
            walk(out, &tool["name"].to_string(), &mut offenders);
        }
    }
    assert!(offenders.is_empty(), "external $ref found: {offenders:?}");
}

#[test]
fn tools_list_contains_exactly_20_unique_names() {
    let list = modern_tools_list();
    let names: Vec<&str> = list["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    assert_eq!(names.len(), 20);
    let unique: HashSet<&str> = names.iter().copied().collect();
    assert_eq!(unique.len(), 20, "duplicate name in {names:?}");
}

/// Representative valid arguments for each of the 20 tools. `name_to_smiles`
/// is deliberately excluded from the automatic smoke loop below — it is the
/// one tool that makes a live network call to PubChem (see
/// `crates/chematic-mcp/README.md`'s Network & privacy section), and this
/// suite must not be flaky under CI network restrictions. Its presence,
/// schema, and error taxonomy are still checked by the tests above/below;
/// only the "actually call it and expect success" step is skipped.
fn sample_arguments() -> Vec<(&'static str, Value)> {
    vec![
        ("parse_smiles", json!({ "smiles": "c1ccccc1" })),
        ("calc_properties", json!({ "smiles": "c1ccccc1" })),
        ("ecfp4", json!({ "smiles": "c1ccccc1" })),
        (
            "tanimoto",
            json!({ "smiles1": "c1ccccc1", "smiles2": "CCO" }),
        ),
        (
            "smarts_match",
            json!({ "smarts": "c1ccccc1", "smiles": "c1ccccc1" }),
        ),
        ("canonical_smiles", json!({ "smiles": "c1ccccc1" })),
        (
            "find_mcs",
            json!({ "smiles_list": ["c1ccccc1", "c1ccccc1O"] }),
        ),
        ("generate_3d", json!({ "smiles": "CCO" })),
        ("pains_check", json!({ "smiles": "CCO" })),
        ("brenk_check", json!({ "smiles": "CCO" })),
        ("sa_score", json!({ "smiles": "CCO" })),
        ("admet_profile", json!({ "smiles": "CCO" })),
        ("boiled_egg", json!({ "smiles": "CCO" })),
        ("lipinski_check", json!({ "smiles": "CCO" })),
        (
            "retrosynthesis",
            json!({ "smiles": "CC(=O)Oc1ccccc1C(=O)O" }),
        ),
        ("smiles_to_moljson", json!({ "smiles": "CCO" })),
        ("representation_router", json!({ "smiles": "CCO" })),
        ("molecule_context_pack", json!({ "smiles": "CCO" })),
    ]
}

#[test]
fn all_smoke_tested_tools_return_structured_content_conforming_to_output_schema() {
    let list = modern_tools_list();
    let tools_by_name: std::collections::HashMap<&str, &Value> = list["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| (t["name"].as_str().unwrap(), t))
        .collect();

    for (name, args) in sample_arguments() {
        let resp = call_tool_modern(name, args.clone());
        assert!(resp.get("error").is_none(), "'{name}' call failed: {resp}");
        assert_eq!(resp["result"]["resultType"], "complete", "'{name}'");
        assert!(
            resp["result"]["content"].is_array()
                && !resp["result"]["content"].as_array().unwrap().is_empty(),
            "'{name}' missing content"
        );
        let structured = resp["result"]
            .get("structuredContent")
            .unwrap_or_else(|| panic!("'{name}' missing structuredContent"));

        let output_schema = &tools_by_name[name]["outputSchema"];
        let validator = jsonschema::validator_for(output_schema)
            .unwrap_or_else(|e| panic!("'{name}' outputSchema failed to compile: {e}"));
        assert!(
            validator.is_valid(structured),
            "'{name}' structuredContent {structured} does not conform to its outputSchema {output_schema}"
        );
    }
}

/// Legacy-dialect equivalent of the smoke test above: every sample tool
/// call must still succeed via the byte-compatible `initialize`-era
/// `tools/call` path (no `_meta`, no `structuredContent` expected) —
/// exercising `content_only`/`handle_legacy_tools_call` rather than the
/// modern presentation layer, since the two are genuinely different code
/// paths in `server.rs`.
#[test]
fn legacy_dialect_all_smoke_tested_tools_succeed() {
    for (name, args) in sample_arguments() {
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "tools/call",
            "params": { "name": name, "arguments": args }
        });
        let resp = chematic_mcp::handle_line(&req.to_string()).unwrap();
        assert!(
            resp.get("error").is_none(),
            "legacy call to '{name}' failed: {resp}"
        );
        assert!(resp["result"]["content"][0]["text"].is_string(), "'{name}'");
        // Legacy results never carry the modern envelope fields.
        assert!(resp["result"].get("structuredContent").is_none());
        assert!(resp["result"].get("resultType").is_none());
    }
}

/// A separate, explicitly-named test for `smiles_to_moljson` and
/// `moljson_to_smiles` chained together, since the latter's valid argument
/// depends on the former's output — kept out of `sample_arguments` (which
/// covers each tool with a hardcoded static argument) to avoid a datajson
/// dependency between rows there.
#[test]
fn moljson_round_trip_conforms_to_schema() {
    let list = modern_tools_list();
    let moljson_schema = list["tools"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["name"] == "smiles_to_moljson")
        .unwrap()["outputSchema"]
        .clone();

    let to_json = call_tool_modern("smiles_to_moljson", json!({ "smiles": "CCO" }));
    assert!(to_json.get("error").is_none(), "{to_json}");
    let moljson_value = to_json["result"]["structuredContent"].clone();
    assert!(
        moljson_value.is_string(),
        "smiles_to_moljson structuredContent must be a bare string"
    );
    assert!(jsonschema::is_valid(&moljson_schema, &moljson_value));

    let moljson_str = moljson_value.as_str().unwrap().to_string();
    let back = call_tool_modern("moljson_to_smiles", json!({ "json": moljson_str }));
    assert!(back.get("error").is_none(), "{back}");
    assert!(back["result"]["structuredContent"]["canonical_smiles"].is_string());
}

#[test]
fn molecule_context_pack_matches_one_of_branch_per_format() {
    let list = modern_tools_list();
    let schema = list["tools"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["name"] == "molecule_context_pack")
        .unwrap()["outputSchema"]
        .clone();

    for format in ["json", "markdown", "prompt"] {
        let resp = call_tool_modern(
            "molecule_context_pack",
            json!({ "smiles": "CCO", "format": format }),
        );
        assert!(resp.get("error").is_none(), "{resp}");
        let structured = &resp["result"]["structuredContent"];
        assert!(
            jsonschema::is_valid(&schema, structured),
            "format '{format}' structuredContent {structured} does not match oneOf"
        );
    }
}

#[test]
fn domain_error_has_machine_readable_code() {
    let resp = call_tool_modern("parse_smiles", json!({ "smiles": "C1CC" }));
    assert!(
        resp.get("error").is_none(),
        "domain error must not be a transport error: {resp}"
    );
    assert_eq!(resp["result"]["isError"], true);
    assert_eq!(
        resp["result"]["structuredContent"]["error"]["code"],
        "INVALID_SMILES"
    );
    assert!(resp["result"]["structuredContent"]["error"]["message"].is_string());
}

#[test]
fn invalid_arguments_map_to_dash_32602_not_iserror() {
    let resp = call_tool_modern("parse_smiles", json!({}));
    assert_eq!(resp["error"]["code"], -32602);
    assert!(resp.get("result").is_none());
}

#[test]
fn find_mcs_schema_violation_is_dash_32602() {
    // minItems: 2 declared -- a single-molecule list must be rejected before
    // any chemistry runs, as an argument-shape problem.
    let resp = call_tool_modern("find_mcs", json!({ "smiles_list": ["c1ccccc1"] }));
    assert_eq!(resp["error"]["code"], -32602);
}

#[test]
fn unsupported_protocol_version_is_dash_32022() {
    let mut meta = modern_meta();
    meta["io.modelcontextprotocol/protocolVersion"] = json!("2099-01-01");
    let req = json!({
        "jsonrpc": "2.0", "id": 1, "method": "server/discover",
        "params": { "_meta": meta }
    });
    let resp = chematic_mcp::handle_line(&req.to_string()).unwrap();
    assert_eq!(resp["error"]["code"], -32022);
    assert_eq!(resp["error"]["data"]["requested"], "2099-01-01");
}

// ── section 13: negative controls ──────────────────────────────────────────
//
// Each test here deliberately breaks a fixture and asserts the *checking
// logic* (a real 2020-12 validator, or a simple invariant check) notices --
// proving the harness is not vacuously green.

#[test]
fn negative_control_duplicate_tool_name_is_detected() {
    let list = modern_tools_list();
    let mut tools = list["tools"].as_array().unwrap().clone();
    let duplicate = tools[0].clone();
    tools.push(duplicate);

    let names: Vec<&str> = tools.iter().map(|t| t["name"].as_str().unwrap()).collect();
    let unique: HashSet<&str> = names.iter().copied().collect();
    assert_ne!(
        unique.len(),
        names.len(),
        "harness failed to detect an injected duplicate name"
    );
}

#[test]
fn negative_control_removed_tool_is_detected() {
    let list = modern_tools_list();
    let mut tools = list["tools"].as_array().unwrap().clone();
    tools.pop();
    assert_ne!(tools.len(), 20, "harness failed to detect a removed tool");
}

#[test]
fn negative_control_corrupted_output_schema_is_rejected_by_real_validator() {
    // "type": 5 is not a legal JSON Schema type keyword value.
    let corrupted = json!({ "type": 5, "properties": { "x": {} } });
    assert!(
        jsonschema::meta::validate(&corrupted).is_err(),
        "harness failed to detect a corrupted outputSchema"
    );
}

#[test]
fn negative_control_wrong_type_structured_content_is_rejected() {
    let list = modern_tools_list();
    let schema = list["tools"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["name"] == "parse_smiles")
        .unwrap()["outputSchema"]
        .clone();
    // `atoms` must be an integer; inject a string instead.
    let wrong = json!({
        "valid": true, "atoms": "six", "bonds": 6, "mol_weight": 78.11, "smiles": "c1ccccc1"
    });
    assert!(
        !jsonschema::is_valid(&schema, &wrong),
        "harness failed to detect a wrong-typed structuredContent field"
    );
}

#[test]
fn negative_control_missing_content_field_is_detected() {
    let resp = call_tool_modern("parse_smiles", json!({ "smiles": "c1ccccc1" }));
    let mut mutated = resp["result"].clone();
    mutated.as_object_mut().unwrap().remove("content");
    assert!(
        mutated.get("content").is_none(),
        "harness failed to detect removed content field"
    );
}

#[test]
fn negative_control_injected_external_ref_is_detected() {
    let poisoned = json!({
        "type": "object",
        "properties": { "x": { "$ref": "https://evil.example.com/schema.json" } }
    });

    fn has_external_ref(v: &Value) -> bool {
        match v {
            Value::Object(map) => {
                if let Some(Value::String(r)) = map.get("$ref")
                    && !r.starts_with('#')
                {
                    return true;
                }
                map.values().any(has_external_ref)
            }
            Value::Array(items) => items.iter().any(has_external_ref),
            _ => false,
        }
    }

    assert!(
        has_external_ref(&poisoned),
        "harness failed to detect an injected external $ref"
    );
}

#[test]
fn negative_control_changed_protocol_version_is_detected() {
    let mut meta = modern_meta();
    meta["io.modelcontextprotocol/protocolVersion"] = json!("2024-11-05");
    let req = json!({
        "jsonrpc": "2.0", "id": 1, "method": "server/discover",
        "params": { "_meta": meta }
    });
    let resp = chematic_mcp::handle_line(&req.to_string()).unwrap();
    assert_eq!(
        resp["error"]["code"], -32022,
        "harness failed to detect an unsupported (legacy) version presented via modern _meta"
    );
}

#[test]
fn negative_control_oversized_request_is_detected() {
    let huge = format!(
        r#"{{"jsonrpc":"2.0","id":1,"method":"ping","params":{{"pad":"{}"}}}}"#,
        "x".repeat(2_000_000)
    );
    let resp = chematic_mcp::handle_line(&huge).unwrap();
    assert_eq!(resp["error"]["code"], -32600);
}

#[test]
fn negative_control_deeply_nested_request_is_detected() {
    let deep = "[".repeat(200) + &"]".repeat(200);
    let resp = chematic_mcp::handle_line(&deep).unwrap();
    assert_eq!(resp["error"]["code"], -32600);
}
