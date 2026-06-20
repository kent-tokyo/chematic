//! chematic-mcp — MCP (Model Context Protocol) server library.
//!
//! Exposes chematic's cheminformatics capabilities as MCP tools
//! callable by AI agents (Claude, GPT-4, etc.) via JSON-RPC 2.0 over stdio.
//!
//! ## Usage (as binary)
//!
//! ```text
//! cargo run -p chematic-mcp
//! ```
//!
//! Register in `.claude/mcp_settings.json`:
//!
//! ```json
//! {
//!   "mcpServers": {
//!     "chematic": {
//!       "command": "chematic-mcp"
//!     }
//!   }
//! }
//! ```
//!
//! ## Available tools
//!
//! | Tool | Description |
//! |------|-------------|
//! | `parse_smiles` | Parse SMILES → atoms, bonds, MW |
//! | `calc_properties` | LogP, TPSA, MW, HBD, HBA, QED |
//! | `ecfp4` | ECFP4 fingerprint (2048-bit hex) |
//! | `tanimoto` | ECFP4 Tanimoto similarity |
//! | `smarts_match` | Substructure search |
//! | `canonical_smiles` | Canonicalize SMILES |
//! | `find_mcs` | Maximum common substructure |
//! | `generate_3d` | 3D coordinates (XYZ) |
//! | `retrosynthesis` | One-step BRICS disconnection, ranked by SA Score |

#![forbid(unsafe_code)]

mod tools;

use serde_json::{Value, json};

/// Process a single JSON-RPC 2.0 request line.
///
/// Returns `Some(response)` when a reply must be written to the client.
/// Returns `None` for notifications (no `id` field) — no reply is sent.
pub fn handle_line(line: &str) -> Option<Value> {
    let req: Value =
        serde_json::from_str(line).unwrap_or_else(|e| json!({ "_parse_error": e.to_string() }));

    // JSON parse failure → -32700
    if let Some(err) = req.get("_parse_error") {
        return Some(json!({
            "jsonrpc": "2.0",
            "id": null,
            "error": { "code": -32700, "message": err }
        }));
    }

    let id = req.get("id");

    // Notifications (no id field): consume but do not respond.
    let id = id?;

    let method = match req.get("method").and_then(|v| v.as_str()) {
        Some(m) => m,
        None => {
            return Some(json!({
                "jsonrpc": "2.0",
                "id": id,
                "error": { "code": -32600, "message": "Invalid request: missing method" }
            }));
        }
    };

    let (result, is_err) = dispatch(method, &req);

    Some(if is_err {
        json!({ "jsonrpc": "2.0", "id": id, "error": result })
    } else {
        json!({ "jsonrpc": "2.0", "id": id, "result": result })
    })
}

fn dispatch(method: &str, req: &Value) -> (Value, bool) {
    match method {
        "initialize" => (
            json!({
                "protocolVersion": "2024-11-05",
                "capabilities": { "tools": {} },
                "serverInfo": {
                    "name": "chematic-mcp",
                    "version": env!("CARGO_PKG_VERSION")
                }
            }),
            false,
        ),

        "tools/list" => (tools::list_tools(), false),

        "tools/call" => {
            let params = req.get("params").unwrap_or(&Value::Null);
            let name = match params.get("name").and_then(|v| v.as_str()) {
                Some(n) => n,
                None => {
                    return (
                        json!({ "code": -32602, "message": "Missing tool name in params" }),
                        true,
                    );
                }
            };
            let arguments = params.get("arguments").unwrap_or(&Value::Null);
            match tools::call_tool(name, arguments) {
                Ok(result) => (result, false),
                Err(msg) => (json!({ "code": -32000, "message": msg }), true),
            }
        }

        "ping" => (json!({}), false),

        _ => (
            json!({ "code": -32601, "message": format!("Method not found: {method}") }),
            true,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialize() {
        let req = r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{}}}"#;
        let resp = handle_line(req).unwrap();
        assert_eq!(resp["jsonrpc"], "2.0");
        assert_eq!(resp["id"], 1);
        assert_eq!(resp["result"]["protocolVersion"], "2024-11-05");
        assert!(resp.get("error").is_none());
    }

    #[test]
    fn test_tools_list() {
        let req = r#"{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp["result"]["tools"].is_array());
        let count = resp["result"]["tools"].as_array().unwrap().len();
        assert_eq!(count, 16);
    }

    #[test]
    fn test_tools_call_parse_smiles() {
        let req = r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp.get("error").is_none(), "unexpected error: {resp}");
        let text = resp["result"]["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["atoms"], 6);
    }

    #[test]
    fn test_tools_call_calc_properties() {
        let req = r#"{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"calc_properties","arguments":{"smiles":"c1ccccc1"}}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp.get("error").is_none());
        let text = resp["result"]["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert!(v["mw"].as_f64().unwrap() > 78.0);
        assert_eq!(v["hbd"], 0);
    }

    #[test]
    fn test_notification_no_response() {
        let req = r#"{"jsonrpc":"2.0","method":"initialized","params":{}}"#;
        let resp = handle_line(req);
        assert!(
            resp.is_none(),
            "notifications should not produce a response"
        );
    }

    #[test]
    fn test_invalid_json() {
        let req = "not valid json{{{{";
        let resp = handle_line(req).unwrap();
        let code = resp["error"]["code"].as_i64().unwrap();
        assert_eq!(code, -32700);
    }

    #[test]
    fn test_unknown_method() {
        let req = r#"{"jsonrpc":"2.0","id":9,"method":"nonexistent/method","params":{}}"#;
        let resp = handle_line(req).unwrap();
        let code = resp["error"]["code"].as_i64().unwrap();
        assert_eq!(code, -32601);
    }

    #[test]
    fn test_ping() {
        let req = r#"{"jsonrpc":"2.0","id":5,"method":"ping","params":{}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp.get("error").is_none());
        assert_eq!(resp["result"], json!({}));
    }
}
