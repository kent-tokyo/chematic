//! chematic-mcp — MCP (Model Context Protocol) server library.
//!
//! Exposes chematic's cheminformatics capabilities as MCP tools callable by
//! AI agents (Claude, GPT-4, etc.) via JSON-RPC 2.0 over stdio.
//!
//! ## Two protocol eras, one stateless tool registry
//!
//! This server speaks two dialects of MCP on the same stdio transport:
//!
//! - **Legacy** (`2024-11-05`-style): `initialize` handshake,
//!   `tools/list`/`tools/call`/`ping`. Byte-compatible with the
//!   pre-2026-07-28 server — see `crates/chematic-mcp/README.md`.
//! - **Modern** (`2026-07-28` stateless tools-only): per-request `_meta`
//!   metadata (no `initialize`), `server/discover`, `tools/list` with
//!   cache hints, `tools/call` with `structuredContent`.
//!
//! A single stdio connection pins to whichever dialect its first
//! id-bearing request speaks (`transport::Connection`); see
//! `docs/mcp/2026-07-28-implementation-rfc.md` for the full design and
//! primary-source citations.
//!
//! ## Available tools
//!
//! 20 tools total; see `crates/chematic-mcp/README.md` for the full,
//! categorized list with descriptions (kept there, not duplicated here, so
//! there is exactly one place to update when a tool is added or changed).
//! `tools::list_tools()` is the actual runtime source of truth.

#![forbid(unsafe_code)]
// `tools::list_tools()` is one large `serde_json::json!` literal (20 tools'
// worth of input/output JSON Schema documents) — the macro's expansion
// depth exceeds rustc's default recursion limit well before it exceeds any
// runtime resource; this is the standard fix for a large static `json!`
// literal, not evidence of unbounded/runaway recursion in actual logic.
#![recursion_limit = "256"]

mod protocol;
mod schema;
mod server;
mod tools;
mod transport;

pub use protocol::{
    ClientMeta, HEADER_MISMATCH, INTERNAL_ERROR, INVALID_PARAMS, INVALID_REQUEST, JsonRpcRequest,
    METHOD_NOT_FOUND, MISSING_REQUIRED_CLIENT_CAPABILITY, PARSE_ERROR, ProtocolEra, RequestContext,
    UNSUPPORTED_PROTOCOL_VERSION,
};
pub use server::McpServer;
pub use tools::{TOOL_COUNT, ToolCallError};
pub use transport::{Connection, run_stdio};

use serde_json::Value;

/// Process a single JSON-RPC 2.0 request line with no memory of any
/// previous call — dialect (legacy or modern) is classified fresh from
/// this one line, same as the pre-2026-07-28 implementation classified
/// every line independently.
///
/// Kept for backward compatibility with existing callers of this exact
/// function signature. Internally, this is equivalent to creating a fresh
/// [`Connection`] and calling [`Connection::handle_line`] once — since a
/// fresh connection has no era pinned yet, the outcome for any single,
/// independent line is unchanged from the pre-2026-07-28 implementation.
/// A caller that wants era *pinning* across multiple lines on one
/// connection (rejecting a mid-session dialect switch, `server/discover`
/// state) must use [`Connection`] directly — that's what the
/// `chematic-mcp` binary does. Without pinning, a modern request still
/// works correctly through this function (it just isn't protected from
/// a differently-dialected line that follows it).
///
/// Returns `Some(response)` when a reply must be written to the client.
/// Returns `None` for notifications (no `id` field) — no reply is sent.
pub fn handle_line(line: &str) -> Option<Value> {
    Connection::new().handle_line(line)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    // ── Legacy fixtures (section 1/13: frozen, byte-compatible wire behavior) ──

    #[test]
    fn test_initialize() {
        let req = r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{}}}"#;
        let resp = handle_line(req).unwrap();
        assert_eq!(resp["jsonrpc"], "2.0");
        assert_eq!(resp["id"], 1);
        assert_eq!(resp["result"]["protocolVersion"], "2024-11-05");
        assert!(resp.get("error").is_none());
        // Legacy results never carry the modern `resultType` field.
        assert!(resp["result"].get("resultType").is_none());
    }

    #[test]
    fn test_tools_list() {
        let req = r#"{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp["result"]["tools"].is_array());
        let count = resp["result"]["tools"].as_array().unwrap().len();
        assert_eq!(count, 20);
        // Legacy tools/list never carries the modern cache envelope fields.
        assert!(resp["result"].get("ttlMs").is_none());
        assert!(resp["result"].get("cacheScope").is_none());
        assert!(resp["result"].get("resultType").is_none());
    }

    #[test]
    fn test_tools_call_parse_smiles() {
        let req = r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}"#;
        let resp = handle_line(req).unwrap();
        assert!(resp.get("error").is_none(), "unexpected error: {resp}");
        let text = resp["result"]["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["atoms"], 6);
        // Legacy tool results never carry structuredContent/resultType.
        assert!(resp["result"].get("structuredContent").is_none());
        assert!(resp["result"].get("resultType").is_none());
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

    // ── Modern fixtures (section 13) ────────────────────────────────────────

    fn modern_meta() -> Value {
        json!({
            "io.modelcontextprotocol/protocolVersion": "2026-07-28",
            "io.modelcontextprotocol/clientInfo": { "name": "test-client", "version": "1.0.0" },
            "io.modelcontextprotocol/clientCapabilities": {}
        })
    }

    #[test]
    fn test_server_discover_success_no_initialize_required() {
        // No `initialize` call precedes this — modern era must not require it.
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "server/discover",
            "params": { "_meta": modern_meta() }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        assert!(resp.get("error").is_none(), "{resp}");
        assert_eq!(resp["result"]["resultType"], "complete");
        assert_eq!(resp["result"]["supportedVersions"], json!(["2026-07-28"]));
        assert_eq!(resp["result"]["serverInfo"]["name"], "chematic-mcp");
        assert_eq!(resp["result"]["capabilities"]["tools"], json!({}));
        assert!(resp["result"]["capabilities"].get("resources").is_none());
        assert!(resp["result"]["capabilities"].get("prompts").is_none());
        assert!(resp["result"]["capabilities"].get("sampling").is_none());
        assert!(resp["result"]["capabilities"].get("logging").is_none());
    }

    #[test]
    fn test_modern_tools_list_all_20_unique_names_with_cache_hints() {
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "tools/list",
            "params": { "_meta": modern_meta() }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        assert_eq!(resp["result"]["ttlMs"], 86_400_000i64);
        assert_eq!(resp["result"]["cacheScope"], "public");
        let names: Vec<&str> = resp["result"]["tools"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| t["name"].as_str().unwrap())
            .collect();
        assert_eq!(names.len(), 20);
        let unique: std::collections::HashSet<&str> = names.iter().copied().collect();
        assert_eq!(unique.len(), 20, "duplicate tool name found: {names:?}");
    }

    #[test]
    fn test_modern_all_tools_have_output_schema() {
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "tools/list",
            "params": { "_meta": modern_meta() }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        for tool in resp["result"]["tools"].as_array().unwrap() {
            assert!(
                tool.get("outputSchema").is_some(),
                "tool '{}' is missing outputSchema",
                tool["name"]
            );
            assert_eq!(tool["inputSchema"]["type"], "object");
        }
    }

    #[test]
    fn test_modern_tools_call_returns_structured_content() {
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "tools/call",
            "params": {
                "_meta": modern_meta(),
                "name": "parse_smiles",
                "arguments": { "smiles": "c1ccccc1" }
            }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        assert!(resp.get("error").is_none(), "{resp}");
        assert_eq!(resp["result"]["resultType"], "complete");
        assert!(resp["result"]["content"].is_array());
        assert!(!resp["result"]["content"].as_array().unwrap().is_empty());
        assert_eq!(resp["result"]["structuredContent"]["atoms"], 6);
    }

    #[test]
    fn test_modern_unsupported_protocol_version() {
        let mut meta = modern_meta();
        meta["io.modelcontextprotocol/protocolVersion"] = json!("1900-01-01");
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "server/discover",
            "params": { "_meta": meta }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        assert_eq!(resp["error"]["code"], -32022);
    }

    #[test]
    fn test_modern_missing_required_metadata_is_typed_error() {
        // _meta present but missing clientInfo
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "server/discover",
            "params": {
                "_meta": {
                    "io.modelcontextprotocol/protocolVersion": "2026-07-28",
                    "io.modelcontextprotocol/clientCapabilities": {}
                }
            }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        assert_eq!(resp["error"]["code"], -32602);
    }

    #[test]
    fn test_legacy_after_modern_pin_rejected() {
        let mut conn = Connection::new();
        let discover = json!({
            "jsonrpc": "2.0", "id": 1, "method": "server/discover",
            "params": { "_meta": modern_meta() }
        });
        assert!(
            conn.handle_line(&discover.to_string())
                .unwrap()
                .get("error")
                .is_none()
        );

        let legacy = json!({ "jsonrpc": "2.0", "id": 2, "method": "initialize", "params": {} });
        let resp = conn.handle_line(&legacy.to_string()).unwrap();
        assert!(resp.get("error").is_some());
    }

    #[test]
    fn test_modern_after_legacy_pin_rejected() {
        let mut conn = Connection::new();
        let legacy = json!({ "jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {} });
        assert!(
            conn.handle_line(&legacy.to_string())
                .unwrap()
                .get("error")
                .is_none()
        );

        let discover = json!({
            "jsonrpc": "2.0", "id": 2, "method": "server/discover",
            "params": { "_meta": modern_meta() }
        });
        let resp = conn.handle_line(&discover.to_string()).unwrap();
        assert!(resp.get("error").is_some());
    }

    #[test]
    fn test_stdout_purity_response_is_single_json_value_per_line() {
        // Every response this server produces must serialize to exactly one
        // JSON value with no embedded newlines (section 12/stdio framing).
        let req = json!({
            "jsonrpc": "2.0", "id": 1, "method": "tools/call",
            "params": {
                "_meta": modern_meta(),
                "name": "molecule_context_pack",
                "arguments": { "smiles": "CCO", "format": "markdown" }
            }
        });
        let resp = handle_line(&req.to_string()).unwrap();
        let serialized = serde_json::to_string(&resp).unwrap();
        assert!(
            !serialized.contains('\n'),
            "serialized response must not contain a literal newline"
        );
    }
}
