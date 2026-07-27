//! Server core: dispatches `server/discover` / `tools/list` / `tools/call` /
//! `ping` / `initialize`, and shapes the response per protocol era.
//!
//! This module is the only place that knows both "which era" and "what the
//! wire shape looks like" — tool implementations (`tools.rs`) never see an
//! era, and the transport layer (`transport.rs`) never shapes a response
//! body, only decides which era a request belongs to.

use std::panic::{self, AssertUnwindSafe};

use serde_json::{Value, json};

use crate::protocol::{self, ClientMeta, JsonRpcRequest, ProtocolEra, RequestContext};
use crate::schema;
use crate::tools::{self, ToolCallError};

/// Implementation-defined (JSON-RPC `-32000..-32019` sub-range) code used
/// by the legacy era for tool-call failures. This is the exact code the
/// pre-refactor server used for every `tools/call` error — kept unchanged
/// for byte-compatibility. Not part of the 2026-07-28 vocabulary (that
/// revision distinguishes `-32602` argument errors from `isError: true`
/// domain errors instead — see `handle_modern_tools_call`).
const LEGACY_TOOL_ERROR: i64 = -32000;

/// Stateless MCP server core. Holds no connection/session state — the only
/// state in this crate's process is the transport layer's era pin
/// (`transport::Connection`), which this type never sees or touches.
#[derive(Default)]
pub struct McpServer;

impl McpServer {
    pub fn new() -> Self {
        McpServer
    }

    /// Handle one JSON-RPC request, already classified into `ctx.era`.
    pub fn handle_request(&self, request: JsonRpcRequest<'_>, ctx: &RequestContext) -> Value {
        let JsonRpcRequest { id, method, params } = request;
        match ctx.era {
            ProtocolEra::Legacy => handle_legacy(id, method, params),
            ProtocolEra::Modern20260728 => handle_modern(id, method, params, ctx),
        }
    }
}

// ── legacy era (byte-compatible with the pre-2026-07-28 server) ──────────

fn handle_legacy(id: &Value, method: &str, params: &Value) -> Value {
    match method {
        "initialize" => protocol::result_response(
            id,
            json!({
                "protocolVersion": protocol::LEGACY_PROTOCOL_VERSION,
                "capabilities": { "tools": {} },
                "serverInfo": {
                    "name": "chematic-mcp",
                    "version": env!("CARGO_PKG_VERSION")
                }
            }),
        ),
        "tools/list" => protocol::result_response(id, tools::list_tools()),
        "tools/call" => handle_legacy_tools_call(id, params),
        "ping" => protocol::result_response(id, json!({})),
        _ => protocol::error_response(
            id,
            protocol::METHOD_NOT_FOUND,
            format!("Method not found: {method}"),
            None,
        ),
    }
}

fn handle_legacy_tools_call(id: &Value, params: &Value) -> Value {
    let name = match params.get("name").and_then(|v| v.as_str()) {
        Some(n) => n,
        None => {
            return protocol::error_response(
                id,
                protocol::INVALID_PARAMS,
                "Missing tool name in params",
                None,
            );
        }
    };
    let arguments = params.get("arguments").cloned().unwrap_or(Value::Null);

    match run_tool_catching_panics(name, &arguments) {
        Ok(Ok(payload)) => protocol::result_response(id, content_only(&payload)),
        Ok(Err(e)) => protocol::error_response(id, LEGACY_TOOL_ERROR, e.legacy_message(), None),
        Err(()) => protocol::error_response(id, protocol::INTERNAL_ERROR, "internal error", None),
    }
}

/// Legacy-era `tools/call` success envelope: unchanged from the
/// pre-refactor server — a single text-content block whose `text` is the
/// JSON-stringified payload.
fn content_only(payload: &Value) -> Value {
    json!({ "content": [{ "type": "text", "text": payload.to_string() }] })
}

// ── modern (2026-07-28) era ────────────────────────────────────────────────

fn handle_modern(id: &Value, method: &str, params: &Value, ctx: &RequestContext) -> Value {
    let meta = match &ctx.client_meta {
        Some(m) => m,
        None => {
            return protocol::error_response(
                id,
                protocol::INVALID_PARAMS,
                "modern request missing required _meta triple (protocolVersion, clientInfo, clientCapabilities)",
                None,
            );
        }
    };

    if let Some(err) = check_protocol_version(id, meta) {
        return err;
    }

    match method {
        "server/discover" => protocol::result_response(id, discover_result()),
        "tools/list" => protocol::result_response(id, modern_tools_list()),
        "tools/call" => handle_modern_tools_call(id, params),
        _ => protocol::error_response(
            id,
            protocol::METHOD_NOT_FOUND,
            format!("Method not found: {method}"),
            None,
        ),
    }
}

fn check_protocol_version(id: &Value, meta: &ClientMeta) -> Option<Value> {
    if meta.protocol_version == protocol::MODERN_PROTOCOL_VERSION {
        return None;
    }
    Some(protocol::error_response(
        id,
        protocol::UNSUPPORTED_PROTOCOL_VERSION,
        "Unsupported protocol version",
        Some(json!({
            "supported": protocol::SUPPORTED_MODERN_VERSIONS,
            "requested": meta.protocol_version
        })),
    ))
}

/// `server/discover` result. `serverInfo` is kept as a top-level field
/// (matching the pinned RC tag's `DiscoverResult` shape and the real
/// `rmcp` v3.0.0-beta.2 SDK's struct) even though the untagged post-RC
/// `main` branch of the spec deletes it in favor of a `_meta.serverInfo`
/// convention — see `docs/mcp/2026-07-28-implementation-rfc.md` for the
/// three-way comparison this decision is based on. `ttlMs`/`cacheScope`
/// (from `CacheableResult`) are included per that same main-branch change,
/// which `rmcp` also adopted.
///
/// Declares only the `tools` capability, with an empty `extensions` map —
/// resources/prompts/sampling/roots/logging/tasks/MCP Apps/subscriptions
/// are all genuinely unimplemented and must not be advertised (section 5).
fn discover_result() -> Value {
    json!({
        "resultType": "complete",
        "supportedVersions": protocol::SUPPORTED_MODERN_VERSIONS,
        "capabilities": {
            "tools": {},
            "extensions": {}
        },
        "serverInfo": {
            "name": "chematic-mcp",
            "version": env!("CARGO_PKG_VERSION")
        },
        "ttlMs": TOOLS_LIST_TTL_MS,
        "cacheScope": "public"
    })
}

/// TTL for the (static, per-process) tool registry: 24h. If the registry
/// ever changes, it changes as part of a new `chematic-mcp` release — a
/// running server process serves one fixed registry for its entire
/// lifetime, so this TTL can never go stale mid-process; a client that
/// restarts the server (a new process/version) gets a fresh, uncached
/// response regardless of a previously-cached TTL window. See
/// `crates/chematic-mcp/README.md`'s "Tool registry caching" section.
const TOOLS_LIST_TTL_MS: i64 = 86_400_000; // 24 hours

fn modern_tools_list() -> Value {
    let mut result = tools::list_tools();
    result["resultType"] = json!("complete");
    result["ttlMs"] = json!(TOOLS_LIST_TTL_MS);
    result["cacheScope"] = json!("public");
    result
}

fn handle_modern_tools_call(id: &Value, params: &Value) -> Value {
    let name = match params.get("name").and_then(|v| v.as_str()) {
        Some(n) => n,
        None => {
            return protocol::error_response(
                id,
                protocol::INVALID_PARAMS,
                "Missing tool name in params",
                None,
            );
        }
    };
    let arguments = params.get("arguments").cloned().unwrap_or(json!({}));

    let input_schema = match find_input_schema(name) {
        Some(s) => s,
        None => {
            return protocol::error_response(
                id,
                protocol::INVALID_PARAMS,
                format!("Unknown tool: {name}"),
                None,
            );
        }
    };
    if let Err(msg) = schema::validate(&input_schema, &arguments) {
        return protocol::error_response(
            id,
            protocol::INVALID_PARAMS,
            format!("invalid arguments for tool '{name}': {msg}"),
            None,
        );
    }

    match run_tool_catching_panics(name, &arguments) {
        Ok(Ok(payload)) => protocol::result_response(id, success_result(&payload)),
        Ok(Err(ToolCallError::InvalidArgs(msg))) => {
            protocol::error_response(id, protocol::INVALID_PARAMS, msg, None)
        }
        Ok(Err(e @ ToolCallError::Domain { .. })) => {
            protocol::result_response(id, domain_error_result(&e))
        }
        Err(()) => protocol::error_response(id, protocol::INTERNAL_ERROR, "internal error", None),
    }
}

fn find_input_schema(name: &str) -> Option<Value> {
    tools::list_tools()["tools"]
        .as_array()?
        .iter()
        .find_map(|t| {
            if t.get("name").and_then(|v| v.as_str()) == Some(name) {
                t.get("inputSchema").cloned()
            } else {
                None
            }
        })
}

/// Modern-era successful `CallToolResult`: the same JSON-stringified text
/// content as the legacy era (matching the RC schema's own
/// `result-with-structured-content.json` example fixture, which uses the
/// stringified payload as `content[0].text` rather than separate prose),
/// plus the machine-readable `structuredContent` and `resultType`.
fn success_result(payload: &Value) -> Value {
    json!({
        "resultType": "complete",
        "content": [{ "type": "text", "text": payload.to_string() }],
        "structuredContent": payload
    })
}

/// Modern-era domain-error `CallToolResult`: `isError: true`, human-readable
/// message in `content`, and a machine-readable `structuredContent.error`
/// object (section 9). This is a *successful* JSON-RPC response — the
/// error is reported inside the tool result, not as a transport error.
fn domain_error_result(err: &ToolCallError) -> Value {
    json!({
        "resultType": "complete",
        "content": [{ "type": "text", "text": err.legacy_message() }],
        "structuredContent": { "error": err.to_structured_error() },
        "isError": true
    })
}

/// Run a tool call, converting a panic (section 9/12: "server
/// implementation panic equivalent -> -32603") into `Err(())`. Requires the
/// default `panic = "unwind"` strategy — verified not overridden by any
/// `[profile]` in the workspace root `Cargo.toml`.
fn run_tool_catching_panics(
    name: &str,
    arguments: &Value,
) -> Result<Result<Value, ToolCallError>, ()> {
    panic::catch_unwind(AssertUnwindSafe(|| tools::call_tool(name, arguments))).map_err(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protocol::{ClientMeta, ProtocolEra, RequestContext};

    fn modern_ctx() -> RequestContext {
        RequestContext {
            era: ProtocolEra::Modern20260728,
            client_meta: Some(ClientMeta {
                protocol_version: protocol::MODERN_PROTOCOL_VERSION.to_string(),
                client_info: json!({ "name": "test", "version": "1.0" }),
                client_capabilities: json!({}),
                log_level: None,
            }),
        }
    }

    fn legacy_ctx() -> RequestContext {
        RequestContext {
            era: ProtocolEra::Legacy,
            client_meta: None,
        }
    }

    fn req<'a>(id: &'a Value, method: &'a str, params: &'a Value) -> JsonRpcRequest<'a> {
        JsonRpcRequest { id, method, params }
    }

    #[test]
    fn legacy_initialize_unchanged() {
        let server = McpServer::new();
        let resp = server.handle_request(req(&json!(1), "initialize", &json!({})), &legacy_ctx());
        assert_eq!(resp["result"]["protocolVersion"], "2024-11-05");
        assert!(resp.get("error").is_none());
    }

    #[test]
    fn legacy_tool_error_uses_dash_32000() {
        let server = McpServer::new();
        let params = json!({ "name": "parse_smiles", "arguments": { "smiles": "C1CC" } });
        let resp = server.handle_request(req(&json!(1), "tools/call", &params), &legacy_ctx());
        assert_eq!(resp["error"]["code"], LEGACY_TOOL_ERROR);
    }

    /// Same as above but for the *other* `ToolCallError` variant
    /// (`InvalidArgs`, e.g. a missing required argument rather than a bad
    /// SMILES string). The modern era splits these two cases into
    /// `-32602`/`isError:true` respectively (see `handle_modern_tools_call`),
    /// but the legacy wire shape predates that taxonomy and must keep
    /// collapsing both into a single `-32000`, unchanged from the
    /// pre-refactor server.
    #[test]
    fn legacy_invalid_args_also_uses_dash_32000_not_dash_32602() {
        let server = McpServer::new();
        let params = json!({ "name": "parse_smiles", "arguments": {} });
        let resp = server.handle_request(req(&json!(1), "tools/call", &params), &legacy_ctx());
        assert_eq!(resp["error"]["code"], LEGACY_TOOL_ERROR);
    }

    #[test]
    fn modern_discover_reports_tools_only() {
        let server = McpServer::new();
        let resp =
            server.handle_request(req(&json!(1), "server/discover", &json!({})), &modern_ctx());
        assert_eq!(resp["result"]["resultType"], "complete");
        assert_eq!(resp["result"]["capabilities"]["tools"], json!({}));
        assert!(resp["result"]["capabilities"].get("resources").is_none());
        assert!(resp["result"]["capabilities"].get("prompts").is_none());
        assert_eq!(resp["result"]["supportedVersions"], json!(["2026-07-28"]));
    }

    #[test]
    fn modern_tools_list_has_cache_hints() {
        let server = McpServer::new();
        let resp = server.handle_request(req(&json!(1), "tools/list", &json!({})), &modern_ctx());
        assert_eq!(resp["result"]["ttlMs"], TOOLS_LIST_TTL_MS);
        assert_eq!(resp["result"]["cacheScope"], "public");
        assert_eq!(resp["result"]["resultType"], "complete");
        assert_eq!(resp["result"]["tools"].as_array().unwrap().len(), 20);
    }

    #[test]
    fn modern_domain_error_is_successful_result_with_is_error() {
        let server = McpServer::new();
        let params = json!({ "name": "parse_smiles", "arguments": { "smiles": "C1CC" } });
        let resp = server.handle_request(req(&json!(1), "tools/call", &params), &modern_ctx());
        assert!(
            resp.get("error").is_none(),
            "domain errors must not be transport errors: {resp}"
        );
        assert_eq!(resp["result"]["isError"], true);
        assert_eq!(
            resp["result"]["structuredContent"]["error"]["code"],
            "INVALID_SMILES"
        );
    }

    #[test]
    fn modern_missing_argument_is_invalid_params() {
        let server = McpServer::new();
        let params = json!({ "name": "parse_smiles", "arguments": {} });
        let resp = server.handle_request(req(&json!(1), "tools/call", &params), &modern_ctx());
        assert_eq!(resp["error"]["code"], protocol::INVALID_PARAMS);
    }

    #[test]
    fn modern_unsupported_version_reports_supported_list() {
        let server = McpServer::new();
        let mut ctx = modern_ctx();
        if let Some(meta) = ctx.client_meta.as_mut() {
            meta.protocol_version = "1900-01-01".to_string();
        }
        let resp = server.handle_request(req(&json!(1), "server/discover", &json!({})), &ctx);
        assert_eq!(
            resp["error"]["code"],
            protocol::UNSUPPORTED_PROTOCOL_VERSION
        );
        assert_eq!(resp["error"]["data"]["requested"], "1900-01-01");
        assert_eq!(resp["error"]["data"]["supported"], json!(["2026-07-28"]));
    }

    #[test]
    fn modern_ping_is_method_not_found() {
        // `ping` was removed in 2026-07-28 (changelog.mdx "Major changes" #5
        // at the pinned RC commit); it is not a dialect violation, just an
        // unrecognized method under the modern vocabulary.
        let server = McpServer::new();
        let resp = server.handle_request(req(&json!(1), "ping", &json!({})), &modern_ctx());
        assert_eq!(resp["error"]["code"], protocol::METHOD_NOT_FOUND);
    }

    #[test]
    fn modern_success_structured_content_matches_content_text() {
        let server = McpServer::new();
        let params = json!({ "name": "parse_smiles", "arguments": { "smiles": "c1ccccc1" } });
        let resp = server.handle_request(req(&json!(1), "tools/call", &params), &modern_ctx());
        assert_eq!(resp["result"]["structuredContent"]["atoms"], 6);
        let text = resp["result"]["content"][0]["text"].as_str().unwrap();
        let reparsed: Value = serde_json::from_str(text).unwrap();
        assert_eq!(reparsed["atoms"], 6);
    }

    #[test]
    fn header_mismatch_and_missing_capability_codes_are_serializable() {
        // HEADER_MISMATCH and MISSING_REQUIRED_CLIENT_CAPABILITY are never
        // emitted by this stdio-only, zero-required-capability server (see
        // `protocol.rs` doc comments) — but the codec vocabulary must still
        // be able to *shape* these errors correctly for a future HTTP
        // adapter / a server that does require a capability. This is a
        // codec-capability check, not a claim that either code is reachable
        // from this server's actual request handling.
        let header_mismatch = protocol::error_response(
            &json!(1),
            protocol::HEADER_MISMATCH,
            "header mismatch",
            None,
        );
        assert_eq!(header_mismatch["error"]["code"], -32020);

        let missing_cap = protocol::error_response(
            &json!(1),
            protocol::MISSING_REQUIRED_CLIENT_CAPABILITY,
            "missing required capability",
            Some(json!({ "requiredCapabilities": { "elicitation": {} } })),
        );
        assert_eq!(missing_cap["error"]["code"], -32021);
        assert_eq!(
            missing_cap["error"]["data"]["requiredCapabilities"]["elicitation"],
            json!({})
        );
    }
}
