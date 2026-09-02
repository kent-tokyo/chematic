//! Protocol codec: JSON-RPC envelope handling, dual-era classification,
//! `_meta` parsing, error-code vocabulary, and adversarial-input limits.
//!
//! This module knows nothing about chemistry or individual tools — it only
//! knows how to read a JSON-RPC request, classify which MCP "era" it speaks,
//! and shape a JSON-RPC response. See `crates/chematic-mcp/README.md` and
//! `docs/mcp/2026-07-28-implementation-rfc.md` for the primary-source
//! citations behind every constant and shape below.

use serde_json::{Map, Value, json};

/// The MCP protocol version this server implements for the modern era.
///
/// Sourced from `LATEST_PROTOCOL_VERSION` in
/// `schema/draft/schema.ts` at the pinned RC commit
/// (`9d700ed62dcf86cb77475c9b81930611a9182f46`, tag `2026-07-28-RC`).
pub const MODERN_PROTOCOL_VERSION: &str = "2026-07-28";

/// The protocol version the legacy (pre-stateless) era speaks. Unchanged
/// from the server's pre-existing behavior.
pub const LEGACY_PROTOCOL_VERSION: &str = "2024-11-05";

// ── JSON-RPC standard error codes (unchanged across every MCP revision) ──

pub const PARSE_ERROR: i64 = -32700;
pub const INVALID_REQUEST: i64 = -32600;
pub const METHOD_NOT_FOUND: i64 = -32601;
pub const INVALID_PARAMS: i64 = -32602;
pub const INTERNAL_ERROR: i64 = -32603;

/// Error code for an HTTP header / body mismatch.
///
/// **Defined but never emitted by this server.** This is a stdio-only
/// transport; there are no HTTP headers to validate. Reserved here per the
/// task brief's explicit request to keep a typed placeholder for a future
/// HTTP adapter PR, and because the value is real (not invented) — see the
/// provenance note below.
///
/// Value provenance: the RC tag pinned in section 0 (commit `9d700ed`) does
/// **not** define this code (it defines `MISSING_REQUIRED_CLIENT_CAPABILITY
/// = -32003` and `UNSUPPORTED_PROTOCOL_VERSION = -32004` instead, with no
/// header-mismatch code at all). This value, and the two below, are read
/// from the spec repo's untagged `main` branch
/// (`7634684382c3d14cf7e9f14073fe40a2d8ace3fa`, 2026-07-23 — no
/// `2026-07-28` final tag exists yet), corroborated independently by the
/// official `rmcp` Rust SDK v3.0.0-beta.2's `ErrorCode` constants and by
/// this task brief's own section 9, which names these exact three numbers.
/// See `docs/mcp/2026-07-28-implementation-rfc.md` for the full
/// three-way comparison table.
pub const HEADER_MISMATCH: i64 = -32020;

/// Error code for a request that requires a client capability the client
/// did not declare.
///
/// **Defined but never emitted by this server.** All 20 tools are plain
/// request/response with no elicitation, sampling, or roots dependency —
/// this server never requires any client capability, so there is no
/// legitimate trigger for this code. A negative-control test exists for the
/// *codec's* ability to serialize this error shape, but the conformance
/// matrix marks the wire scenario `not_applicable` rather than wiring a fake
/// trigger to turn it green.
pub const MISSING_REQUIRED_CLIENT_CAPABILITY: i64 = -32021;

/// Error code for a request whose protocol version this server does not
/// support. Actively emitted — see [`classify_and_validate`].
pub const UNSUPPORTED_PROTOCOL_VERSION: i64 = -32022;

/// Every protocol version this server can negotiate. `2024-11-05` is only
/// reachable via the legacy `initialize` handshake (never via a `_meta`
/// triple — the modern era did not exist in that revision), so it is
/// listed here for `UNSUPPORTED_PROTOCOL_VERSION`'s `data.supported` array
/// but is never a valid `io.modelcontextprotocol/protocolVersion` value.
pub const SUPPORTED_MODERN_VERSIONS: &[&str] = &[MODERN_PROTOCOL_VERSION];

// ── Reserved `_meta` keys (schema.ts `RequestMetaObject`, RC commit) ──────

pub const META_PROTOCOL_VERSION: &str = "io.modelcontextprotocol/protocolVersion";
pub const META_CLIENT_INFO: &str = "io.modelcontextprotocol/clientInfo";
pub const META_CLIENT_CAPABILITIES: &str = "io.modelcontextprotocol/clientCapabilities";
pub const META_LOG_LEVEL: &str = "io.modelcontextprotocol/logLevel";

/// Which dialect of MCP a connection (or a single request) is speaking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProtocolEra {
    /// `initialize` handshake, `tools/list`/`tools/call`/`ping` without a
    /// modern `_meta` triple. Byte-compatible with the pre-2026-07-28 server.
    Legacy,
    /// Per-request `_meta.io.modelcontextprotocol/*` triple, `server/discover`,
    /// no `initialize`/`ping`.
    Modern20260728,
}

/// Parsed reserved `_meta` fields from a modern-era request.
#[derive(Debug, Clone)]
pub struct ClientMeta {
    pub protocol_version: String,
    pub client_info: Value,
    pub client_capabilities: Value,
    pub log_level: Option<Value>,
}

/// Per-request context threaded from the transport layer into the
/// protocol-agnostic server core. Tool implementations never see this —
/// only the presentation layer (legacy vs. modern response shaping) does.
#[derive(Debug, Clone)]
pub struct RequestContext {
    pub era: ProtocolEra,
    pub client_meta: Option<ClientMeta>,
}

/// Classification of a single incoming request's *dialect* — independent of
/// any previously-pinned era. The transport layer compares this against the
/// connection's pin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequestDialect {
    Legacy,
    Modern,
}

/// A minimal JSON-RPC request view used by the protocol/server layers.
pub struct JsonRpcRequest<'a> {
    pub id: &'a Value,
    pub method: &'a str,
    pub params: &'a Value,
}

/// Build a JSON-RPC error response object.
pub fn error_response(
    id: &Value,
    code: i64,
    message: impl Into<String>,
    data: Option<Value>,
) -> Value {
    let mut error = Map::new();
    error.insert("code".to_string(), json!(code));
    error.insert("message".to_string(), json!(message.into()));
    if let Some(d) = data {
        error.insert("data".to_string(), d);
    }
    json!({ "jsonrpc": "2.0", "id": id, "error": Value::Object(error) })
}

/// Build a JSON-RPC success response object.
pub fn result_response(id: &Value, result: Value) -> Value {
    json!({ "jsonrpc": "2.0", "id": id, "result": result })
}

/// Extract the reserved `_meta` triple from a request's `params`, if present.
///
/// Returns:
/// - `Ok(Some(meta))` if a full, well-typed modern `_meta` triple is present.
/// - `Ok(None)` if `_meta` (or the `protocolVersion` key specifically) is
///   entirely absent — this request is not modern-shaped.
/// - `Err(message)` if `_meta.protocolVersion` is present (so the caller is
///   clearly attempting the modern dialect) but the triple is malformed —
///   e.g. `clientInfo`/`clientCapabilities` missing or wrong-typed. This is
///   an argument-shape problem, mapped to `-32602` by the caller.
pub fn parse_client_meta(params: &Value) -> Result<Option<ClientMeta>, String> {
    let meta = match params.get("_meta") {
        Some(m) if m.is_object() => m,
        Some(_) => return Err("_meta must be an object".to_string()),
        None => return Ok(None),
    };

    let protocol_version = match meta.get(META_PROTOCOL_VERSION) {
        Some(v) => v,
        None => return Ok(None), // no protocolVersion => not modern-shaped at all
    };
    let protocol_version = protocol_version
        .as_str()
        .ok_or_else(|| format!("{META_PROTOCOL_VERSION} must be a string"))?
        .to_string();

    let client_info = meta
        .get(META_CLIENT_INFO)
        .ok_or_else(|| format!("missing required _meta key: {META_CLIENT_INFO}"))?;
    if !client_info.is_object()
        || client_info.get("name").and_then(|v| v.as_str()).is_none()
        || client_info
            .get("version")
            .and_then(|v| v.as_str())
            .is_none()
    {
        return Err(format!(
            "{META_CLIENT_INFO} must be an object with string `name` and `version`"
        ));
    }

    let client_capabilities = meta
        .get(META_CLIENT_CAPABILITIES)
        .ok_or_else(|| format!("missing required _meta key: {META_CLIENT_CAPABILITIES}"))?;
    if !client_capabilities.is_object() {
        return Err(format!("{META_CLIENT_CAPABILITIES} must be an object"));
    }

    let log_level = meta.get(META_LOG_LEVEL).cloned();

    Ok(Some(ClientMeta {
        protocol_version,
        client_info: client_info.clone(),
        client_capabilities: client_capabilities.clone(),
        log_level,
    }))
}

/// Classify which dialect a single request (by method + `_meta` shape) is
/// speaking, without reference to any connection pin.
///
/// The presence of a modern `_meta.protocolVersion` key is the primary
/// signal, checked before any method-name special-casing: a genuine legacy
/// client never sends `_meta` at all, so its presence is a stronger signal
/// than the method name for a method that happens to exist in both eras'
/// vocabularies (or, as with `ping`/`initialize` below, in neither once a
/// modern envelope is attached to them).
///
/// - `server/discover` is modern-only (it did not exist before this
///   revision — legacy clients never send it) regardless of `_meta`.
/// - `initialize`/`ping` were both removed in the modern era (see
///   `docs/specification/draft/changelog.mdx` "Major changes" item 5 at
///   commit `9d700ed`). Sent *without* `_meta`, they're the ordinary legacy
///   methods. Sent *with* a modern `_meta` triple attached (a malformed or
///   confused client), the request is still modern-shaped — it just names
///   a method that doesn't exist in that vocabulary, so it falls through to
///   `-32601 Method not found` inside the modern dispatcher rather than
///   being treated as a dialect violation.
/// - `tools/list` / `tools/call` exist in both eras; disambiguated purely
///   by `_meta` presence, same as everything else.
pub fn classify_dialect(method: &str, params: &Value) -> Result<RequestDialect, String> {
    if parse_client_meta(params)?.is_some() {
        return Ok(RequestDialect::Modern);
    }
    match method {
        "server/discover" => Ok(RequestDialect::Modern),
        _ => Ok(RequestDialect::Legacy),
    }
}

// ── Adversarial-input limits (section 12) ─────────────────────────────────

/// Maximum accepted raw request line length, in bytes. Checked *before*
/// attempting to parse JSON, so an oversized line never reaches the parser.
pub const MAX_REQUEST_BYTES: usize = 1 << 20; // 1 MiB

/// Maximum serialized JSON-RPC response size written to stdio.
pub const MAX_RESPONSE_BYTES: usize = 1 << 20; // 1 MiB

/// Maximum accepted JSON nesting depth (objects and arrays combined).
/// Checked on the raw byte stream before parsing (so an adversarially deep
/// document is rejected before `serde_json` ever recurses into it), and
/// again on the parsed `Value` tree as a second, independent check.
pub const MAX_JSON_DEPTH: usize = 64;

/// Maximum accepted elements in any single JSON array in a request.
pub const MAX_ARRAY_LEN: usize = 10_000;

/// Maximum accepted length (in UTF-8 bytes) of any single JSON string in a
/// request.
pub const MAX_STRING_LEN: usize = 1 << 18; // 256 KiB

/// Scan a raw JSON-RPC line for nesting depth and byte-size problems
/// *before* handing it to `serde_json`. Deliberately simple: a single pass
/// over bytes tracking `{`/`[`/`}`/`]` nesting while skipping over string
/// literals (respecting `\"` escapes), so brackets inside a string payload
/// are not miscounted as structural nesting.
///
/// Returns `Err(message)` on the first violation found.
pub fn prescan_raw_request(line: &str) -> Result<(), String> {
    if line.len() > MAX_REQUEST_BYTES {
        return Err(format!(
            "request exceeds maximum size of {MAX_REQUEST_BYTES} bytes"
        ));
    }

    let mut depth: usize = 0;
    let mut in_string = false;
    let mut escaped = false;

    for b in line.bytes() {
        if in_string {
            if escaped {
                escaped = false;
            } else if b == b'\\' {
                escaped = true;
            } else if b == b'"' {
                in_string = false;
            }
            continue;
        }
        match b {
            b'"' => in_string = true,
            b'{' | b'[' => {
                depth += 1;
                if depth > MAX_JSON_DEPTH {
                    return Err(format!(
                        "request exceeds maximum nesting depth of {MAX_JSON_DEPTH}"
                    ));
                }
            }
            b'}' | b']' => depth = depth.saturating_sub(1),
            _ => {}
        }
    }
    Ok(())
}

/// Recursively validate a parsed JSON `Value` against array/string/depth
/// limits. This is a second, independent check after `prescan_raw_request`
/// (defense in depth — the two checks use different mechanisms), and the
/// only check capable of enforcing per-array and per-string limits (the
/// byte-level prescan cannot distinguish "one huge array" from "one huge
/// string" cheaply).
pub fn validate_value_limits(value: &Value) -> Result<(), String> {
    fn walk(value: &Value, depth: usize) -> Result<(), String> {
        if depth > MAX_JSON_DEPTH {
            return Err(format!(
                "request exceeds maximum nesting depth of {MAX_JSON_DEPTH}"
            ));
        }
        match value {
            Value::String(s) => {
                if s.len() > MAX_STRING_LEN {
                    return Err(format!(
                        "a string value exceeds maximum length of {MAX_STRING_LEN} bytes"
                    ));
                }
                Ok(())
            }
            Value::Array(items) => {
                if items.len() > MAX_ARRAY_LEN {
                    return Err(format!(
                        "an array exceeds maximum length of {MAX_ARRAY_LEN} elements"
                    ));
                }
                for item in items {
                    walk(item, depth + 1)?;
                }
                Ok(())
            }
            Value::Object(map) => {
                for v in map.values() {
                    walk(v, depth + 1)?;
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }
    walk(value, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn modern_meta_json() -> Value {
        json!({
            "_meta": {
                META_PROTOCOL_VERSION: MODERN_PROTOCOL_VERSION,
                META_CLIENT_INFO: { "name": "test", "version": "1.0" },
                META_CLIENT_CAPABILITIES: {}
            }
        })
    }

    #[test]
    fn classify_discover_is_always_modern() {
        assert_eq!(
            classify_dialect("server/discover", &json!({})).unwrap(),
            RequestDialect::Modern
        );
    }

    #[test]
    fn classify_initialize_is_always_legacy() {
        assert_eq!(
            classify_dialect("initialize", &json!({})).unwrap(),
            RequestDialect::Legacy
        );
    }

    #[test]
    fn classify_ping_is_always_legacy() {
        assert_eq!(
            classify_dialect("ping", &json!({})).unwrap(),
            RequestDialect::Legacy
        );
    }

    #[test]
    fn classify_tools_call_without_meta_is_legacy() {
        let params = json!({ "name": "parse_smiles", "arguments": { "smiles": "C" } });
        assert_eq!(
            classify_dialect("tools/call", &params).unwrap(),
            RequestDialect::Legacy
        );
    }

    #[test]
    fn classify_tools_call_with_meta_is_modern() {
        let mut params = modern_meta_json();
        params["name"] = json!("parse_smiles");
        assert_eq!(
            classify_dialect("tools/call", &params).unwrap(),
            RequestDialect::Modern
        );
    }

    #[test]
    fn parse_client_meta_full_triple() {
        let meta = parse_client_meta(&modern_meta_json()).unwrap().unwrap();
        assert_eq!(meta.protocol_version, MODERN_PROTOCOL_VERSION);
    }

    #[test]
    fn parse_client_meta_missing_client_info_errors() {
        let params = json!({
            "_meta": {
                META_PROTOCOL_VERSION: MODERN_PROTOCOL_VERSION,
                META_CLIENT_CAPABILITIES: {}
            }
        });
        assert!(parse_client_meta(&params).is_err());
    }

    #[test]
    fn parse_client_meta_absent_is_none() {
        assert!(parse_client_meta(&json!({})).unwrap().is_none());
    }

    #[test]
    fn prescan_rejects_oversized_request() {
        let huge = "x".repeat(MAX_REQUEST_BYTES + 1);
        assert!(prescan_raw_request(&huge).is_err());
    }

    #[test]
    fn prescan_rejects_deep_nesting() {
        let deep = "[".repeat(MAX_JSON_DEPTH + 1) + &"]".repeat(MAX_JSON_DEPTH + 1);
        assert!(prescan_raw_request(&deep).is_err());
    }

    #[test]
    fn prescan_accepts_brackets_inside_strings() {
        // A string literal full of brackets must not count toward nesting depth.
        let line = format!(r#"{{"a":"{}"}}"#, "[".repeat(MAX_JSON_DEPTH + 5));
        assert!(prescan_raw_request(&line).is_ok());
    }

    #[test]
    fn validate_value_limits_rejects_huge_array() {
        let arr = Value::Array(vec![Value::Null; MAX_ARRAY_LEN + 1]);
        assert!(validate_value_limits(&json!({ "a": arr })).is_err());
    }

    #[test]
    fn validate_value_limits_rejects_huge_string() {
        let s = "x".repeat(MAX_STRING_LEN + 1);
        assert!(validate_value_limits(&json!({ "a": s })).is_err());
    }

    #[test]
    fn validate_value_limits_accepts_normal_request() {
        let v = json!({ "smiles": "c1ccccc1", "list": [1, 2, 3] });
        assert!(validate_value_limits(&v).is_ok());
    }
}
