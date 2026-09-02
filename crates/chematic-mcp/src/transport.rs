//! stdio transport: newline-delimited JSON-RPC framing, plus the
//! connection-pinned protocol-era state that section 10 explicitly
//! whitelists as the *only* allowed transport state (everything else about
//! this server is stateless request/response — see
//! `docs/specification/draft/basic/lifecycle.mdx` at the pinned RC commit).
//!
//! `Connection` owns the era pin. `McpServer` (in `server.rs`) never sees
//! it — it only ever receives an already-classified `RequestContext`.

use std::io::{self, BufRead, Write};

use serde_json::{Value, json};

use crate::protocol::{self, JsonRpcRequest, ProtocolEra, RequestContext};
use crate::server::McpServer;

/// A single stdio connection's lifetime. Pins the protocol era (legacy vs.
/// modern) to whichever dialect the first id-bearing request speaks;
/// notifications never pin (they carry no `id`, so there's nothing to
/// reply to or classify against). A request whose dialect disagrees with
/// an already-set pin is rejected with a typed protocol error rather than
/// silently switching — per the task brief's explicit requirement — while
/// every individual modern request still supplies and is validated against
/// its own full `_meta` triple; the pin never substitutes a cached value
/// for a later request's own metadata (that would violate the protocol's
/// statelessness, not just the brief).
pub struct Connection {
    server: McpServer,
    pinned_era: Option<ProtocolEra>,
}

/// Read one newline-delimited frame without allocating beyond the protocol
/// request limit. A frame that exceeds the limit is fatal to the stdio loop;
/// stopping there avoids treating its remainder as a new JSON-RPC request.
fn read_bounded_line<R: BufRead>(reader: &mut R, max_bytes: usize) -> io::Result<Option<Vec<u8>>> {
    let mut line = Vec::with_capacity(max_bytes.min(8192));
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return if line.is_empty() {
                Ok(None)
            } else {
                Ok(Some(line))
            };
        }
        let newline = buffer.iter().position(|byte| *byte == b'\n');
        let available = newline.map_or(buffer.len(), |index| index + 1);
        if line.len().saturating_add(available) > max_bytes.saturating_add(1) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("request exceeds maximum size of {max_bytes} bytes"),
            ));
        }
        line.extend_from_slice(&buffer[..available]);
        reader.consume(available);
        if newline.is_some() {
            return Ok(Some(line));
        }
    }
}

impl Default for Connection {
    fn default() -> Self {
        Self::new()
    }
}

impl Connection {
    pub fn new() -> Self {
        Connection {
            server: McpServer::new(),
            pinned_era: None,
        }
    }

    /// Process one raw JSON-RPC line. Returns `Some(response)` when a reply
    /// must be written back; `None` for notifications and blank lines.
    pub fn handle_line(&mut self, line: &str) -> Option<Value> {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            return None;
        }

        if let Err(msg) = protocol::prescan_raw_request(trimmed) {
            return Some(protocol::error_response(
                &Value::Null,
                protocol::INVALID_REQUEST,
                msg,
                None,
            ));
        }

        let parsed: Value = match serde_json::from_str(trimmed) {
            Ok(v) => v,
            Err(e) => {
                return Some(protocol::error_response(
                    &Value::Null,
                    protocol::PARSE_ERROR,
                    e.to_string(),
                    None,
                ));
            }
        };

        if let Err(msg) = protocol::validate_value_limits(&parsed) {
            return Some(protocol::error_response(
                &Value::Null,
                protocol::INVALID_REQUEST,
                msg,
                None,
            ));
        }

        // No `id` => a notification. Consumed, no response, no era pin.
        let id = parsed.get("id")?.clone();

        let method = match parsed.get("method").and_then(|v| v.as_str()) {
            Some(m) => m.to_string(),
            None => {
                return Some(protocol::error_response(
                    &id,
                    protocol::INVALID_REQUEST,
                    "Invalid request: missing method",
                    None,
                ));
            }
        };

        let params = parsed.get("params").cloned().unwrap_or_else(|| json!({}));

        let dialect = match protocol::classify_dialect(&method, &params) {
            Ok(d) => d,
            Err(msg) => {
                return Some(protocol::error_response(
                    &id,
                    protocol::INVALID_PARAMS,
                    msg,
                    None,
                ));
            }
        };
        let era = match dialect {
            protocol::RequestDialect::Legacy => ProtocolEra::Legacy,
            protocol::RequestDialect::Modern => ProtocolEra::Modern20260728,
        };

        match self.pinned_era {
            None => self.pinned_era = Some(era),
            Some(pinned) if pinned != era => {
                return Some(protocol::error_response(
                    &id,
                    protocol::INVALID_REQUEST,
                    format!(
                        "protocol era mismatch: this connection is pinned to {pinned:?} \
                         (set by an earlier request), but method '{method}' speaks {era:?}. \
                         A single stdio connection may not switch dialects mid-session."
                    ),
                    Some(
                        json!({ "pinnedEra": format!("{pinned:?}"), "requestEra": format!("{era:?}") }),
                    ),
                ));
            }
            Some(_) => {}
        }

        let client_meta = if era == ProtocolEra::Modern20260728 {
            match protocol::parse_client_meta(&params) {
                Ok(m) => m,
                Err(msg) => {
                    return Some(protocol::error_response(
                        &id,
                        protocol::INVALID_PARAMS,
                        msg,
                        None,
                    ));
                }
            }
        } else {
            None
        };

        let ctx = RequestContext { era, client_meta };
        let request = JsonRpcRequest {
            id: &id,
            method: &method,
            params: &params,
        };
        let response = self.server.handle_request(request, &ctx);
        match serde_json::to_vec(&response) {
            Ok(bytes) if bytes.len() <= protocol::MAX_RESPONSE_BYTES => Some(response),
            Ok(_) => Some(protocol::error_response(
                &id,
                protocol::INTERNAL_ERROR,
                format!(
                    "response exceeds maximum size of {} bytes",
                    protocol::MAX_RESPONSE_BYTES
                ),
                None,
            )),
            Err(_) => Some(protocol::error_response(
                &id,
                protocol::INTERNAL_ERROR,
                "response serialization failed",
                None,
            )),
        }
    }
}

/// Run the stdio server loop: read newline-delimited JSON-RPC from `stdin`,
/// write responses to `stdout`, one message per line. Diagnostics only ever
/// go to `stderr` (never `stdout` — section 12's stdout-purity requirement).
pub fn run_stdio() {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let mut connection = Connection::new();

    let mut input = stdin.lock();
    loop {
        let line = match read_bounded_line(&mut input, protocol::MAX_REQUEST_BYTES) {
            Ok(Some(bytes)) => match String::from_utf8(bytes) {
                Ok(line) => line,
                Err(e) => {
                    eprintln!("chematic-mcp: invalid UTF-8 request: {e}");
                    break;
                }
            },
            Ok(None) => break,
            Err(e) => {
                eprintln!("chematic-mcp: read error: {e}");
                break;
            }
        };

        if let Some(response) = connection.handle_line(&line) {
            let serialized = match serde_json::to_string(&response) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("chematic-mcp: serialization error: {e}");
                    continue;
                }
            };
            if let Err(e) = writeln!(out, "{serialized}") {
                eprintln!("chematic-mcp: write error: {e}");
                break;
            }
            if let Err(e) = out.flush() {
                eprintln!("chematic-mcp: flush error: {e}");
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_reader_rejects_an_oversized_frame_before_full_allocation() {
        let mut input = std::io::Cursor::new(b"123456789\nnext\n".to_vec());
        let error = read_bounded_line(&mut input, 8).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("maximum size of 8 bytes"));
    }

    #[test]
    fn legacy_then_legacy_is_fine() {
        let mut conn = Connection::new();
        let r1 = conn
            .handle_line(r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}"#)
            .unwrap();
        assert!(r1.get("error").is_none());
        let r2 = conn
            .handle_line(r#"{"jsonrpc":"2.0","id":2,"method":"ping","params":{}}"#)
            .unwrap();
        assert!(r2.get("error").is_none());
    }

    #[test]
    fn modern_then_legacy_is_rejected() {
        let mut conn = Connection::new();
        let discover = r#"{"jsonrpc":"2.0","id":1,"method":"server/discover","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientInfo":{"name":"t","version":"1"},"io.modelcontextprotocol/clientCapabilities":{}}}}"#;
        let r1 = conn.handle_line(discover).unwrap();
        assert!(r1.get("error").is_none(), "{r1}");

        let legacy_init = r#"{"jsonrpc":"2.0","id":2,"method":"initialize","params":{}}"#;
        let r2 = conn.handle_line(legacy_init).unwrap();
        assert_eq!(r2["error"]["code"], protocol::INVALID_REQUEST);
    }

    #[test]
    fn legacy_then_modern_is_rejected() {
        let mut conn = Connection::new();
        let r1 = conn
            .handle_line(r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}"#)
            .unwrap();
        assert!(r1.get("error").is_none());

        let discover = r#"{"jsonrpc":"2.0","id":2,"method":"server/discover","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientInfo":{"name":"t","version":"1"},"io.modelcontextprotocol/clientCapabilities":{}}}}"#;
        let r2 = conn.handle_line(discover).unwrap();
        assert_eq!(r2["error"]["code"], protocol::INVALID_REQUEST);
    }

    #[test]
    fn stdio_backward_compat_probe_gets_method_not_found_before_any_pin() {
        // A dual-era client probes with server/discover first; per
        // lifecycle.mdx, a purely-legacy server would answer -32601, which
        // is this server's cue that it must be legacy-only. This server
        // supports both, so it must NOT return -32601 here.
        let mut conn = Connection::new();
        let discover = r#"{"jsonrpc":"2.0","id":1,"method":"server/discover","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientInfo":{"name":"t","version":"1"},"io.modelcontextprotocol/clientCapabilities":{}}}}"#;
        let r1 = conn.handle_line(discover).unwrap();
        assert_ne!(r1["error"]["code"], protocol::METHOD_NOT_FOUND);
        assert!(r1.get("error").is_none());
    }

    #[test]
    fn notifications_do_not_pin_era() {
        let mut conn = Connection::new();
        // A notification (no `id`) carrying legacy shape must not pin.
        let notif = r#"{"jsonrpc":"2.0","method":"notifications/initialized","params":{}}"#;
        assert!(conn.handle_line(notif).is_none());

        // The connection should still be free to pin modern on the first
        // real (id-bearing) request.
        let discover = r#"{"jsonrpc":"2.0","id":1,"method":"server/discover","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientInfo":{"name":"t","version":"1"},"io.modelcontextprotocol/clientCapabilities":{}}}}"#;
        let r = conn.handle_line(discover).unwrap();
        assert!(r.get("error").is_none(), "{r}");
    }

    #[test]
    fn oversized_request_rejected() {
        let mut conn = Connection::new();
        let huge = "x".repeat(protocol::MAX_REQUEST_BYTES + 10);
        let r = conn.handle_line(&huge).unwrap();
        assert_eq!(r["error"]["code"], protocol::INVALID_REQUEST);
    }

    #[test]
    fn deeply_nested_request_rejected() {
        let mut conn = Connection::new();
        let deep =
            "[".repeat(protocol::MAX_JSON_DEPTH + 5) + &"]".repeat(protocol::MAX_JSON_DEPTH + 5);
        let r = conn.handle_line(&deep).unwrap();
        assert_eq!(r["error"]["code"], protocol::INVALID_REQUEST);
    }
}
