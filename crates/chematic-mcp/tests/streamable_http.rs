//! The opt-in Streamable HTTP transport (#779), started on an ephemeral
//! loopback port. Every answer that reaches the tool layer is compared with
//! the stdio answer to the same JSON-RPC body, so the two transports cannot
//! drift apart.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{SocketAddr, TcpStream};
use std::process::{Command, Stdio};
use std::time::Duration;

use chematic_mcp::{HttpConfig, HttpServer};
use serde_json::{Value, json};

fn start(config: HttpConfig) -> SocketAddr {
    let server = HttpServer::bind(HttpConfig { port: 0, ..config }).unwrap();
    let addr = server.local_addr().unwrap();
    std::thread::spawn(move || server.serve());
    addr
}

fn default_server() -> SocketAddr {
    start(HttpConfig::default())
}

struct Reply {
    status: u16,
    headers: Vec<(String, String)>,
    body: Vec<u8>,
}

impl Reply {
    fn json(&self) -> Value {
        serde_json::from_slice(&self.body)
            .unwrap_or_else(|e| panic!("{e}: {}", String::from_utf8_lossy(&self.body)))
    }
    fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(n, _)| n.eq_ignore_ascii_case(name))
            .map(|(_, v)| v.as_str())
    }
}

fn raw(addr: SocketAddr, request: &[u8]) -> Reply {
    let mut stream = TcpStream::connect(addr).unwrap();
    stream
        .set_read_timeout(Some(Duration::from_secs(20)))
        .unwrap();
    stream.write_all(request).unwrap();
    let mut reader = BufReader::new(stream);
    let mut status_line = String::new();
    reader.read_line(&mut status_line).unwrap();
    let status: u16 = status_line.split(' ').nth(1).unwrap().parse().unwrap();
    let mut headers = Vec::new();
    let mut length = 0;
    loop {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        let line = line.trim_end();
        if line.is_empty() {
            break;
        }
        let (n, v) = line.split_once(':').unwrap();
        if n.eq_ignore_ascii_case("content-length") {
            length = v.trim().parse().unwrap();
        }
        headers.push((n.to_string(), v.trim().to_string()));
    }
    let mut body = vec![0; length];
    reader.read_exact(&mut body).unwrap();
    Reply {
        status,
        headers,
        body,
    }
}

fn post(addr: SocketAddr, headers: &[(&str, &str)], body: &str) -> Reply {
    let mut request = format!(
        "POST /mcp HTTP/1.1\r\nHost: 127.0.0.1:{}\r\nContent-Length: {}\r\n",
        addr.port(),
        body.len()
    );
    for (n, v) in headers {
        request.push_str(&format!("{n}: {v}\r\n"));
    }
    request.push_str("\r\n");
    request.push_str(body);
    raw(addr, request.as_bytes())
}

fn meta() -> Value {
    json!({
        "io.modelcontextprotocol/protocolVersion": "2026-07-28",
        "io.modelcontextprotocol/clientInfo": { "name": "http-test", "version": "1.0.0" },
        "io.modelcontextprotocol/clientCapabilities": {}
    })
}

fn body(id: i64, method: &str, mut params: Value) -> String {
    params
        .as_object_mut()
        .unwrap()
        .insert("_meta".into(), meta());
    json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params }).to_string()
}

/// A conforming request: every mirrored header set from the body.
fn call(addr: SocketAddr, method: &str, params: Value) -> (Reply, String) {
    let name = params
        .get("name")
        .and_then(Value::as_str)
        .map(str::to_string);
    let text = body(1, method, params);
    let mut headers = vec![
        ("Content-Type", "application/json"),
        ("Accept", "application/json, text/event-stream"),
        ("MCP-Protocol-Version", "2026-07-28"),
        ("Mcp-Method", method),
    ];
    if let Some(n) = name.as_deref() {
        headers.push(("Mcp-Name", n));
    }
    (post(addr, &headers, &text), text)
}

fn stdio_answer(text: &str) -> Value {
    chematic_mcp::handle_line(text).unwrap()
}

#[test]
fn discover_list_and_call_answer_as_stdio_does() {
    let addr = default_server();

    let (reply, text) = call(addr, "server/discover", json!({}));
    assert_eq!(reply.status, 200);
    assert_eq!(reply.header("content-type"), Some("application/json"));
    assert_eq!(reply.json(), stdio_answer(&text));
    assert!(reply.json()["result"].is_object());

    let (reply, text) = call(addr, "tools/list", json!({}));
    assert_eq!(reply.status, 200);
    assert_eq!(reply.json(), stdio_answer(&text));
    assert_eq!(
        reply.json()["result"]["tools"].as_array().unwrap().len(),
        chematic_mcp::TOOL_COUNT
    );

    let (reply, text) = call(
        addr,
        "tools/call",
        json!({ "name": "calc_properties", "arguments": { "smiles": "CC(=O)Oc1ccccc1C(=O)O" } }),
    );
    assert_eq!(reply.status, 200);
    let answer = reply.json();
    assert_eq!(answer, stdio_answer(&text));
    assert!(
        answer["result"]["structuredContent"].is_object(),
        "{answer}"
    );
    assert_ne!(answer["result"]["isError"], true);
}

#[test]
fn invalid_arguments_and_chemistry_errors_are_the_stdio_errors() {
    let addr = default_server();
    // Arguments that fail the input schema.
    let (reply, text) = call(
        addr,
        "tools/call",
        json!({ "name": "parse_smiles", "arguments": {} }),
    );
    assert_eq!(reply.json(), stdio_answer(&text));
    assert_eq!(reply.status, 400);
    assert_eq!(reply.json()["error"]["code"], chematic_mcp::INVALID_PARAMS);
    // A typed chemistry error (an unparsable SMILES).
    let (reply, text) = call(
        addr,
        "tools/call",
        json!({ "name": "parse_smiles", "arguments": { "smiles": "C1CC(" } }),
    );
    assert_eq!(reply.status, 200);
    assert_eq!(reply.json(), stdio_answer(&text));
    assert_eq!(reply.json()["result"]["isError"], true);
    assert_eq!(
        reply.json()["result"]["structuredContent"]["error"]["code"],
        "INVALID_SMILES"
    );
}

#[test]
fn mirrored_headers_must_match_the_body() {
    let addr = default_server();
    let text = body(
        7,
        "tools/call",
        json!({ "name": "parse_smiles", "arguments": { "smiles": "C" } }),
    );
    let base = [("Content-Type", "application/json")];
    let cases: Vec<(Vec<(&str, &str)>, i64)> = vec![
        // No MCP-Protocol-Version.
        (
            vec![("Mcp-Method", "tools/call"), ("Mcp-Name", "parse_smiles")],
            -32020,
        ),
        // Mcp-Method differs.
        (
            vec![
                ("MCP-Protocol-Version", "2026-07-28"),
                ("Mcp-Method", "tools/list"),
                ("Mcp-Name", "parse_smiles"),
            ],
            -32020,
        ),
        // Mcp-Name missing for tools/call.
        (
            vec![
                ("MCP-Protocol-Version", "2026-07-28"),
                ("Mcp-Method", "tools/call"),
            ],
            -32020,
        ),
        // Mcp-Name differs.
        (
            vec![
                ("MCP-Protocol-Version", "2026-07-28"),
                ("Mcp-Method", "tools/call"),
                ("Mcp-Name", "ecfp4"),
            ],
            -32020,
        ),
        // A version this server does not speak.
        (
            vec![
                ("MCP-Protocol-Version", "2025-06-18"),
                ("Mcp-Method", "tools/call"),
                ("Mcp-Name", "parse_smiles"),
            ],
            -32022,
        ),
    ];
    for (headers, code) in cases {
        let all: Vec<(&str, &str)> = base.iter().copied().chain(headers.clone()).collect();
        let reply = post(addr, &all, &text);
        assert_eq!(reply.status, 400, "{headers:?}");
        assert_eq!(reply.json()["error"]["code"], code, "{headers:?}");
        assert_eq!(reply.json()["id"], 7);
    }
    // The Base64 sentinel form of Mcp-Name is decoded before comparing.
    let reply = post(
        addr,
        &[
            ("Content-Type", "application/json"),
            ("MCP-Protocol-Version", "2026-07-28"),
            ("Mcp-Method", "tools/call"),
            ("Mcp-Name", "=?base64?cGFyc2Vfc21pbGVz?="),
        ],
        &text,
    );
    assert_eq!(
        reply.status,
        200,
        "{}",
        String::from_utf8_lossy(&reply.body)
    );
    // The legacy initialize handshake is not served over HTTP.
    let legacy =
        json!({ "jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {} }).to_string();
    let reply = post(
        addr,
        &[
            ("Content-Type", "application/json"),
            ("MCP-Protocol-Version", "2026-07-28"),
            ("Mcp-Method", "initialize"),
        ],
        &legacy,
    );
    assert_eq!(reply.status, 400);
    assert_eq!(reply.json()["error"]["code"], chematic_mcp::INVALID_PARAMS);
    // So is any request whose `_meta` lacks the protocol version, whatever
    // the headers say (missing required `_meta` field: -32602, 400).
    for params in [
        json!({}),
        json!({ "_meta": { "io.modelcontextprotocol/clientCapabilities": {} } }),
    ] {
        let text = json!({ "jsonrpc": "2.0", "id": 9, "method": "tools/list", "params": params })
            .to_string();
        let reply = post(
            addr,
            &[
                ("Content-Type", "application/json"),
                ("MCP-Protocol-Version", "2026-07-28"),
                ("Mcp-Method", "tools/list"),
            ],
            &text,
        );
        assert_eq!(reply.status, 400, "{text}");
        assert_eq!(reply.json()["error"]["code"], chematic_mcp::INVALID_PARAMS);
        assert_eq!(reply.json()["id"], 9);
    }
    // Optional whitespace around a header value is not part of it
    // (RFC 9110 section 5.5).
    let reply = post(
        addr,
        &[
            ("Content-Type", "application/json"),
            ("MCP-Protocol-Version", " 2026-07-28 "),
            ("Mcp-Method", "\ttools/call "),
            ("Mcp-Name", "  parse_smiles  "),
        ],
        &text,
    );
    assert_eq!(
        reply.status,
        200,
        "{}",
        String::from_utf8_lossy(&reply.body)
    );
}

#[test]
fn client_info_is_optional() {
    // `io.modelcontextprotocol/clientInfo` is a SHOULD for clients; a
    // request without it is served, over HTTP and stdio alike.
    let addr = default_server();
    let text = json!({
        "jsonrpc": "2.0", "id": 3, "method": "tools/list",
        "params": { "_meta": {
            "io.modelcontextprotocol/protocolVersion": "2026-07-28",
            "io.modelcontextprotocol/clientCapabilities": {}
        } }
    })
    .to_string();
    let reply = post(
        addr,
        &[
            ("Content-Type", "application/json"),
            ("MCP-Protocol-Version", "2026-07-28"),
            ("Mcp-Method", "tools/list"),
        ],
        &text,
    );
    assert_eq!(reply.status, 200);
    assert_eq!(reply.json(), stdio_answer(&text));
    assert!(reply.json()["result"]["tools"].is_array());
}

#[test]
fn unknown_method_is_404_with_the_jsonrpc_error() {
    let addr = default_server();
    let (reply, text) = call(addr, "resources/list", json!({}));
    assert_eq!(reply.status, 404);
    assert_eq!(
        reply.json()["error"]["code"],
        chematic_mcp::METHOD_NOT_FOUND
    );
    assert_eq!(reply.json(), stdio_answer(&text));
}

#[test]
fn request_limits_and_malformed_requests_are_typed() {
    let addr = default_server();
    // A body larger than the protocol limit is refused before it is read.
    let request = format!(
        "POST /mcp HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\nContent-Length: {}\r\n\r\n",
        (1 << 20) + 1
    );
    let reply = raw(addr, request.as_bytes());
    assert_eq!(reply.status, 413);
    assert_eq!(reply.json()["error"]["code"], chematic_mcp::INVALID_REQUEST);
    // Too deeply nested.
    let deep = "[".repeat(100) + &"]".repeat(100);
    let reply = post(addr, &[("Content-Type", "application/json")], &deep);
    assert_eq!(reply.status, 400);
    assert_eq!(reply.json()["error"]["code"], chematic_mcp::INVALID_REQUEST);
    // Not JSON.
    let reply = post(addr, &[("Content-Type", "application/json")], "{not json");
    assert_eq!(reply.status, 400);
    assert_eq!(reply.json()["error"]["code"], chematic_mcp::PARSE_ERROR);
    // A batch.
    let reply = post(addr, &[("Content-Type", "application/json")], "[{}]");
    assert_eq!(reply.status, 400);
    // No Content-Length.
    let reply = raw(
        addr,
        b"POST /mcp HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n\r\n",
    );
    assert_eq!(reply.status, 411);
    // Wrong media type, unacceptable Accept.
    let (_, text) = call(addr, "tools/list", json!({}));
    assert_eq!(
        post(addr, &[("Content-Type", "text/plain")], &text).status,
        415
    );
    let reply = post(
        addr,
        &[
            ("Content-Type", "application/json"),
            ("Accept", "text/html"),
        ],
        &text,
    );
    assert_eq!(reply.status, 406);
}

#[test]
fn only_post_at_the_endpoint_and_no_cors() {
    let addr = default_server();
    for method in ["GET", "DELETE", "OPTIONS"] {
        let request = format!("{method} /mcp HTTP/1.1\r\nHost: 127.0.0.1\r\n\r\n");
        let reply = raw(addr, request.as_bytes());
        assert_eq!(reply.status, 405, "{method}");
        assert_eq!(reply.header("allow"), Some("POST"));
        assert!(reply.header("access-control-allow-origin").is_none());
    }
    let reply = raw(
        addr,
        b"POST /other HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Length: 0\r\n\r\n",
    );
    assert_eq!(reply.status, 404);
    // A notification is accepted with no body.
    let notification =
        json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }).to_string();
    let reply = post(addr, &[("Content-Type", "application/json")], &notification);
    assert_eq!(reply.status, 202);
    assert!(reply.body.is_empty());
}

#[test]
fn origin_and_host_are_checked_against_dns_rebinding() {
    let addr = start(HttpConfig {
        allowed_origins: vec!["https://app.example".into()],
        ..HttpConfig::default()
    });
    let (_, text) = call(addr, "server/discover", json!({}));
    let headers = |origin: &'static str| {
        vec![
            ("Content-Type", "application/json"),
            ("MCP-Protocol-Version", "2026-07-28"),
            ("Mcp-Method", "server/discover"),
            ("Origin", origin),
        ]
    };
    let reply = post(addr, &headers("http://evil.example"), &text);
    assert_eq!(reply.status, 403);
    assert!(reply.header("access-control-allow-origin").is_none());
    let reply = post(addr, &headers("null"), &text);
    assert_eq!(reply.status, 403);
    // Loopback origins pass on a loopback bind; other origins only when
    // listed with --allowed-origin.
    for origin in [
        "http://localhost:5173",
        "http://127.0.0.1:3000",
        "https://app.example",
    ] {
        let reply = post(addr, &headers(origin), &text);
        assert_eq!(reply.status, 200, "{origin}");
        assert!(reply.header("access-control-allow-origin").is_none());
    }
    let reply = post(addr, &headers("https://other.example"), &text);
    assert_eq!(reply.status, 403);
    // A rebound host name reaching the loopback port.
    let request = format!(
        "POST /mcp HTTP/1.1\r\nHost: evil.example:{}\r\nContent-Type: application/json\r\nMCP-Protocol-Version: 2026-07-28\r\nMcp-Method: server/discover\r\nContent-Length: {}\r\n\r\n{text}",
        addr.port(),
        text.len()
    );
    assert_eq!(raw(addr, request.as_bytes()).status, 403);
}

#[test]
fn concurrency_is_bounded() {
    let addr = start(HttpConfig {
        max_concurrency: 1,
        io_timeout: Duration::from_secs(3),
        ..HttpConfig::default()
    });
    // Hold the only worker with a connection that sends nothing.
    let held = TcpStream::connect(addr).unwrap();
    std::thread::sleep(Duration::from_millis(200));
    let (reply, _) = call(addr, "tools/list", json!({}));
    assert_eq!(reply.status, 503);
    assert_eq!(reply.header("retry-after"), Some("1"));
    // Once the held connection times out, requests are served again.
    std::thread::sleep(Duration::from_secs(4));
    drop(held);
    let (reply, _) = call(addr, "tools/list", json!({}));
    assert_eq!(reply.status, 200);
}

#[test]
fn non_loopback_bind_requires_the_flag() {
    let refused = HttpServer::bind(HttpConfig {
        bind: "0.0.0.0".parse().unwrap(),
        port: 0,
        ..HttpConfig::default()
    });
    assert!(refused.is_err());
}

#[test]
fn the_binary_serves_http_on_request_and_stdio_by_default() {
    let exe = env!("CARGO_BIN_EXE_chematic-mcp");
    let mut child = Command::new(exe)
        .args(["--transport", "streamable-http", "--port", "0"])
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(child.stderr.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    let addr: SocketAddr = line
        .split("http://")
        .nth(1)
        .and_then(|rest| rest.split('/').next())
        .unwrap()
        .parse()
        .unwrap();
    let (reply, text) = call(addr, "tools/list", json!({}));
    let _ = child.kill();
    let _ = child.wait();
    assert_eq!(reply.status, 200);
    assert_eq!(reply.json(), stdio_answer(&text));

    // Refuses a non-loopback address without the flag.
    let status = Command::new(exe)
        .args([
            "--transport",
            "streamable-http",
            "--bind",
            "0.0.0.0",
            "--port",
            "0",
        ])
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(!status.success());
    // HTTP options without the HTTP transport are a usage error.
    let status = Command::new(exe)
        .args(["--port", "3000"])
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert_eq!(status.code(), Some(2));

    // No arguments: the stdio server, unchanged.
    let mut child = Command::new(exe)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let request = body(1, "tools/list", json!({}));
    {
        let mut stdin = child.stdin.take().unwrap();
        writeln!(stdin, "{request}").unwrap();
    }
    let mut out = String::new();
    BufReader::new(child.stdout.take().unwrap())
        .read_line(&mut out)
        .unwrap();
    let _ = child.wait();
    let answer: Value = serde_json::from_str(&out).unwrap();
    assert_eq!(answer, stdio_answer(&request));
}
