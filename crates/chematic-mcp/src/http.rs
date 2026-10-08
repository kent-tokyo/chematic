//! Opt-in Streamable HTTP transport (MCP 2026-07-28 stateless dialect, #779).
//!
//! A self-hosted adapter, not a hosted service: one MCP endpoint (`/mcp`)
//! that accepts POST, one JSON-RPC request per HTTP request, answered with a
//! single `application/json` object. There are no sessions, no GET stream
//! and no SSE (this server sends no request-scoped notifications), as the
//! 2026-07-28 revision of the transport allows (`basic/transports/
//! streamable-http.mdx` in the specification repository).
//!
//! The adapter reuses the stdio stack: after HTTP-level checks it hands the
//! body to a fresh [`Connection`] (so every request is classified and
//! validated on its own, as the stateless dialect requires) and the same
//! [`crate::McpServer`] and tool registry answer it. Nothing about tool
//! behaviour lives here.
//!
//! HTTP-level rules, in the order they are checked:
//!
//! 1. Only `/mcp` exists (404 otherwise); only POST is allowed (405 with
//!    `Allow: POST` for GET, DELETE and anything else).
//! 2. Bound to a loopback address, the `Host` header must name a loopback
//!    host, and an `Origin` header must be a loopback origin
//!    (`http(s)://localhost|127.0.0.1|[::1][:port]`) or one of
//!    `--allowed-origin` (403 otherwise; DNS rebinding). Bound elsewhere,
//!    only `--allowed-origin` values pass. No CORS headers are ever sent.
//! 3. `Content-Type` must be `application/json` (415), `Accept` (when sent)
//!    must admit `application/json` (406), the body must carry a
//!    `Content-Length` (411; chunked bodies are not read) of at most
//!    [`crate::protocol::MAX_REQUEST_BYTES`] (413).
//! 4. The body must be one JSON-RPC object (batches are not part of the
//!    transport; 400).
//! 5. For a request (an `id`-bearing message), the mirrored headers must
//!    match the body: `MCP-Protocol-Version` against `_meta`'s
//!    `io.modelcontextprotocol/protocolVersion`, `Mcp-Method` against
//!    `method`, and `Mcp-Name` (Base64 sentinel `=?base64?…?=` decoded)
//!    against `params.name` for `tools/call`. A missing, malformed or
//!    different header is `400` with JSON-RPC error `-32020`
//!    (`HeaderMismatch`); a protocol version this server does not speak is
//!    `400` with `-32022` and the supported list. A request without
//!    `_meta.io.modelcontextprotocol/protocolVersion` is `400` with `-32602`
//!    first: the legacy `initialize` dialect is not served over HTTP.
//! 6. The JSON-RPC response's status: `200` for a result (a tool's typed
//!    chemistry error is a result with `isError`, exactly as on stdio),
//!    `404` for `-32601`, `400` for parse, invalid-request, invalid-params
//!    (arguments that fail the tool's input schema) and the three
//!    protocol-defined errors, `500` for an internal error. The body is the
//!    same JSON-RPC response stdio writes. A notification gets `202` and no
//!    body.
//!
//! Limits: at most `max_concurrency` connections are served at once (more
//! get `503` with `Retry-After`), request headers are bounded, and reads and
//! writes time out, so a slow client cannot hold a worker forever. Every
//! response closes its connection.

use std::io::{self, BufRead, BufReader, Write};
use std::net::{IpAddr, Ipv4Addr, SocketAddr, TcpListener, TcpStream};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

use serde_json::{Value, json};

use crate::protocol;
use crate::transport::Connection;

/// The MCP endpoint path.
pub const MCP_ENDPOINT: &str = "/mcp";
/// Largest request line accepted.
const MAX_REQUEST_LINE_BYTES: usize = 8 * 1024;
/// Largest header block accepted (all header lines together).
const MAX_HEADER_BYTES: usize = 32 * 1024;
/// Most header fields accepted.
const MAX_HEADERS: usize = 100;

/// Options of the Streamable HTTP adapter.
#[derive(Debug, Clone)]
pub struct HttpConfig {
    /// Address to bind. Loopback by default.
    pub bind: IpAddr,
    /// Port to bind; `0` picks a free one ([`HttpServer::local_addr`]).
    pub port: u16,
    /// Required to bind a non-loopback address: this first slice has no
    /// authentication, so anyone who can reach the port can call the tools.
    pub allow_non_loopback: bool,
    /// Extra browser origins allowed to call the endpoint (exact
    /// `scheme://host[:port]` strings). Empty by default: on a loopback bind
    /// only loopback origins pass, on any other bind every request carrying
    /// an `Origin` header is refused.
    pub allowed_origins: Vec<String>,
    /// Connections served at once; further connections get `503`.
    pub max_concurrency: usize,
    /// Read and write timeout per connection.
    pub io_timeout: Duration,
}

impl Default for HttpConfig {
    fn default() -> Self {
        HttpConfig {
            bind: IpAddr::V4(Ipv4Addr::LOCALHOST),
            port: 3000,
            allow_non_loopback: false,
            allowed_origins: Vec::new(),
            max_concurrency: 16,
            io_timeout: Duration::from_secs(30),
        }
    }
}

/// A bound Streamable HTTP adapter; [`HttpServer::serve`] runs it.
pub struct HttpServer {
    listener: TcpListener,
    config: Arc<HttpConfig>,
}

impl HttpServer {
    /// Bind the adapter. Refuses a non-loopback address unless
    /// `allow_non_loopback` is set, and a `max_concurrency` of 0.
    pub fn bind(config: HttpConfig) -> io::Result<Self> {
        if !config.bind.is_loopback() && !config.allow_non_loopback {
            return Err(io::Error::new(
                io::ErrorKind::PermissionDenied,
                format!(
                    "refusing to bind {}: the HTTP transport has no authentication; \
                     pass --allow-non-loopback to expose it beyond this host",
                    config.bind
                ),
            ));
        }
        if config.max_concurrency == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "max_concurrency must be at least 1",
            ));
        }
        let listener = TcpListener::bind(SocketAddr::new(config.bind, config.port))?;
        Ok(HttpServer {
            listener,
            config: Arc::new(config),
        })
    }

    /// The address actually bound (the port when `port` was 0).
    pub fn local_addr(&self) -> io::Result<SocketAddr> {
        self.listener.local_addr()
    }

    /// Accept connections until the listener fails. Each connection is
    /// served on its own thread, at most `max_concurrency` at once.
    pub fn serve(self) -> io::Result<()> {
        let active = Arc::new(AtomicUsize::new(0));
        for stream in self.listener.incoming() {
            let mut stream = match stream {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("chematic-mcp: accept error: {e}");
                    continue;
                }
            };
            let _ = stream.set_read_timeout(Some(self.config.io_timeout));
            let _ = stream.set_write_timeout(Some(self.config.io_timeout));
            if active.fetch_add(1, Ordering::SeqCst) >= self.config.max_concurrency {
                active.fetch_sub(1, Ordering::SeqCst);
                let response = HttpResponse::json(
                    503,
                    &protocol::error_response(
                        &Value::Null,
                        protocol::INTERNAL_ERROR,
                        "server busy: too many concurrent requests",
                        Some(json!({ "maxConcurrency": self.config.max_concurrency })),
                    ),
                )
                .with_header("Retry-After", "1");
                let _ = response.write_to(&mut stream);
                continue;
            }
            let config = Arc::clone(&self.config);
            let active = Arc::clone(&active);
            std::thread::spawn(move || {
                struct Release(Arc<AtomicUsize>);
                impl Drop for Release {
                    fn drop(&mut self) {
                        self.0.fetch_sub(1, Ordering::SeqCst);
                    }
                }
                let _release = Release(active);
                serve_connection(stream, &config);
            });
        }
        Ok(())
    }
}

/// Bind and serve with `config`, printing the endpoint to stderr.
pub fn run_http(config: HttpConfig) -> io::Result<()> {
    let server = HttpServer::bind(config)?;
    let addr = server.local_addr()?;
    eprintln!(
        "chematic-mcp: Streamable HTTP (MCP 2026-07-28) listening on http://{addr}{MCP_ENDPOINT} \
         (no authentication)"
    );
    server.serve()
}

/// A parsed HTTP request: method, path, headers and body.
#[derive(Debug, Clone, Default)]
pub struct HttpRequest {
    pub method: String,
    pub path: String,
    /// Header names lower-cased; values as received.
    pub headers: Vec<(String, String)>,
    pub body: Vec<u8>,
}

impl HttpRequest {
    /// The first header named `name` (case-insensitive).
    pub fn header(&self, name: &str) -> Option<&str> {
        let name = name.to_ascii_lowercase();
        self.headers
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
    }

    fn header_count(&self, name: &str) -> usize {
        let name = name.to_ascii_lowercase();
        self.headers.iter().filter(|(n, _)| *n == name).count()
    }
}

/// An HTTP response the adapter writes.
#[derive(Debug, Clone)]
pub struct HttpResponse {
    pub status: u16,
    pub headers: Vec<(String, String)>,
    pub body: Vec<u8>,
}

impl HttpResponse {
    fn json(status: u16, body: &Value) -> Self {
        HttpResponse {
            status,
            headers: vec![("Content-Type".into(), "application/json".into())],
            body: serde_json::to_vec(body).unwrap_or_default(),
        }
    }

    fn empty(status: u16) -> Self {
        HttpResponse {
            status,
            headers: Vec::new(),
            body: Vec::new(),
        }
    }

    fn with_header(mut self, name: &str, value: &str) -> Self {
        self.headers.push((name.into(), value.into()));
        self
    }

    /// The response body parsed as JSON, if it is JSON.
    pub fn json_body(&self) -> Option<Value> {
        serde_json::from_slice(&self.body).ok()
    }

    fn write_to<W: Write>(&self, out: &mut W) -> io::Result<()> {
        let mut head = format!("HTTP/1.1 {} {}\r\n", self.status, reason(self.status));
        for (name, value) in &self.headers {
            head.push_str(&format!("{name}: {value}\r\n"));
        }
        head.push_str(&format!("Content-Length: {}\r\n", self.body.len()));
        head.push_str("Cache-Control: no-store\r\n");
        head.push_str("X-Content-Type-Options: nosniff\r\n");
        head.push_str("Connection: close\r\n\r\n");
        out.write_all(head.as_bytes())?;
        out.write_all(&self.body)?;
        out.flush()
    }
}

fn reason(status: u16) -> &'static str {
    match status {
        200 => "OK",
        202 => "Accepted",
        400 => "Bad Request",
        403 => "Forbidden",
        404 => "Not Found",
        405 => "Method Not Allowed",
        406 => "Not Acceptable",
        408 => "Request Timeout",
        411 => "Length Required",
        413 => "Content Too Large",
        415 => "Unsupported Media Type",
        431 => "Request Header Fields Too Large",
        501 => "Not Implemented",
        503 => "Service Unavailable",
        _ => "Error",
    }
}

/// A JSON-RPC error with no `id`, for failures before the body is read.
fn http_error(status: u16, code: i64, message: impl Into<String>) -> HttpResponse {
    HttpResponse::json(
        status,
        &protocol::error_response(&Value::Null, code, message, None),
    )
}

fn serve_connection(stream: TcpStream, config: &HttpConfig) {
    let mut writer = match stream.try_clone() {
        Ok(w) => w,
        Err(_) => return,
    };
    let mut reader = BufReader::new(stream);
    let response = match read_request(&mut reader) {
        Ok(request) => handle_http_request(&request, config),
        Err(response) => response,
    };
    let _ = response.write_to(&mut writer);
}

/// Read one bounded line ending in `\n`; `None` at end of input.
fn read_line<R: BufRead>(reader: &mut R, max: usize) -> Result<Option<Vec<u8>>, HttpResponse> {
    let mut line = Vec::new();
    loop {
        let buffer = match reader.fill_buf() {
            Ok(b) => b,
            Err(e)
                if matches!(
                    e.kind(),
                    io::ErrorKind::WouldBlock | io::ErrorKind::TimedOut
                ) =>
            {
                return Err(http_error(
                    408,
                    protocol::INVALID_REQUEST,
                    "request timed out",
                ));
            }
            Err(_) => return Err(http_error(400, protocol::INVALID_REQUEST, "read error")),
        };
        if buffer.is_empty() {
            return Ok(if line.is_empty() { None } else { Some(line) });
        }
        let newline = buffer.iter().position(|b| *b == b'\n');
        let take = newline.map_or(buffer.len(), |i| i + 1);
        if line.len() + take > max {
            return Err(http_error(
                431,
                protocol::INVALID_REQUEST,
                "request line or headers too large",
            ));
        }
        line.extend_from_slice(&buffer[..take]);
        reader.consume(take);
        if newline.is_some() {
            while matches!(line.last(), Some(b'\n' | b'\r')) {
                line.pop();
            }
            return Ok(Some(line));
        }
    }
}

/// Read and bound one HTTP/1.1 request (request line, headers, a
/// `Content-Length` body). An `Err` is the response to send instead.
pub fn read_request<R: BufRead>(reader: &mut R) -> Result<HttpRequest, HttpResponse> {
    let line = read_line(reader, MAX_REQUEST_LINE_BYTES)?
        .ok_or_else(|| http_error(400, protocol::INVALID_REQUEST, "empty request"))?;
    let line = String::from_utf8(line)
        .map_err(|_| http_error(400, protocol::INVALID_REQUEST, "request line is not UTF-8"))?;
    let mut parts = line.split(' ');
    let (Some(method), Some(target), Some(version), None) =
        (parts.next(), parts.next(), parts.next(), parts.next())
    else {
        return Err(http_error(
            400,
            protocol::INVALID_REQUEST,
            "malformed request line",
        ));
    };
    if !version.starts_with("HTTP/1.") {
        return Err(http_error(
            400,
            protocol::INVALID_REQUEST,
            "only HTTP/1.x is served",
        ));
    }
    let path = target.split('?').next().unwrap_or("").to_string();

    let mut headers = Vec::new();
    let mut header_bytes = 0usize;
    loop {
        let line = read_line(reader, MAX_HEADER_BYTES)?
            .ok_or_else(|| http_error(400, protocol::INVALID_REQUEST, "headers ended early"))?;
        if line.is_empty() {
            break;
        }
        header_bytes += line.len();
        if header_bytes > MAX_HEADER_BYTES || headers.len() >= MAX_HEADERS {
            return Err(http_error(
                431,
                protocol::INVALID_REQUEST,
                "request headers too large",
            ));
        }
        let line = String::from_utf8(line)
            .map_err(|_| http_error(400, protocol::INVALID_REQUEST, "header is not UTF-8"))?;
        let Some((name, value)) = line.split_once(':') else {
            return Err(http_error(
                400,
                protocol::INVALID_REQUEST,
                "malformed header",
            ));
        };
        if name.is_empty() || !name.bytes().all(is_tchar) {
            return Err(http_error(
                400,
                protocol::INVALID_REQUEST,
                "malformed header name",
            ));
        }
        headers.push((
            name.to_ascii_lowercase(),
            value.trim_matches(|c| c == ' ' || c == '\t').to_string(),
        ));
    }

    let mut request = HttpRequest {
        method: method.to_string(),
        path,
        headers,
        body: Vec::new(),
    };
    if request.method != "POST" {
        return Ok(request);
    }
    if request.header("transfer-encoding").is_some() {
        return Err(http_error(
            411,
            protocol::INVALID_REQUEST,
            "send the body with Content-Length; chunked bodies are not read",
        ));
    }
    if request.header_count("content-length") > 1 {
        return Err(http_error(
            400,
            protocol::INVALID_REQUEST,
            "more than one Content-Length",
        ));
    }
    let length: usize = match request.header("content-length") {
        None => {
            return Err(http_error(
                411,
                protocol::INVALID_REQUEST,
                "Content-Length is required",
            ));
        }
        Some(v) => v
            .parse()
            .map_err(|_| http_error(400, protocol::INVALID_REQUEST, "malformed Content-Length"))?,
    };
    if length > protocol::MAX_REQUEST_BYTES {
        return Err(http_error(
            413,
            protocol::INVALID_REQUEST,
            format!(
                "request exceeds maximum size of {} bytes",
                protocol::MAX_REQUEST_BYTES
            ),
        ));
    }
    let mut body = vec![0u8; length];
    reader.read_exact(&mut body).map_err(|e| match e.kind() {
        io::ErrorKind::WouldBlock | io::ErrorKind::TimedOut => {
            http_error(408, protocol::INVALID_REQUEST, "request timed out")
        }
        _ => http_error(
            400,
            protocol::INVALID_REQUEST,
            "body shorter than Content-Length",
        ),
    })?;
    request.body = body;
    Ok(request)
}

/// RFC 9110 `tchar`.
fn is_tchar(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b"!#$%&'*+-.^_`|~".contains(&b)
}

/// Visible ASCII, space and tab only (RFC 9110 field values).
fn safe_header_value(value: &str) -> bool {
    value
        .bytes()
        .all(|b| b == b'\t' || (0x20..=0x7e).contains(&b))
}

/// Whether `host` (a `Host` header value) names a loopback host.
fn loopback_host(host: &str) -> bool {
    let name = if let Some(rest) = host.strip_prefix('[') {
        rest.split(']').next().unwrap_or("")
    } else {
        host.rsplit_once(':').map_or(host, |(h, port)| {
            if port.bytes().all(|b| b.is_ascii_digit()) {
                h
            } else {
                host
            }
        })
    };
    name.eq_ignore_ascii_case("localhost")
        || name.parse::<IpAddr>().is_ok_and(|ip| ip.is_loopback())
}

/// Whether `origin` (an `Origin` header value) is `http` or `https` on a
/// loopback host, with an optional port. A page served from this machine is
/// as trusted as the local process that can already reach the port; a
/// rebinding page (`http://evil.example`) is not.
fn loopback_origin(origin: &str) -> bool {
    let rest = origin
        .strip_prefix("http://")
        .or_else(|| origin.strip_prefix("https://"));
    rest.is_some_and(|host| !host.is_empty() && !host.contains(['/', '@']) && loopback_host(host))
}

/// Decode the `=?base64?…?=` sentinel form of a mirrored header value.
fn decode_header_value(value: &str) -> Option<String> {
    match value
        .strip_prefix("=?base64?")
        .and_then(|v| v.strip_suffix("?="))
    {
        Some(encoded) => String::from_utf8(base64_decode(encoded)?).ok(),
        None => Some(value.to_string()),
    }
}

fn base64_decode(input: &str) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(input.len() * 3 / 4);
    let mut buf = 0u32;
    let mut bits = 0u32;
    let mut padding = 0;
    for c in input.bytes() {
        let v = match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => {
                padding += 1;
                continue;
            }
            _ => return None,
        };
        if padding > 0 {
            return None;
        }
        buf = (buf << 6) | u32::from(v);
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buf >> bits) as u8);
            buf &= (1 << bits) - 1;
        }
    }
    (padding <= 2 && input.len().is_multiple_of(4)).then_some(out)
}

fn header_mismatch(id: &Value, message: impl Into<String>) -> HttpResponse {
    HttpResponse::json(
        400,
        &protocol::error_response(id, protocol::HEADER_MISMATCH, message, None),
    )
}

/// Answer one parsed request (everything after reading it off the socket).
pub fn handle_http_request(request: &HttpRequest, config: &HttpConfig) -> HttpResponse {
    if request.path != MCP_ENDPOINT {
        return http_error(
            404,
            protocol::INVALID_REQUEST,
            format!("not found: the MCP endpoint is {MCP_ENDPOINT}"),
        );
    }
    if request.method != "POST" {
        // GET (the pre-2026-07-28 stream) and DELETE (session end) are not
        // part of this revision; OPTIONS gets no CORS preflight answer.
        return http_error(
            405,
            protocol::INVALID_REQUEST,
            "only POST is served at the MCP endpoint",
        )
        .with_header("Allow", "POST");
    }
    if config.bind.is_loopback() {
        match request.header("host") {
            Some(host) if loopback_host(host) => {}
            _ => {
                return http_error(
                    403,
                    protocol::INVALID_REQUEST,
                    "Host must name this loopback server",
                );
            }
        }
    }
    if let Some(origin) = request.header("origin")
        && !(config.bind.is_loopback() && loopback_origin(origin))
        && !config.allowed_origins.iter().any(|o| o == origin)
    {
        return http_error(403, protocol::INVALID_REQUEST, "Origin not allowed");
    }
    let media = |v: &str| {
        v.split(';')
            .next()
            .unwrap_or("")
            .trim()
            .to_ascii_lowercase()
    };
    if request.header("content-type").map(media).as_deref() != Some("application/json") {
        return http_error(
            415,
            protocol::INVALID_REQUEST,
            "Content-Type must be application/json",
        );
    }
    if let Some(accept) = request.header("accept") {
        let admits = accept
            .split(',')
            .map(media)
            .any(|m| m == "application/json" || m == "application/*" || m == "*/*");
        if !admits {
            return http_error(
                406,
                protocol::INVALID_REQUEST,
                "Accept must admit application/json",
            );
        }
    }

    let Ok(body) = std::str::from_utf8(&request.body) else {
        return http_error(400, protocol::PARSE_ERROR, "body is not UTF-8");
    };
    if let Err(msg) = protocol::prescan_raw_request(body) {
        return http_error(400, protocol::INVALID_REQUEST, msg);
    }
    let parsed: Value = match serde_json::from_str(body) {
        Ok(v) => v,
        Err(e) => return http_error(400, protocol::PARSE_ERROR, e.to_string()),
    };
    if !parsed.is_object() {
        return http_error(
            400,
            protocol::INVALID_REQUEST,
            "the body must be one JSON-RPC message (batches are not supported)",
        );
    }

    if let Some(id) = parsed.get("id")
        && let Some(response) = validate_mirrored_headers(request, &parsed, id)
    {
        return response;
    }

    match Connection::new().handle_line(body) {
        None => HttpResponse::empty(202),
        Some(response) => {
            let status = match response["error"]["code"].as_i64() {
                Some(protocol::METHOD_NOT_FOUND) => 404,
                Some(
                    protocol::PARSE_ERROR
                    | protocol::INVALID_REQUEST
                    | protocol::INVALID_PARAMS
                    | protocol::HEADER_MISMATCH
                    | protocol::MISSING_REQUIRED_CLIENT_CAPABILITY
                    | protocol::UNSUPPORTED_PROTOCOL_VERSION,
                ) => 400,
                Some(_) => 500,
                None => 200,
            };
            HttpResponse::json(status, &response)
        }
    }
}

/// The request-metadata header checks of the transport (see the module
/// doc); `Some(response)` rejects the request.
fn validate_mirrored_headers(
    request: &HttpRequest,
    body: &Value,
    id: &Value,
) -> Option<HttpResponse> {
    for name in ["mcp-protocol-version", "mcp-method", "mcp-name"] {
        if request.header_count(name) > 1 {
            return Some(header_mismatch(id, format!("more than one {name} header")));
        }
        if request.header(name).is_some_and(|v| !safe_header_value(v)) {
            return Some(header_mismatch(
                id,
                format!("{name} header has invalid characters"),
            ));
        }
    }
    // The stateless dialect requires `_meta.protocolVersion` on every
    // request; a body without it is malformed (`-32602`, 400) whatever the
    // headers say. This is also where a legacy `initialize` body, which
    // carries no `_meta`, is turned away: HTTP serves 2026-07-28 only.
    let params = body.get("params").cloned().unwrap_or_else(|| json!({}));
    let body_version = match params.get("_meta") {
        Some(meta) if meta.is_object() => meta.get(protocol::META_PROTOCOL_VERSION),
        _ => None,
    };
    let Some(body_version) = body_version.and_then(Value::as_str) else {
        return Some(HttpResponse::json(
            400,
            &protocol::error_response(
                id,
                protocol::INVALID_PARAMS,
                format!(
                    "missing required _meta key: {} (Streamable HTTP serves only the \
                     2026-07-28 stateless dialect)",
                    protocol::META_PROTOCOL_VERSION
                ),
                None,
            ),
        ));
    };
    let Some(version) = request.header("mcp-protocol-version") else {
        return Some(header_mismatch(id, "missing MCP-Protocol-Version header"));
    };
    if !protocol::SUPPORTED_MODERN_VERSIONS.contains(&version) {
        return Some(HttpResponse::json(
            400,
            &protocol::error_response(
                id,
                protocol::UNSUPPORTED_PROTOCOL_VERSION,
                "Unsupported protocol version",
                Some(json!({
                    "supported": protocol::SUPPORTED_MODERN_VERSIONS,
                    "requested": version
                })),
            ),
        ));
    }
    if body_version != version {
        return Some(header_mismatch(
            id,
            format!(
                "Header mismatch: MCP-Protocol-Version header value '{version}' does not match \
                 body value '{body_version}'"
            ),
        ));
    }
    let method = body.get("method").and_then(Value::as_str).unwrap_or("");
    match request.header("mcp-method") {
        None => return Some(header_mismatch(id, "missing Mcp-Method header")),
        Some(m) if m != method => {
            return Some(header_mismatch(
                id,
                format!(
                    "Header mismatch: Mcp-Method header value '{m}' does not match body value \
                     '{method}'"
                ),
            ));
        }
        Some(_) => {}
    }
    let name_field = match method {
        "tools/call" | "prompts/get" => Some("name"),
        "resources/read" => Some("uri"),
        _ => None,
    };
    if let Some(field) = name_field {
        let body_name = params.get(field).and_then(Value::as_str).unwrap_or("");
        match request.header("mcp-name").map(decode_header_value) {
            None => return Some(header_mismatch(id, "missing Mcp-Name header")),
            Some(None) => {
                return Some(header_mismatch(id, "Mcp-Name header is not valid Base64"));
            }
            Some(Some(n)) if n != body_name => {
                return Some(header_mismatch(
                    id,
                    format!(
                        "Header mismatch: Mcp-Name header value '{n}' does not match body value \
                         '{body_name}'"
                    ),
                ));
            }
            Some(Some(_)) => {}
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loopback_hosts() {
        for host in [
            "localhost",
            "localhost:3000",
            "127.0.0.1:8080",
            "[::1]:3000",
            "127.0.0.2",
        ] {
            assert!(loopback_host(host), "{host}");
        }
        for host in [
            "example.com",
            "evil.com:3000",
            "10.0.0.1",
            "localhost.evil.com",
        ] {
            assert!(!loopback_host(host), "{host}");
        }
    }

    #[test]
    fn loopback_origins() {
        for origin in [
            "http://localhost",
            "http://localhost:3000",
            "https://127.0.0.1:8443",
            "http://[::1]:3000",
        ] {
            assert!(loopback_origin(origin), "{origin}");
        }
        for origin in [
            "null",
            "http://evil.com",
            "http://localhost.evil.com",
            "http://localhost@evil.com",
            "file://localhost",
            "http://",
            "http://localhost/path",
        ] {
            assert!(!loopback_origin(origin), "{origin}");
        }
    }

    #[test]
    fn base64_sentinel() {
        assert_eq!(
            decode_header_value("parse_smiles").as_deref(),
            Some("parse_smiles")
        );
        assert_eq!(
            decode_header_value("=?base64?SGVsbG8sIOS4lueVjA==?=").as_deref(),
            Some("Hello, 世界")
        );
        assert_eq!(decode_header_value("=?base64?***?="), None);
    }

    #[test]
    fn non_loopback_bind_needs_the_flag() {
        let config = HttpConfig {
            bind: IpAddr::V4(Ipv4Addr::UNSPECIFIED),
            port: 0,
            ..HttpConfig::default()
        };
        let error = HttpServer::bind(config).err().unwrap();
        assert_eq!(error.kind(), io::ErrorKind::PermissionDenied);
    }
}
