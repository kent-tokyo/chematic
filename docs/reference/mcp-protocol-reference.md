# MCP protocol and transport reference

`chematic-mcp` exposes chematic tools to MCP clients. It is a tools-only,
stateless server: stdio is the default transport, and self-hosted Streamable
HTTP is opt-in. Start with the
[`chematic-mcp` README](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mcp/README.md)
for a minimal setup.

## Supported surface

| Capability | Status |
|---|---|
| Legacy `2024-11-05`-style stdio | Supported |
| 2026-07-28 stateless stdio | Supported |
| 2026-07-28 stateless HTTP | Supported with `--transport streamable-http` |
| HTTP sessions, GET stream, or SSE | Not supported |
| Authentication, OAuth, or hosted endpoint | Not provided |
| Resources, prompts, sampling, roots, and tasks | Not advertised or implemented |

The server exposes 20 tools. Nineteen run locally. `name_to_smiles` is the only
tool that sends input to PubChem. The registry in
[`crates/chematic-mcp/src/tools.rs`](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mcp/src/tools.rs)
is the authoritative inventory and schema source.

## Install and run over stdio

```bash
cargo install chematic-mcp --version 1.2.1 --locked
```

Configure an MCP client to launch the installed binary:

```json
{
  "mcpServers": {
    "chematic": {
      "command": "chematic-mcp",
      "args": []
    }
  }
}
```

Contributors can run `cargo run -p chematic-mcp --release`. Requests and
responses are newline-delimited JSON-RPC 2.0. The first request pins the stdio
connection to one protocol dialect; switching dialects later returns a typed
protocol error.

## Protocol dialects

### Legacy stdio

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{}}}
{"jsonrpc":"2.0","method":"notifications/initialized","params":{}}
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}
```

### Modern stateless dialect

Modern requests use `server/discover`, `tools/list`, and `tools/call`. Every
request includes:

- `io.modelcontextprotocol/protocolVersion` — required;
- `io.modelcontextprotocol/clientCapabilities` — required;
- `io.modelcontextprotocol/clientInfo` — optional.

```json
{
  "jsonrpc": "2.0",
  "id": 1,
  "method": "tools/call",
  "params": {
    "_meta": {
      "io.modelcontextprotocol/protocolVersion": "2026-07-28",
      "io.modelcontextprotocol/clientCapabilities": {}
    },
    "name": "parse_smiles",
    "arguments": {"smiles": "c1ccccc1"}
  }
}
```

Modern results add `resultType` and `structuredContent` while retaining MCP
`content`. Tool schemas use JSON Schema 2020-12 and do not resolve external
`$ref` values.

## Streamable HTTP

```bash
chematic-mcp --transport streamable-http
chematic-mcp --transport streamable-http --port 8080
chematic-mcp --transport streamable-http --port 0
```

The endpoint is `POST /mcp`, bound to `127.0.0.1:3000` by default. Each request
gets one `application/json` response; there are no sessions or server-sent
events.

| Option | Default | Meaning |
|---|---|---|
| `--transport stdio\|streamable-http` | `stdio` | Transport |
| `--bind <ip>` | `127.0.0.1` | Listen address |
| `--port <n>` | `3000` | Listen port; `0` selects a free port |
| `--allow-non-loopback` | off | Required for a non-loopback bind |
| `--allowed-origin <origin>` | none | Additional exact browser origin; repeatable |
| `--max-concurrency <n>` | `16` | Concurrent connection limit |

Requests require `Content-Type: application/json`, an `Accept` value that
admits JSON, and mirrored metadata headers:

- `MCP-Protocol-Version: 2026-07-28`;
- `Mcp-Method: <method>`;
- `Mcp-Name: <tool>` for `tools/call`.

```bash
curl -s http://127.0.0.1:3000/mcp \
  -H 'Content-Type: application/json' \
  -H 'Accept: application/json' \
  -H 'MCP-Protocol-Version: 2026-07-28' \
  -H 'Mcp-Method: tools/call' \
  -H 'Mcp-Name: parse_smiles' \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientCapabilities":{}},"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}'
```

## HTTP errors and limits

| Status | Meaning |
|---|---|
| `200` | JSON-RPC result, including typed chemistry errors with `isError: true` |
| `202` | Notification accepted |
| `400` | Invalid JSON-RPC, metadata, schema, or mirrored headers |
| `403` | Rejected `Host` or `Origin` |
| `404` | Unknown path or unsupported RPC method |
| `405` | Method other than POST |
| `406` / `415` | Unsupported `Accept` / `Content-Type` |
| `408` | Request timeout |
| `411` / `413` | Missing length / body over 1 MiB |
| `431` | Oversized headers |
| `503` | Concurrency limit reached |

Request bodies are limited to 1 MiB. Headers are limited to an 8 KiB request
line, 32 KiB total, and 100 fields. Connections have a 30-second read/write
timeout. The same depth and payload checks apply to stdio.

## Security and privacy

- The default listener is loopback-only. Non-loopback binds require
  `--allow-non-loopback`.
- Loopback mode checks `Host` and permits only loopback or explicitly allowed
  browser origins. No CORS headers are emitted.
- The server provides no TLS, authentication, OAuth, or authorization. Put a
  non-loopback deployment behind controls you operate.
- Nineteen tools keep their inputs in the local process. `name_to_smiles` sends
  the supplied name to the public PubChem HTTPS API with a 10-second timeout.

## Conformance and implementation

The official `@modelcontextprotocol/conformance` 0.2.0-alpha.12 suite was run
against the 2026-07-28 HTTP adapter. Applicable tools-only checks pass; skipped
or failing suite cases require fixture tools or capabilities the server does
not advertise. The machine-readable record is
`validation/results/mcp_2026_07_28_http_conformance.json`.

Both transports reuse the same JSON-RPC codec, server, tool registry, schemas,
and typed errors. Dispatch panics are converted to `-32603 Internal Error`.
The crate forbids unsafe code. Browser applications should use
`chematic-wasm`; `chematic-mcp` requires an OS process for stdio or TCP.
