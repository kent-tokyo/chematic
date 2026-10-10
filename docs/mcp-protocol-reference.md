# MCP protocol and transport reference

This is the detailed reference for the `chematic-mcp` protocol, transports,
limits, security model, tool schemas, and conformance evidence. For installation
and a minimal client configuration, start with the
[`chematic-mcp` README](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mcp/README.md).

**MCP 2026-07-28 tools-only stateless server** for chematic — call
cheminformatics tools from AI agents. stdio by default; an opt-in,
self-hosted Streamable HTTP transport is available with
`--transport streamable-http` (see "Streamable HTTP" below). Over stdio it
also speaks the legacy (`2024-11-05`-style) dialect on the same connection,
byte-compatible with earlier `chematic-mcp` releases.

## Overview

`chematic-mcp` exposes 20 cheminformatics tools via JSON-RPC 2.0 over stdio (or, opt-in, Streamable HTTP),
making them directly callable by Claude and other MCP-compatible AI agents.

**Transport status**: stdio is the default, and running `chematic-mcp`
with no arguments behaves exactly as before. The server then runs as a local
OS process reading newline-delimited JSON-RPC 2.0 from stdin and writing
responses to stdout. `chematic-mcp --transport streamable-http` instead
serves the same protocol, server and tool registry over HTTP at `/mcp`,
bound to `127.0.0.1:3000` by default. That adapter is something you run
yourself: **there is no hosted Remote MCP endpoint, no authentication or
OAuth, and no service SLA.** Apart from the HTTP listener you start
yourself, nothing here is reachable over the network except the one tool
noted under "Network & privacy".

| Capability | Status |
|---|---|
| Legacy stdio (`2024-11-05`-style `initialize` handshake) | **Supported**, byte-compatible |
| 2026-07-28 stateless stdio (`server/discover`, per-request `_meta`) | **Supported** |
| Streamable HTTP, 2026-07-28 stateless dialect (`--transport streamable-http`, self-hosted, loopback by default) | **Supported**, opt-in — no hosted endpoint |
| Legacy dialect or sessions over HTTP (`initialize`, `Mcp-Session-Id`, GET stream) | **Unsupported** — HTTP serves 2026-07-28 only |
| Authentication / OAuth | **Unsupported** |
| Tasks extension | **Unsupported** |
| MCP Apps | **Unsupported** |
| Resources / Prompts / Sampling / Roots / Logging / Subscriptions | **Unsupported** — not advertised in `server/discover`'s capabilities |

See this README for the protocol behavior and compatibility notes; the tool
registry in `src/tools.rs` is the source of truth for the available tools.

## Quick start

Install the published binary from crates.io. Pinning the version keeps the
server and the documented protocol surface in sync:

```bash
cargo install chematic-mcp --version 1.1.0 --locked
```

Make sure Cargo's binary directory is on `PATH` (`$CARGO_HOME/bin`, normally
`~/.cargo/bin`), then configure an MCP client to run the installed binary over
the default stdio transport:

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

No repository checkout is required for this path. Streamable HTTP remains an
opt-in, self-hosted transport; this package does not provide a hosted endpoint,
authentication, OAuth, or a service SLA.

## Running from a source checkout

Contributors working in this repository can run the workspace binary directly:

```bash
cargo run -p chematic-mcp --release
```

The server reads newline-delimited JSON-RPC 2.0 requests from stdin and writes
responses to stdout, and auto-detects which protocol era the client speaks
from the *first* request on the connection (see "Protocol eras" below). A
checkout-based client configuration is:

```json
{
  "mcpServers": {
    "chematic": {
      "command": "cargo",
      "args": ["run", "-p", "chematic-mcp", "--release", "--quiet"]
    }
  }
}
```

Workspace users that link the protocol/server crate directly can keep the
published version and local path explicit:

```toml
[dependencies]
chematic-mcp = { version = "1.1.0", path = "../chematic-mcp" }
```

## Protocol eras

A single stdio connection pins to whichever dialect its first request
speaks and stays there for the rest of the connection (a request that tries
to switch dialects mid-connection gets a typed protocol error, not a
silent dialect change).

### Legacy (`2024-11-05`-style)

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{}}}
{"jsonrpc":"2.0","method":"notifications/initialized","params":{}}
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}
```

Byte-compatible with pre-2026-07-28 releases: `tools/call` results are
`{"content":[{"type":"text","text":"<json>"}]}`, and tool failures (both
argument-shape and chemistry-domain) surface as a JSON-RPC error with the
implementation-defined code `-32000`.

### Modern (2026-07-28 stateless)

No `initialize` handshake — every request carries its protocol version,
client identity, and capabilities inline in a reserved `_meta` triple:

```json
{
  "jsonrpc": "2.0", "id": 1, "method": "server/discover",
  "params": {
    "_meta": {
      "io.modelcontextprotocol/protocolVersion": "2026-07-28",
      "io.modelcontextprotocol/clientInfo": {"name": "my-client", "version": "1.0.0"},
      "io.modelcontextprotocol/clientCapabilities": {}
    }
  }
}
```

```json
{
  "jsonrpc": "2.0", "id": 2, "method": "tools/call",
  "params": {
    "_meta": {
      "io.modelcontextprotocol/protocolVersion": "2026-07-28",
      "io.modelcontextprotocol/clientInfo": {"name": "my-client", "version": "1.0.0"},
      "io.modelcontextprotocol/clientCapabilities": {}
    },
    "name": "parse_smiles",
    "arguments": {"smiles": "c1ccccc1"}
  }
}
```

`tools/call` results carry `resultType`, `content`, and `structuredContent`
(validated against each tool's `outputSchema`). A chemistry-domain failure
(e.g. an invalid SMILES string) is a **successful** RPC with
`isError: true` and a machine-readable `structuredContent.error.code`
(e.g. `"INVALID_SMILES"`) — never a JSON-RPC transport error, so an LLM
client can see what went wrong and retry. An argument-shape/schema
violation (missing/wrong-typed argument, unknown tool) *is* a JSON-RPC
error, `-32602 Invalid Params`.

`io.modelcontextprotocol/protocolVersion` and
`io.modelcontextprotocol/clientCapabilities` are required on every request.
`io.modelcontextprotocol/clientInfo` is optional: clients SHOULD send it, but
a request without it is served. A `clientInfo` that is present but malformed
is still rejected with `-32602`.

### Tool registry caching

The modern `tools/list` result carries a cache hint:

```json
{ "ttlMs": 86400000, "cacheScope": "public" }
```

The 20-tool registry is static and identical for every caller within one
running process, so `ttlMs` is 24 hours and `cacheScope` is `"public"`. If
the registry ever changes, it changes as part of a new `chematic-mcp`
release (a new process) — a client restarting the server always gets a
fresh response regardless of a previously-cached TTL window; there is no
runtime path that mutates the registry mid-process.

## Streamable HTTP (opt-in)

```bash
chematic-mcp --transport streamable-http                 # http://127.0.0.1:3000/mcp
chematic-mcp --transport streamable-http --port 8080
chematic-mcp --transport streamable-http --port 0        # any free port, printed to stderr
```

The HTTP transport serves the **2026-07-28 stateless dialect only**. It reuses
the stdio stack: every HTTP request becomes one JSON-RPC message, which the
same codec, `McpServer` and tool registry answer. Tool schemas,
`structuredContent` and the typed errors are therefore exactly the ones stdio
returns, and the integration tests compare the two byte for byte. There are
no sessions, no GET stream and no SSE. Each POST gets one `application/json`
answer, and the connection closes after it.

| Option | Default | Meaning |
|---|---|---|
| `--transport stdio\|streamable-http` | `stdio` | Which transport to serve |
| `--bind <ip>` | `127.0.0.1` | Address to listen on |
| `--port <n>` | `3000` | Port (`0` picks a free one) |
| `--allow-non-loopback` | off | Required to bind any non-loopback address |
| `--allowed-origin <origin>` | none | Extra browser `Origin` allowed (repeatable, exact `scheme://host[:port]`) |
| `--max-concurrency <n>` | `16` | Connections served at once; more get `503` with `Retry-After: 1` |

HTTP-only options given together with `--transport stdio` are rejected, as
are unknown options (exit status 2).

**Requests.** Send `POST /mcp` with `Content-Type: application/json` and an
`Accept` that admits `application/json`. Each request also carries the
mirrored metadata headers the transport requires:

- `MCP-Protocol-Version: 2026-07-28`
- `Mcp-Method: <method>`
- for `tools/call` only, `Mcp-Name: <tool name>` (the `=?base64?…?=` form is
  decoded)

```bash
curl -s http://127.0.0.1:3000/mcp \
  -H 'Content-Type: application/json' -H 'Accept: application/json, text/event-stream' \
  -H 'MCP-Protocol-Version: 2026-07-28' -H 'Mcp-Method: tools/call' -H 'Mcp-Name: parse_smiles' \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientCapabilities":{}},"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}'
```

**Status codes.**

| Status | When |
|---|---|
| `200` | A result. This includes a tool's typed chemistry error, which is a result with `isError: true`, as on stdio. |
| `202` | A notification. The body is empty. |
| `400` | One of the following: malformed JSON (`-32700`); not one JSON-RPC object, including batches (`-32600`); `_meta` without `protocolVersion`, or arguments that fail the tool's input schema (`-32602`); a mirrored header that is missing or does not match the body (`-32020`); an unsupported protocol version (`-32022`, with the supported list). |
| `403` | A `Host` that is not loopback on a loopback bind, or an `Origin` that is not allowed. |
| `404` | Any path other than `/mcp`, or an unknown method (`-32601`). This includes `initialize`, `ping`, `resources/*` and `prompts/*`. |
| `405` | Any method other than POST. The response carries `Allow: POST`. |
| `406` | An `Accept` header that does not admit `application/json`. |
| `408` | The client was too slow to send the request. |
| `411` | No `Content-Length`. Chunked bodies are not read. |
| `413` | A body over 1 MiB. This is the same limit as stdio's `MAX_REQUEST_BYTES`. |
| `415` | A `Content-Type` other than `application/json`. |
| `431` | Request headers that are too large. |
| `503` | The concurrency limit is reached. |

**Limits.**

- Bodies are bounded at 1 MiB, and the same depth and size checks as stdio
  apply.
- Request headers are bounded at 8 KiB for the request line, 32 KiB in
  total, and 100 fields.
- Every connection has a 30 s read/write timeout.
- At most `--max-concurrency` connections are served at once.

**Security.**

- **Authentication.** There is none. Anyone who can reach the port can call
  every tool, including `name_to_smiles`, which makes outbound requests to
  PubChem.
- **Loopback by default.** The server binds `127.0.0.1` unless told
  otherwise. On a loopback bind it requires a loopback `Host` header, and
  accepts an `Origin` only if it is loopback
  (`http(s)://localhost|127.0.0.1|[::1][:port]`) or listed with
  `--allowed-origin`. Everything else gets `403`, which blocks DNS rebinding.
  No CORS headers are ever sent, and `OPTIONS` gets `405`, so a browser page
  on another origin cannot read responses.
- **Binding outside loopback.** Binding to a non-loopback address (for
  example `--bind 0.0.0.0`) is refused unless you also pass
  `--allow-non-loopback`. If you do, put the server behind your own TLS and
  authentication (a reverse proxy, a VPN, a firewall rule). This crate
  provides neither, and there is no OAuth flow.
- **No hosted service.** This is not a hosted chematic service and carries
  no SLA.

**Conformance.** The official MCP conformance suite
(`@modelcontextprotocol/conformance` 0.2.0-alpha.12) was run against this
adapter for every 2026-07-28 server scenario; see
`validation/results/mcp_2026_07_28_http_conformance.json` for the results.
Every check that applies to a tools-only server passes. The checks that do
not pass fall into two groups:

- checks that need the suite's own fixture tools (`test_*`, a
  `json_schema_2020_12_tool`, tools with `x-mcp-header`);
- checks that exercise capabilities this server does not implement or
  advertise (resources, prompts, completion, input-required results).

## Available tools (20)

### Name resolution (requires internet)

This is the **only** tool of the 20 that makes a network call. The other 19
are pure local computation — see [Network & privacy](#network-privacy) below.

| Tool | Description |
|------|-------------|
| `name_to_smiles` | Convert a chemical name (IUPAC, common, or trade name) to SMILES via PubChem |

### Parsing & basic info

| Tool | Description |
|------|-------------|
| `parse_smiles` | Parse SMILES → atom count, bond count, MW |
| `canonical_smiles` | Canonicalize a SMILES string |

### Molecular properties

| Tool | Description |
|------|-------------|
| `calc_properties` | MW, exact mass, LogP, TPSA, HBD, HBA, rotatable bonds, QED |
| `lipinski_check` | Lipinski Rule-of-Five with per-rule breakdown |
| `sa_score` | Synthetic accessibility score (1 = easy, 10 = hard; < 6 = synthesizable) |

### Drug-likeness & safety filters

| Tool | Description |
|------|-------------|
| `pains_check` | PAINS structural alerts (HTS false-positive filter) |
| `brenk_check` | Brenk toxicity / instability alerts |

### ADMET / pharmacokinetics

| Tool | Description |
|------|-------------|
| `admet_profile` | Full ADMET: BBB, Caco-2, hERG, CYP3A4, AMES, PPB, hepatic clearance |
| `boiled_egg` | BOILED-Egg method (Daina & Zoete 2016) — GI absorption + BBB zone prediction |

### Similarity & substructure

| Tool | Description |
|------|-------------|
| `ecfp4` | ECFP4 fingerprint as 2048-bit hex + popcount |
| `tanimoto` | Tanimoto similarity (ECFP4) between two molecules |
| `smarts_match` | SMARTS substructure search — match count + atom maps |
| `find_mcs` | Maximum common substructure across a list of molecules |

### 3D

| Tool | Description |
|------|-------------|
| `generate_3d` | 3D coordinates via rule-based placement + DREIDING minimization (XYZ) |

### Retrosynthesis

| Tool | Description |
|------|-------------|
| `retrosynthesis` | One-step BRICS disconnection, all breakable bonds cut individually, ranked by max fragment SA Score |

### Format conversion & LLM integration

| Tool | Description |
|------|-------------|
| `smiles_to_moljson` | SMILES → MolJSON (explicit atom/bond JSON representation for LLM consumption) |
| `moljson_to_smiles` | MolJSON → canonical SMILES |
| `representation_router` | Route SMILES to the molecular text representation best suited to a given LLM task (MolJSON/CML/InChI/canonical SMILES) |
| `molecule_context_pack` | Assemble identifiers, properties, drug-likeness, ADMET, and MolJSON into a single LLM/RAG context object (json/markdown/prompt output) |

## Example call

Legacy era:

```json
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"boiled_egg","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}
```

Response:

```json
{"jsonrpc":"2.0","id":1,"result":{"content":[{"type":"text","text":"{\"gi_absorbed\":true,\"bbb_penetrant\":false,\"logp\":1.316,\"tpsa\":63.6,\"method\":\"BOILED-Egg (Daina & Zoete 2016)\"}"}]}}
```

Modern era (same tool, with the `_meta` triple — see "Protocol eras" above):

```json
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientInfo":{"name":"my-client","version":"1.0.0"},"io.modelcontextprotocol/clientCapabilities":{}},"name":"boiled_egg","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}
```

Response — same `content`, plus `resultType`/`structuredContent`:

```json
{"jsonrpc":"2.0","id":1,"result":{"resultType":"complete","content":[{"type":"text","text":"{\"gi_absorbed\":true,\"bbb_penetrant\":false,\"logp\":1.316,\"tpsa\":63.6,\"method\":\"BOILED-Egg (Daina & Zoete 2016)\"}"}],"structuredContent":{"gi_absorbed":true,"bbb_penetrant":false,"logp":1.316,"tpsa":63.6,"method":"BOILED-Egg (Daina & Zoete 2016)"}}}
```

## Network & privacy

- **19 of 20 tools are pure local computation** — no network I/O, no external
  service dependency, nothing leaves the process.
- **`name_to_smiles` is the one exception.** The chemical name string you pass
  it is sent, URL-encoded, to the public PubChem REST API
  (`pubchem.ncbi.nlm.nih.gov`) over HTTPS with a 10-second timeout. If PubChem
  is unreachable, slow, or returns an unexpected response, the tool call
  fails — it does not fall back to local computation, since there is none for
  name resolution.
- If you're working with proprietary or privacy-sensitive compound names,
  be aware `name_to_smiles` is the only tool where your input crosses the
  network; the other 19 never do.

## Design

- **No unsafe code** — `#![forbid(unsafe_code)]` enforced.
- **WASM-incompatible** — both transports need an OS process (stdio, or a
  TCP listener for the opt-in HTTP adapter); use `chematic-wasm` for browser.
- **No new dependencies for HTTP** — the Streamable HTTP adapter
  (`src/http.rs`) uses only `std::net`; it is a thin framing layer in front
  of the same `Connection`/`McpServer` the stdio transport uses.
- **Layered core**: `transport` (stdio framing + connection-pinned protocol
  era; `http` adds HTTP framing in front of a fresh connection per request) → `protocol` (JSON-RPC codec, `_meta` parsing, error vocabulary,
  adversarial-input limits) → `server` (method dispatch + per-era response
  shaping) → `tools` (chemistry, protocol-agnostic). Every tool computes its
  result exactly once; the presentation layer decides how to wrap it per
  era — see the protocol sections above.
- **JSON Schema 2020-12**: every tool has both `inputSchema` and
  `outputSchema`; no external `$ref` is ever resolved. Runtime argument
  validation uses a small hand-rolled subset validator (no new runtime
  dependency); the real `jsonschema` crate is a `dev-dependency` used only
  in tests to independently confirm every schema is itself valid 2020-12.
- **Panics never reach the wire** — every tool dispatch is wrapped in
  `catch_unwind` and mapped to `-32603 Internal Error`.
