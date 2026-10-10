# chematic-mcp

Local, tools-only MCP server for chematic. It exposes 20 cheminformatics tools
over stdio by default, with an opt-in self-hosted Streamable HTTP transport.

There is no hosted endpoint, authentication, OAuth, or service SLA. Nineteen
tools run locally; `name_to_smiles` is the only tool that calls an external
service (PubChem).

## Install and configure

```bash
cargo install chematic-mcp --version 1.1.0 --locked
```

Ensure Cargo's binary directory is on `PATH` (`$CARGO_HOME/bin`, normally
`~/.cargo/bin`), then configure your MCP client:

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

No repository checkout is required. Running without arguments uses stdio and
auto-detects the protocol dialect from the first request.

Contributors can run the workspace binary directly:

```bash
cargo run -p chematic-mcp --release
```

```toml
[dependencies]
chematic-mcp = { version = "1.1.0", path = "../chematic-mcp" }
```

## Transports and protocol

| Capability | Status |
|---|---|
| Legacy stdio (`2024-11-05`-style initialize) | Supported; byte-compatible |
| 2026-07-28 stateless stdio | Supported |
| 2026-07-28 stateless Streamable HTTP | Supported; opt-in and self-hosted |
| HTTP sessions, GET stream, or SSE | Unsupported |
| Authentication or OAuth | Unsupported |
| Resources, prompts, sampling, roots, logging | Unsupported and not advertised |

The modern dialect returns `structuredContent` validated against each tool's
schema. Chemistry-domain failures are successful RPC results with
`isError: true`; malformed arguments remain JSON-RPC `-32602` errors.

The [full protocol reference](https://kent-tokyo.github.io/chematic/mcp-protocol-reference/)
documents request shapes, headers, status codes, limits, legacy behavior, and
conformance results.

## Streamable HTTP

```bash
chematic-mcp --transport streamable-http
chematic-mcp --transport streamable-http --port 8080
```

The listener binds to `127.0.0.1:3000` by default and serves `POST /mcp` only.
Binding outside loopback requires `--allow-non-loopback`; put such a listener
behind your own TLS, authentication, and network controls.

Important defaults:

- request body: 1 MiB maximum;
- request headers: 32 KiB total and 100 fields;
- connection timeout: 30 seconds;
- concurrency: 16 connections;
- loopback `Host`/`Origin` checks; no CORS response headers.

Example call:

```bash
curl -s http://127.0.0.1:3000/mcp \
  -H 'Content-Type: application/json' \
  -H 'Accept: application/json' \
  -H 'MCP-Protocol-Version: 2026-07-28' \
  -H 'Mcp-Method: tools/call' \
  -H 'Mcp-Name: parse_smiles' \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"_meta":{"io.modelcontextprotocol/protocolVersion":"2026-07-28","io.modelcontextprotocol/clientCapabilities":{}},"name":"parse_smiles","arguments":{"smiles":"c1ccccc1"}}}'
```

## Tools

| Area | Tools |
|---|---|
| Name resolution | `name_to_smiles` |
| Parsing | `parse_smiles`, `canonical_smiles` |
| Properties | `calc_properties`, `lipinski_check`, `sa_score` |
| Filters | `pains_check`, `brenk_check` |
| ADMET | `admet_profile`, `boiled_egg` |
| Search | `ecfp4`, `tanimoto`, `smarts_match`, `find_mcs` |
| 3D | `generate_3d` |
| Reactions | `retrosynthesis` |
| Formats and LLM context | `smiles_to_moljson`, `moljson_to_smiles`, `representation_router`, `molecule_context_pack` |

The registry in `src/tools.rs` is the source of truth.

## Network and privacy

`name_to_smiles` sends the supplied chemical name to the public PubChem REST
API over HTTPS with a 10-second timeout. The other 19 tools do not perform
network I/O. Do not use `name_to_smiles` for a proprietary name unless that
external disclosure is acceptable.

## Design and verification

- `#![forbid(unsafe_code)]` is enforced.
- Both transports use the same codec, `McpServer`, schemas, and tool registry.
- Requests are depth/size bounded; panics are converted to internal errors.
- Every tool publishes JSON Schema 2020-12 input and output schemas.
- The official MCP conformance suite was run for the supported tools-only
  2026-07-28 surface; unsupported capabilities are not advertised.

Evidence: [MCP HTTP record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-09-issues-769-779.md)
and [machine-readable result](https://github.com/kent-tokyo/chematic/blob/main/validation/results/mcp_2026_07_28_http_conformance.json).

## Development

```bash
cargo test -p chematic-mcp --all-targets --locked
cargo clippy -p chematic-mcp --all-targets --locked -- -D warnings
```
