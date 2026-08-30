# chematic-mcp

**MCP 2026-07-28 tools-only stateless stdio server** for chematic — call
cheminformatics tools from AI agents. Also speaks the legacy
(`2024-11-05`-style) stdio dialect on the same connection, byte-compatible
with earlier `chematic-mcp` releases.

## Overview

`chematic-mcp` exposes 20 cheminformatics tools via JSON-RPC 2.0 over stdio,
making them directly callable by Claude and other MCP-compatible AI agents.

**Transport status**: stdio only. The server runs as a local OS process reading
newline-delimited JSON-RPC 2.0 from stdin and writing responses to stdout.
The transport/protocol-codec/server-core/tool-registry layering exists so a
Streamable HTTP adapter *could* be added later without rewriting tool logic,
but no HTTP code exists in this crate today — there is no hosted Remote MCP
endpoint, no authentication, and no public service SLA. Nothing here is
reachable over the network except the one tool noted below.

| Capability | Status |
|---|---|
| Legacy stdio (`2024-11-05`-style `initialize` handshake) | **Supported**, byte-compatible |
| 2026-07-28 stateless stdio (`server/discover`, per-request `_meta`) | **Supported** |
| Remote HTTP (Streamable HTTP) | **Unsupported** |
| Authentication / OAuth | **Unsupported** |
| Tasks extension | **Unsupported** |
| MCP Apps | **Unsupported** |
| Resources / Prompts / Sampling / Roots / Logging / Subscriptions | **Unsupported** — not advertised in `server/discover`'s capabilities |

See this README for the protocol behavior and compatibility notes; the tool
registry in `src/tools.rs` is the source of truth for the available tools.

```toml
[dependencies]
chematic-mcp = { version = "0.4", path = "../chematic-mcp" }
```

## Running the server

```bash
cargo run -p chematic-mcp --release
```

The server reads newline-delimited JSON-RPC 2.0 requests from stdin and writes
responses to stdout, and auto-detects which protocol era the client speaks
from the *first* request on the connection (see "Protocol eras" below). Add
it to your Claude Desktop or Claude Code MCP config:

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

## Available tools (20)

### Name resolution (requires internet)

This is the **only** tool of the 20 that makes a network call. The other 19
are pure local computation — see [Network & privacy](#network--privacy) below.

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
- **WASM-incompatible** — stdio transport requires OS process; use `chematic-wasm` for browser.
- **Layered core**: `transport` (stdio framing + connection-pinned protocol
  era) → `protocol` (JSON-RPC codec, `_meta` parsing, error vocabulary,
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
