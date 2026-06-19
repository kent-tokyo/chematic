# chematic-mcp

MCP (Model Context Protocol) server for chematic — call cheminformatics tools from AI agents.

## Overview

`chematic-mcp` exposes 15 cheminformatics tools via JSON-RPC 2.0 over stdio,
making them directly callable by Claude and other MCP-compatible AI agents.

```toml
[dependencies]
chematic-mcp = { version = "0.4", path = "../chematic-mcp" }
```

## Running the server

```bash
cargo run -p chematic-mcp --release
```

The server reads newline-delimited JSON-RPC 2.0 requests from stdin and writes
responses to stdout. Add it to your Claude Desktop or Claude Code MCP config:

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

## Available tools (15)

### Name resolution (requires internet)

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

## Example call

```json
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"boiled_egg","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}
```

Response:

```json
{"jsonrpc":"2.0","id":1,"result":{"content":[{"type":"text","text":"{\"gi_absorbed\":true,\"bbb_penetrant\":false,\"logp\":1.316,\"tpsa\":63.6,\"method\":\"BOILED-Egg (Daina & Zoete 2016)\"}"}]}}
```

## Design

- **Mostly local** — 14 of 15 tools are pure computation with no network I/O. `name_to_smiles` is the exception (PubChem REST, requires internet).
- **No unsafe code** — `#![forbid(unsafe_code)]` enforced.
- **WASM-incompatible** — stdio transport requires OS process; use `chematic-wasm` for browser.
