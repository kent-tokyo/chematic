# chematic-mcp

MCP (Model Context Protocol) server for chematic — call cheminformatics tools from AI agents.

## Overview

`chematic-mcp` exposes 20 cheminformatics tools via JSON-RPC 2.0 over stdio,
making them directly callable by Claude and other MCP-compatible AI agents.

**Transport status**: stdio only. The server runs as a local OS process reading
newline-delimited JSON-RPC 2.0 from stdin and writing responses to stdout —
there is no hosted Remote MCP endpoint, no authentication, and no public
service SLA. A remote-ready refactor (transport-neutral protocol handling, so
a Streamable HTTP adapter could be added without rewriting tool logic) is
under consideration but not implemented; nothing here is reachable over the
network except the one tool noted below.

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

```json
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"boiled_egg","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}
```

Response:

```json
{"jsonrpc":"2.0","id":1,"result":{"content":[{"type":"text","text":"{\"gi_absorbed\":true,\"bbb_penetrant\":false,\"logp\":1.316,\"tpsa\":63.6,\"method\":\"BOILED-Egg (Daina & Zoete 2016)\"}"}]}}
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
