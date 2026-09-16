# AI-assisted molecular analysis with MCP

`chematic-mcp` is a local stdio server that lets an MCP-compatible agent call
named cheminformatics operations on SMILES. The agent supplies language and
workflow logic; chematic returns structured chemical results. It is not a
remote service, a general reaction planner, or a substitute for experimental
or clinical judgement.

## Start locally

Build the server from this checkout:

```bash
cargo build --release -p chematic-mcp
```

Then point the MCP client at the resulting local binary:

```json
{
  "mcpServers": {
    "chematic": {
      "command": "/absolute/path/to/target/release/chematic-mcp"
    }
  }
}
```

The server uses newline-delimited JSON-RPC over stdio. It has no HTTP
endpoint, authentication layer, or hosted SLA. See the
[`chematic-mcp` README](https://github.com/kent-tokyo/chematic/blob/main/crates/chematic-mcp/README.md) for the legacy
and modern protocol envelopes.

## What an agent can call

The current server exposes 20 tools. The useful groups are:

| Need | Tools |
| --- | --- |
| Parse and identify a molecule | `parse_smiles`, `canonical_smiles`, `smiles_to_moljson`, `moljson_to_smiles` |
| Calculate screening properties | `calc_properties`, `lipinski_check`, `sa_score`, `pains_check`, `brenk_check`, `admet_profile`, `boiled_egg` |
| Compare structures | `ecfp4`, pairwise `tanimoto`, `smarts_match`, `find_mcs` |
| Produce bounded follow-up representations | `representation_router`, `molecule_context_pack`, `generate_3d`, `retrosynthesis` |
| Resolve a public chemical name | `name_to_smiles` |

There is no built-in indexed `similarity_search` or `standardize` MCP tool.
For a caller-owned collection, the agent or application must iterate pairwise
`tanimoto` calls or use a binding-level index, while retaining input IDs and
all rejected records. `retrosynthesis` is a bounded one-step BRICS
disconnection, not a synthetic route recommendation. `generate_3d` and the
ADMET-related tools have their documented bounded or experimental scope.

## Example: inspect one input

For a SMILES supplied by the user, an agent can call the property and alert
tools independently, then cite the returned fields instead of inventing a
chemical interpretation. A modern `tools/call` request looks like this:

```json
{
  "jsonrpc": "2.0",
  "id": 1,
  "method": "tools/call",
  "params": {
    "_meta": {
      "io.modelcontextprotocol/protocolVersion": "2026-07-28",
      "io.modelcontextprotocol/clientInfo": {"name": "my-agent", "version": "1.0"},
      "io.modelcontextprotocol/clientCapabilities": {}
    },
    "name": "calc_properties",
    "arguments": {"smiles": "CC(=O)Nc1ccc(O)cc1"}
  }
}
```

Successful modern calls include machine-readable `structuredContent`.
Invalid chemistry input instead returns a successful RPC with `isError: true`
and a structured error code such as `INVALID_SMILES`; invalid argument shape
is JSON-RPC `-32602`. Treat either outcome as data in an agent workflow rather
than silently retrying with altered structures.

## Local-data boundary

Nineteen tools operate locally. `name_to_smiles` is the exception: it sends
the supplied chemical name to PubChem over HTTPS. Do not use that tool for
proprietary names. The stdio server currently has bounded request and tool
inputs, but it does not implement JSON-RPC cancellation; avoid presenting
long-running work as interruptible. Exact limits and failure behaviour are in
[Error and resource limits](../error-and-limits.md).

For a Python-created comparison prompt, see
[`examples/ai_agent_demo.py`](https://github.com/kent-tokyo/chematic/blob/main/examples/ai_agent_demo.py). For
version-pinned chemistry compatibility and validation scope, start with
[Validation](../validation.md) and the [RDKit migration guide](../rdkit-migration.md).
