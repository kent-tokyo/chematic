# Use case: AI-assisted molecular analysis with MCP

## Problem

You want an AI agent (Claude, GPT-4, etc.) to reason about molecular structures — but LLMs have no native chemistry tools. You need a bridge that lets the AI call real cheminformatics functions on actual SMILES strings.

## Solution

chematic ships a built-in MCP server (`chematic-mcp`) that exposes chemistry tools to any MCP-compatible AI agent. Wire it to Claude Desktop once; the AI can then evaluate, search, and synthesise molecules in natural conversation — no Python environment on the client side.

## Setup

```bash
# Install chematic Python package (no C/C++ toolchain needed)
pip install chematic

# Or build the MCP server binary from source
cargo build -p chematic-mcp --release
```

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "chematic": {
      "command": "/path/to/chematic-mcp"
    }
  }
}
```

## What the agent can do

Once connected, Claude Desktop can call these tools:

| Tool | What it does |
|------|-------------|
| `evaluate_molecule` | MW, LogP, TPSA, ADMET, PAINS, Lipinski |
| `similarity_search` | Find similar structures in a SMILES list |
| `substructure_search` | SMARTS-based filter across a set |
| `generate_3d` | ETKDG 3D coordinates + MMFF94 minimization |
| `retro_disconnect` | Retrosynthesis (60 reaction templates) |
| `standardize` | Salt stripping, neutralization, tautomer canon. |

## Example conversation

> **User:** I have a hit compound `CC(=O)Nc1ccc(O)cc1`. Is it drug-like?

```
chematic evaluates the molecule and returns:

Molecular weight 151.2 Da, formula C8H9NO2.
LogP 0.46 (mildly lipophilic), TPSA 49.3 Å².
HBD 2, HBA 3, 2 rotatable bonds, 1 aromatic ring.
Drug-likeness: no Lipinski rule-of-5 violations. Likely orally bioavailable.
QED 0.67. No structural alerts (PAINS / Brenk clean).
```

> **User:** Find me 5 analogs from this SMILES list that are most similar.

The agent calls `similarity_search` and returns ranked results with Tanimoto scores.

## Scripting with Python

```python
import chematic

mol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol

# Natural-language summary (same as what the MCP tool returns)
print(mol.describe())

# Full ADMET profile
print(mol.admet())

# Retrosynthesis
routes = mol.retro_disconnect(max_results=5)
for r in routes:
    print(r)
```

## Typical agentic workflow

```
User query
    ↓
Claude reads SMILES / name
    ↓
chematic MCP: standardize → evaluate → similarity_search
    ↓
Claude interprets results and suggests next experiment
    ↓
chematic MCP: retro_disconnect (if synthesis needed)
    ↓
Claude writes experimental plan
```

chematic is the chemistry engine; the AI handles language, reasoning, and decision-making.
No Python environment is needed on the client — the MCP server binary handles everything.

## Related APIs

- [`mol.describe()`](../api/chematic.md) — natural-language property summary (the same text the MCP server returns)
- [`mol.admet()`](../api/chematic.md) — BBB, Caco-2, hERG, CYP3A4 in one call
- [`mol.retro_disconnect()`](../api/chematic.md) — retrosynthesis with 60 reaction templates
- `chematic-mcp` binary — run `cargo build -p chematic-mcp --release`
