"""
LLM / AI agent context builder
================================
Shows how mol.describe() and mol.diff() produce structured, natural-language
context that you can drop directly into an LLM prompt or MCP tool response.

Demonstrates:
  - mol.describe()          → one-paragraph property summary
  - mol.diff(other)         → element-level + descriptor delta between two molecules
  - Building a comparison prompt for Claude / GPT-4

Run:
    python examples/ai_agent_demo.py
    python examples/ai_agent_demo.py "CC(=O)Oc1ccccc1C(=O)O" "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"

Dependencies:
    pip install chematic
"""
import sys
import chematic

ASPIRIN    = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN  = "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"


def build_comparison_prompt(smiles1: str, smiles2: str, name1: str = "Compound A", name2: str = "Compound B") -> str:
    mol1 = chematic.from_smiles(smiles1)
    mol2 = chematic.from_smiles(smiles2)
    if mol1 is None or mol2 is None:
        return "Error: invalid SMILES"

    d = mol1.diff(mol2)

    return f"""You are a medicinal chemistry assistant.

### {name1}
{mol1.describe()}

### {name2}
{mol2.describe()}

### Structural difference ({name1} → {name2})
{d['summary']}
Common scaffold: {d['common_atoms']} heavy atoms (MCS).
Element changes: {d['delta_elements']}

### Question
Compare these two compounds from a drug-likeness and ADMET perspective. Which would you prioritise for oral administration and why?"""


def main() -> None:
    if len(sys.argv) == 3:
        smi1, smi2 = sys.argv[1], sys.argv[2]
        name1, name2 = "Compound A", "Compound B"
    else:
        smi1, smi2 = ASPIRIN, IBUPROFEN
        name1, name2 = "Aspirin", "Ibuprofen"

    print("=" * 60)
    print("SINGLE-MOLECULE DESCRIBE OUTPUT")
    print("=" * 60)
    mol = chematic.from_smiles(smi1)
    print(f"SMILES: {smi1}")
    print(mol.describe())

    print()
    print("=" * 60)
    print("STRUCTURAL DIFF")
    print("=" * 60)
    d = chematic.from_smiles(smi1).diff(chematic.from_smiles(smi2))
    print(f"{name1} → {name2}:")
    print(f"  Summary:    {d['summary']}")
    print(f"  ΔMW:        {d['delta_mw']:+.1f} Da")
    print(f"  ΔLogP:      {d['delta_logp']:+.2f}")
    print(f"  ΔTPSA:      {d['delta_tpsa']:+.1f} Å²")
    print(f"  ΔHBD:       {d['delta_hbd']:+d}")
    print(f"  Common MCS: {d['common_atoms']} atoms")

    print()
    print("=" * 60)
    print("EXAMPLE LLM PROMPT (paste into Claude / GPT-4)")
    print("=" * 60)
    print(build_comparison_prompt(smi1, smi2, name1, name2))


if __name__ == "__main__":
    main()
