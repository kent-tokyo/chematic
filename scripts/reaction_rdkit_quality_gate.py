#!/usr/bin/env python3
"""Compare the bounded reaction-SMARTS presence contract with RDKit.

This is intentionally a small semantic gate, not a reaction-yield or product
prediction benchmark. Both engines are evaluated on the same fixed cases and
the report retains every case plus the exact tool/version boundary.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "validation/reaction_rdkit_quality_cases.json"
REPORT = ROOT / "validation/results/reaction-rdkit-quality-v1.0.13.json"
CLI = ROOT / "target/debug/chematic"


def split_section(value: str) -> list[list[str]]:
    return [group.split(".") for group in value.split("|") if group]


def molecule_matches(pattern: str, smiles: str) -> bool:
    query = Chem.MolFromSmarts(pattern)
    molecule = Chem.MolFromSmiles(smiles)
    return query is not None and molecule is not None and molecule.HasSubstructMatch(query)


def groups_match(groups: list[list[str]], molecules: list[str]) -> bool:
    def components_match(components: list[str]) -> bool:
        def assign(component_index: int, assigned: set[int]) -> bool:
            if component_index == len(components):
                return True
            return any(
                index not in assigned and molecule_matches(pattern, molecules[index])
                and assign(component_index + 1, assigned | {index})
                for index in range(len(molecules))
                for pattern in [components[component_index]]
            )

        return assign(0, set())

    # `|` separates alternatives; `.` inside one group requires all
    # components to be assigned injectively.
    return not groups or any(components_match(group) for group in groups)


def rdkit_match(reaction: str, query: str) -> bool:
    reaction_parts = reaction.split(">")
    query_parts = query.split(">")
    if len(reaction_parts) != 3 or len(query_parts) != 3:
        raise ValueError("only three-section reaction strings are admitted")
    reactants, agents, products = ([part.split(".") if part else [] for part in reaction_parts])
    q_reactants, q_agents, q_products = query_parts
    return (
        groups_match(split_section(q_reactants), reactants)
        and groups_match(split_section(q_agents), agents)
        and groups_match(split_section(q_products), products)
    )


def chematic_match(reaction: str, query: str) -> bool:
    result = subprocess.run(
        [str(CLI), "reaction-match", reaction, query],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return bool(json.loads(result.stdout)["matched"])


def main() -> int:
    cases = json.loads(CASES.read_text(encoding="utf-8"))
    rows = []
    for case in cases:
        rdkit = rdkit_match(case["reaction"], case["query"])
        chematic = chematic_match(case["reaction"], case["query"])
        rows.append({**case, "rdkit": rdkit, "chematic": chematic, "agree": rdkit == chematic})
    report = {
        "schema_version": 1,
        "target_version": "1.0.13",
        "scope": "bounded reaction-SMARTS presence matching",
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_binary": str(CLI.relative_to(ROOT)),
        "cases": rows,
        "agreement": {"matched": sum(row["agree"] for row in rows), "total": len(rows)},
        "not_claimed": ["reaction yield", "product prediction", "precision/recall beyond this fixed presence corpus", "Indigo equivalence"],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not all(row["agree"] and row["expected"] == row["rdkit"] for row in rows):
        raise SystemExit("reaction quality gate failed")
    print(f"reaction RDKit quality gate: {len(rows)}/{len(rows)} cases agree")


if __name__ == "__main__":
    main()
