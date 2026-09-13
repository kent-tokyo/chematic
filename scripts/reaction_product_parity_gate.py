#!/usr/bin/env python3
"""Compare bounded SMIRKS product sets with RDKit on fixed local cases.

This gate compares mapped canonical product structures only.  It does not
claim reaction yield, selectivity, or broad reaction-model precision/recall.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]


def canonical_product_set(products) -> list[list[str]]:
    # RunReactants returns one product tuple per embedding.  A product *set*
    # must not count repeated embeddings as distinct outcomes; chematic's
    # transform layer already removes those duplicates.  Deduplicate both
    # molecules within a tuple and identical tuples across embeddings.
    sets = set()
    for product_set in products:
        values = []
        for mol in product_set:
            if mol is None:
                continue
            values.append(Chem.MolToSmiles(mol, canonical=True))
        sets.add(tuple(sorted(values)))
    return [list(values) for values in sorted(sets)]


def rdkit_products(smirks: str, reactants: list[str]) -> list[list[str]]:
    reaction = AllChem.ReactionFromSmarts(smirks)
    if reaction is None:
        raise ValueError("RDKit could not parse SMIRKS")
    molecules = [Chem.MolFromSmiles(value) for value in reactants]
    if any(mol is None for mol in molecules):
        raise ValueError("RDKit could not parse a reactant")
    return canonical_product_set(reaction.RunReactants(tuple(molecules)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import chematic

    cases = json.loads(args.cases.read_text(encoding="utf-8"))
    rows = []
    for case in cases:
        rdkit = rdkit_products(case["smirks"], case["reactants"])
        raw = chematic.run_smirks(
            case["smirks"], [chematic.from_smiles(value) for value in case["reactants"]]
        )
        chematic_value = set()
        for product_set in raw:
            molecules = [Chem.MolFromSmiles(product.smiles) for product in product_set]
            if any(mol is None for mol in molecules):
                raise SystemExit(f"{case['id']}: generated product is not RDKit-readable")
            chematic_value.add(tuple(sorted(Chem.MolToSmiles(mol, canonical=True) for mol in molecules)))
        chematic_value = [list(values) for values in sorted(chematic_value)]
        rows.append({**case, "rdkit": rdkit, "chematic": chematic_value, "agree": rdkit == chematic_value})

    report = {
        "schema_version": 1,
        "profile": "bounded_reaction_product_structure_parity",
        "rdkit_version": rdBase.rdkitVersion,
        "case_file_sha256": hashlib.sha256(args.cases.read_bytes()).hexdigest(),
        "cases": rows,
        "agreement": {"matched": sum(row["agree"] for row in rows), "total": len(rows)},
        "not_claimed": ["reaction yield", "selectivity", "broad reaction precision/recall", "Indigo equivalence"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"reaction product parity: {report['agreement']['matched']}/{report['agreement']['total']} cases")
    return 0 if all(row["agree"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
