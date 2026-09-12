#!/usr/bin/env python3
"""Verify fail-closed reaction behavior for incompatible stereo inputs.

This is intentionally separate from the RDKit product-set parity gate.  The
pinned RDKit reaction runner currently emits a product for the opposite
enantiomer because it does not enforce reactant-template chirality during
``RunReactants``.  Schematic's safer contract rejects that input.  The
behavior is recorded, not counted as RDKit-equivalent parity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases",
        type=Path,
        default=ROOT / "validation/reaction_product_stereo_safety_cases.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import chematic

    cases = json.loads(args.cases.read_text(encoding="utf-8"))
    rows = []
    for case in cases:
        raw = chematic.run_smirks(
            case["smirks"], [chematic.from_smiles(value) for value in case["reactants"]]
        )
        schematic_product_sets = [
            sorted(Chem.MolToSmiles(Chem.MolFromSmiles(product.smiles), canonical=True)
                   for product in product_set)
            for product_set in raw
        ]
        schematic_product_sets = sorted(set(tuple(values) for values in schematic_product_sets))
        schematic_product_sets = [list(values) for values in schematic_product_sets]

        reaction = AllChem.ReactionFromSmarts(case["smirks"])
        reactants = tuple(Chem.MolFromSmiles(value) for value in case["reactants"])
        rdkit_product_count = len(reaction.RunReactants(reactants))
        expected = case["expected_chematic_product_sets"]
        safe_refusal = schematic_product_sets == expected
        rows.append({
            **case,
            "chematic_product_sets": schematic_product_sets,
            "rdkit_run_reactants_product_set_count": rdkit_product_count,
            "safe_refusal": safe_refusal,
        })

    report = {
        "schema_version": 1,
        "profile": "reaction_product_stereo_fail_closed",
        "rdkit_version": rdBase.rdkitVersion,
        "case_file_sha256": hashlib.sha256(args.cases.read_bytes()).hexdigest(),
        "cases": rows,
        "agreement": {
            "safe_refusals": sum(row["safe_refusal"] for row in rows),
            "total": len(rows),
        },
        "not_claimed": ["RDKit reaction-runner stereo parity", "reaction yield", "selectivity"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "reaction stereo safety: "
        f"{report['agreement']['safe_refusals']}/{report['agreement']['total']} safe refusals"
    )
    return 0 if all(row["safe_refusal"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
