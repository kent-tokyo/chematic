#!/usr/bin/env python3
"""Compare bounded reaction, standardization, and R-group workflows with RDKit.

Accuracy is measured first. Timing is then restricted to rows where both
engines produce the same semantic result, so a faster wrong answer cannot win.
The two engines alternate first position across blocks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
REACTIONS = (
    ROOT / "validation/reaction_product_parity_cases.json",
    ROOT / "validation/reaction_product_parity_strata_v2.json",
)
STANDARDIZATION = (
    ROOT / "validation/standardization_phase1_fixtures.jsonl",
    ROOT / "validation/standardization_phase1_holdout.jsonl",
)
RGROUPS = ROOT / "validation/rdkit/rgroup_labeled_cases.json"
EXPECTED_RDKIT = "2026.03.6"
JOINT_INVALID_REACTION_IDS = {
    "v2_ester_missing_second_reactant",
    "v2_invalid_unclosed_reactant",
    "v2_invalid_smirks_branch",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def canonical(smiles: str) -> str:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse {smiles!r}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def payload_without_dummy(smiles: str) -> str:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse R-group payload {smiles!r}")
    editable = Chem.RWMol(mol)
    for atom in sorted(
        (atom for atom in editable.GetAtoms() if atom.GetAtomicNum() == 0),
        key=lambda atom: atom.GetIdx(),
        reverse=True,
    ):
        editable.RemoveAtom(atom.GetIdx())
    payload = editable.GetMol()
    if payload.GetNumAtoms() == 0:
        return ""
    Chem.SanitizeMol(payload)
    return Chem.MolToSmiles(payload, canonical=True, isomericSmiles=True)


def ratio_summary(rdkit: list[int], chematic: list[int]) -> dict[str, Any]:
    ratios = [left / right for left, right in zip(rdkit, chematic)]
    rng = random.Random(17)
    medians = [
        statistics.median(ratios[rng.randrange(len(ratios))] for _ in ratios)
        for _ in range(10_000)
    ]
    medians.sort()
    return {
        "rdkit_over_chematic_median": statistics.median(ratios),
        "ci95": [medians[249], medians[9749]],
        "winning_blocks": sum(value > 1.0 for value in ratios),
        "blocks": len(ratios),
    }


def timed_blocks(
    chematic_fn: Callable[[], Any], rdkit_fn: Callable[[], Any], blocks: int
) -> dict[str, Any]:
    samples = {"chematic": [], "rdkit": []}
    chematic_fn()
    rdkit_fn()
    functions = {"chematic": chematic_fn, "rdkit": rdkit_fn}
    for block in range(blocks):
        order = ["chematic", "rdkit"] if block % 2 == 0 else ["rdkit", "chematic"]
        for engine in order:
            started = time.perf_counter_ns()
            functions[engine]()
            samples[engine].append(time.perf_counter_ns() - started)
    return {
        "samples_ns": samples,
        "comparison": ratio_summary(samples["rdkit"], samples["chematic"]),
    }


def reaction_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions

    rows = []
    for path in REACTIONS:
        payload = json.loads(path.read_text())
        rows.extend(payload if isinstance(payload, list) else payload["cases"])

    def rdkit_one(row: dict[str, Any]) -> list[str]:
        try:
            reaction = rdChemReactions.ReactionFromSmarts(row["smirks"])
            reactants = tuple(Chem.MolFromSmiles(text) for text in row["reactants"])
        except Exception:
            return []
        if reaction is None or any(mol is None for mol in reactants):
            return []
        products = set()
        try:
            product_sets = reaction.RunReactants(reactants)
        except Exception:
            return []
        for product_set in product_sets:
            try:
                products.add(
                    ".".join(
                        sorted(
                            Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
                            for product in product_set
                        )
                    )
                )
            except Exception:
                continue
        return sorted(products)

    def chematic_one(row: dict[str, Any]) -> list[str]:
        try:
            reactants = [chematic.from_smiles(text) for text in row["reactants"]]
            result = chematic.run_smirks_checked(
                row["smirks"], reactants, rdkit_compat=True
            )
        except Exception:
            return []
        return sorted(
            {
                ".".join(sorted(product.rdkit_smiles for product in product_set))
                for product_set in result["products"]
            }
        )

    observed = {
        "rdkit": {row["id"]: rdkit_one(row) for row in rows},
        "chematic": {row["id"]: chematic_one(row) for row in rows},
    }
    exact = [
        row
        for row in rows
        if observed["rdkit"][row["id"]] == observed["chematic"][row["id"]]
    ]

    timed_rows = [
        row for row in exact if row["id"] not in JOINT_INVALID_REACTION_IDS
    ]

    def chematic_values() -> list[list[str]]:
        return [chematic_one(row) for row in timed_rows]

    def rdkit_values() -> list[list[str]]:
        return [rdkit_one(row) for row in timed_rows]

    end_to_end = timed_blocks(chematic_values, rdkit_values, blocks)
    end_to_end.update(
        rows=len(timed_rows),
        selection="valid rows with identical product sets",
        mode="template lookup/cache plus application; reactant parsing included",
    )
    prepared_rows = timed_rows
    chematic_preparsed = [
        (row, [chematic.from_smiles(text) for text in row["reactants"]])
        for row in prepared_rows
    ]
    rdkit_prepared = [
        (
            row,
            rdChemReactions.ReactionFromSmarts(row["smirks"]),
            tuple(Chem.MolFromSmiles(text) for text in row["reactants"]),
        )
        for row in prepared_rows
    ]

    def chematic_prepared_values() -> list[list[str]]:
        values = []
        for row, reactants in chematic_preparsed:
            result = chematic.run_smirks_checked(
                row["smirks"], reactants, rdkit_compat=True
            )
            values.append(
                sorted(
                    {
                        ".".join(
                            sorted(product.rdkit_smiles for product in product_set)
                        )
                        for product_set in result["products"]
                    }
                )
            )
        return values

    def rdkit_prepared_values() -> list[list[str]]:
        values = []
        for _, reaction, reactants in rdkit_prepared:
            products = set()
            for product_set in reaction.RunReactants(reactants):
                try:
                    products.add(
                        ".".join(
                            sorted(
                                Chem.MolToSmiles(
                                    product, canonical=True, isomericSmiles=True
                                )
                                for product in product_set
                            )
                        )
                    )
                except Exception:
                    continue
            values.append(sorted(products))
        return values

    prepared = timed_blocks(
        chematic_prepared_values, rdkit_prepared_values, blocks
    )
    prepared.update(
        rows=len(prepared_rows),
        selection="valid rows with identical product sets",
        mode="preparsed reactants and precompiled/warm-cached templates",
    )
    return {
        "rows": len(rows),
        "semantic_exact": len(exact) - len(JOINT_INVALID_REACTION_IDS),
        "joint_invalid": len(JOINT_INVALID_REACTION_IDS),
        "mismatches": [
            {
                "id": row["id"],
                "rdkit": observed["rdkit"][row["id"]],
                "chematic": observed["chematic"][row["id"]],
            }
            for row in rows
            if row not in exact
        ],
        "timing": {"end_to_end": end_to_end, "prepared": prepared},
    }


def standardization_rows() -> list[dict[str, str]]:
    rows = []
    for path in STANDARDIZATION:
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if row.get("input") and row.get("expected_kept_fragment"):
                rows.append(row)
    return rows


def standardization_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    rows = standardization_rows()
    expected = [canonical(row["expected_kept_fragment"]) for row in rows]

    def chematic_one(row: dict[str, str]) -> str:
        return canonical(chematic.from_smiles(row["input"]).largest_fragment().smiles)

    def rdkit_one(row: dict[str, str]) -> str:
        mol = Chem.MolFromSmiles(row["input"])
        if mol is None:
            return ""
        return Chem.MolToSmiles(
            rdMolStandardize.LargestFragmentChooser().choose(mol),
            canonical=True,
            isomericSmiles=True,
        )

    observed = {
        "rdkit": [rdkit_one(row) for row in rows],
        "chematic": [chematic_one(row) for row in rows],
    }
    accuracy = {
        engine: {
            "match": sum(value == target for value, target in zip(values, expected)),
            "total": len(rows),
        }
        for engine, values in observed.items()
    }
    exact_indices = [
        i
        for i, target in enumerate(expected)
        if observed["rdkit"][i] == target and observed["chematic"][i] == target
    ]
    exact = [rows[i] for i in exact_indices]

    def chematic_values() -> list[str]:
        return [chematic_one(row) for row in exact]

    def rdkit_values() -> list[str]:
        return [rdkit_one(row) for row in exact]

    timing = timed_blocks(chematic_values, rdkit_values, blocks)
    timing.update(
        rows=len(exact),
        selection="intersection matching the fixture exactly",
        mode="parse plus largest-fragment selection plus canonical write",
    )
    return {"rows": len(rows), "accuracy": accuracy, "timing": timing}


def rgroup_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from rdkit import Chem
    from rdkit.Chem import rdRGroupDecomposition

    cases = json.loads(RGROUPS.read_text())["cases"]

    def chematic_case(case: dict[str, Any]) -> list[dict[str, str]]:
        molecules = [chematic.from_smiles(text) for text in case["molecules"]]
        results = chematic.rgroup_decompose_labeled(case["core"], molecules)
        return [
            {}
            if result is None
            else {
                key: payload_without_dummy(value)
                for key, value in result.items()
                if key.startswith("R")
            }
            for result in results
        ]

    def rdkit_case(case: dict[str, Any]) -> list[dict[str, str]]:
        core = Chem.MolFromSmarts(case["core"])
        molecules = [Chem.MolFromSmiles(text) for text in case["molecules"]]
        if core is None or any(mol is None for mol in molecules):
            return [{} for _ in molecules]
        rows, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [core], molecules, asSmiles=True, asRows=True
        )
        by_input: list[list[str]] = []
        matched = iter(rows)
        unmatched_set = set(unmatched)
        for index in range(len(molecules)):
            if index in unmatched_set:
                by_input.append({})
                continue
            row = next(matched)
            by_input.append(
                {
                    key: payload_without_dummy(value)
                    for key, value in row.items()
                    if key.startswith("R")
                }
            )
        return by_input

    observed = {
        "rdkit": [rdkit_case(case) for case in cases],
        "chematic": [chematic_case(case) for case in cases],
    }
    expected = [case["expected_labeled"] for case in cases]
    accuracy = {
        engine: {
            "match": sum(
                actual_row == expected_row
                for actual_case, expected_case in zip(actual, expected)
                for actual_row, expected_row in zip(actual_case, expected_case)
            ),
            "total": sum(len(case) for case in expected),
        }
        for engine, actual in observed.items()
    }
    exact_cases = [
        case
        for index, case in enumerate(cases)
        if observed["rdkit"][index] == expected[index]
        and observed["chematic"][index] == expected[index]
    ]

    def chematic_values() -> list[list[dict[str, str]]]:
        return [chematic_case(case) for case in exact_cases]

    def rdkit_values() -> list[list[dict[str, str]]]:
        return [rdkit_case(case) for case in exact_cases]

    timing = timed_blocks(chematic_values, rdkit_values, blocks)
    timing.update(
        rows=sum(len(case["molecules"]) for case in exact_cases),
        selection="whole cases where both engines matched every labelled row",
        mode="parse plus decomposition plus payload canonicalization",
    )
    return {"rows": sum(len(case["molecules"]) for case in cases), "accuracy": accuracy, "timing": timing}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=21)
    args = parser.parse_args()
    if args.blocks < 3:
        parser.error("--blocks must be at least 3")

    import chematic
    import rdkit

    if rdkit.__version__ != EXPECTED_RDKIT:
        raise SystemExit(
            f"RDKit mismatch: expected {EXPECTED_RDKIT}, found {rdkit.__version__}"
        )
    report = {
        "schema_version": 1,
        "source_commit": git_head(),
        "engines": {"chematic": chematic.__version__, "rdkit": rdkit.__version__},
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
        },
        "corpora": {
            "reaction": [
                {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
                for path in REACTIONS
            ],
            "standardization": [
                {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
                for path in STANDARDIZATION
            ],
            "rgroup": {"path": str(RGROUPS.relative_to(ROOT)), "sha256": sha256(RGROUPS)},
        },
        "results": {
            "reaction": reaction_gate(args.blocks),
            "standardization": standardization_gate(args.blocks),
            "rgroup": rgroup_gate(args.blocks),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["results"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
