#!/usr/bin/env python3
"""Bounded accuracy and speed comparison for workflow APIs versus Indigo."""

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

ROOT = Path(__file__).resolve().parents[2]
REACTION_CORPUS = ROOT / "validation/reaction_product_parity_strata_v2.json"
STANDARDIZATION_CORPORA = [
    ROOT / "validation/standardization_phase1_fixtures.jsonl",
    ROOT / "validation/standardization_phase1_holdout.jsonl",
]
RGROUP_CORPUS = Path(__file__).with_name("rgroup_cases.json")
EXPECTED_CHEMATIC = "1.1.0"
EXPECTED_INDIGO = "1.46.0"
EXPECTED_RDKIT = "2026.03.6"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_commit() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def canonical_rdkit(smiles: str) -> str:
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


def ratio_summary(indigo: list[int], chematic: list[int]) -> dict[str, Any]:
    ratios = [left / right for left, right in zip(indigo, chematic)]
    rng = random.Random(17)
    medians = [
        statistics.median(ratios[rng.randrange(len(ratios))] for _ in ratios)
        for _ in range(10_000)
    ]
    medians.sort()
    return {
        "indigo_over_chematic_median": statistics.median(ratios),
        "ci95": [medians[249], medians[9749]],
        "winning_blocks": sum(value > 1.0 for value in ratios),
        "blocks": len(ratios),
    }


def timed_blocks(
    chematic_fn: Callable[[], Any],
    indigo_fn: Callable[[], Any],
    blocks: int,
) -> dict[str, Any]:
    samples = {"chematic": [], "indigo": []}
    chematic_fn()
    indigo_fn()
    functions = {"chematic": chematic_fn, "indigo": indigo_fn}
    for block in range(blocks):
        order = ["chematic", "indigo"] if block % 2 == 0 else ["indigo", "chematic"]
        for engine in order:
            start = time.perf_counter_ns()
            functions[engine]()
            samples[engine].append(time.perf_counter_ns() - start)
    return {"samples_ns": samples, "comparison": ratio_summary(samples["indigo"], samples["chematic"])}


def standardization_rows() -> list[dict[str, str]]:
    rows = []
    for path in STANDARDIZATION_CORPORA:
        for line in path.read_text().splitlines():
            row = json.loads(line)
            expected = row.get("expected_kept_fragment")
            if row.get("input") and expected:
                rows.append({"id": row["id"], "input": row["input"], "expected": expected})
    return rows


def standardization_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from indigo import Indigo

    rows = standardization_rows()
    indigo = Indigo()
    indigo.setOption("standardize-keep-largest", True)

    expected = [canonical_rdkit(row["expected"]) for row in rows]

    def chematic_values(selected: list[dict[str, str]] = rows) -> list[str]:
        return [
            canonical_rdkit(chematic.from_smiles(row["input"]).largest_fragment().smiles)
            for row in selected
        ]

    def indigo_values(selected: list[dict[str, str]] = rows) -> list[str]:
        values = []
        for row in selected:
            mol = indigo.loadMolecule(row["input"])
            mol.standardize()
            values.append(canonical_rdkit(mol.smiles()))
        return values

    observed = {"chematic": chematic_values(), "indigo": indigo_values()}
    accuracy = {
        engine: {
            "match": sum(value == target for value, target in zip(values, expected)),
            "total": len(rows),
            "mismatches": [
                {"id": row["id"], "expected": target, "actual": value}
                for row, value, target in zip(rows, values, expected)
                if value != target
            ],
        }
        for engine, values in observed.items()
    }
    timed_rows = [
        row
        for index, row in enumerate(rows)
        if all(observed[engine][index] == expected[index] for engine in observed)
    ]
    timing = timed_blocks(
        lambda: chematic_values(timed_rows),
        lambda: indigo_values(timed_rows),
        blocks,
    )
    timing["rows"] = len(timed_rows)
    timing["selection"] = "intersection of rows both engines matched exactly"
    return {"rows": len(rows), "accuracy": accuracy, "timing": timing}


def rdkit_reaction_products(smirks: str, reactants: list[str]) -> list[str]:
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions

    reaction = rdChemReactions.ReactionFromSmarts(smirks)
    molecules = [Chem.MolFromSmiles(text) for text in reactants]
    if reaction is None or any(mol is None for mol in molecules):
        raise ValueError("RDKit could not parse reaction row")
    products = set()
    for product_set in reaction.RunReactants(tuple(molecules)):
        try:
            products.add(".".join(sorted(Chem.MolToSmiles(mol, True) for mol in product_set)))
        except Exception:
            continue
    return sorted(products)


def reaction_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from indigo import Indigo

    all_rows = json.loads(REACTION_CORPUS.read_text())["cases"]
    rows = [row for row in all_rows if len(row["reactants"]) == 1]
    expected = {}
    supported = []
    for row in rows:
        try:
            products = rdkit_reaction_products(row["smirks"], row["reactants"])
        except Exception:
            continue
        # Indigo's transform() surface returns one transformed molecule rather
        # than RDKit/CheMatic's list of product sets. Keep this lane semantic:
        # only rows with one single-product answer are comparable.
        if len(products) == 1 and "." not in products[0]:
            expected[row["id"]] = products
            supported.append(row)
    rows = supported

    def chematic_values(selected: list[dict[str, Any]] = rows) -> dict[str, list[str]]:
        values = {}
        for row in selected:
            reactant = chematic.from_smiles(row["reactants"][0])
            result = chematic.run_smirks_checked(
                row["smirks"], [reactant], rdkit_compat=True
            )
            product_sets = result["products"]
            values[row["id"]] = sorted(
                {".".join(sorted(product.rdkit_smiles for product in products)) for products in product_sets}
            )
        return values

    def indigo_values(selected: list[dict[str, Any]] = rows) -> dict[str, list[str]]:
        indigo = Indigo()
        values = {}
        for row in selected:
            try:
                reaction = indigo.loadQueryReaction(row["smirks"])
                molecule = indigo.loadMolecule(row["reactants"][0])
                mapping = indigo.transform(reaction, molecule)
                if mapping is None:
                    values[row["id"]] = []
                else:
                    values[row["id"]] = [canonical_rdkit(molecule.smiles())]
            except Exception:
                values[row["id"]] = []
        return values

    observed = {"chematic": chematic_values(), "indigo": indigo_values()}
    accuracy = {
        engine: {
            "match": sum(values[row["id"]] == expected[row["id"]] for row in rows),
            "total": len(rows),
            "mismatches": [
                {"id": row["id"], "expected": expected[row["id"]], "actual": values[row["id"]]}
                for row in rows
                if values[row["id"]] != expected[row["id"]]
            ],
        }
        for engine, values in observed.items()
    }
    timed_rows = [
        row
        for row in rows
        if all(observed[engine][row["id"]] == expected[row["id"]] for engine in observed)
    ]
    timing = timed_blocks(
        lambda: chematic_values(timed_rows),
        lambda: indigo_values(timed_rows),
        blocks,
    )
    timing["rows"] = len(timed_rows)
    timing["selection"] = "intersection of rows both engines matched exactly"
    return {"rows": len(rows), "accuracy": accuracy, "timing": timing}


def indigo_payloads(indigo, core: str, molecules: list[str]) -> list[list[str]]:
    query = indigo.loadQueryMolecule(core)
    decomposer = indigo.createDecomposer(query)
    rows = []
    for text in molecules:
        molecule = indigo.loadMolecule(text)
        try:
            item = decomposer.decomposeMolecule(molecule)
            decompositions = list(item.iterateDecompositions())
        except Exception:
            rows.append([])
            continue
        if not decompositions:
            rows.append([])
            continue
        result = decompositions[0].decomposedMoleculeWithRGroups()
        payloads = []
        for rgroup in result.iterateRGroups():
            for fragment in rgroup.iterateRGroupFragments():
                from rdkit import Chem

                mol = Chem.MolFromMolBlock(fragment.molfile(), removeHs=False)
                if mol is None:
                    raise ValueError("RDKit could not parse Indigo R-group fragment")
                payloads.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
        rows.append(sorted(payloads))
    return rows


def rgroup_gate(blocks: int) -> dict[str, Any]:
    import chematic
    from indigo import Indigo

    cases = json.loads(RGROUP_CORPUS.read_text())["cases"]
    expected = [sorted(row) for case in cases for row in case["expected"]]

    def chematic_values() -> list[list[str]]:
        rows = []
        for case in cases:
            molecules = [chematic.from_smiles(text) for text in case["molecules"]]
            results = chematic.rgroup_decompose_labeled(case["core"], molecules)
            for result in results:
                if result is None:
                    rows.append([])
                else:
                    rows.append(sorted(payload_without_dummy(value) for key, value in result.items() if key.startswith("R")))
        return rows

    def indigo_values() -> list[list[str]]:
        indigo = Indigo()
        rows = []
        for case in cases:
            rows.extend(indigo_payloads(indigo, case.get("indigo_core", case["core"]), case["molecules"]))
        return rows

    observed = {"chematic": chematic_values(), "indigo": indigo_values()}
    accuracy = {
        engine: {
            "match": sum(value == target for value, target in zip(values, expected)),
            "total": len(expected),
            "mismatches": [
                {"row": index, "expected": target, "actual": value}
                for index, (value, target) in enumerate(zip(values, expected))
                if value != target
            ],
        }
        for engine, values in observed.items()
    }
    timing = None
    if all(result["match"] == len(expected) for result in accuracy.values()):
        timing = timed_blocks(chematic_values, indigo_values, blocks)
        timing["rows"] = len(expected)
        timing["selection"] = "all rows; both engines matched exactly"
    return {"rows": len(expected), "accuracy": accuracy, "timing": timing}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=21)
    args = parser.parse_args()
    if args.blocks < 3:
        parser.error("--blocks must be at least 3")

    import chematic
    import indigo
    import rdkit

    versions = {
        "chematic": chematic.__version__,
        "indigo": indigo.__version__,
        "rdkit": rdkit.__version__,
    }
    required = {
        "chematic": EXPECTED_CHEMATIC,
        "indigo": EXPECTED_INDIGO,
        "rdkit": EXPECTED_RDKIT,
    }
    if versions != required:
        raise SystemExit(f"comparator version mismatch: expected {required}, found {versions}")

    report = {
        "schema_version": 1,
        "source_commit": source_commit(),
        "engines": versions,
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
        },
        "corpora": {
            "reaction": {"path": str(REACTION_CORPUS.relative_to(ROOT)), "sha256": file_hash(REACTION_CORPUS)},
            "standardization": [
                {"path": str(path.relative_to(ROOT)), "sha256": file_hash(path)}
                for path in STANDARDIZATION_CORPORA
            ],
            "rgroup": {"path": str(RGROUP_CORPUS.relative_to(ROOT)), "sha256": file_hash(RGROUP_CORPUS)},
        },
        "blocks": args.blocks,
        "workflows": {
            "standardization": standardization_gate(args.blocks),
            "reaction": reaction_gate(args.blocks),
            "rgroup": rgroup_gate(args.blocks),
        },
        "claim_boundary": "Results apply only to the named versions, corpora, observers, source commit, and host.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
