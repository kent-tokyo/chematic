#!/usr/bin/env python3
"""Paired A6 conformer-generation throughput: chematic versus RDKit.

Parsing is excluded. Both engines generate the requested number of candidates,
run MMFF94 with a 200-iteration budget, use all available CPU cores within one
molecule, and process molecules in the same fixed A6 order. Repeats alternate
engine order to reduce thermal/order bias.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

import chematic
from public_package_3d_chematic import (
    RANDOM_SEED,
    load_rows,
    pipeline_config,
    run_lowest_of_n,
)


RDLogger.DisableLog("rdApp.*")


def run_rdkit(molecules: list[Chem.Mol], count: int) -> tuple[int, int]:
    successes = 0
    conformers = 0
    for source in molecules:
        mol = Chem.AddHs(Chem.Mol(source))
        params = AllChem.ETKDGv3()
        params.randomSeed = RANDOM_SEED
        params.numThreads = 0
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=count, params=params))
        if not ids:
            continue
        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            numThreads=0,
            maxIters=200,
        )
        finite = [energy for status, energy in results if status >= 0]
        if finite:
            successes += 1
            conformers += len(finite)
    return successes, conformers


def run_chematic(molecules: list[object], count: int) -> tuple[int, int]:
    config = pipeline_config(
        chematic,
        "mmff94_bond_angle_strict",
        stereo_safe=True,
    )
    successes = 0
    conformers = 0
    for mol in molecules:
        result = run_lowest_of_n(
            mol,
            chematic,
            config,
            best_of_n=count,
            explicit_h_ensemble=False,
            full_h_refine=True,
        )
        if result.get("status") == "success":
            successes += 1
            conformers += int(result.get("explicit_h_refined", 0))
    return successes, conformers


def timed(callable_) -> tuple[float, tuple[int, int]]:
    started = time.perf_counter()
    result = callable_()
    return time.perf_counter() - started, result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tiers", choices=("A", "B", "AB"), default="AB")
    parser.add_argument("--count", type=int)
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.conformers <= 0 or args.repeats <= 0:
        parser.error("--conformers and --repeats must be positive")

    rows = load_rows(args.tiers)
    if args.count is not None:
        rows = rows[: args.count]
    smiles = [str(row["smiles"]) for row in rows]
    chematic_molecules = [chematic.from_smiles(value) for value in smiles]
    rdkit_molecules = [Chem.MolFromSmiles(value) for value in smiles]
    if any(mol is None for mol in rdkit_molecules):
        raise RuntimeError("fixed A6 corpus contains an RDKit parse failure")

    samples = {"chematic": [], "rdkit": []}
    coverage: dict[str, list[tuple[int, int]]] = {"chematic": [], "rdkit": []}
    runners = {
        "chematic": lambda: run_chematic(chematic_molecules, args.conformers),
        "rdkit": lambda: run_rdkit(rdkit_molecules, args.conformers),
    }
    # Untimed warm-up covers imports, thread pools, and lazy parameter tables.
    run_chematic(chematic_molecules[:3], min(2, args.conformers))
    run_rdkit(rdkit_molecules[:3], min(2, args.conformers))
    order_log: list[list[str]] = []
    for repeat in range(args.repeats):
        order = ["chematic", "rdkit"] if repeat % 2 == 0 else ["rdkit", "chematic"]
        order_log.append(order)
        for engine in order:
            elapsed, result = timed(runners[engine])
            samples[engine].append(elapsed)
            coverage[engine].append(result)

    medians = {engine: statistics.median(values) for engine, values in samples.items()}
    report = {
        "schema": "a6-conformer-speed-vs-rdkit/v1",
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", None),
        "molecules": len(rows),
        "conformers_per_molecule": args.conformers,
        "repeats": args.repeats,
        "parsing_included": False,
        "order": order_log,
        "seconds": samples,
        "median_seconds": medians,
        "chematic_speedup_vs_rdkit": medians["rdkit"] / medians["chematic"],
        "coverage": coverage,
        "boundary": (
            "Same host and fixed corpus. Both engines use all cores within one molecule, "
            "the same candidate count, and a 200-iteration MMFF94 budget. This measures "
            "throughput, not chemical quality or public-registry artifacts."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
