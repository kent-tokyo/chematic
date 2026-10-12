#!/usr/bin/env python3
"""Compare CheMatic and RDKit ETKDG + MMFF94 on the frozen A6 corpus.

Both engines use explicit hydrogens, the same conformer count, seed and MMFF94
iteration budget, all available CPU cores, and alternating execution order.
Parsing is excluded. The report records both throughput and the lowest finite
MMFF94 energy selected for every molecule.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import time
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

from public_package_3d_chematic import RANDOM_SEED, load_rows


RDLogger.DisableLog("rdApp.*")


def run_chematic(
    molecules: list[object],
    count: int,
    max_mmff_iters: int,
    retry_candidates: int,
    retry_max_iters: int,
) -> list[dict[str, object]]:
    rows = []
    for molecule in molecules:
        conformers = molecule.rdkit_embed_multiple_confs_mmff(
            count,
            random_seed=RANDOM_SEED,
            max_mmff_iters=max_mmff_iters,
            retry_candidates=retry_candidates,
            retry_mmff_iters=retry_max_iters,
        )
        finite = [(status, energy) for status, energy, _ in conformers if energy == energy]
        rows.append(
            {
                "generated": len(conformers),
                "converged": sum(status == 0 for status, _ in finite),
                "best_energy": min((energy for _, energy in finite), default=None),
            }
        )
    return rows


def run_rdkit(
    molecules: list[Chem.Mol], count: int, max_mmff_iters: int
) -> list[dict[str, object]]:
    rows = []
    for source in molecules:
        molecule = Chem.AddHs(Chem.Mol(source))
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = RANDOM_SEED
        parameters.numThreads = 0
        identifiers = list(
            AllChem.EmbedMultipleConfs(molecule, numConfs=count, params=parameters)
        )
        optimized = (
            AllChem.MMFFOptimizeMoleculeConfs(
                molecule, numThreads=0, maxIters=max_mmff_iters
            )
            if identifiers
            else []
        )
        finite = [(status, energy) for status, energy in optimized if energy == energy]
        rows.append(
            {
                "generated": len(identifiers),
                "converged": sum(status == 0 for status, _ in finite),
                "best_energy": min((energy for _, energy in finite), default=None),
            }
        )
    return rows


def summarize_quality(
    source_rows: list[dict[str, object]],
    chematic_rows: list[dict[str, object]],
    rdkit_rows: list[dict[str, object]],
) -> dict[str, object]:
    deltas = []
    failures = []
    for source, chematic_row, rdkit_row in zip(source_rows, chematic_rows, rdkit_rows):
        chematic_energy = chematic_row["best_energy"]
        rdkit_energy = rdkit_row["best_energy"]
        if chematic_energy is None or rdkit_energy is None:
            failures.append(
                {
                    "id": source.get("id"),
                    "chematic": chematic_energy,
                    "rdkit": rdkit_energy,
                }
            )
            continue
        deltas.append(float(chematic_energy) - float(rdkit_energy))
    ordered = sorted(deltas)
    return {
        "compared": len(deltas),
        "failures": failures,
        "median_delta_kcal_mol": statistics.median(deltas) if deltas else None,
        "p90_delta_kcal_mol": ordered[int(0.9 * (len(ordered) - 1))] if ordered else None,
        "max_delta_kcal_mol": max(deltas, default=None),
        "at_or_below_rdkit": sum(delta <= 1.0e-9 for delta in deltas),
        "within_1_kcal_mol": sum(delta <= 1.0 for delta in deltas),
        "within_2_kcal_mol": sum(delta <= 2.0 for delta in deltas),
    }


def per_molecule_results(
    source_rows: list[dict[str, object]],
    chematic_rows: list[dict[str, object]],
    rdkit_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    output = []
    for index, (source, chematic_row, rdkit_row) in enumerate(
        zip(source_rows, chematic_rows, rdkit_rows)
    ):
        chematic_energy = chematic_row["best_energy"]
        rdkit_energy = rdkit_row["best_energy"]
        output.append(
            {
                "row_index": index,
                "id": source.get("id"),
                "name": source.get("name"),
                "chematic": chematic_row,
                "rdkit": rdkit_row,
                "delta_kcal_mol": (
                    float(chematic_energy) - float(rdkit_energy)
                    if chematic_energy is not None and rdkit_energy is not None
                    else None
                ),
            }
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tiers", choices=("A", "B", "AB"), default="AB")
    parser.add_argument("--count", type=int)
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--max-mmff-iters", type=int, default=200)
    parser.add_argument("--retry-candidates", type=int, default=0)
    parser.add_argument("--retry-max-iters", type=int, default=2000)
    args = parser.parse_args()
    if (
        args.conformers <= 0
        or args.repeats <= 0
        or args.max_mmff_iters <= 0
        or args.retry_candidates < 0
    ):
        parser.error(
            "conformers/repeats/max MMFF iterations must be positive and retry candidates non-negative"
        )

    source_rows = load_rows(args.tiers)
    if args.count is not None:
        source_rows = source_rows[: args.count]
    smiles = [str(row["smiles"]) for row in source_rows]
    chematic_molecules = [chematic.from_smiles(value).add_hydrogens() for value in smiles]
    rdkit_molecules = [Chem.MolFromSmiles(value) for value in smiles]
    if any(molecule is None for molecule in rdkit_molecules):
        raise RuntimeError("fixed A6 corpus contains an RDKit parse failure")

    runners = {
        "chematic": lambda: run_chematic(
            chematic_molecules,
            args.conformers,
            args.max_mmff_iters,
            args.retry_candidates,
            args.retry_max_iters,
        ),
        "rdkit": lambda: run_rdkit(
            rdkit_molecules, args.conformers, args.max_mmff_iters
        ),
    }
    run_chematic(
        chematic_molecules[:2],
        min(args.conformers, 2),
        args.max_mmff_iters,
        args.retry_candidates,
        args.retry_max_iters,
    )
    run_rdkit(
        rdkit_molecules[:2], min(args.conformers, 2), args.max_mmff_iters
    )
    samples = {"chematic": [], "rdkit": []}
    results: dict[str, list[dict[str, object]]] = {}
    execution_order = []
    for repeat in range(args.repeats):
        order = ["chematic", "rdkit"] if repeat % 2 == 0 else ["rdkit", "chematic"]
        execution_order.append(order)
        for engine in order:
            started = time.perf_counter()
            results[engine] = runners[engine]()
            samples[engine].append(time.perf_counter() - started)

    medians = {engine: statistics.median(values) for engine, values in samples.items()}
    report = {
        "schema": "a6-rdkit-compatible-etkdg-mmff/v1",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.platform(),
        "python": platform.python_version(),
        "chematic": getattr(chematic, "__version__", None),
        "rdkit": rdBase.rdkitVersion,
        "molecules": len(source_rows),
        "conformers_per_molecule": args.conformers,
        "max_mmff_iters": args.max_mmff_iters,
        "adaptive_retry_candidates": args.retry_candidates,
        "adaptive_retry_max_iters": args.retry_max_iters,
        "chematic_etkdg_oracle": "RDKit 2026.09.1 ETKDGv3",
        "repeats": args.repeats,
        "execution_order": execution_order,
        "seconds": samples,
        "median_seconds": medians,
        "speedup_rdkit_over_chematic": medians["rdkit"] / medians["chematic"],
        "quality": summarize_quality(
            source_rows,
            results["chematic"],
            results["rdkit"],
        ),
        "rows": per_molecule_results(
            source_rows,
            results["chematic"],
            results["rdkit"],
        ),
        "coverage": {
            engine: {
                "generated": sum(int(row["generated"]) for row in rows),
                "converged": sum(int(row["converged"]) for row in rows),
            }
            for engine, rows in results.items()
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
