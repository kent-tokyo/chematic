#!/usr/bin/env python3
"""Time the bounded RDKit side of the conformer-ensemble comparison.

This script deliberately records raw samples. It does not decide that the
RDKit and chematic ensembles are chemically equivalent; that belongs to the
separate quality gate described in the dated benchmark record.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import time

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


SMILES = (
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
)


def run(num_threads: int) -> float:
    started = time.perf_counter()
    for source in SMILES:
        molecule = Chem.AddHs(Chem.MolFromSmiles(source))
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 42
        parameters.numThreads = num_threads
        identifiers = AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=10,
            params=parameters,
        )
        if len(identifiers) != 10:
            raise RuntimeError(f"expected 10 conformers for {source}, got {len(identifiers)}")
        results = AllChem.MMFFOptimizeMoleculeConfs(
            molecule,
            numThreads=num_threads,
            maxIters=200,
        )
        if len(results) != 10:
            raise RuntimeError(f"expected 10 MMFF results for {source}, got {len(results)}")
    return (time.perf_counter() - started) * 1_000.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=21)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")

    single_thread_ms: list[float] = []
    all_threads_ms: list[float] = []
    for block in range(args.repeats):
        if block % 2 == 0:
            single_thread_ms.append(run(1))
            all_threads_ms.append(run(0))
        else:
            all_threads_ms.append(run(0))
            single_thread_ms.append(run(1))

    print(
        json.dumps(
            {
                "schema": "rdkit-conformer-ensemble-bounded/v1",
                "python": platform.python_version(),
                "rdkit": rdBase.rdkitVersion,
                "molecules": list(SMILES),
                "conformers_per_molecule": 10,
                "random_seed": 42,
                "mmff_max_iterations": 200,
                "repeats": args.repeats,
                "single_thread_samples_ms": single_thread_ms,
                "all_threads_samples_ms": all_threads_ms,
                "single_thread_median_ms": statistics.median(single_thread_ms),
                "all_threads_median_ms": statistics.median(all_threads_ms),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
