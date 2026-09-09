#!/usr/bin/env python3
"""Measure selected-column scaling for chematic.bulk.descriptors_array.

This is a release-artifact benchmark, not a Rust source benchmark.  It checks
the public Python contract at the same time as measuring wall time:

* the requested column set is preserved (the mapping's insertion order is not
  treated as a contract);
* every requested column has one value per valid input molecule;
* repeated calls have the same SHA-256 digest; and
* Python-visible allocations are recorded with tracemalloc.

Native allocations made by the Rust extension are intentionally not claimed by
tracemalloc.  Use a platform profiler or RSS measurement when native allocation
accounting is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import tracemalloc
from pathlib import Path

import numpy as np

import chematic


THREE = ["mw", "logp", "tpsa"]
EIGHT = THREE + ["hbd", "hba", "qed", "wiener_index", "kappa1"]
ALL = [
    "mw",
    "exact_mass",
    "tpsa",
    "logp",
    "molar_refractivity",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "ring_count",
    "aromatic_ring_count",
    "num_heteroatoms",
    "num_stereocenters",
    "num_spiro_atoms",
    "num_bridgehead_atoms",
    "fsp3",
    "qed",
    "sa_score",
    "formal_charge",
    "labute_asa",
    "bertz_ct",
    "wiener_index",
    "kappa1",
    "kappa2",
    "kappa3",
    "chi0",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi0v",
    "chi1v",
    "chi2v",
    "chi3v",
    "chi4v",
    "num_aromatic_heterocycles",
    "num_aliphatic_heterocycles",
    "num_saturated_rings",
    "num_aliphatic_rings",
    "num_unspecified_stereocenters",
    "sum_estate",
    "max_estate",
    "min_estate",
    "lipinski_passes",
    "veber_passes",
    "egan_passes",
    "ghose_passes",
    "reos_passes",
    "pains_passes",
    "bbb_passes",
    "bbb_score",
    "caco2",
    "herg_risk",
    "cyp3a4_risk",
    "pka_acid",
    "pka_base",
    "schultz_mti",
    "gutman_mti",
    "vabc",
    "gravitational_index",
]


def read_smiles(path: Path) -> list[str]:
    values: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(line.split()[0])
    return values


def digest(result: dict[str, np.ndarray], columns: list[str]) -> str:
    hasher = hashlib.sha256()
    for name in columns:
        array = np.asarray(result[name])
        hasher.update(name.encode("utf-8"))
        hasher.update(str(array.dtype).encode("ascii"))
        hasher.update(array.tobytes(order="C"))
    return hasher.hexdigest()


def measure(smiles: list[str], columns: list[str], repeats: int) -> dict[str, object]:
    result = chematic.bulk.descriptors_array(smiles, columns)
    rows = len(next(iter(result.values()))) if result else 0
    if set(result) != set(columns):
        raise AssertionError(f"column set changed: {list(result)!r}")
    if any(len(np.asarray(result[name])) != rows for name in columns):
        raise AssertionError("descriptor columns have different row counts")
    first_digest = digest(result, columns)

    tracemalloc.start()
    started = time.perf_counter()
    last_digest = first_digest
    for _ in range(repeats):
        result = chematic.bulk.descriptors_array(smiles, columns)
        last_digest = digest(result, columns)
    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    if last_digest != first_digest:
        raise AssertionError("descriptor output is not deterministic")

    return {
        "columns": columns,
        "column_count": len(columns),
        "rows": rows,
        "repeats": repeats,
        "seconds": elapsed,
        "calls_per_second": repeats / elapsed if elapsed else None,
        "digest_sha256": first_digest,
        "python_tracemalloc_peak_bytes": peak,
        "native_allocation_note": "tracemalloc excludes allocations inside the Rust extension",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--limit", type=int, help="benchmark only the first N input SMILES")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")

    smiles = read_smiles(args.corpus)
    if args.limit is not None:
        if args.limit < 1:
            parser.error("--limit must be positive")
        smiles = smiles[: args.limit]
    reports = {
        "3": measure(smiles, THREE, args.repeats),
        "8": measure(smiles, EIGHT, args.repeats),
        "all": measure(smiles, ALL, args.repeats),
    }
    report = {
        "operation": "chematic.bulk.descriptors_array",
        "input": str(args.corpus),
        "input_smiles": len(smiles),
        "reports": reports,
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
