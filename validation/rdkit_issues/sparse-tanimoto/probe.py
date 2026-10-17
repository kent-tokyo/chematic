#!/usr/bin/env python3
"""Original bounded set corpus; run against candidate chematic and live RDKit.

PYTHONPATH=<candidate>:<rdkit> python3 probe.py <output-directory>
No upstream fixtures or implementation are copied. Empty/empty policy is
reported separately from ordinary nonempty-union comparisons.
"""
import gzip
import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

import chematic
import numpy as np
from rdkit import DataStructs, rdBase
from rdkit.DataManip.Metric import rdMetricMatrixCalc as calc


def reference(rows, metric):
    sets = list(map(set, rows))
    values = []
    for i in range(1, len(sets)):
        for j in range(i):
            union = sets[i] | sets[j]
            similarity = len(sets[i] & sets[j]) / len(union) if union else 1.0
            values.append(similarity if metric == "similarity" else 1 - similarity)
    return values


def vectors(rows, dimension, cls):
    result = []
    for row in rows:
        vector = cls(dimension)
        for bit in row:
            vector.SetBit(bit)
        result.append(vector)
    return result


def main():
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    records = []
    errors = []
    native_runs = oracle_runs = compared_values = 0
    max_error = 0.0
    for n, dimension, density, seed, empty_first in itertools.product(
        [2, 3, 4, 8, 16], [16, 128, 2048], [0.1, 0.8], range(4), [False, True]
    ):
        rng = random.Random(seed)
        rows = []
        for i in range(n):
            row = [j for j in range(dimension) if rng.random() < density]
            if not row:
                row = [rng.randrange(dimension)]
            row += row[:3]
            rng.shuffle(row)
            rows.append([] if empty_first and i == 0 else row)
        dense = vectors(rows, dimension, DataStructs.ExplicitBitVect)
        sparse = vectors(rows, dimension, DataStructs.SparseBitVect)
        record = dict(
            n=n,
            dimension=dimension,
            density=density,
            seed=seed,
            empty_first=empty_first,
            rows=rows,
            results={},
        )
        for metric in ["similarity", "distance"]:
            expected = reference(rows, metric)
            functions = (
                calc.GetTanimotoSimMat
                if metric == "similarity"
                else calc.GetTanimotoDistMat
            )
            results = {
                "set_reference": expected,
                "rdkit_explicit": list(functions(dense)),
                "rdkit_sparse": list(functions(sparse)),
            }
            oracle_runs += 2
            representations = {
                "lists": rows,
                "generators": (iter(row) for row in rows),
                "numpy": [np.asarray(row, dtype=np.uint32) for row in rows],
                "explicit_sources": dense,
                "sparse_sources": sparse,
            }
            for name, source in representations.items():
                results[name] = chematic.sparse_tanimoto_matrix_checked(
                    source, metric=metric
                )
                native_runs += 1
            for name, values in results.items():
                if name == "set_reference":
                    continue
                if len(values) != len(expected):
                    errors.append([len(records), metric, name, "length"])
                for a, b in zip(values, expected):
                    delta = abs(a - b)
                    max_error = max(max_error, delta)
                    compared_values += 1
                    if delta > 1e-15:
                        errors.append([len(records), metric, name, a, b])
            reversed_result = chematic.sparse_tanimoto_matrix_checked(
                rows[::-1], metric=metric
            )
            native_runs += 1
            if reversed_result != reference(rows[::-1], metric):
                errors.append([len(records), metric, "row permutation"])
            results["reversed_rows"] = reversed_result
            record["results"][metric] = results
        records.append(record)
    high_rows = [[0, 2**31 - 2], [2**31 - 2], [2**31 - 3]]
    high_sources = vectors(high_rows, 2**32 - 1, DataStructs.SparseBitVect)
    high = {"rows": high_rows}
    for metric in ["similarity", "distance"]:
        native = chematic.sparse_tanimoto_matrix_checked(high_sources, metric=metric)
        oracle = list(
            (
                calc.GetTanimotoSimMat
                if metric == "similarity"
                else calc.GetTanimotoDistMat
            )(high_sources)
        )
        expected = reference(high_rows, metric)
        assert native == oracle == expected
        high[metric] = dict(native=native, rdkit_sparse=oracle, reference=expected)
    wrapped = vectors([[0, 2**32 - 2]], 2**32 - 1, DataStructs.SparseBitVect)
    signed_view = {"rdkit_onbits": list(wrapped[0].GetOnBits())}
    try:
        signed_view["native"] = chematic.sparse_tanimoto_matrix_checked(wrapped)
    except OverflowError as error:
        signed_view["native_error"] = str(error)
    signed_view["interpretation"] = (
        "RDKit GetOnBits exposes -2 for bit 2**32-2 on this build. Native unsigned input refuses the negative view; caller-supplied unsigned IDs remain supported."
    )
    max_rows = [[0, 2**32 - 1], [2**32 - 1], [2**32 - 2]]
    assert chematic.sparse_tanimoto_matrix_checked(
        max_rows, metric="similarity"
    ) == reference(max_rows, "similarity")
    empty = vectors([[], []], 128, DataStructs.SparseBitVect)
    policy = {
        "native_similarity": chematic.sparse_tanimoto_matrix_checked(
            [[], []], metric="similarity"
        ),
        "native_distance": chematic.sparse_tanimoto_matrix_checked([[], []]),
        "rdkit_similarity": list(calc.GetTanimotoSimMat(empty)),
        "rdkit_distance": list(calc.GetTanimotoDistMat(empty)),
        "interpretation": "Native empty/empty similarity is 1, matching existing BitVec2048; current RDKit differs. This diagnostic is excluded from ordinary parity counts.",
    }
    singleton = {"native": chematic.sparse_tanimoto_matrix_checked([[1]])}
    try:
        singleton["rdkit"] = list(
            calc.GetTanimotoDistMat(vectors([[1]], 128, DataStructs.SparseBitVect))
        )
    except RuntimeError as error:
        singleton["rdkit_error"] = str(error)
    artifact = Path(chematic.chematic.__file__)
    summary = dict(
        schema_version=1,
        scope="supplied binary on-bit sets, not fingerprint-generation or full RDKit wrapper parity",
        rdkit_version=rdBase.rdkitVersion,
        numpy_version=np.__version__,
        python_version=sys.version,
        native_binary_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest(),
        datasets=len(records),
        native_runs=native_runs,
        oracle_runs=oracle_runs,
        compared_values=compared_values,
        max_absolute_error=max_error,
        ordinary_comparison_failures=len(errors),
        failures=errors,
        empty_empty_policy=policy,
        singleton_policy=singleton,
        high_sparse_ids=high,
        signed_rdkit_view=signed_view,
        u32_max_native_similarity=reference(max_rows, "similarity"),
    )
    with gzip.open(output / "corpus-results.json.gz", "wt") as stream:
        json.dump(records, stream, separators=(",", ":"))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
