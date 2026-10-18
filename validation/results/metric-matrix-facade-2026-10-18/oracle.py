"""Live RDKit metric namespace and independent math.dist/binary-set checks."""

import hashlib
import json
import math
from pathlib import Path
import random

import numpy as np
import rdkit
from rdkit import DataStructs as D, RDLogger
from rdkit.DataManip.Metric import rdMetricMatrixCalc as R

from chematic.rdkit_compat import ExplicitBitVect, rdMetricMatrixCalc as C

HERE = Path(__file__).parent
rng = random.Random(9601)
corpus = hashlib.sha256()
euclidean_comparisons = binary_comparisons = compared_pairs = 0
max_euclidean_error = max_math_error = max_binary_error = 0.0


def fingerprint(cls, width, bits):
    result = cls(width)
    for bit in bits:
        result.SetBit(bit)
    return result


for case in range(128):
    count, features = 2 + case % 31, 1 + case % 17
    rows = [[rng.uniform(-20, 20) for _ in range(features)] for _ in range(count)]
    corpus.update(json.dumps(rows).encode())
    reference = R.GetEuclideanDistMat(rows)
    independent = [math.dist(rows[i], rows[j]) for i in range(count) for j in range(i)]
    # A strided NumPy source carries exactly the same numeric rows.
    backing = np.empty((count, features * 2))
    backing[:, ::2] = rows
    for source in [
        rows,
        tuple(map(tuple, rows)),
        backing[:, ::2],
        (iter(row) for row in rows),
    ]:
        actual = C.GetEuclideanDistMat(descripMat=source)
        assert actual.dtype == np.float64 and actual.flags.owndata
        np.testing.assert_allclose(actual, reference, atol=1e-12, rtol=1e-14)
        np.testing.assert_allclose(actual, independent, atol=1e-12, rtol=1e-14)
        max_euclidean_error = max(
            max_euclidean_error, float(np.max(np.abs(actual - reference)))
        )
        max_math_error = max(
            max_math_error, float(np.max(np.abs(actual - independent)))
        )
        euclidean_comparisons += 1
        compared_pairs += len(actual)

    width = [32, 64, 128, 257][case % 4]
    sets = [
        set(rng.sample(range(width), rng.randrange(min(width, 12) + 1)))
        for _ in range(count)
    ]
    if case % 3 == 0:
        sets[0] = sets[1] = set()
    corpus.update(json.dumps([width, [sorted(s) for s in sets]]).encode())
    explicit = [fingerprint(D.ExplicitBitVect, width, bits) for bits in sets]
    sparse = [fingerprint(D.SparseBitVect, width, bits) for bits in sets]
    schematic = [fingerprint(ExplicitBitVect, width, bits) for bits in sets]
    for name in ["GetTanimotoDistMat", "GetTanimotoSimMat"]:
        distance = name == "GetTanimotoDistMat"
        reference = getattr(R, name)(explicit)
        sparse_reference = getattr(R, name)(sparse)
        np.testing.assert_array_equal(reference, sparse_reference)
        independent = []
        for i in range(count):
            for j in range(i):
                union = sets[i] | sets[j]
                sim = len(sets[i] & sets[j]) / len(union) if union else 0.0
                independent.append(1 - sim if distance else sim)
        sources = [
            explicit,
            tuple(sparse),
            (explicit[i] if i % 2 else sparse[i] for i in range(count)),
            schematic,
        ]
        for source in sources:
            actual = getattr(C, name)(bitVectList=source)
            np.testing.assert_allclose(actual, reference, atol=1e-15, rtol=0)
            np.testing.assert_array_equal(actual, independent)
            max_binary_error = max(
                max_binary_error, float(np.max(np.abs(actual - reference)))
            )
            binary_comparisons += 1
            compared_pairs += len(actual)

diagnostics = []
for scale in [1e200, 1e-200]:
    rows = [[0.0], [scale]]
    actual = C.GetEuclideanDistMat(rows)[0]
    reference = R.GetEuclideanDistMat(rows)[0]
    assert actual == math.dist(*rows) == scale
    diagnostics.append(
        dict(
            case="extreme_scale",
            scale=scale,
            chematic=float(actual),
            rdkit=str(reference),
            independent=scale,
        )
    )
# Publicly documented narrower/extended profiles are separate from ordinary parity.
for name in ["GetTanimotoSimMat", "GetTanimotoDistMat"]:
    RDLogger.DisableLog("rdApp.error")
    try:
        getattr(R, name)([])
    except RuntimeError:
        diagnostics.append(
            dict(
                case="empty_ranking",
                metric=name,
                rdkit="RuntimeError",
                chematic=getattr(C, name)([]).tolist(),
            )
        )
    else:
        raise AssertionError("reference empty-ranking behavior changed")
    finally:
        RDLogger.EnableLog("rdApp.error")
    mixed = [
        fingerprint(D.ExplicitBitVect, 8, [1]),
        fingerprint(D.SparseBitVect, 8, [1, 2]),
    ]
    try:
        getattr(R, name)(mixed)
    except ValueError:
        diagnostics.append(
            dict(
                case="mixed_source_types",
                metric=name,
                rdkit="ValueError",
                chematic=getattr(C, name)(mixed).tolist(),
            )
        )
    else:
        raise AssertionError("reference mixed-source behavior changed")

high = fingerprint(D.SparseBitVect, 2**32 - 1, [2**32 - 2])
signed_view = list(high.GetOnBits())
assert signed_view == [-2]
try:
    C.GetTanimotoSimMat([high, high])
except OverflowError:
    diagnostics.append(
        dict(case="high_signed_sparse_view", view=signed_view, chematic="OverflowError")
    )
else:
    raise AssertionError("negative signed IDs must not be reinterpreted")
unequal = [D.ExplicitBitVect(8), D.ExplicitBitVect(16)]
reference = R.GetTanimotoSimMat(unequal).tolist()
try:
    C.GetTanimotoSimMat(unequal)
except ValueError:
    diagnostics.append(
        dict(case="unequal_widths", rdkit=reference, chematic="ValueError")
    )
else:
    raise AssertionError("unequal declared widths must refuse")

summary = dict(
    rdkit=rdkit.__version__,
    numpy=np.__version__,
    seed=9601,
    euclidean_comparisons=euclidean_comparisons,
    binary_comparisons=binary_comparisons,
    compared_pairs=compared_pairs,
    max_rdkit_euclidean_error=max_euclidean_error,
    max_math_dist_error=max_math_error,
    max_rdkit_binary_error=max_binary_error,
    diagnostics=len(diagnostics),
    inputs_sha256=corpus.hexdigest(),
    scope="Ordinary namespace/output comparisons; numeric and input-profile differences are recorded separately",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(HERE / "reference-diagnostics.json").write_text(
    json.dumps(diagnostics, indent=2, allow_nan=False) + "\n"
)
print(json.dumps(summary, indent=2))
