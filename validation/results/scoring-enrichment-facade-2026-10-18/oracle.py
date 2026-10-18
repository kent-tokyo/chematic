"""Live RDKit ordinary-cutoff comparison and independent exact count checks."""

from fractions import Fraction
import hashlib
import inspect
import json
import math
from pathlib import Path
import random

import numpy as np
import rdkit
from rdkit.ML.Scoring import Scoring as R
from chematic.rdkit_compat import Scoring as C

HERE = Path(__file__).parent
rng = random.Random(7981)
comparisons = paired_records = 0
max_reference_error = max_count_error = 0.0
inputs_hash = hashlib.sha256()


def expected(labels, fractions):
    active = sum(labels)
    n = len(labels)
    return [
        float(Fraction(sum(labels[: math.ceil(n * f)]) * n, math.ceil(n * f) * active))
        for f in fractions
    ]


for n in range(2, 66):
    for trial in range(16):
        labels = [rng.randrange(2) for _ in range(n)]
        labels[rng.randrange(n)] = 1
        rows = [(float(n - i), label) for i, label in enumerate(labels)]
        selected = rng.randrange(1, n)
        fraction = (selected - 0.25) / n
        column = -1 if trial % 2 else 1
        want = R.CalcEnrichment(rows, column, [fraction])
        count_formula = expected(labels, [fraction])
        assert len(want) == 1
        inputs_hash.update(json.dumps([labels, column, fraction]).encode())
        for representation in [rows, tuple(rows), np.array(rows), (r for r in rows)]:
            actual = C.CalcEnrichment(representation, column, (f for f in [fraction]))
            max_reference_error = max(max_reference_error, abs(actual[0] - want[0]))
            max_count_error = max(max_count_error, abs(actual[0] - count_formula[0]))
            assert abs(actual[0] - want[0]) < 1e-12
            assert abs(actual[0] - count_formula[0]) < 1e-12
            comparisons += 1
            paired_records += n

boundaries = []
for labels, fractions in [
    ([1], [0.01, 0.9, 1.0]),
    ([1, 1, 0, 0, 0], [0.9]),
    ([1, 1, 0, 0, 0], [1.0]),
    ([1, 0, 1, 0, 0], [1.0, 0.2, 0.01, 0.2, 0.5, 0.9]),
    ([1, 0, 1, 0, 0], [0.1, 0.11, 0.12]),
    ([0, 0, 1, 1], [0.5, 1.0]),
]:
    rows = [(float(len(labels) - i), label) for i, label in enumerate(labels)]
    actual = C.CalcEnrichment(rows, 1, fractions)
    formula = expected(labels, fractions)
    assert len(actual) == len(fractions)
    np.testing.assert_allclose(actual, formula, atol=1e-12, rtol=1e-12)
    reference = R.CalcEnrichment(rows, 1, fractions)
    boundaries.append(
        dict(
            labels=labels,
            fractions=fractions,
            chematic=actual,
            independent=formula,
            rdkit=reference,
        )
    )
assert C.CalcEnrichment([(1.0, 1), (0.0, 0)], 1, [1.0]) == [1.0]
summary = dict(
    seed=7981,
    rdkit=rdkit.__version__,
    numpy=np.__version__,
    ordinary_rankings=1024,
    ordinary_representation_comparisons=comparisons,
    paired_records=paired_records,
    max_rdkit_error=max_reference_error,
    max_independent_error=max_count_error,
    corrected_boundary_cases=len(boundaries),
    upstream_boundary_differences=sum(b["chematic"] != b["rdkit"] for b in boundaries),
    inputs_sha256=inputs_hash.hexdigest(),
    rdkit_function_sha256=hashlib.sha256(
        inspect.getsource(R.CalcEnrichment).encode()
    ).hexdigest(),
    scope="Ordinary single-cutoff parity; boundary corrections validated separately against exact count formula",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(HERE / "boundary-comparisons.json").write_text(json.dumps(boundaries, indent=2) + "\n")
print(json.dumps(summary, indent=2))
