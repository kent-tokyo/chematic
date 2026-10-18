"""Live RDKit oracle; PYTHONPATH must include the candidate and reference installs."""

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import rdkit
from rdkit.DataStructs import ExplicitBitVect, SparseBitVect
from rdkit.ML.InfoTheory import rdInfoTheory as R
from chematic.rdkit_compat import rdInfoTheory as C

SEED = 96010204
rng = np.random.default_rng(SEED)
comparisons = scores_checked = votes = ordering_differences = 0
max_score_error = max_independent_error = 0.0
anomalies = []
input_hash = hashlib.sha256()


def independent(hits, totals, preferred, metric):
    fractions = [Fraction(k, n) for k, n in zip(hits, totals)]
    if 0 < len(preferred) < len(totals):
        if max(fractions[c] for c in preferred) < max(
            fractions[c] for c in range(len(totals)) if c not in preferred
        ):
            return 0.0
    table = [hits, [n - k for n, k in zip(totals, hits)]]
    size = sum(totals)
    if metric == 4:
        value = Fraction(0)
        for row in table:
            for c, observed in enumerate(row):
                expected = Fraction(sum(row) * totals[c], size)
                if expected:
                    value += (observed - expected) ** 2 / expected
        return float(value)
    # Direct mutual information, independent of the entropy-difference formula.
    return sum(
        float(Fraction(k, size))
        * math.log2(float(Fraction(k * size, sum(row) * totals[c])))
        for row in table
        for c, k in enumerate(row)
        if k
    )


for classes, trial in itertools.product((2, 3, 4), range(24)):
    width = 9
    totals = [int(x) for x in rng.integers(3, 18, size=classes)]
    hits = [[int(rng.integers(1, n)) for n in totals] for _ in range(width)]
    preferred_sets = [[], [0], list(range(1, classes)), [classes - 1]]
    input_hash.update(json.dumps([totals, hits], separators=(",", ":")).encode())
    for metric, preferred, vector_type in itertools.product(
        (2, 4), preferred_sets, (ExplicitBitVect, SparseBitVect)
    ):
        ours = C.InfoBitRanker(width, classes, metric)
        reference = R.InfoBitRanker(width, classes, R.InfoType.values[metric])
        ours.SetBiasList(iter(preferred))
        reference.SetBiasList(preferred)
        for label, count in enumerate(totals):
            for row in range(count):
                vector = vector_type(width)
                for bit in range(width):
                    if row < hits[bit][label]:
                        vector.SetBit(bit)
                ours.AccumulateVotes(vector, label)
                reference.AccumulateVotes(vector, label)
                votes += 1
        actual = ours.GetTopN(width)
        expected = reference.GetTopN(width)
        for record in actual:
            bit = int(record[0])
            exact = independent(hits[bit], totals, preferred, metric)
            error = abs(float(record[1]) - exact)
            max_independent_error = max(max_independent_error, error)
            assert error < 1e-12
            np.testing.assert_array_equal(record[2:], hits[bit])
            scores_checked += 1
        ids = [int(x) for x in expected[:, 0]]
        if sorted(ids) != list(range(width)):
            anomalies.append(
                dict(
                    classes=classes,
                    trial=trial,
                    metric=metric,
                    vector=vector_type.__name__,
                    preferred=preferred,
                    totals=totals,
                    hits=hits,
                    ours=actual.tolist(),
                    rdkit=expected.tolist(),
                )
            )
            continue
        actual_by_id = actual[np.argsort(actual[:, 0])]
        expected_by_id = expected[np.argsort(expected[:, 0])]
        np.testing.assert_array_equal(actual_by_id[:, 2:], expected_by_id[:, 2:])
        np.testing.assert_allclose(
            actual_by_id[:, 1], expected_by_id[:, 1], atol=1e-12, rtol=1e-12
        )
        max_score_error = max(
            max_score_error,
            float(np.max(np.abs(actual_by_id[:, 1] - expected_by_id[:, 1]))),
        )
        ordering_differences += int(not np.array_equal(actual[:, 0], expected[:, 0]))
        comparisons += 1

summary = dict(
    seed=SEED,
    rdkit=rdkit.__version__,
    numpy=np.__version__,
    datasets=72,
    attempted_rankers=1152,
    complete_reference_comparisons=comparisons,
    anomalous_reference_rankers=len(anomalies),
    paired_votes=votes,
    independently_checked_scores=scores_checked,
    max_rdkit_score_error=max_score_error,
    max_independent_score_error=max_independent_error,
    ordering_differences=ordering_differences,
    inputs_sha256=input_hash.hexdigest(),
    scope="Per-bit scores/counts; deterministic chematic ties can differ from RDKit heap order.",
)
path = Path(__file__).parent
(path / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(path / "reference-anomalies.json").write_text(json.dumps(anomalies, indent=2) + "\n")
print(json.dumps(summary, indent=2))
