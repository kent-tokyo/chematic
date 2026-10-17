"""Independent exact-rational Pearson tables plus live RDKit score/count oracle."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import rdkit
from rdkit import DataStructs
from rdkit.ML.InfoTheory import rdInfoTheory as R
from chematic.rdkit_compat import rdInfoTheory as C

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(20261018)
report = dict(rdkit=rdkit.__version__, numpy=np.__version__, seed=20261018,
              attempted_rankers=0, paired_rankers=0, paired_votes=0, rational_scores=0,
              reference_invalid_ids=[],
              reversed_rankers=0, doubled_rankers=0, degenerate_rankers=0,
              order_differences=0, max_chi_difference=0., max_entropy_difference=0.,
              max_rational_difference=0.)


def vectors(matrix, sparse):
    constructor = DataStructs.SparseBitVect if sparse else DataStructs.ExplicitBitVect
    result = []
    for row in matrix:
        vector = constructor(matrix.shape[1])
        for bit in np.flatnonzero(row):
            vector.SetBit(int(bit))
        result.append(vector)
    return result


def rank(module, items, labels, width, classes, metric):
    ranker = module.InfoBitRanker(width, classes, metric)
    for vector, label in zip(items, labels):
        ranker.AccumulateVotes(vector, int(label))
    return ranker.GetTopN(width)


def indexed(rows):
    return rows[np.argsort(rows[:, 0])]


def rational(matrix, labels, classes):
    totals = [int(np.count_nonzero(labels == c)) for c in range(classes)]
    n = len(labels)
    result = []
    for bit in range(matrix.shape[1]):
        on = [int(np.count_nonzero(matrix[labels == c, bit])) for c in range(classes)]
        score = Fraction(0)
        for observed in (on, [total - value for total, value in zip(totals, on)]):
            for value, total in zip(observed, totals):
                expected = Fraction(sum(observed) * total, n)
                if expected:
                    score += (value - expected) ** 2 / expected
        result.append(float(score))
    return np.array(result)


for width in (1, 3, 8, 31):
    for classes in (2, 3, 5):
        for trial in range(16):
            n = classes + 17 + trial * 3
            labels = np.arange(n) % classes
            rng.shuffle(labels)
            matrix = rng.integers(0, 2, size=(n, width), dtype=np.uint8)
            # Nondegenerate live reference: every bit has on and off observations.
            matrix[0] = 0
            matrix[1] = 1
            exact = rational(matrix, labels, classes)
            report['rational_scores'] += width
            for sparse in (False, True):
                items = vectors(matrix, sparse)
                for metric in (C.InfoType.ENTROPY, C.InfoType.CHISQUARE):
                    candidate = rank(C, items, labels, width, classes, metric)
                    reference = rank(R, items, labels, width, classes,
                                     getattr(R.InfoType, metric.name))
                    a, b = indexed(candidate), indexed(reference)
                    np.testing.assert_array_equal(a[:, 0], np.arange(width))
                    report['attempted_rankers'] += 1
                    report['paired_votes'] += n
                    if np.array_equal(b[:, 0], np.arange(width)):
                        np.testing.assert_array_equal(a[:, 2:], b[:, 2:])
                        np.testing.assert_allclose(a[:, 1], b[:, 1], rtol=1e-11, atol=1e-11)
                        difference = float(np.max(np.abs(a[:, 1] - b[:, 1])))
                        key = 'max_chi_difference' if metric == C.InfoType.CHISQUARE else 'max_entropy_difference'
                        report[key] = max(report[key], difference)
                        report['paired_rankers'] += 1
                        report['order_differences'] += int(not np.array_equal(candidate[:, 0], reference[:, 0]))
                    else:
                        report['reference_invalid_ids'].append(dict(
                            width=width, classes=classes, trial=trial, sparse=sparse,
                            metric=metric.name, reference=reference.tolist(),
                            candidate=candidate.tolist(), matrix=matrix.tolist(),
                            labels=labels.tolist()))
                    if metric != C.InfoType.CHISQUARE:
                        continue
                    np.testing.assert_allclose(a[:, 1], exact, rtol=1e-12, atol=1e-12)
                    report['max_rational_difference'] = max(report['max_rational_difference'], float(np.max(np.abs(a[:, 1] - exact))))
                    reverse = rank(C, items[::-1], labels[::-1], width, classes, metric)
                    np.testing.assert_array_equal(reverse, candidate)
                    twice = rank(C, items * 2, np.tile(labels, 2), width, classes, metric)
                    np.testing.assert_array_equal(twice[:, 0], candidate[:, 0])
                    np.testing.assert_allclose(twice[:, 1:], candidate[:, 1:] * 2, rtol=1e-12, atol=1e-12)
                    report['reversed_rankers'] += 1
                    report['doubled_rankers'] += 1

# Empty classes and constant features are independently specified zero-contribution
# cases, kept outside live parity to avoid RDKit's undefined zero denominators.
for classes in (1, 2, 5):
    labels = np.array([0, 0, 0, classes - 1])
    for matrix in (np.zeros((4, 3), dtype=np.uint8), np.ones((4, 3), dtype=np.uint8),
                   np.array([[0, 1, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0]], dtype=np.uint8)):
        exact = rational(matrix, labels, classes)
        for sparse in (False, True):
            with np.errstate(all='raise'):
                rows = indexed(rank(C, vectors(matrix, sparse), labels, 3, classes, C.InfoType.CHISQUARE))
            np.testing.assert_allclose(rows[:, 1], exact, rtol=1e-12, atol=1e-12)
            assert np.isfinite(rows).all()
            report['degenerate_rankers'] += 1

report['oracle_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(OUT / 'oracle.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: len(value) if isinstance(value, list) else value
                  for key, value in report.items()}, sort_keys=True))
