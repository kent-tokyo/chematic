"""Independent RDKit entropy-ranking checks on original synthetic datasets."""

import hashlib
import json
from pathlib import Path

import numpy as np
import rdkit
from rdkit import DataStructs
from rdkit.ML.InfoTheory import rdInfoTheory as reference
import chematic
from chematic.rdkit_compat import rdInfoTheory as candidate


class Sparse:
    def __init__(self, width, on):
        self.width, self.on = width, on

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return iter(self.on)


def by_id(rows):
    return rows[np.argsort(rows[:, 0])]


def main():
    rng = np.random.default_rng(9601)
    digest = hashlib.sha256()
    datasets = votes = metamorphic = 0
    maximum_difference = 0.0
    tie_order_differences = 0
    for width in (1, 4, 16, 64):
        for classes in (2, 3, 5):
            for size in (1, 3, 11, 40):
                for seed in range(8):
                    bits = rng.random((size, width)) < (0.1 + seed * 0.1)
                    labels = rng.integers(0, classes, size)
                    digest.update(bits.tobytes() + labels.tobytes())
                    for cls in (DataStructs.ExplicitBitVect, DataStructs.SparseBitVect):
                        actual = candidate.InfoBitRanker(width, classes)
                        expected = reference.InfoBitRanker(
                            width, classes, reference.InfoType.ENTROPY
                        )
                        for row, label in zip(bits, labels):
                            vector = cls(width)
                            for bit in np.flatnonzero(row):
                                vector.SetBit(int(bit))
                            actual.AccumulateVotes(vector, int(label))
                            expected.AccumulateVotes(vector, int(label))
                            votes += 1
                        a, b = actual.GetTopN(width), expected.GetTopN(width)
                        ordered_a, ordered_b = by_id(a), by_id(b)
                        np.testing.assert_allclose(
                            a[:, 1], b[:, 1], rtol=1e-11, atol=1e-12
                        )
                        np.testing.assert_array_equal(ordered_a[:, 0], np.arange(width))
                        np.testing.assert_array_equal(
                            ordered_a[:, 2:], ordered_b[:, 2:]
                        )
                        np.testing.assert_allclose(
                            ordered_a[:, 1], ordered_b[:, 1], rtol=1e-11, atol=1e-12
                        )
                        maximum_difference = max(
                            maximum_difference,
                            float(np.abs(ordered_a[:, 1] - ordered_b[:, 1]).max()),
                        )
                        tie_order_differences += int(
                            not np.array_equal(a[:, 0], b[:, 0])
                        )
                        np.testing.assert_array_equal(
                            actual.GetTopN(min(5, width)), a[: min(5, width)]
                        )
                        assert np.all(a[:-1, 1] >= a[1:, 1])
                        reverse = candidate.InfoBitRanker(width, classes)
                        doubled = candidate.InfoBitRanker(width, classes)
                        for row, label in zip(bits[::-1], labels[::-1]):
                            vector = Sparse(width, np.flatnonzero(row)[::-1].tolist())
                            reverse.AccumulateVotes(vector, int(label))
                            doubled.AccumulateVotes(vector, int(label))
                            doubled.AccumulateVotes(vector, int(label))
                        np.testing.assert_array_equal(reverse.GetTopN(width), a)
                        d = by_id(doubled.GetTopN(width))
                        np.testing.assert_allclose(d[:, 1], ordered_a[:, 1], atol=1e-14)
                        np.testing.assert_array_equal(d[:, 2:], 2 * ordered_a[:, 2:])
                        metamorphic += 2
                        datasets += 1
    print(
        json.dumps(
            {
                "scope": "ENTROPY information gain and per-class on-counts; ties use candidate bit-ID order",
                "rdkit_version": rdkit.__version__,
                "numpy_version": np.__version__,
                "ranker_comparisons": datasets,
                "paired_votes": votes,
                "metamorphic_checks": metamorphic,
                "maximum_absolute_information_gain_difference": maximum_difference,
                "rank_order_difference_diagnostic": tie_order_differences,
                "corpus_sha256": digest.hexdigest(),
                "native_fixture_sha256": hashlib.sha256(
                    Path(chematic.chematic.__file__).read_bytes()
                ).hexdigest(),
                "native_fixture_source_commit": "a56bc1353bd70eb5381dce38f818c0b78755d027",
                "native_fixture_scope": "reused import fixture; new implementation calls Python/NumPy only; no native rebuild",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
