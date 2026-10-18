"""RDKit co-occurrence matrices and independent integer pair-membership sums."""

import hashlib
import json
from pathlib import Path
import random

import numpy as np
import rdkit
from rdkit import DataStructs as D
from rdkit.ML.InfoTheory import rdInfoTheory as R

from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as C

HERE = Path(__file__).parent
rng = random.Random(960101)
inputs = hashlib.sha256()
comparisons = pair_values = prefix_checks = metamorphic_checks = paired_votes = 0


class Lazy:
    def __init__(self, width, bits):
        self.width, self.bits = width, bits

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return (np.uint32(bit) for bit in reversed(sorted(self.bits)))


def vector(cls, width, bits):
    if cls is Lazy:
        return Lazy(width, bits)
    result = cls(width)
    for bit in bits:
        result.SetBit(bit)
    return result


def exact(selected, votes):
    return [
        sum(selected[i] in vote and selected[j] in vote for vote in votes)
        for i in range(len(selected))
        for j in range(i)
    ]


for case in range(128):
    width = 16 + case % 49
    selected = rng.sample(range(width), 2 + case % 15)
    votes = [
        set(rng.sample(range(width), rng.randrange(width // 2 + 1)))
        for _ in range(8 + case % 33)
    ]
    votes[0] = set()
    inputs.update(
        json.dumps([width, selected, [sorted(vote) for vote in votes]]).encode()
    )
    expected = exact(selected, votes)
    reference = R.BitCorrMatGenerator()
    reference.SetBitList(selected)
    sparse_reference = R.BitCorrMatGenerator()
    sparse_reference.SetBitList(selected)
    for vote in votes:
        reference.CollectVotes(vector(D.ExplicitBitVect, width, vote))
        sparse_reference.CollectVotes(vector(D.SparseBitVect, width, vote))
    np.testing.assert_array_equal(reference.GetCorrMatrix(), expected)
    np.testing.assert_array_equal(sparse_reference.GetCorrMatrix(), expected)
    for cls in [D.ExplicitBitVect, D.SparseBitVect, Lazy, ExplicitBitVect]:
        generator = C.BitCorrMatGenerator()
        generator.SetBitList((np.int64(bit) for bit in selected))
        for index, vote in enumerate(votes):
            generator.CollectVotes(vector(cls, width, vote))
            if index in {0, len(votes) // 2, len(votes) - 1}:
                np.testing.assert_array_equal(
                    generator.GetCorrMatrix(), exact(selected, votes[: index + 1])
                )
                prefix_checks += 1
        actual = generator.GetCorrMatrix()
        assert actual.dtype == np.float64 and actual.flags.owndata
        np.testing.assert_array_equal(actual, reference.GetCorrMatrix())
        comparisons += 1
        pair_values += len(actual)
        paired_votes += len(votes)
        for transformed, factor in [(list(reversed(votes)), 1), (votes + votes, 2)]:
            generator.SetBitList(selected)
            assert np.all(generator.GetCorrMatrix() == 0)
            for vote in transformed:
                generator.CollectVotes(vector(cls, width, vote))
            np.testing.assert_array_equal(
                generator.GetCorrMatrix(), np.array(expected) * factor
            )
            metamorphic_checks += 1
        # Reset clears votes and permits a new width, as confirmed upstream.
        generator.SetBitList(list(reversed(selected)))
        reference.SetBitList(list(reversed(selected)))
        assert np.all(reference.GetCorrMatrix() == 0)
        generator.CollectVotes(vector(cls, width + 3, votes[-1]))
        reference.CollectVotes(vector(D.ExplicitBitVect, width + 3, votes[-1]))
        np.testing.assert_array_equal(
            generator.GetCorrMatrix(), reference.GetCorrMatrix()
        )
        metamorphic_checks += 1
        # Restore the shared reference for the next representation.
        reference.SetBitList(selected)
        for vote in votes:
            reference.CollectVotes(vector(D.ExplicitBitVect, width, vote))

diagnostics = []
c = C.BitCorrMatGenerator()
c.SetBitList([1, 3])
r = R.BitCorrMatGenerator()
r.SetBitList([1, 3])
for width in [8, 16]:
    fp = vector(D.ExplicitBitVect, width, {1, 3})
    r.CollectVotes(fp)
    try:
        c.CollectVotes(fp)
    except ValueError:
        diagnostics.append(
            dict(
                case="width_change_without_reset",
                rdkit=r.GetCorrMatrix().tolist(),
                chematic=c.GetCorrMatrix().tolist(),
                chematic_error="ValueError",
            )
        )
high = vector(D.SparseBitVect, 2**32 - 1, {0, 2**32 - 2})
view = list(high.GetOnBits())
assert set(view) == {0, -2}
c.SetBitList([0, 2**32 - 2])
try:
    c.CollectVotes(high)
except ValueError:
    diagnostics.append(
        dict(
            case="high_signed_sparse_view",
            view=view,
            chematic_error="ValueError",
            counts=c.GetCorrMatrix().tolist(),
        )
    )
else:
    raise AssertionError("negative signed ID must not be reinterpreted")
assert len(diagnostics) == 2
summary = dict(
    rdkit=rdkit.__version__,
    numpy=np.__version__,
    seed=960101,
    datasets=128,
    final_comparisons=comparisons,
    pair_values=pair_values,
    paired_votes=paired_votes,
    prefix_checks=prefix_checks,
    metamorphic_checks=metamorphic_checks,
    max_count_error=0,
    diagnostics=len(diagnostics),
    inputs_sha256=inputs.hexdigest(),
    scope="Exact selected-pair counts and lifecycle; source-profile differences separate from ordinary parity",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(HERE / "reference-diagnostics.json").write_text(
    json.dumps(diagnostics, indent=2) + "\n"
)
print(json.dumps(summary, indent=2))
