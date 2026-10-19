#!/usr/bin/env python3
"""Compare InfoTheory file formatting and saved-ranking lifecycle with RDKit."""

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
from rdkit import DataStructs, rdBase
from rdkit.ML.InfoTheory import rdInfoTheory as reference

from chematic.rdkit_compat import rdInfoTheory as candidate
import chematic._info_theory as implementation


def expected_table(rows, classes):
    header = f"{'Bit':>12}{'InfoContent':>12}"
    header += "".join(f"{'class':>10}{i}" for i in range(classes))
    lines = [header + "\n"]
    for row in rows:
        text = f"{int(row[0]):>12}{float(row[1]):>12.5g}"
        text += "".join(f"{int(v):>10}" for v in row[2:])
        lines.append(text + "\n")
    return "".join(lines).encode("ascii")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(9601)
    reports = []
    with tempfile.TemporaryDirectory(prefix="chematic-ranker-export-") as directory:
        left_path = Path(directory) / "candidate.txt"
        right_path = Path(directory) / "reference.txt"
        for case in range(64):
            width = int(rng.integers(2, 25))
            classes = int(rng.integers(2, 5))
            records = int(rng.integers(classes * 4, classes * 16))
            labels = np.arange(records) % classes
            matrix = rng.integers(0, 2, size=(records, width))
            for mode in range(1, 5):
                for sparse in (False, True):
                    left = candidate.InfoBitRanker(width, classes, mode)
                    right = reference.InfoBitRanker(
                        width, classes, reference.InfoType.values[mode]
                    )
                    if mode in (2, 4):
                        left.SetBiasList([0])
                        right.SetBiasList([0])
                    for values, label in zip(matrix, labels):
                        vector = (
                            DataStructs.SparseBitVect(width)
                            if sparse
                            else DataStructs.ExplicitBitVect(width)
                        )
                        for bit in np.flatnonzero(values):
                            vector.SetBit(int(bit))
                        left.AccumulateVotes(vector, int(label))
                        right.AccumulateVotes(vector, int(label))
                    left_rows = left.GetTopN(width)
                    right_rows = right.GetTopN(width)
                    left.WriteTopBitsToFile(left_path)
                    right.WriteTopBitsToFile(str(right_path))
                    left_bytes = left_path.read_bytes()
                    right_bytes = right_path.read_bytes()
                    assert left_bytes == expected_table(left_rows, classes)
                    assert right_bytes == expected_table(right_rows, classes)
                    a = left_rows[np.argsort(left_rows[:, 0])]
                    b = right_rows[np.argsort(right_rows[:, 0])]
                    np.testing.assert_array_equal(a[:, 0], np.arange(width))
                    comparable = np.array_equal(a[:, 0], b[:, 0])
                    score_fields_differ = None
                    max_score_difference = None
                    if comparable:
                        np.testing.assert_array_equal(a[:, 2:], b[:, 2:])
                        np.testing.assert_allclose(
                            a[:, 1], b[:, 1], rtol=1e-12, atol=1e-12
                        )
                        score_fields_differ = sum(
                            format(x, ".5g") != format(y, ".5g")
                            for x, y in zip(a[:, 1], b[:, 1])
                        )
                        max_score_difference = float(np.max(np.abs(a[:, 1] - b[:, 1])))
                    # Both runtimes keep the last ranking even after the returned
                    # array is overwritten and an additional vote is counted.
                    left_rows[:] = -123
                    right_rows[:] = -123
                    empty = DataStructs.ExplicitBitVect(width)
                    left.AccumulateVotes(empty, 0)
                    right.AccumulateVotes(empty, 0)
                    left.WriteTopBitsToFile(left_path)
                    right.WriteTopBitsToFile(str(right_path))
                    assert left_path.read_bytes() == left_bytes
                    assert right_path.read_bytes() == right_bytes
                    reports.append(
                        {
                            "case": case,
                            "mode": mode,
                            "sparse": sparse,
                            "bits": width,
                            "classes": classes,
                            "records": records,
                            "comparison_status": (
                                "comparable"
                                if comparable
                                else "reference_bit_ids_not_permutation"
                            ),
                            "reference_rows": None if comparable else b.tolist(),
                            "identical_file": left_bytes == right_bytes,
                            "score_fields_differ": score_fields_differ,
                            "max_score_difference": max_score_difference,
                            "candidate_sha256": hashlib.sha256(left_bytes).hexdigest(),
                            "reference_sha256": hashlib.sha256(right_bytes).hexdigest(),
                        }
                    )
    result = {
        "scope": "Python file-format and ranking-snapshot comparison, not native wheel integration",
        "reference_version": rdBase.rdkitVersion,
        "seed": 9601,
        "implementation_sha256": hashlib.sha256(
            Path(implementation.__file__).read_bytes()
        ).hexdigest(),
        "comparisons": len(reports),
        "format_checks": len(reports) * 2,
        "snapshot_checks": len(reports) * 2,
        "identical_files": sum(row["identical_file"] for row in reports),
        "comparable_rankings": sum(
            row["comparison_status"] == "comparable" for row in reports
        ),
        "reference_nonpermutations": sum(
            row["comparison_status"] != "comparable" for row in reports
        ),
        "score_fields_differ": sum(row["score_fields_differ"] or 0 for row in reports),
        "maximum_score_difference": max(
            row["max_score_difference"] or 0 for row in reports
        ),
        "cases": reports,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}))


if __name__ == "__main__":
    main()
