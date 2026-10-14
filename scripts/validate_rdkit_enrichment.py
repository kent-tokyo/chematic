#!/usr/bin/env python3
"""Compare local ranked EF to integer counts and live RDKit on ordinary cutoffs."""
import argparse
import hashlib
import inspect
import json
import math
from pathlib import Path
import platform
import subprocess

import chematic
import rdkit
from rdkit.ML.Scoring.Scoring import CalcEnrichment

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference(labels, fractions):
    n, actives = len(labels), sum(labels)
    return [n * sum(labels[:math.ceil(n*f)]) / (math.ceil(n*f)*actives) for f in fractions]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    comparisons = 0
    max_error = 0.
    for n in range(2, 9):
        for mask in range(1, 1 << n):
            labels = [bool(mask & (1 << i)) for i in range(n)]
            for k in range(1, n):
                fraction = k/n
                if math.ceil(n*fraction) == n:
                    continue
                actual = chematic.enrichment_factors(labels, [fraction])
                oracle = CalcEnrichment([(0., active) for active in labels], 1, [fraction])
                expected = reference(labels, [fraction])
                assert len(oracle) == len(actual) == 1
                error = max(abs(actual[0]-oracle[0]), abs(actual[0]-expected[0]))
                assert error < 1e-12, (labels, fraction, actual, oracle, expected)
                max_error = max(error, max_error)
                comparisons += 1
    fixture_path = ROOT / "validation/rdkit_issues/screening/enrichment-7981.json"
    fixture = json.loads(fixture_path.read_text())
    boundaries = []
    for labels, fractions in [
        (fixture["ranked_actives"], [.9]),
        (fixture["ranked_actives"], [1.]),
        ([True], [1.]),
        ([True, False, True, False, False], [1., .2, .01, .2, .5, .9]),
    ]:
        actual = chematic.enrichment_factors(labels, fractions)
        expected = reference(labels, fractions)
        assert len(actual) == len(fractions)
        assert all(abs(a-b) < 1e-12 for a,b in zip(actual, expected))
        boundaries.append(dict(ranked_actives=labels, fractions=fractions,
                               chematic=actual, count_reference=expected,
                               rdkit=CalcEnrichment([(0., a) for a in labels], 1, fractions)))
    undefined = []
    for labels, fractions in [([], [1.]), ([True], []), ([False], [1.]), ([True], [0.]), ([True], [math.nan])]:
        try:
            chematic.enrichment_factors(labels, fractions)
        except ValueError as exc:
            undefined.append(dict(ranked_actives=labels, fractions=[repr(f) for f in fractions],
                                  error=str(exc)))
        else:
            raise AssertionError((labels, fractions))
    module = Path(chematic.__file__).resolve()
    report = dict(
        issue=fixture["issue"],
        artifact_lane="local source extension; uncommitted working tree; not a release",
        source_baseline=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        python=platform.python_version(), rdkit_version=rdkit.__version__,
        rdkit_calc_enrichment_sha256=hashlib.sha256(inspect.getsource(CalcEnrichment).encode()).hexdigest(),
        module_path=str(module), module_sha256=sha(module),
        source_sha256={str(path.relative_to(ROOT)): sha(path) for path in [
            ROOT / "crates/chematic-chem/src/screening.rs", ROOT / "crates/chematic-py/src/similarity.rs"]},
        fixture_sha256=sha(fixture_path),
        ordinary_single_cutoff_comparisons=comparisons, max_absolute_error=max_error,
        boundary_comparisons=boundaries, undefined_input_rejections=undefined,
        scope="ranked binary labels; ceil(N*f); no sorting, ties, weights, ROC, RIE or BEDROC",
    )
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(comparisons=comparisons, max_absolute_error=max_error,
                          boundary_cases=len(boundaries), undefined_cases=len(undefined))))


if __name__ == "__main__":
    main()
