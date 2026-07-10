#!/usr/bin/env python3
"""
Controls for ringinfo_parity.py.

ringinfo_parity.py already has a positive control (self_test(), on hand-built
metric dicts). It has no negative control on the real pipeline: nothing
proves that calling chematic_ring_metrics() on the same molecule twice
reports "stable" rather than some artifact of how the values are typed and
compared (see test_ring_collateral_damage.py for why that gap matters --
that exact bug class shipped once in a sibling harness). This file adds that
negative control against the real function, plus re-runs the existing
positive control for completeness.

Usage: python scripts/tests/test_ringinfo_parity.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ringinfo_parity import chematic_ring_metrics, self_test  # noqa: E402


def negative_control(n=10):
    """The same molecule, reparsed N times, must report an identical metrics dict."""
    import chematic

    smi = "c1ccc2ccccc2c1"  # naphthalene: exercises every ring-bearing field
    results = []
    for _ in range(n):
        mol = chematic.from_smiles(smi)
        results.append(chematic_ring_metrics(mol))

    first = results[0]
    mismatched = [k for k in first if len(set(r[k] for r in results)) > 1]
    assert not mismatched, (
        f"NEGATIVE CONTROL FAILED: identical input reported as unstable for {mismatched} "
        "-- comparison is not value-based or chematic_ring_metrics() is non-deterministic."
    )
    print(f"negative control: PASS (identical input stable across all {len(first)} fields, n={n})")


if __name__ == "__main__":
    self_test()
    negative_control()
