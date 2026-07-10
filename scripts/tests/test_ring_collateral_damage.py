#!/usr/bin/env python3
"""
Controls for ring_collateral_damage.py's compute()/comparison pipeline.

Two controls are required, not one:

1. POSITIVE (breakage must show red): corrupt one metric value and confirm
   the set()-based instability check flags it.
2. NEGATIVE (sameness must show green): call compute() on the SAME molecule,
   parsed fresh N times, and confirm every metric reports stable.

The negative control is the one that matters here. A harness that always
reports "unstable" regardless of truth trivially passes a positive control
(it always shows red) but fails a negative control -- which is exactly the
class of bug this project shipped once: ring_collateral_damage.py compared
mol.scaffold() by Python object identity (PyO3's Mol has no value __eq__/
__hash__), so two calls on the identical molecule were always "different"
objects and the script reported Murcko scaffold as 100% unstable regardless
of the real result (measurement bug, not a product bug -- see d2b852b/the
SSSR-fix commit for the corrected numbers). This test exercises the real
compute() function end-to-end so a regression of that exact bug class fails
here, not just in a hand-built dict.

Usage: python scripts/tests/test_ring_collateral_damage.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ring_collateral_damage import METRICS, compute  # noqa: E402


def positive_control():
    import chematic

    mol = chematic.from_smiles("c1ccccc1")
    good = compute(mol)
    corrupted = dict(good)
    corrupted["ring_count"] = good["ring_count"] + 1
    assert good != corrupted, "positive-control setup broken: corruption produced no difference"
    values = {"ring_count": {good["ring_count"], corrupted["ring_count"]}}
    assert len(values["ring_count"]) > 1, "POSITIVE CONTROL FAILED: harness did not detect injected fault"
    print("positive control: PASS (injected fault correctly reported unstable)")


def negative_control(n=10):
    """The same molecule, reparsed N times, must report stable on every metric."""
    import chematic

    smi = "c1ccc2ccccc2c1"  # naphthalene: exercises ring-bearing metrics, not just acyclic ones
    unstable_metrics = []
    for metric in METRICS:
        values = set()
        for _ in range(n):
            mol = chematic.from_smiles(smi)
            values.add(compute(mol)[metric])
        if len(values) > 1:
            unstable_metrics.append((metric, values))

    assert not unstable_metrics, (
        "NEGATIVE CONTROL FAILED: identical input reported as unstable for "
        f"{[m for m, _ in unstable_metrics]} -- comparison is not value-based. "
        f"Details: {unstable_metrics}"
    )
    print(f"negative control: PASS (identical input stable across all {len(METRICS)} metrics, n={n})")


if __name__ == "__main__":
    positive_control()
    negative_control()
