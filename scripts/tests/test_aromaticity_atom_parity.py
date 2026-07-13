#!/usr/bin/env python3
"""
Controls for aromaticity_atom_parity.py's per_atom_arom() and its
worst-of-N comparison approach.

aromaticity_atom_parity.py has no extracted self_test(); this file adds
both controls directly against the real per_atom_arom() function.

Usage: python scripts/tests/test_aromaticity_atom_parity.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from aromaticity_atom_parity import per_atom_arom  # noqa: E402


def positive_control():
    """A kekulized (non-aromatic-flagged) parse must differ from the aromatic-flagged parse."""
    import chematic

    aromatic = per_atom_arom(chematic.from_smiles("c1ccccc1"))
    kekulized = per_atom_arom(chematic.from_smiles("C1=CC=CC=C1"))
    assert aromatic != kekulized, (
        "POSITIVE CONTROL FAILED: aromatic-flagged and Kekulized benzene compared equal "
        "-- per_atom_arom() is not sensitive to the aromatic flag"
    )
    print("positive control: PASS (aromatic-flag difference correctly detected)")


def negative_control(n=10):
    """The same molecule, reparsed N times, must report the identical per-atom flag list."""
    import chematic

    smi = "c1cnc2[nH]cnc2n1"  # purine: mixed aromatic heteroatoms, exercises real flag logic
    results = [tuple(per_atom_arom(chematic.from_smiles(smi))) for _ in range(n)]
    assert len(set(results)) == 1, (
        "NEGATIVE CONTROL FAILED: identical input produced different per-atom aromaticity "
        f"flags across {n} reparses -- comparison is not value-based or parsing is "
        f"non-deterministic. Distinct results: {set(results)}"
    )
    print(f"negative control: PASS (identical input stable across {n} reparses)")


if __name__ == "__main__":
    positive_control()
    negative_control()
