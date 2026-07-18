#!/usr/bin/env python3
"""
Generate the complete RDKit isotope-mass-delta table used by
`crates/chematic-fp/src/ecfp.rs`'s `EcfpInvariantMode::RdkitMorgan` atom
invariant, from RDKit's own `PeriodicTable` -- not a hand-picked subset of
"common organic elements", which left a real gap (e.g. carbon-11 truncated
to a different delta than RDKit's own value, since it fell back to the
`mass_number as f64` approximation).

Formula, matching RDKit's actual `getConnectivityInvariants` C++ source:
    deltaMass = static_cast<int>(atom->getMass() - getAtomicWeight(Z))
i.e. the atom's mass (see below) minus the element's average (standard)
atomic weight, truncated toward zero.

`atom->getMass()` itself has two cases, both reproduced here:
  - A real, RDKit-recognized isotope: `GetMassForIsotope(Z, A)`'s exact
    physical mass (verified: `GetMassForIsotope` returns exactly `0.0` for
    any (Z, A) RDKit doesn't recognize -- used as the "skip" sentinel when
    building RDKIT_ISOTOPE_DELTA_TABLE below, not a real 0 u isotope).
  - No isotope specified, OR an explicit isotope RDKit does NOT recognize:
    confirmed directly via `Atom.GetMass()` on parsed molecules (not
    inferred) -- unspecified returns the element's own average atomic
    weight (`C` -> `12.011`, not the monoisotopic `12.000`), and an
    unrecognized explicit isotope (e.g. `[500CH4]`) returns the raw mass
    number itself (`500.0`), not 0 and not a nearby real isotope's mass.
    RDKIT_ATOMIC_WEIGHTS below covers this case: chematic-fp's own
    rdkit_isotope_delta() computes `mass_number as f64 - average` directly
    for isotopes missing from RDKIT_ISOTOPE_DELTA_TABLE, using this table
    for `average` (NOT chematic-core's `Element::atomic_mass()`, which is a
    different, monoisotopic quantity -- using it here silently reproduces
    the exact bug this table exists to close).

Usage:
    .venv/bin/python scripts/gen_rdkit_isotope_delta_table.py > \\
        crates/chematic-fp/src/rdkit_isotope_delta_table.rs
"""

import sys

from rdkit import Chem, rdBase

MAX_MASS_NUMBER = 511  # covers every isotope RDKit's PeriodicTable knows


def main():
    pt = Chem.GetPeriodicTable()
    max_z = pt.GetMaxAtomicNumber()

    rows = []
    weights = []
    for atomic_number in range(1, max_z + 1):
        average = pt.GetAtomicWeight(atomic_number)
        weights.append(average)
        for mass_number in range(1, MAX_MASS_NUMBER + 1):
            mass = pt.GetMassForIsotope(atomic_number, mass_number)
            if mass <= 0.0:
                continue
            delta = int(mass - average)  # truncates toward zero, matching C++ static_cast<int>
            rows.append((atomic_number, mass_number, delta))

    print("//! RDKit isotope-mass-delta table for `EcfpInvariantMode::RdkitMorgan`.")
    print("//!")
    print("//! Generated from RDKit's own `PeriodicTable`, not a hand-picked subset --")
    print("//! covers every isotope RDKit's `GetMassForIsotope` recognizes, for every")
    print("//! element, plus every element's own average atomic weight (needed as the")
    print("//! subtrahend for isotopes RDKit itself doesn't recognize -- see")
    print("//! rdkit_isotope_delta() in ecfp.rs).")
    print("//!")
    print("//! Formula: `int(GetMassForIsotope(Z, A) - GetAtomicWeight(Z))`")
    print("//! (RDKit's `getConnectivityInvariants` deltaMass computation, truncated")
    print("//! toward zero, matching Rust's `f64 as i32` cast semantics.)")
    print("//!")
    print("//! DO NOT EDIT MANUALLY.")
    print(f"//! Generator: scripts/gen_rdkit_isotope_delta_table.py, RDKit {rdBase.rdkitVersion}")
    print()
    print("/// `(atomic_number, mass_number, delta)`, sorted by `(atomic_number, mass_number)`")
    print("/// for binary search.")
    print("#[allow(clippy::large_const_arrays)]")
    print(f"pub(crate) const RDKIT_ISOTOPE_DELTA_TABLE: [(u8, u16, i16); {len(rows)}] = [")
    for z, a, d in rows:
        print(f"    ({z}, {a}, {d}),")
    print("];")
    print()
    print("/// `GetAtomicWeight(atomic_number)`, indexed by atomic number directly")
    print("/// (index 0 is an unused placeholder -- valid atomic numbers start at 1).")
    print(f"pub(crate) const RDKIT_ATOMIC_WEIGHTS: [f64; {max_z + 1}] = [")
    print("    0.0, // index 0 unused")
    for w in weights:
        print(f"    {w},")
    print("];")

    print(f"-- generated {len(rows)} isotope rows, {len(weights)} atomic weights", file=sys.stderr)


if __name__ == "__main__":
    main()
