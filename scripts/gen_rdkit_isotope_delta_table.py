#!/usr/bin/env python3
"""
Generate the complete RDKit isotope-mass-delta table used by
`crates/chematic-fp/src/ecfp.rs`'s `EcfpInvariantMode::RdkitMorgan` atom
invariant, from RDKit's own `PeriodicTable` -- not a hand-picked subset of
"common organic elements", which left a real gap (e.g. carbon-11 truncated
to a different delta than RDKit's own value, since it fell back to the
`mass_number as f64` approximation).

Formula, matching RDKit's actual `getConnectivityInvariants` C++ source:
    deltaMass = static_cast<int>(getMassForIsotope(Z, A) - getAtomicWeight(Z))
i.e. the isotope's exact mass minus the element's average (standard) atomic
weight, truncated toward zero. `GetMassForIsotope` returns exactly `0.0` for
any (Z, A) that isn't a real, RDKit-known isotope -- used here as the
"skip" sentinel, not a real 0 u isotope.

Usage:
    .venv/bin/python scripts/gen_rdkit_isotope_delta_table.py > \\
        crates/chematic-fp/src/rdkit_isotope_delta_table.rs
"""

from rdkit import Chem, rdBase

MAX_MASS_NUMBER = 511  # covers every isotope RDKit's PeriodicTable knows


def main():
    pt = Chem.GetPeriodicTable()
    rows = []
    for atomic_number in range(1, pt.GetMaxAtomicNumber() + 1):
        average = pt.GetAtomicWeight(atomic_number)
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
    print("//! element.")
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

    import sys

    print(f"-- generated {len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
