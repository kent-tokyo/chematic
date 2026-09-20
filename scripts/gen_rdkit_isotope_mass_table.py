#!/usr/bin/env python3
"""Generate the full RDKit isotope-mass table for descriptor compatibility.

``Descriptors.MolWt`` uses ``PeriodicTable.GetMassForIsotope`` for every
recognized explicit isotope.  Keeping a short hand-written list makes an
otherwise valid labelled atom silently fall back to its integer mass number.
This generator records every recognized RDKit isotope and leaves that integer
fallback to the Rust lookup for labels RDKit does not recognize.

Usage:
    python scripts/gen_rdkit_isotope_mass_table.py > \
        crates/chematic-chem/src/rdkit_isotope_mass_table.rs
"""

from rdkit import Chem, rdBase


MAX_MASS_NUMBER = 511


def rust_float(value: float) -> str:
    """Render an unambiguous Rust `f64` literal without losing precision."""
    rendered = format(value, ".12g")
    return rendered if any(marker in rendered for marker in (".", "e", "E")) else f"{rendered}.0"


def main() -> None:
    periodic_table = Chem.GetPeriodicTable()
    rows: list[tuple[int, int, float]] = []
    for atomic_number in range(1, periodic_table.GetMaxAtomicNumber() + 1):
        for mass_number in range(1, MAX_MASS_NUMBER + 1):
            mass = periodic_table.GetMassForIsotope(atomic_number, mass_number)
            if mass > 0.0:
                rows.append((atomic_number, mass_number, mass))

    print("//! RDKit explicit-isotope mass table for `rdkit_molecular_weight`.")
    print("//!")
    print("//! Generated from RDKit's `PeriodicTable.GetMassForIsotope` for every")
    print("//! recognized isotope. Unknown syntactically valid labels deliberately")
    print("//! retain RDKit's integer mass-number fallback in the caller.")
    print("//!")
    print("//! DO NOT EDIT MANUALLY.")
    print(f"//! Generator: scripts/gen_rdkit_isotope_mass_table.py, RDKit {rdBase.rdkitVersion}")
    print()
    print("/// `(atomic_number, mass_number, nuclide_mass_da)`, sorted for binary search.")
    print("#[allow(clippy::large_const_arrays)]")
    print(f"pub(crate) const RDKIT_ISOTOPE_MASS_TABLE: [(u8, u16, f64); {len(rows)}] = [")
    for atomic_number, mass_number, mass in rows:
        print(f"    ({atomic_number}, {mass_number}, {rust_float(mass)}),")
    print("];" )


if __name__ == "__main__":
    main()
