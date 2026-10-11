//! Aromaticity-A1-1b-0: a faithful, independent reproduction of RDKit's
//! *default* aromaticity model (`AROMATICITY_RDKIT`/`AROMATICITY_DEFAULT`),
//! ported directly from RDKit's own source
//! (`Code/GraphMol/Aromaticity.cpp`, functions `getAtomDonorTypeArom`,
//! `countAtomElec`, `isAtomCandForArom`, `applyHuckel`, `applyHuckelToFused`,
//! `aromaticityHelper`'s `includeFused` branch — the exact path
//! `setAromaticity(mol, AROMATICITY_RDKIT, ...)` calls).
//!
//! `AromaticityAlgorithm::RdkitLike` now routes through this engine when the
//! input can be kekulized; explicit, self-consistent aromatic input is preserved
//! unless re-perception strictly extends it without removing any supplied
//! aromatic atom or bond. The historical per-ring
//! implementation remains an infallible fallback for inputs that cannot satisfy
//! either representation. `Huckel` remains the default and is unchanged. This module
//! also backs a separate, explicitly opt-in, fallible production API:
//! [`assign_aromaticity_rdkit_parity_experimental`] and
//! [`apply_aromaticity_rdkit_parity_experimental`], re-exported from the
//! crate root. Every other item in this module (the low-level donor-type/
//! Hückel machinery) is crate-private; the only way to reach it from
//! outside the crate is through those two functions, or, for diagnostics,
//! `diagnostics::rdkit_parity_aromaticity` behind the `diagnostics` feature.
//! See `docs/rfcs/aromaticity_a1_rfc.md`'s "A1-1b-0"/"A1-1b-1" sections for the
//! full design writeup, the calibration battery, and the corpus gate.
//!
//! Ported from RDKit commit `e89c9f656a694fab4105139844cba88d2e013354`, an
//! ancestor of release tag `Release_2026_03_4` (which resolves to
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`, independently verified via
//! the GitHub tags API during Morgan M4-A0). `Code/GraphMol/Aromaticity.cpp`
//! is byte-identical between the two commits (independently diffed during
//! M4-A0's provenance audit — the 130 commits between them never touch this
//! file), so the 5 functions cited above are unaffected either way. See
//! `THIRD_PARTY_NOTICES.md` at the repo root for the required BSD 3-Clause
//! attribution and license text.
//!
//! Unlike this crate's own `ring_pi_electrons`/`evaluate_atom_pi_contribution`
//! (which evaluate an atom's contribution *per candidate ring/component*),
//! RDKit computes each atom's [`ElectronDonorType`] **once, globally, per
//! molecule** — whether a multiple bond "counts" for aromaticity purposes
//! depends on whether that bond is part of *any* SSSR ring in the whole
//! molecule (`RingInfo::numBondRings(bond) > 0`), not on whether it's inside
//! the *specific* candidate ring currently being evaluated. This is the
//! precise, source-verified point where this crate's own `ring_pi_electrons`
//! diverges from RDKit for the SMARTS-A0/PR #86 false-positive family: its
//! `CarbonExocyclicHeteroatomDouble` rule checks "is the double-bond partner
//! outside *this ring's* atom set" where RDKit checks "is this bond outside
//! *every* ring in the molecule" — an exocyclic-to-the-candidate-ring double
//! bond whose partner is itself a *different* ring's atom (e.g. this crate's
//! reproducer's atom 8, `C=N` where the N is in a second fused ring) still
//! counts as a normal one-electron donor under RDKit's rule, not a
//! zero-electron "spent on the exocyclic bond" donor.
//!
//! Prefers pre-kekulized input (no `BondOrder::Aromatic`), matching RDKit's own
//! pipeline (`Kekulize` normally runs before `setAromaticity`), but retains a
//! validated explicit-aromatic representation when that conversion is impossible.

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::aromaticity::AromaticityModel;

include!("rdkit_parity/typing.rs");
include!("rdkit_parity/view.rs");

#[cfg(test)]
include!("rdkit_parity_tests.rs");
