//! Hückel aromaticity perception with antiaromaticity detection.
//!
//! Works on kekulized molecules (no `Aromatic` bond orders) **or** on molecules
//! that retain `Aromatic` bond orders from the SMILES parser (pre-kekulization).
//! Call `kekulize` + `apply_kekule` from `chematic-core` before calling
//! `assign_aromaticity` if you need the explicit double-bond form.
//!
//! Algorithm:
//! 1. Find all SSSR rings via `find_sssr`.
//! 2. **Pass 1**: evaluate each ring independently using Hückel electron counting.
//!    Aromatic (`BondOrder::Aromatic`) bonds are treated equivalently to double bonds
//!    so that pre-kekulization input is handled correctly.
//!    A special "bridgehead N" rule covers fused-ring N atoms whose entire valence
//!    is satisfied by single σ-bonds (like indolizine's junction nitrogen).
//! 3. **Pass 2**: iterative propagation. Rings that were `NonAromatic` or
//!    indeterminate in Pass 1 are re-evaluated using the already-aromatic atom set
//!    as context: confirmed-aromatic atoms contribute 1π unconditionally, allowing
//!    fused rings to be recognised bottom-up (e.g. the 6-ring of indolizine).
//! 4. Classify rings by electron count:
//!    - 4n+2 electrons (n >= 0): aromatic (favorable)
//!    - 4n electrons (n > 0): antiaromatic (unfavorable, strongly disfavored)
//!    - Other: non-aromatic
//! 5. Record all aromatic atoms, bonds, and antiaromatic rings in an `AromaticityModel`.

include!("aromaticity/model.rs");
include!("aromaticity/rings.rs");
include!("aromaticity/trace.rs");

#[cfg(test)]
include!("aromaticity_tests.rs");
