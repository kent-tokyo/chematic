//! Tautomer enumeration and canonical tautomer selection.
//!
//! Provides:
//! - [`canonical_tautomer`]: return the canonical (preferred) tautomer form.
//! - [`enumerate_tautomers`]: enumerate all reachable tautomers up to a cap.
//! - [`tautomer_parent`]: like `canonical_tautomer`, but budget-limited via
//!   [`TautomerLimits`] and returns an explainable [`TautomerAuditRecord`]
//!   instead of a bare `Molecule` -- see
//!   `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` section 4.
//!
//! Pattern matching is implemented directly on the molecule graph (no SMARTS).

#![forbid(unsafe_code)]

use std::collections::HashSet;
use std::time::Instant;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, implicit_hcount,
};

use crate::parent::{InvalidInputReason, ParentAudit, ParentComputationStatus, ParentResult};
use crate::standardize::MoleculeSnapshot;

include!("tautomer/rules.rs");
include!("tautomer/aromatic.rs");
include!("tautomer/transforms.rs");
include!("tautomer/search.rs");
include!("tautomer/audit.rs");

#[cfg(test)]
include!("tautomer_tests.rs");
