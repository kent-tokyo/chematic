//! Canonical SMILES generation via the Morgan (extended connectivity) algorithm.
//!
//! A canonical SMILES is a unique string representation of a molecule:
//! two molecules that are identical (same graph, same atom properties) will
//! always produce the same canonical SMILES string.
//!
//! Algorithm:
//! 1. Assign initial invariants to each atom (atomic number, degree, charge, …).
//! 2. Iteratively update ranks using Morgan-style neighbor aggregation until
//!    the number of distinct ranks stabilises.
//! 3. Use the resulting ranks to impose a canonical DFS traversal order.
//!    Critically, both the ring-closure discovery DFS and the write DFS
//!    use the *same* canonical traversal order so the output is stable.
//! 4. Tie-breaking when two atoms have equal Morgan rank is resolved by
//!    atomic_number → isotope → charge → aromaticity → degree.
//!
//! Reference: Weininger, D. (1988) J. Chem. Inf. Comput. Sci. 28, 31-36.

use std::collections::{HashMap, HashSet};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL, STEREO_VACANT_SENTINEL,
    implicit_hcount, remap_square_planar_tag, remap_tetrahedral_parity, valence_inferred_hcount,
};
use smallvec::SmallVec;

use crate::parser::parse;
use crate::writer::{
    bond_token_from, emit_bracket_hydrogens, square_planar_token, suppress_standalone_wedge,
};

include!("canonical/core.rs");
include!("canonical/ez_planning.rs");
include!("canonical/writer.rs");
include!("canonical/stereo.rs");

#[cfg(test)]
include!("canonical_tests.rs");
