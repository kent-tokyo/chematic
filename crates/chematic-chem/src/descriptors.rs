//! Molecular descriptor functions for drug-likeness and physical property estimation.
//!
//! All functions accept a `&Molecule` reference.  Molecules with aromatic bonds
//! (SMILES lowercase notation) are kekulized internally where hydrogen counts
//! are required; the caller's molecule is never mutated.

use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::OnceLock;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule, implicit_hcount};
use chematic_perception::{all_ring_list, aromatic_ring_list, find_ring_families, find_sssr};
use chematic_smarts::parse_smarts;

include!("descriptors/mass_and_basic.rs");
include!("descriptors/tpsa.rs");
include!("descriptors/crippen.rs");
include!("descriptors/rings_and_filters.rs");
include!("descriptors/mqn.rs");
include!("descriptors/topology_and_shape.rs");
include!("descriptors/composition.rs");
include!("descriptors/information_and_bcut.rs");

#[cfg(test)]
include!("descriptors_tests.rs");

#[cfg(test)]
#[path = "descriptors/rotor_order_tests.rs"]
mod rotor_order_tests;
