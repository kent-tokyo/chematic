//! Smallest Set of Smallest Rings (SSSR) via Horton's algorithm.
//!
//! Algorithm overview:
//! 1. Compute the cycle rank r = E - V + C (Euler characteristic),
//!    where C is the number of connected components.
//! 2. For every vertex v (as a candidate root) and every ring-eligible edge
//!    (x, y), form the candidate cycle SP(v, x) + edge(x, y) + SP(y, v),
//!    where SP is the shortest path in v's BFS tree. Keep it only if it's a
//!    genuine simple cycle (the two paths share no vertex other than v).
//!    This produces O(V*E) candidates and is guaranteed (Horton, 1987) to
//!    contain a minimum-weight cycle basis — unlike a single spanning tree's
//!    fundamental-cycle set (exactly r candidates, no redundancy), which can
//!    only ever report *a* valid basis, never guaranteed to be minimal.
//! 3. Represent each cycle as a set of bond indices (for GF(2) XOR independence
//!    testing) and sort candidates by (length, canonical tie-break) — the
//!    tie-break uses a local Weisfeiler-Leman-style atom ranking (see
//!    `canonical_atom_ranks`) so ring *selection* doesn't depend on input
//!    atom-numbering (i.e. SMILES parse/traversal order), only on molecular
//!    graph structure.
//! 4. Use Gaussian elimination over GF(2) to greedily build an independent
//!    basis of r cycles from that sorted candidate list.
//! 5. Convert the chosen bond-sets back to ordered atom sequences for the public API.

use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::collections::VecDeque;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder};

include!("sssr/basis.rs");
include!("sssr/candidates.rs");
include!("sssr/symmetrized.rs");
include!("sssr/canonical.rs");

#[cfg(test)]
include!("sssr_tests.rs");
