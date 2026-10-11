//! RDKit's canonical atom ranking and the Kekulé structure it picks with it.
//!
//! A port of RDKit 2026.03.6 `new_canon.{h,cpp}` (`rankFragmentAtoms` with
//! every atom and bond in play, the ranking `Kekulize(mol, canonical=true)`
//! uses) and of `Kekulize.cpp`'s traversal (`kekulizeWorker`). When an
//! aromatic system has several Kekulé structures RDKit takes the one its
//! walk in canonical-rank order reaches first; algorithms that read bond
//! orders off that structure (MMFF94's aromaticity pass, whose fullerene
//! cages are decided by which hexagons hold three double bonds) need the
//! same one, not just a valid one.
//!
//! The ranking is exact for the invariants chematic carries (degree, element,
//! isotope, hydrogen count, charge, atom map, tetrahedral parity, ring
//! membership in RDKit's SSSR). Double-bond stereo enters only through
//! [`rdkit_canonical_atom_ranks_with_bond_stereo`], whose caller supplies
//! RDKit's `STEREOE`/`STEREOZ` labels; RDKit's ring-stereo atom property is
//! not modelled there. [`rdkit_rank_mol_atoms`] takes every invariant in
//! RDKit's own terms (including `_ringStereoAtoms` and cis/trans stereo
//! atoms) for callers that model RDKit's molecule themselves.
//!
//! Which atoms must take a double bond is read from a valid Kekulé structure
//! supplied by the caller (`chematic_core::kekulize`): every aromatic atom
//! that holds a double bond there is one of RDKit's double-bond candidates.

use smallvec::SmallVec;
use std::collections::VecDeque;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, KekuleResult, Molecule, implicit_hcount,
};

include!("rdkit_canon/ranker.rs");
include!("rdkit_canon/api.rs");

#[cfg(test)]
include!("rdkit_canon_tests.rs");
