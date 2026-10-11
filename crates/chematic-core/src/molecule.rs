//! Molecule graph: atoms, bonds, and adjacency list.

use crate::atom::{Atom, RGroupLabel};
use crate::bond::{BondEntry, BondOrder};
use crate::element::Element;
use crate::stereo_group::StereoGroup;

/// Newtype index for an atom in a Molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomIdx(pub u32);

/// Newtype index for a bond in a Molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BondIdx(pub u32);

/// Error types for molecule construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolError {
    /// Atom index out of range.
    InvalidAtomIdx(AtomIdx),
    /// Duplicate bond between the same pair of atoms.
    DuplicateBond(AtomIdx, AtomIdx),
}

impl core::fmt::Display for MolError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidAtomIdx(idx) => write!(f, "invalid atom index: {}", idx.0),
            Self::DuplicateBond(a, b) => {
                write!(f, "duplicate bond between atoms {} and {}", a.0, b.0)
            }
        }
    }
}

impl std::error::Error for MolError {}

/// An immutable molecular graph built via [`MoleculeBuilder`].
///
/// Representation: atom list + bond list + per-atom adjacency list.
/// No external graph library is used; all graph traversal is domain-aware.
/// Sentinel used in `stereo_neighbor_order` to represent the implicit H in a bracket atom.
pub const STEREO_H_SENTINEL: u32 = u32::MAX;

#[derive(Clone)]
pub struct Molecule {
    atoms: Vec<Atom>,
    /// Non-chemical caller tags, indexed by atom. Empty until a tag is set.
    /// Kept outside `Atom` to preserve its public struct and equality contract.
    atom_tags: Vec<Option<core::num::NonZeroU16>>,
    bonds: Vec<BondEntry>,
    /// adjacency[atom_idx] = list of (neighbor_atom_idx, bond_idx)
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
    /// Lossless display/interchange identity for `R`/`R<n>` wildcard atoms.
    /// Kept outside [`Atom`] so chemistry primitives and the public atom
    /// struct remain stable and element-focused.
    r_groups: std::collections::HashMap<u32, RGroupLabel>,
    /// Enhanced stereo groups (ChemDraw V3000 Absolute / Or / And).
    stereo_groups: Vec<StereoGroup>,
    /// SMILES-text-order neighbor sequence for chiral atoms.
    ///
    /// Keyed by atom index.  Each value lists the atom indices of neighbors in
    /// the order they appeared in the SMILES string (including ring-closure
    /// partners), with [`STEREO_H_SENTINEL`] (`u32::MAX`) standing in for the
    /// implicit bracket H.  Populated by the SMILES parser; absent for atoms
    /// not parsed from SMILES or without recorded stereo.
    stereo_neighbor_order: std::collections::HashMap<u32, Vec<u32>>,
    /// Directional (`/`, `\`) bond marker, stashed for bonds whose `order` was
    /// overwritten to `Aromatic` (e.g. an exocyclic C=N adjacent to an
    /// aromatic ring atom — `order` must stay `Aromatic` for SMARTS `:a`
    /// matching, but the E/Z direction would otherwise be lost). Keyed by
    /// bond index; value is `BondOrder::Up` or `BondOrder::Down`. Absent for
    /// bonds whose direction is already carried directly by `order`.
    bond_directions: std::collections::HashMap<u32, BondOrder>,
    /// Parser-side anchor endpoint for a directional marker stashed on an
    /// aromatic bond. This preserves which endpoint introduced the marker
    /// when atom/bond indices are rebuilt during canonical search.
    bond_direction_anchors: std::collections::HashMap<u32, AtomIdx>,
    /// SMILES ring-closure bonds, keyed by bond index: the ring label in the
    /// high 32 bits and the label's occurrence in the low 32. Recorded by the
    /// SMILES parser so [`Molecule::rdkit_bond_order`] can list bonds the way
    /// RDKit's parser creates them; empty for other inputs and cleared by any
    /// edit that renumbers bonds.
    smiles_ring_closure_keys: std::collections::HashMap<u32, u64>,
    /// Memoized derived perception data (SSSR, aromatic view, ...). Cleared
    /// by every mutation and empty after `clone()`.
    derived: crate::derived_cache::DerivedCache,
}

include!("molecule/access.rs");
include!("molecule/graph.rs");
include!("molecule/metadata.rs");
include!("molecule/cache.rs");
include!("molecule/builder.rs");

#[cfg(test)]
include!("molecule_tests.rs");
