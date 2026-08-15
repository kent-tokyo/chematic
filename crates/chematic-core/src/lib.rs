#![forbid(unsafe_code)]
//! `chematic-core` — core cheminformatics types for the chematic ecosystem.
//!
//! Provides:
//! - [`Element`]: periodic-table element indexed by atomic number.
//! - [`Atom`]: a single atom with element, charge, isotope, aromaticity, chirality.
//! - [`BondOrder`] / [`BondEntry`]: bond types and graph edges.
//! - [`Molecule`] / [`MoleculeBuilder`]: undirected molecular graph (no external graph lib).
//! - [`implicit_hcount`]: valence-based implicit hydrogen count for organic-subset atoms.
//! - [`kekulize`] / [`apply_kekule`]: assign alternating single/double bonds to aromatic systems.
//!
//! # Design principles
//! - Pure Rust, no C/C++ FFI, no unsafe.
//! - Zero external dependencies: compiles to wasm32-unknown-unknown without modification.
//! - Domain-aware abstractions: graph types carry chemical semantics, not just topology.

pub mod atom;
pub mod bond;
pub mod coords3d;
pub mod element;
pub mod kekulization;
pub mod molecule;
pub mod stereo_geometry;
pub mod stereo_group;
pub mod valence;

// Re-export the most commonly used types at crate root.
pub use atom::{Atom, Chirality, CipCode, SquarePlanarPermutation};
pub use bond::{BondEntry, BondOrder};
pub use coords3d::{Coords3D, Point3};
pub use element::Element;
pub use kekulization::{KekuleError, KekuleResult, apply_kekule, kekulize};
pub use molecule::{AtomIdx, BondIdx, MolError, Molecule, MoleculeBuilder, STEREO_H_SENTINEL};
// `StereoConfiguration`/`CanonicalStereoConfiguration`/`canonicalize_configuration`/
// `equivalent_under_rotation` are deliberately `pub(crate)`, not re-exported
// here: all four are hardcoded to `[u32; 4]`, which only fits today's two
// 4-coordinate geometries. Only the two arity-fixed-by-their-own-tag-semantics
// bridge functions are public; see `stereo_geometry.rs`'s own doc comments
// (on `StereoConfiguration` and each bridge function) for the full reasoning.
pub use stereo_geometry::{
    StereoGeometry, StereoGeometryError, remap_square_planar_tag, remap_tetrahedral_parity,
};
pub use stereo_group::{StereoGroup, StereoGroupKind};
#[allow(deprecated)]
pub use valence::total_hcount;
pub use valence::{
    ValenceError, bond_order_sum, implicit_hcount, valence_inferred_hcount, validate_valence,
};
