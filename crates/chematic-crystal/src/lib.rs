#![forbid(unsafe_code)]
//! `chematic-crystal` — periodic (crystal) structure representation and
//! geometry for the chematic ecosystem.
//!
//! This crate is a pure structural/geometric foundation: [`Lattice`],
//! fractional/Cartesian coordinates, [`PeriodicSite`]/[`PeriodicStructure`],
//! periodic-boundary-condition displacement and distance, periodic neighbor
//! enumeration, and diagonal supercells. It is deliberately **not** an
//! extension of `chematic_core::Molecule` -- see the crate README and
//! `docs/rfcs/chematic_crystal_foundation.md` for why a periodic crystal
//! and a molecular bond graph are kept as distinct first-class types.
//!
//! # Scope
//!
//! No symmetry (space groups, Wyckoff positions, Niggli reduction), no CIF
//! parser, no XRD, no DFT, no materials-property prediction. See
//! `docs/crystal_scope.md` for the full non-goal list.
//!
//! # Design principles
//!
//! - Pure Rust, `#![forbid(unsafe_code)]`, zero required dependencies
//!   beyond `chematic-core` (optional `serde`).
//! - Compiles to `wasm32-unknown-unknown` without modification.
//! - Public constructors reject `NaN`/`Infinity` rather than propagate them.
//! - Deterministic output ordering everywhere (no float-keyed sorts).

pub mod error;
pub mod lattice;
pub mod neighbor;
pub mod periodic;
pub mod poscar;
pub mod site;
pub mod structure;
pub mod supercell;
pub mod validation;

pub use error::CrystalError;
pub use lattice::Lattice;
pub use neighbor::PeriodicNeighbor;
pub use periodic::{PeriodicDisplacement, minimum_image};
pub use poscar::{
    PoscarDocument, PoscarError, PredictorCorrector, parse_contcar, parse_poscar, write_poscar,
};
pub use site::{CartesianCoord, FractionalCoord, Occupancy, PeriodicSite, SiteSpecies};
pub use structure::PeriodicStructure;
