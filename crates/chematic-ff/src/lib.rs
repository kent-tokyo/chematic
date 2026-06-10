#![forbid(unsafe_code)]
//! `chematic-ff` — Force field atom typing and parameters.
//!
//! Provides pure-Rust force field support for molecular mechanics calculations:
//! - **DREIDING**: general-purpose force field (existing)
//! - **MMFF94**: Merck Molecular Force Field (new, for small molecules)
//!
//! Includes atom type enumerations, assignment functions, and parameter lookups.

pub mod dreiding;
pub mod mmff94;
pub mod params;

pub use dreiding::{DREIDINGType, assign_dreiding_types};
pub use mmff94::{MMFF94Type, assign_mmff94_types, AssignError};
pub use params::{dreiding_angle, dreiding_bond_len, dreiding_torsion_barrier, dreiding_vdw};
