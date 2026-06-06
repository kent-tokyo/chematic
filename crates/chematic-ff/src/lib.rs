#![forbid(unsafe_code)]
//! `chematic-ff` — DREIDING force field atom typing and parameters.
//!
//! Provides pure-Rust DREIDING force field support for molecular mechanics calculations.
//! Includes:
//! - [`DREIDINGType`]: enumeration of DREIDING atom types
//! - [`assign_dreiding_types`]: automatic atom type assignment
//! - VDW, bond length, and angle parameters via module `params`

pub mod dreiding;
pub mod params;

pub use dreiding::{assign_dreiding_types, DREIDINGType};
pub use params::{dreiding_angle, dreiding_bond_len, dreiding_torsion_barrier, dreiding_vdw};
