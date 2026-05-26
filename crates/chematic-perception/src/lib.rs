//! `chematic-perception` — molecular perception algorithms.
//!
//! Provides:
//! - [`sssr`]: Smallest Set of Smallest Rings (SSSR) via Balducci-Pearlman algorithm.
//! - [`aromaticity`]: Hückel aromaticity perception for kekulized molecules.

#![forbid(unsafe_code)]

pub mod aromaticity;
pub mod sssr;

pub use aromaticity::{AromaticityModel, assign_aromaticity};
pub use sssr::{RingSet, find_sssr};
