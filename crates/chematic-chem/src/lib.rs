//! `chematic-chem` — molecular descriptors for chematic.
//!
//! All descriptor functions take a `&Molecule` reference.
//! Values are approximate; calibrated against RDKit defaults.

#![forbid(unsafe_code)]

pub mod cip;
pub mod descriptors;
pub mod scaffold;
pub mod standardize;
pub mod tautomer;

pub use cip::{assign_cip, CipAssignment};
pub use descriptors::{
    exact_mass, heavy_atom_count, hba_count, hbd_count, lipinski_passes,
    logp_crippen, molecular_weight, rotatable_bond_count, tpsa,
};

pub use scaffold::{generic_murcko_scaffold, murcko_scaffold};
pub use standardize::{largest_fragment, neutralize_charges};
pub use tautomer::{canonical_tautomer, enumerate_tautomers};
