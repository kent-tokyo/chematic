//! `chematic-chem` — molecular descriptors for chematic.
//!
//! All descriptor functions take a `&Molecule` reference.
//! Values are approximate; calibrated against RDKit defaults.

#![forbid(unsafe_code)]

pub mod alerts;
pub mod brics;
pub mod cip;
pub mod descriptors;
pub mod qed;
pub mod scaffold;
pub mod standardize;
pub mod tautomer;

pub use cip::{assign_cip, CipAssignment};
pub use descriptors::{
    aromatic_ring_count, egan_passes, exact_mass, formal_charge_sum, fsp3,
    ghose_passes, heavy_atom_count, hba_count, hbd_count, lipinski_passes,
    logp_crippen, molar_refractivity, molecular_weight, num_aliphatic_heterocycles,
    num_aliphatic_rings, num_aromatic_heterocycles, num_bridgehead_atoms,
    num_heteroatoms, num_saturated_heterocycles, num_saturated_rings,
    num_spiro_atoms, num_stereocenters, num_unspecified_stereocenters,
    reos_passes, ring_count, rotatable_bond_count, tpsa, veber_passes,
};

pub use alerts::{pains_matches, pains_passes};
pub use brics::{brics_bonds, brics_fragments};
pub use scaffold::{generic_murcko_scaffold, murcko_scaffold};
pub use standardize::{largest_fragment, neutralize_charges};
pub use tautomer::{canonical_tautomer, enumerate_tautomers};
pub use qed::qed;
