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
pub mod mmff94_advanced;
pub mod mmff94_bci;
pub mod mmff94_energy;
pub mod mmff94_numeric;
pub mod mmff94_params;
pub mod params;

pub use dreiding::{DREIDINGType, assign_dreiding_types};
pub use mmff94::{MMFF94Type, assign_mmff94_types, AssignError, mmff94_charges_3d};
pub use mmff94_bci::{bci, mmff94_charges_bci, mmff94_formal_charge};
pub use mmff94_advanced::{
    ElectrostaticMatrix, MMFF94BatchProperties,
};
pub use mmff94_params::{
    mmff94_angle_params, mmff94_bond_dipole, mmff94_bond_params, mmff94_charge_params,
    mmff94_electrostatic_scaling_1_4, mmff94_torsion_params, mmff94_vdw_params,
    AngleParams, BondDipoleParams, BondParams, ChargeParams, ElectrostaticScalingParams,
    MMFF94MoleculeProperties, TorsionParams, VdWParams,
};
pub use mmff94_energy::{
    mmff94_angle_energy, mmff94_bond_energy, mmff94_torsion_energy, mmff94_vdw_combined,
    mmff94_vdw_energy, AngleEnergyParams, BondEnergyParams, TorsionEnergyParams, VdwEnergyParams,
};
pub use mmff94_numeric::{
    assign_mmff94_numeric_types, mmff94_charges_numeric, pbci_for, NumericTypeError,
};
pub use params::{dreiding_angle, dreiding_bond_len, dreiding_torsion_barrier, dreiding_vdw};
