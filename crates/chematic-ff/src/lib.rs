#![forbid(unsafe_code)]
//! `chematic-ff` — Force field atom typing and parameters.
//!
//! Provides pure-Rust force field support for molecular mechanics calculations:
//! - **DREIDING**: general-purpose force field
//! - **MMFF94**: Merck Molecular Force Field (small organic molecules)
//! - **UFF**: Universal Force Field (all elements, metal complexes)
//!
//! Includes atom type enumerations, assignment functions, and parameter lookups.

pub mod dreiding;
pub mod mmff94;
pub mod mmff94_advanced;
pub mod mmff94_bci;
pub mod mmff94_energy;
pub mod mmff94_minimizer;
pub mod mmff94_numeric;
pub mod mmff94_numeric_type_registry;
pub mod mmff94_params;
pub mod params;
pub mod uff;

pub use dreiding::{DREIDINGType, assign_dreiding_types};
pub use mmff94::{AssignError, MMFF94Type, assign_mmff94_types, mmff94_charges_3d};
pub use mmff94_advanced::{ElectrostaticMatrix, MMFF94BatchProperties};
pub use mmff94_bci::{bci, mmff94_charges_bci, mmff94_formal_charge};
pub use mmff94_energy::{
    AngleEnergyParams, BondEnergyParams, Mmff94Resolution, TorsionEnergyParams, VdwEnergyParams,
    mmff94_angle_energy, mmff94_angle_energy_resolved, mmff94_bond_energy,
    mmff94_bond_energy_resolved, mmff94_oop, mmff94_stbn, mmff94_stbn_type_only,
    mmff94_torsion_energy, mmff94_vdw_combined, mmff94_vdw_energy,
};
pub use mmff94_minimizer::{
    EnergyBreakdown, MLTB_TYPES, MinimizeResult, MinimizerError, OOP_SP2_TYPES, angle_type_for,
    bond_type_for, is_angle_in_ring_of_size_3_or_4, minimize_mmff94_full, minimize_mmff94_lbfgs,
    mmff94_energy_breakdown, mmff94_torsion_scan, mmff94_total_energy, stretch_bend_type_for,
    torsion_no_term_by_design, torsion_type_for,
};
pub use mmff94_numeric::{
    NumericTypeError, assign_mmff94_numeric_types, assign_mmff94_numeric_types_with_view,
    mmff94_charges_numeric, pbci_for,
};
pub use mmff94_numeric_type_registry::{Mmff94NumericTypeInfo, mmff94_numeric_type_info};
pub use mmff94_params::{
    AngleParams, BondDipoleParams, BondParams, ChargeParams, ElectrostaticScalingParams,
    MMFF94MoleculeProperties, TorsionParams, VdWParams, mmff94_angle_params, mmff94_bond_dipole,
    mmff94_bond_params, mmff94_charge_params, mmff94_electrostatic_scaling_1_4,
    mmff94_torsion_params, mmff94_vdw_params,
};
pub use params::{dreiding_angle, dreiding_bond_len, dreiding_torsion_barrier, dreiding_vdw};
pub use uff::{UffMinimizeResult, UffType, assign_uff_types, minimize_uff, uff_total_energy};
