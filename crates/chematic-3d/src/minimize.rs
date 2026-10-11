//! Simplified force-field geometry minimization for molecular structures.
//!
//! Uses bounded force-field minimizers over bond, angle, torsion, out-of-plane,
//! van der Waals, and electrostatic terms. The production MMFF94 bridge uses
//! the prepared analytic energy/gradient pair; reference and legacy paths
//! retain finite-difference or gradient-descent implementations.

#[allow(unused_imports)]
use crate::fmath::DetMath;
use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_ff::{
    EnergyBreakdown, MinimizerError, Mmff94Convergence, Mmff94EnergyModel, Mmff94TerminationReason,
    NumericTypeError, OOP_SP2_TYPES, UffType, angle_type_for,
    assign_mmff94_numeric_types_with_view, assign_uff_types, bond_type_for,
    is_angle_in_ring_of_size_3_or_4, minimize_uff as ff_minimize_uff,
    minimize_uff_with_constraint as ff_minimize_uff_with_constraint, mmff94_angle_energy_resolved,
    mmff94_bond_energy_resolved, mmff94_oop, mmff94_stbn, mmff94_torsion_term_params,
    stretch_bend_type_for, uff_total_energy,
};
use chematic_ff::{
    assign_dreiding_types, assign_mmff94_types, dreiding_angle, dreiding_bond_len, dreiding_vdw,
    mmff94_angle_params, mmff94_bond_params, mmff94_charges_3d, mmff94_vdw_params,
};
use chematic_perception::find_sssr;

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Force field parameters
// ---------------------------------------------------------------------------

/// Bond stretching spring constant (kcal/mol/Ų).
/// Used in both DREIDING and generic force fields.
const BOND_SPRING_CONSTANT: f64 = 700.0;

/// Angle bending spring constant (kcal/mol/rad²).
/// Used in both DREIDING and generic force fields.
const ANGLE_SPRING_CONSTANT: f64 = 100.0;

/// Van der Waals interaction cutoff distance (Ångströms).
/// Interactions beyond this distance are ignored.
const VDW_CUTOFF: f64 = 8.0;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

include!("minimize/legacy.rs");
include!("minimize/policy.rs");
include!("minimize/bridge.rs");

#[cfg(test)]
#[path = "minimize_tests.rs"]
mod embedded_tests;
