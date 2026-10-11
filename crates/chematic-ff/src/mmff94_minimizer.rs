//! MMFF94 geometry minimizer using complete Halgren 1996 parameters.
//!
//! Provides full MMFF94 energy evaluation (bond + angle + torsion + vdW + electrostatic)
//! and geometry optimization with two algorithms:
//! - **Steepest descent** (`minimize_mmff94_full`) — robust, simple
//! - **L-BFGS** (`minimize_mmff94_lbfgs`) — faster convergence, quasi-Newton
//!
//! ## Energy terms
//! - **Bond**: cubic-corrected harmonic (Halgren MMFF.II eq. 1)
//! - **Angle**: cubic-corrected harmonic (Halgren MMFF.III eq. 2)
//! - **Torsion**: three-term Fourier (Halgren MMFF.IV)
//! - **vdW**: buffered 14-7 potential with Slater-Kirkwood combining rule (Halgren MMFF.I eq. 2)
//! - **Electrostatic**: Coulomb with δ buffer (Halgren MMFF.V eq. 14)

#[allow(unused_imports)]
use crate::fmath::DetMath;
use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;
use rayon::prelude::*;

use crate::mmff94_energy::{
    AngleEnergyParams, BondEnergyParams, TorsionEnergyParams, mmff94_angle_energy_resolved,
    mmff94_bond_energy_resolved, mmff94_oop, mmff94_stbn_oriented, mmff94_vdw_combined,
};
use crate::mmff94_numeric::{
    NumericTypeError, assign_mmff94_numeric_types_with_view, mmff94_charges_numeric,
};

type CoordVec = Vec<[f64; 3]>;
type LbfgsHistory = VecDeque<(CoordVec, CoordVec, f64)>;
const LBFGS_HISTORY_SIZE: usize = 80;

/// RDKit's `MDYNE_A_TO_KCAL_MOL` (mdyn Å → kcal/mol).
const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;
const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
/// Angle and out-of-plane prefactor, `MDYNE_A_TO_KCAL_MOL * DEG2RAD^2`
/// as RDKit computes it (0.0438438…; Halgren rounds it to 0.043844).
const MMFF_ANGLE_CONV: f64 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
/// Magnitude of the cubic angle-bend constant, −0.4 rad⁻¹ in degrees as
/// RDKit writes it (`cb = -0.006981317`; Halgren rounds it to 0.007).
const MMFF_ANGLE_CUBIC: f64 = 0.006981317;
/// Stretch-bend prefactor, `MDYNE_A_TO_KCAL_MOL * DEG2RAD` (2.512076…).
const MMFF_STBN_CONV: f64 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
/// RDKit's `nonBondedThresh` default: nonbonded pairs farther apart than
/// this (Å) contribute nothing.
const NONBONDED_THRESH: f64 = 100.0;
type VdwPairs = Vec<PreparedVdwPair>;
type ElectrostaticPairs = Vec<PreparedElectrostaticPair>;

/// Below this size a linear pass over prepared vdW pairs is faster than
/// rebuilding spatial hash bins for every line-search energy/gradient probe.
/// The cutoff and pair set are identical; only candidate enumeration changes.
const DIRECT_VDW_ATOM_THRESHOLD: usize = 64;

include!("mmff94_minimizer/model.rs");
include!("mmff94_minimizer/optimize.rs");
include!("mmff94_minimizer/terms.rs");
include!("mmff94_minimizer/energy.rs");
include!("mmff94_minimizer/typing.rs");

#[cfg(test)]
#[path = "mmff94_minimizer_tests.rs"]
mod embedded_tests;
