//! MMFF94 angle bending parameters (Halgren 1996 Table V).

use super::{AngleEnergyParams, Mmff94Resolution};
#[allow(unused_imports)]
use crate::fmath::DetMath;

mod table;

pub use table::MMFF94_ANGLE_ENERGY;
use table::MMFF94_EQ_LEVEL;

/// `stage` in `0..4`, matching RDKit's Level 2/3/4/5. Falls back to `atom_type`
/// itself (identity, i.e. equivalent to an exact-match retry) for any type
/// not in [`MMFF94_EQ_LEVEL`] -- see that table's own doc for why this is safe.
fn eq_level(atom_type: u8, stage: usize) -> u8 {
    MMFF94_EQ_LEVEL
        .binary_search_by_key(&atom_type, |&(t, _)| t)
        .map(|idx| MMFF94_EQ_LEVEL[idx].1[stage])
        .ok()
        .or_else(|| {
            crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(atom_type)
                .map(|info| info.equivalence_levels[stage])
        })
        .unwrap_or(atom_type)
}

fn has_eq_level_row(atom_type: u8) -> bool {
    MMFF94_EQ_LEVEL
        .binary_search_by_key(&atom_type, |&(t, _)| t)
        .is_ok()
        || crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(atom_type).is_some()
}

/// Outcome of the full table search chain (RDKit's real eqLevel ladder plus
/// chematic's own generic-`angle_type` net), before any empirical fallback.
enum TableSearch {
    /// A real, usable row (`ka != 0.0`).
    Hit(AngleEnergyParams, Mmff94Resolution),
    /// A row was found but `ka == 0.0` -- RDKit's `isDoubleZero` placeholder
    /// (issue #227 Stage C: 97 `type_i == 0 && type_k == 0` rows restored to
    /// [`MMFF94_ANGLE_ENERGY`], present in RDKit's real
    /// `defaultMMFFAngleData` but absent from chematic's first port). Only
    /// `theta0` is meaningful; `ka` must come from the empirical rule.
    ZeroKa { theta0: f64 },
    /// No row at any stage, and every stage that mattered used a
    /// well-defined `eqLevel` substitution (or none at all).
    Exhausted,
    /// The ladder needed to substitute a numeric type absent from
    /// [`MMFF94_EQ_LEVEL`] (`type_i`/`type_k` outside the 1-55 range) to
    /// reach the point of giving up. RDKit's real
    /// `MMFFDefCollection::operator()` returns `nullptr` for such a type and
    /// `MMFFAngleCollection::operator()` dereferences it with no null
    /// check -- undefined behavior in RDKit's own C++, confirmed by direct
    /// source read at the pinned commit. There is no well-defined "real"
    /// answer to reproduce here, so this is never treated as resolvable --
    /// see `scripts/mmff94_provenance/PROVENANCE.md`'s Stage C entry for the
    /// oracle investigation this is based on (the `(0, 43, 18, 63)` tuple).
    UndefinedSubstitution,
}

fn search_table(angle_type: u8, type_i: u8, type_j: u8, type_k: u8) -> TableSearch {
    let search = |at: u8, ti: u8, tk: u8, kind: Mmff94Resolution| {
        MMFF94_ANGLE_ENERGY
            .binary_search_by_key(&(at, ti, type_j, tk), |&(a0, a, b, c, _, _)| (a0, a, b, c))
            .ok()
            .map(|idx| {
                let (_, _, _, _, ka, theta0) = MMFF94_ANGLE_ENERGY[idx];
                (ka, theta0, kind)
            })
    };
    let eq_level_search = |stage: usize, level: u8| {
        let (ci, ck) = (eq_level(type_i, stage), eq_level(type_k, stage));
        let (lo, hi) = if ci <= ck { (ci, ck) } else { (ck, ci) };
        search(
            angle_type,
            lo,
            hi,
            Mmff94Resolution::EquivalentType { level },
        )
    };
    let hit = search(angle_type, type_i, type_k, Mmff94Resolution::DirectTable)
        .or_else(|| search(angle_type, type_k, type_i, Mmff94Resolution::DirectTable))
        .or_else(|| eq_level_search(1, 3))
        .or_else(|| eq_level_search(2, 4))
        // RDKit's `MMFFAngleCollection` stops here: no retry at angle type 0;
        // an angle the four levels miss takes the empirical rule (a 4-ring
        // 3-40-20 angle: theta0 90, not the open-chain type-0 row's 112).
        .or_else(|| eq_level_search(3, 5));
    match hit {
        Some((ka, theta0, kind)) if ka != 0.0 => {
            TableSearch::Hit(AngleEnergyParams { ka, theta0 }, kind)
        }
        Some((_, theta0, _)) => TableSearch::ZeroKa { theta0 },
        None if has_eq_level_row(type_i) && has_eq_level_row(type_k) => TableSearch::Exhausted,
        None => TableSearch::UndefinedSubstitution,
    }
}

/// Look up angle bending parameters by MMFF94 numeric atom types.
///
/// `type_j` is the central atom. Both orderings (ti,tj,tk) and (tk,tj,ti) are
/// tried at the exact atom types first, then via RDKit's real `eqLevel`
/// canonical-type-substitution ladder (`MMFFAngleCollection::operator()`,
/// `Code/ForceField/MMFF/Params.h` at the pinned commit -- issue #227 Stage
/// B): `type_i`/`type_k` are substituted through 4 equivalence-class stages
/// (`type_j` and `angle_type` stay fixed throughout -- RDKit's real
/// algorithm never substitutes the center or changes the angle type). Only
/// once the full eqLevel ladder is exhausted does this fall back to the
/// generic `angle_type=0` row for the same triple — a chematic-specific
/// safety net RDKit's own algorithm doesn't have, kept from before this
/// change: without it, a *correct* ring/bond-type classification for a
/// triple the (much smaller) specialized table doesn't happen to cover
/// would silently drop the angle term entirely, which is worse than the
/// un-classified behavior it replaces.
///
/// A `ka == 0.0` placeholder row (issue #227 Stage C, see [`TableSearch`])
/// is never surfaced here -- it is not a usable parameter on its own, only
/// [`mmff94_angle_energy_resolved`] can turn it into one via the empirical
/// rule. This keeps this function's contract identical to before Stage C
/// restored those rows to the table: every input that used to return
/// `Some`/`None` still does.
pub fn mmff94_angle_energy(
    angle_type: u8,
    type_i: u8,
    type_j: u8,
    type_k: u8,
) -> Option<AngleEnergyParams> {
    match search_table(angle_type, type_i, type_j, type_k) {
        TableSearch::Hit(params, _) => Some(params),
        TableSearch::ZeroKa { .. }
        | TableSearch::Exhausted
        | TableSearch::UndefinedSubstitution => None,
    }
}

static MMFF94_ANGLE_Z: &[(u8, f64)] = &[
    (1, 1.395),
    (6, 2.494),
    (7, 2.711),
    (8, 3.045),
    (9, 2.847),
    (14, 2.350),
    (15, 2.350),
    (16, 2.980),
    (17, 2.909),
    (35, 3.017),
    (53, 3.086),
];

static MMFF94_ANGLE_C: &[(u8, f64)] = &[
    (6, 1.016),
    (7, 1.113),
    (8, 1.337),
    (14, 0.811),
    (15, 1.068),
    (16, 1.249),
    (17, 1.078),
];

fn angle_z(atomic_number: u8) -> Option<f64> {
    MMFF94_ANGLE_Z
        .binary_search_by_key(&atomic_number, |&(z, _)| z)
        .ok()
        .map(|idx| MMFF94_ANGLE_Z[idx].1)
}

fn angle_c(atomic_number: u8) -> Option<f64> {
    MMFF94_ANGLE_C
        .binary_search_by_key(&atomic_number, |&(z, _)| z)
        .ok()
        .map(|idx| MMFF94_ANGLE_C[idx].1)
}

/// Halgren MMFF.V eq. 20's `theta0` sub-rule (page 627-628, immediately
/// above eq. 20 itself), run only for [`TableSearch::Exhausted`] -- when a
/// `ka == 0.0` placeholder row IS found, its own `theta0` is reused verbatim
/// instead and this sub-rule never runs at all
/// (`getMMFFAngleBendEmpiricalRuleParams`, `AtomTyper.cpp`). The ring-size
/// 60°/90° override below is specific to this from-scratch branch and is
/// NOT applied in the reuse case -- confirmed by direct source read at the
/// pinned commit (`AtomTyper.cpp` lines ~2746-2785: the override lives
/// inside the `if (!oldMMFFAngleParams)` branch only).
fn angle_empirical_theta0(
    atomic_number_j: u8,
    crd: u8,
    val: u8,
    mltb: u8,
    linear: bool,
    ring_size: u8,
) -> f64 {
    let mut theta0 = 120.0;
    match crd {
        4 => theta0 = 109.45,
        2 => {
            if atomic_number_j == 8 {
                theta0 = 105.0;
            } else if linear {
                theta0 = 180.0;
            }
        }
        3 if val == 3 && mltb == 0 => {
            theta0 = if atomic_number_j == 7 { 107.0 } else { 92.0 };
        }
        _ => {}
    }
    match ring_size {
        3 => theta0 = 60.0,
        4 => theta0 = 90.0,
        _ => {}
    }
    theta0
}

/// Halgren MMFF.V eq. 20 empirical angle-bend force constant (page 628,
/// Table VI's per-element `Z`/`C` constants). Returns `None` if any of the
/// three atoms' element has no `Z`/`C` entry -- RDKit would silently use
/// `0.0` there, a physically-invalid zero force constant; chematic fails
/// closed instead (issue #227 Stage C acceptance gate).
fn angle_empirical_ka(
    atomic_number_i: u8,
    atomic_number_j: u8,
    atomic_number_k: u8,
    r0_ij: f64,
    r0_jk: f64,
    ring_size: u8,
    theta0: f64,
) -> Option<f64> {
    let zi = angle_z(atomic_number_i)?;
    let zk = angle_z(atomic_number_k)?;
    let cj = angle_c(atomic_number_j)?;
    let mut beta = 1.75;
    match ring_size {
        4 => beta *= 0.85,
        3 => beta *= 0.05,
        _ => {}
    }
    let d = (r0_ij - r0_jk).powi(2) / (r0_ij + r0_jk).powi(2);
    let theta0_rad = theta0.to_radians();
    Some(beta * zi * cj * zk / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * (2.0 * d).dexp()))
}

/// Look up angle bending parameters, falling back to Halgren's MMFF.V eq. 20
/// empirical rule when the table (including the eqLevel ladder) has no
/// usable force constant, and reporting which mechanism actually produced
/// the result (issue #227 Stage C).
///
/// `r0_ij`/`r0_jk` are the two flanking bonds' equilibrium lengths (Å),
/// supplied by the caller rather than re-resolved here -- a caller that
/// already computed them (e.g. for the adjacent stretch-bend term) does not
/// pay for, or risk diverging from, a second bond resolution. `ring_size` is
/// `3`/`4` if the triple is in a 3- or 4-membered ring
/// (`isAngleInRingOfSize3or4`, local bond adjacency, NOT SSSR -- see
/// `crate::mmff94_minimizer::is_angle_in_ring_of_size_3_or_4`), else `0`.
///
/// Returns `None` when the triple cannot be resolved at all: the eqLevel
/// ladder needed an undefined substitution ([`TableSearch::UndefinedSubstitution`]),
/// an outer/central atom type has no [`crate::mmff94_numeric_type_registry`]
/// entry, or [`angle_empirical_ka`] has no `Z`/`C` constant for one of the
/// three elements.
pub fn mmff94_angle_energy_resolved(
    angle_type: u8,
    type_i: u8,
    type_j: u8,
    type_k: u8,
    r0_ij: f64,
    r0_jk: f64,
    ring_size: u8,
) -> Option<(AngleEnergyParams, Mmff94Resolution)> {
    use crate::mmff94_numeric_type_registry::mmff94_numeric_type_info;

    match search_table(angle_type, type_i, type_j, type_k) {
        TableSearch::Hit(params, kind) => Some((params, kind)),
        TableSearch::ZeroKa { theta0 } => {
            let info_i = mmff94_numeric_type_info(type_i)?;
            let info_j = mmff94_numeric_type_info(type_j)?;
            let info_k = mmff94_numeric_type_info(type_k)?;
            let ka = angle_empirical_ka(
                info_i.atomic_number,
                info_j.atomic_number,
                info_k.atomic_number,
                r0_ij,
                r0_jk,
                ring_size,
                theta0,
            )?;
            Some((
                AngleEnergyParams { ka, theta0 },
                Mmff94Resolution::EmpiricalAngle,
            ))
        }
        TableSearch::Exhausted => {
            let info_i = mmff94_numeric_type_info(type_i)?;
            let info_j = mmff94_numeric_type_info(type_j)?;
            let info_k = mmff94_numeric_type_info(type_k)?;
            let theta0 = angle_empirical_theta0(
                info_j.atomic_number,
                info_j.coordination,
                info_j.valence,
                info_j.multiple_bond_count,
                info_j.linear,
                ring_size,
            );
            let ka = angle_empirical_ka(
                info_i.atomic_number,
                info_j.atomic_number,
                info_k.atomic_number,
                r0_ij,
                r0_jk,
                ring_size,
                theta0,
            )?;
            Some((
                AngleEnergyParams { ka, theta0 },
                Mmff94Resolution::EmpiricalAngle,
            ))
        }
        TableSearch::UndefinedSubstitution => None,
    }
}

#[cfg(test)]
#[path = "angle_tests.rs"]
mod embedded_tests;
