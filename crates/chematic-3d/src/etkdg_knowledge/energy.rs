//! Torsion energy evaluation and torsion-space optimization (Wave 2 spec §7).
//!
//! # Why torsion-space, not Cartesian, optimization
//!
//! `optimize_torsions` rotates rigid molecular fragments about acyclic
//! (bridge) single bonds only -- it never does Cartesian gradient descent
//! on atom positions directly. This is a deliberate design choice, not a
//! missing feature:
//!
//! - A rotation about a genuine bridge bond `B-C` is, by construction, a
//!   **rigid, proper (determinant +1) rotation** of everything on one side
//!   of that bond. Rigid rotations preserve every pairwise distance within
//!   the rotated fragment exactly (to floating-point precision) -- so bond
//!   lengths and ring-closure distances *within* either side are
//!   structurally guaranteed to be unchanged, not merely checked
//!   after the fact. Cartesian gradient descent on torsion energy has no
//!   such guarantee: it can and typically will stretch bonds to "cheat" a
//!   better torsion angle, exactly the failure mode spec §7 forbids.
//! - Ring and macrocycle bonds are never rotated (rotating a bond that is
//!   *not* a bridge -- i.e. is part of a ring -- has no well-defined single
//!   rigid two-side split, since both "sides" stay connected via the rest
//!   of the ring). Their [`TorsionPotential`]s are still **scored** by
//!   [`evaluate_torsion_energy`] (so the report reflects the true total
//!   torsion energy, ring contributions included), just never
//!   mechanically adjusted by [`optimize_torsions`] -- an honest,
//!   documented limitation, not a silent gap.
//! - A proper rotation cannot change configuration (R/S at a stereocenter,
//!   E/Z across a double bond) -- only conformation. A chiral center's own
//!   four substituents are either entirely on the fixed side or entirely on
//!   the rotated side of any bridge bond that does not pass through the
//!   center itself (rotating the bond does not touch the *relative*
//!   arrangement of the stereocenter's own neighbors), and this
//!   implementation never rotates about a double bond in the first place
//!   (`classify_bond` excludes double/triple bonds from torsion matching
//!   entirely, so no [`TorsionPotential`] this crate produces ever has a
//!   double bond as its central bond). This is verified empirically in this
//!   module's tests (`rotation_about_bridge_bond_preserves_chirality_sign`:
//!   signed chiral volume before/after, on a real, non-degenerate 4-distinct-
//!   atom quadruple, with an explicit non-zero-movement assertion so the
//!   test cannot silently degrade into comparing unmoved coordinates) --
//!   not merely asserted. An earlier version of this test used a
//!   duplicate-atom-index quadruple (the exact defect this crate's own audit
//!   flags in the *legacy* test suite) that made the rotation a geometric
//!   no-op, so the assertion passed without checking anything; fixed after
//!   a later independent review pass caught it, see that test's own comment
//!   for the full account.

use chematic_core::{AtomIdx, Molecule};
use std::collections::{HashMap, HashSet};

use crate::coords::{Coords3D, Point3};

use super::types::{TorsionKnowledgeError, TorsionPotential};

// ---------------------------------------------------------------------------
// Dihedral geometry
// ---------------------------------------------------------------------------

/// Dihedral angle (degrees, in `(-180, 180]`) of the A-B-C-D quadruple given
/// their 3D coordinates. Standard atan2-based formula (numerically stable
/// near 0/180, unlike an acos-based formula); returns a finite value (0.0)
/// rather than panicking for degenerate (collinear) geometry -- `atan2(0,0)`
/// is well-defined in Rust as `0.0`.
fn dihedral_deg(coords: &Coords3D, atoms: [AtomIdx; 4]) -> f64 {
    let p0 = coords.get(atoms[0]);
    let p1 = coords.get(atoms[1]);
    let p2 = coords.get(atoms[2]);
    let p3 = coords.get(atoms[3]);

    let b1 = p1.sub(&p0);
    let b2 = p2.sub(&p1);
    let b3 = p3.sub(&p2);

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let b2_unit = b2.try_normalize().unwrap_or(Point3::zero());
    let m1 = n1.cross(&b2_unit);

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    y.atan2(x).to_degrees()
}

// ---------------------------------------------------------------------------
// Energy report
// ---------------------------------------------------------------------------

/// Result of scoring every potential in `potentials` against one geometry.
#[derive(Clone, Debug, Default)]
pub struct TorsionEnergyReport {
    pub total_energy: f64,
    /// `(rule_id, energy)` for every potential, in the same order as the
    /// input slice.
    pub per_potential_energy: Vec<(String, f64)>,
    pub max_abs_gradient_deg: f64,
    /// Count of potentials whose energy or gradient came out non-finite
    /// (excluded from `total_energy`/`max_abs_gradient_deg`, never silently
    /// folded in as if they were zero).
    pub n_non_finite: usize,
}

/// Evaluate every potential in `potentials` against `coords`. Pure/read-only
/// (never mutates `coords`).
pub fn evaluate_torsion_energy(
    mol: &Molecule,
    coords: &Coords3D,
    potentials: &[TorsionPotential],
) -> Result<TorsionEnergyReport, TorsionKnowledgeError> {
    if coords.atom_count() != mol.atom_count() {
        return Err(TorsionKnowledgeError::CoordsAtomCountMismatch);
    }
    let mut report = TorsionEnergyReport::default();
    for pot in potentials {
        for &a in &pot.atoms {
            if a.0 as usize >= coords.atom_count() {
                return Err(TorsionKnowledgeError::InvalidTopology);
            }
        }
        let phi = dihedral_deg(coords, pot.atoms);
        let e = pot.energy(phi);
        let g = pot.d_energy_d_phi_deg(phi);
        report.per_potential_energy.push((pot.rule_id.clone(), e));
        if !e.is_finite() || !g.is_finite() {
            report.n_non_finite += 1;
            continue;
        }
        report.total_energy += e;
        if g.abs() > report.max_abs_gradient_deg {
            report.max_abs_gradient_deg = g.abs();
        }
    }
    Ok(report)
}

// ---------------------------------------------------------------------------
// Optimization
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsionOptimizationConfig {
    pub max_iterations: usize,
    /// Fixed rotation step (degrees) applied per accepted move -- a
    /// fixed-*magnitude* step (halved on backtracking) rather than a
    /// gradient-scaled one, so behavior does not depend on the arbitrary
    /// energy-unit scale of `amplitude`. Direction is NOT taken from the
    /// analytic gradient's sign: `rotate_fragment`'s Rodrigues-rotation angle
    /// sign and `dihedral_deg`'s phi-sign convention are two independently
    /// derived formulas with no proven relationship between them, so the
    /// line search tries both +step and -step at each magnitude and keeps
    /// whichever improves energy (see `optimize_torsions`'s inner loop).
    /// An earlier version of this doc claimed "direction from the gradient's
    /// sign", which stopped being true once the both-directions search was
    /// added to fix 3 non-convergence failures (step-size persistence was
    /// the other fix) -- corrected after independent review caught the
    /// comment/code mismatch.
    pub step_deg: f64,
    /// Converged when every rotatable bond's `|dE/dphi|` (summed over its
    /// potentials) is below this threshold.
    pub convergence_grad_deg: f64,
    /// Maximum halvings of `step_deg` tried per bond per iteration before
    /// giving up on improving that bond this iteration.
    pub max_line_search_steps: usize,
}

impl Default for TorsionOptimizationConfig {
    fn default() -> Self {
        Self {
            max_iterations: 200,
            step_deg: 5.0,
            convergence_grad_deg: 0.01,
            max_line_search_steps: 8,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct TorsionOptimizationReport {
    pub energy_before: f64,
    pub energy_after: f64,
    pub iterations_used: usize,
    pub converged: bool,
    /// Maximum change in any real bond's length, before vs. after (should
    /// be ~0 by construction -- rigid rotation about a bridge bond cannot
    /// change any bond length; measured and reported rather than assumed,
    /// per this program's standing "verify, don't assume" practice).
    pub max_bond_length_delta: f64,
    /// Maximum change in any ring-membership atom pair's distance (any two
    /// atoms that co-occur in the same ring), before vs. after -- same
    /// structural guarantee as bond length, measured independently.
    pub max_ring_closure_delta: f64,
    pub rotated_bond_count: usize,
}

/// Rotate the rigid fragment on `mol`'s "moving" side of bridge bond
/// `(b, c)` (the connected component containing `c` once edge `(b,c)` is
/// removed) by `delta_deg` about the `b`->`c` axis, in place on `coords`.
fn rotate_fragment(
    mol: &Molecule,
    coords: &mut Coords3D,
    b: AtomIdx,
    c: AtomIdx,
    delta_deg: f64,
) -> usize {
    let moving = component_excluding_edge(mol, c, b, c);
    let origin = coords.get(b);
    let axis = coords.get(c).sub(&origin);
    let Some(axis_unit) = axis.try_normalize() else {
        return 0; // degenerate (coincident b/c) -- nothing sensible to rotate
    };
    let theta = delta_deg.to_radians();
    let (sin_t, cos_t) = (theta.sin(), theta.cos());
    for &atom in &moving {
        let p = coords.get(atom).sub(&origin);
        // Rodrigues' rotation formula about `axis_unit` through the origin.
        let rotated = p
            .scale(cos_t)
            .add(&axis_unit.cross(&p).scale(sin_t))
            .add(&axis_unit.scale(axis_unit.dot(&p) * (1.0 - cos_t)));
        coords.set(atom, rotated.add(&origin));
    }
    moving.len()
}

/// Every atom reachable from `start` without crossing the edge
/// `(avoid_a, avoid_b)` (in either direction).
fn component_excluding_edge(
    mol: &Molecule,
    start: AtomIdx,
    avoid_a: AtomIdx,
    avoid_b: AtomIdx,
) -> HashSet<AtomIdx> {
    let mut seen = HashSet::new();
    let mut stack = vec![start];
    seen.insert(start);
    while let Some(cur) = stack.pop() {
        for (nbr, _) in mol.neighbors(cur) {
            if (cur == avoid_a && nbr == avoid_b) || (cur == avoid_b && nbr == avoid_a) {
                continue;
            }
            if seen.insert(nbr) {
                stack.push(nbr);
            }
        }
    }
    seen
}

/// `true` if bond `(a, b)` is a bridge (removing it disconnects `a` from
/// `b`) -- i.e. it is not part of any ring. Only bridge bonds are ever
/// mechanically rotated by [`optimize_torsions`].
///
/// `pub(crate)`: this is the ONE definition of "is this central bond
/// mechanically rotatable" in the crate. Wave 2/3 Coordinator integration
/// (`pipeline_v2.rs`) reuses it directly (via `crate::etkdg_knowledge::is_bridge_bond`,
/// re-exported `pub(crate)` from the parent module) to decide, per
/// [`TorsionPotential`], whether it was actually applied to geometry vs. scored-only
/// -- rather than a second, independently-derived predicate (e.g. `ring_size.is_some()`,
/// which is wrong: `rules_basic::flat_ring` potentials carry `ring_size: None` despite
/// targeting a ring bond -- or re-deriving `classify_bond(..).ring == NotInRing`,
/// which is a genuine second definition of the same fact that could silently drift
/// from this one). One predicate, no drift.
pub(crate) fn is_bridge_bond(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    !component_excluding_edge(mol, a, a, b).contains(&b)
}

fn max_bond_length_delta(mol: &Molecule, before: &Coords3D, after: &Coords3D) -> f64 {
    mol.bonds()
        .map(|(_, bond)| {
            let d0 = before.get(bond.atom1).distance(&before.get(bond.atom2));
            let d1 = after.get(bond.atom1).distance(&after.get(bond.atom2));
            (d1 - d0).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn max_ring_closure_delta(rings: &[Vec<AtomIdx>], before: &Coords3D, after: &Coords3D) -> f64 {
    let mut worst = 0.0_f64;
    for ring in rings {
        let n = ring.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let d0 = before.get(ring[i]).distance(&before.get(ring[j]));
                let d1 = after.get(ring[i]).distance(&after.get(ring[j]));
                worst = worst.max((d1 - d0).abs());
            }
        }
    }
    worst
}

/// Optimize dihedral angles to reduce total torsion energy (from
/// `potentials`), touching only genuinely rotatable (bridge, single, non-
/// ring) central bonds. Never changes topology, never rotates a ring bond,
/// and is a rigid rotation so it structurally cannot stretch a bond or
/// break a ring (measured and reported in
/// [`TorsionOptimizationReport::max_bond_length_delta`]/
/// [`max_ring_closure_delta`] as a verification, not merely assumed).
///
/// Deterministic given identical `(mol, coords, potentials, config)` -- pure
/// steepest-descent-direction coordinate descent, no randomness.
///
/// # Errors
/// - [`TorsionKnowledgeError::CoordsAtomCountMismatch`] /
///   [`InvalidTopology`](TorsionKnowledgeError::InvalidTopology): as
///   [`evaluate_torsion_energy`].
/// - [`TorsionKnowledgeError::NonFiniteEnergy`]: any potential ever produces
///   a non-finite energy/gradient during optimization (never silently
///   skipped mid-run).
/// - [`TorsionKnowledgeError::RingIntegrityViolated`]: the post-hoc bond-
///   length or ring-closure self-check (see module docs -- these should be
///   structurally impossible given rigid-rotation-only moves, but are
///   checked, not assumed) found a change beyond floating-point tolerance.
/// - [`TorsionKnowledgeError::NonConvergence`]: `max_iterations` was reached
///   without every rotatable bond's gradient falling below
///   `config.convergence_grad_deg`. A typed failure, not a silently-flagged
///   `converged: false` buried in an `Ok` -- matching spec §7's explicit
///   "non-convergence is a typed failure, not silent."
pub fn optimize_torsions(
    mol: &Molecule,
    coords: &Coords3D,
    potentials: &[TorsionPotential],
    config: &TorsionOptimizationConfig,
) -> Result<(Coords3D, TorsionOptimizationReport), TorsionKnowledgeError> {
    if coords.atom_count() != mol.atom_count() {
        return Err(TorsionKnowledgeError::CoordsAtomCountMismatch);
    }

    let before_report = evaluate_torsion_energy(mol, coords, potentials)?;
    if before_report.n_non_finite > 0 {
        return Err(TorsionKnowledgeError::NonFiniteEnergy);
    }

    let mut working = coords.clone();

    // Partition potentials: only those whose central bond is a genuine
    // bridge (acyclic, single, non-ring) are ever mechanically rotated.
    let rotatable: Vec<&TorsionPotential> = potentials
        .iter()
        .filter(|p| {
            let (b, c) = p.central_bond;
            mol.bond_between(b, c).is_some() && is_bridge_bond(mol, b, c)
        })
        .collect();

    let rings =
        chematic_perception::augmented_ring_set(mol, chematic_perception::find_sssr(mol).rings());

    let mut iterations_used = 0;
    let mut converged = false;

    // Per-bond step size, persisted *across* outer iterations (not reset to
    // `config.step_deg` every iteration): a plain fixed-step steepest
    // descent oscillates forever around a minimum without ever getting the
    // gradient below a tight tolerance. Backtracking line search already
    // halves the step within one iteration when a move doesn't improve
    // energy; carrying that shrunk value into the *next* iteration too
    // (rather than resetting) gives the standard diminishing-step-size
    // convergence guarantee for a unimodal (or locally unimodal) objective,
    // deterministically and without introducing randomness.
    let mut steps: HashMap<(AtomIdx, AtomIdx), f64> = rotatable
        .iter()
        .map(|p| (p.central_bond, config.step_deg))
        .collect();

    for _iter in 0..config.max_iterations.max(1) {
        iterations_used += 1;
        let mut max_grad = 0.0_f64;

        for pot in &rotatable {
            let (b, c) = pot.central_bond;
            let phi = dihedral_deg(&working, pot.atoms);
            let grad = pot.d_energy_d_phi_deg(phi);
            if !grad.is_finite() {
                return Err(TorsionKnowledgeError::NonFiniteEnergy);
            }
            if grad.abs() > max_grad {
                max_grad = grad.abs();
            }
            if grad.abs() < config.convergence_grad_deg {
                continue;
            }

            let e_before = pot.energy(phi);
            // Try both rotation directions at each step size, rather than
            // committing to `-grad.signum()` alone: `rotate_fragment`'s
            // rotation-angle sign and `dihedral_deg`'s phi-sign convention
            // are two independently-derived formulas, and this module does
            // not assume (without checking) that a positive rotation angle
            // maps to a positive phi change -- trying both sides makes the
            // descent correct regardless of that relationship, at the cost
            // of at most 2x the trials.
            let step_slot = steps.entry(pot.central_bond).or_insert(config.step_deg);
            let mut step = *step_slot;
            let mut accepted = false;
            for _ls in 0..config.max_line_search_steps.max(1) {
                let mut best: Option<(Coords3D, f64)> = None;
                for &sign in &[1.0_f64, -1.0] {
                    let mut trial = working.clone();
                    rotate_fragment(mol, &mut trial, b, c, sign * step);
                    let phi_trial = dihedral_deg(&trial, pot.atoms);
                    let e_trial = pot.energy(phi_trial);
                    if !e_trial.is_finite() {
                        return Err(TorsionKnowledgeError::NonFiniteEnergy);
                    }
                    if best.as_ref().is_none_or(|(_, e_best)| e_trial < *e_best) {
                        best = Some((trial, e_trial));
                    }
                }
                let (trial, e_trial) = best.expect("both directions were tried");
                if e_trial < e_before {
                    working = trial;
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            let _ = accepted; // no action needed either way; loop continues to next bond
            *step_slot = step;
        }

        if max_grad < config.convergence_grad_deg {
            converged = true;
            break;
        }
    }

    let bond_delta = max_bond_length_delta(mol, coords, &working);
    let ring_delta = max_ring_closure_delta(rings.as_slice(), coords, &working);
    const GEOMETRY_INTEGRITY_TOLERANCE: f64 = 1e-6;
    if bond_delta > GEOMETRY_INTEGRITY_TOLERANCE || ring_delta > GEOMETRY_INTEGRITY_TOLERANCE {
        return Err(TorsionKnowledgeError::RingIntegrityViolated);
    }

    if !converged {
        return Err(TorsionKnowledgeError::NonConvergence);
    }

    let after_report = evaluate_torsion_energy(mol, &working, potentials)?;

    let report = TorsionOptimizationReport {
        energy_before: before_report.total_energy,
        energy_after: after_report.total_energy,
        iterations_used,
        converged,
        max_bond_length_delta: bond_delta,
        max_ring_closure_delta: ring_delta,
        rotated_bond_count: rotatable.len(),
    };

    Ok((working, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::Coords3D;
    use crate::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};
    use crate::etkdg_knowledge::{FourierTorsionTerm, TorsionKnowledgeSource};
    use chematic_smiles::parse;

    fn butane_coords() -> (chematic_core::Molecule, Coords3D) {
        let mol = parse("CCCC").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        (mol, coords)
    }

    fn butane_potential() -> TorsionPotential {
        TorsionPotential {
            atoms: [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
            central_bond: (AtomIdx(1), AtomIdx(2)),
            source: TorsionKnowledgeSource::StandardExperimental,
            rule_id: "test:butane_anti".to_string(),
            terms: vec![FourierTorsionTerm::from_rdkit(1, 1, 5.0)], // prefers 180 deg
            ring_size: None,
        }
    }

    #[test]
    fn dihedral_is_periodic_and_bounded() {
        let (_, coords) = butane_coords();
        let phi = dihedral_deg(&coords, [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)]);
        assert!(phi.is_finite());
        assert!((-180.0..=180.0).contains(&phi), "{phi}");
    }

    /// `rotate_fragment` and `dihedral_deg` are two independently derived
    /// formulas (Rodrigues' rotation vs. an atan2-based dihedral) with no
    /// *proven* relationship between their angle-sign conventions --
    /// `optimize_torsions` works around that by trying both rotation
    /// directions in its line search (see `TorsionOptimizationConfig::
    /// step_deg`'s doc). That workaround makes the optimizer correct either
    /// way, but it would also silently paper over a genuinely wrong "moving"
    /// atom set or a degenerate rotation axis producing a phi change that
    /// isn't `delta_deg` at all (e.g. a change of ~0 or some unrelated
    /// magnitude). This test checks the magnitude directly, independent of
    /// sign: after a single `rotate_fragment(delta_deg)` about bond (1,2),
    /// `|phi_after - phi_before|` must equal `|delta_deg|` (mod wraparound).
    #[test]
    fn rotate_fragment_changes_dihedral_by_exactly_delta_magnitude() {
        let (mol, mut coords) = butane_coords();
        let atoms = [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)];
        let phi_before = dihedral_deg(&coords, atoms);
        let delta = 17.0;
        let moved = rotate_fragment(&mol, &mut coords, AtomIdx(1), AtomIdx(2), delta);
        assert!(moved > 0, "rotate_fragment must move at least one atom");
        let phi_after = dihedral_deg(&coords, atoms);
        let mut raw_diff = (phi_after - phi_before).abs();
        if raw_diff > 180.0 {
            raw_diff = 360.0 - raw_diff; // wraparound at the +/-180 boundary
        }
        assert!(
            (raw_diff - delta.abs()).abs() < 1e-6,
            "expected |phi change| == {}, got {} (before={}, after={})",
            delta.abs(),
            raw_diff,
            phi_before,
            phi_after
        );
    }

    #[test]
    fn evaluate_energy_rejects_mismatched_atom_count() {
        let (mol, _coords) = butane_coords();
        let wrong = Coords3D::new_zeroed(1);
        let pot = butane_potential();
        let err = evaluate_torsion_energy(&mol, &wrong, &[pot]).unwrap_err();
        assert_eq!(err, TorsionKnowledgeError::CoordsAtomCountMismatch);
    }

    #[test]
    fn evaluate_energy_is_finite_for_a_real_geometry() {
        let (mol, coords) = butane_coords();
        let pot = butane_potential();
        let report = evaluate_torsion_energy(&mol, &coords, &[pot]).unwrap();
        assert!(report.total_energy.is_finite());
        assert_eq!(report.n_non_finite, 0);
    }

    #[test]
    fn analytic_gradient_matches_finite_difference() {
        let term = FourierTorsionTerm::from_rdkit(2, -1, 8.0);
        for phi in [-170.0, -90.0, -1.0, 0.0, 37.0, 91.5, 179.0] {
            let h = 1e-4;
            let fd = (term.energy(phi + h) - term.energy(phi - h)) / (2.0 * h);
            let analytic = term.d_energy_d_phi_deg(phi);
            assert!(
                (fd - analytic).abs() < 1e-4,
                "phi={phi}: fd={fd}, analytic={analytic}"
            );
        }
    }

    #[test]
    fn optimize_torsions_reduces_or_holds_energy() {
        let (mol, coords) = butane_coords();
        let pot = butane_potential();
        let config = TorsionOptimizationConfig::default();
        let (_new_coords, report) = optimize_torsions(&mol, &coords, &[pot], &config).unwrap();
        assert!(
            report.energy_after <= report.energy_before + 1e-9,
            "{report:?}"
        );
        assert!(report.converged);
        assert!(report.max_bond_length_delta < 1e-6, "{report:?}");
        assert!(report.max_ring_closure_delta < 1e-6, "{report:?}");
    }

    #[test]
    fn optimize_torsions_never_stretches_a_bond() {
        let (mol, coords) = butane_coords();
        let pot = butane_potential();
        let config = TorsionOptimizationConfig::default();
        let (new_coords, _report) = optimize_torsions(&mol, &coords, &[pot], &config).unwrap();
        for (_, bond) in mol.bonds() {
            let before = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
            let after = new_coords
                .get(bond.atom1)
                .distance(&new_coords.get(bond.atom2));
            assert!(
                (after - before).abs() < 1e-6,
                "bond stretched: {before} -> {after}"
            );
        }
    }

    #[test]
    fn optimize_torsions_never_breaks_a_ring() {
        let mol = parse("C1CCCCC1").unwrap(); // cyclohexane
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        // A potential whose central bond IS a ring bond: must be scored but
        // never mechanically rotated (no bridge exists for a ring bond).
        let pot = TorsionPotential {
            atoms: [AtomIdx(5), AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            central_bond: (AtomIdx(0), AtomIdx(1)),
            source: TorsionKnowledgeSource::SmallRingExperimental,
            rule_id: "test:ring_bond".to_string(),
            terms: vec![FourierTorsionTerm::from_rdkit(3, 1, 10.0)],
            ring_size: Some(6),
        };
        let config = TorsionOptimizationConfig::default();
        let (new_coords, report) = optimize_torsions(&mol, &coords, &[pot], &config).unwrap();
        assert_eq!(
            report.rotated_bond_count, 0,
            "ring bond must not be rotated"
        );
        // Ring geometry must be essentially untouched.
        for i in 0..mol.atom_count() {
            let p0 = coords.get(AtomIdx(i as u32));
            let p1 = new_coords.get(AtomIdx(i as u32));
            assert!(p0.distance(&p1) < 1e-6, "ring atom {i} moved");
        }
    }

    #[test]
    fn optimize_torsions_typed_non_convergence() {
        let (mol, coords) = butane_coords();
        let pot = butane_potential();
        let config = TorsionOptimizationConfig {
            max_iterations: 1,
            step_deg: 0.0001,
            convergence_grad_deg: 1e-12, // unreachable in 1 tiny-step iteration
            max_line_search_steps: 1,
        };
        let err = optimize_torsions(&mol, &coords, &[pot], &config).unwrap_err();
        assert_eq!(err, TorsionKnowledgeError::NonConvergence);
    }

    #[test]
    fn optimize_torsions_rejects_nonfinite_energy_upfront() {
        let (mol, coords) = butane_coords();
        let mut pot = butane_potential();
        pot.terms = vec![FourierTorsionTerm::new(1, 0.0, f64::NAN)];
        let config = TorsionOptimizationConfig::default();
        let err = optimize_torsions(&mol, &coords, &[pot], &config).unwrap_err();
        assert_eq!(err, TorsionKnowledgeError::NonFiniteEnergy);
    }

    #[test]
    fn deterministic_given_identical_input() {
        let (mol, coords) = butane_coords();
        let pot = butane_potential();
        let config = TorsionOptimizationConfig::default();
        let (c1, r1) =
            optimize_torsions(&mol, &coords, std::slice::from_ref(&pot), &config).unwrap();
        let (c2, r2) =
            optimize_torsions(&mol, &coords, std::slice::from_ref(&pot), &config).unwrap();
        for i in 0..mol.atom_count() {
            assert_eq!(c1.get(AtomIdx(i as u32)), c2.get(AtomIdx(i as u32)));
        }
        assert_eq!(r1.iterations_used, r2.iterations_used);
    }

    #[test]
    fn rotation_about_bridge_bond_preserves_chirality_sign() {
        // 2-pentanol: a stereocenter (C1) with a rotatable, genuinely
        // 4-heavy-atom C3-C4 bond further down the chain (atom 4 has its own
        // further neighbor, atom 5 -- NOT a duplicate-atom-index placeholder;
        // see the fix note below). Rotating that distal bond must never flip
        // the stereocenter's own signed (chiral) volume.
        //
        // An earlier version of this test used 2-butanol with atoms
        // `[1,3,4,4]` (atom 4 duplicated, since 2-butanol's terminal ethyl
        // carbon has no 4th atom beyond it) -- exactly the "duplicate atom
        // index" defect this crate's own audit (`docs/3d_torsion_knowledge_
        // audit.md` §4) flags in the *legacy* test suite, found here in this
        // PR's own new tests by a later independent review pass. With a
        // duplicate index, `dihedral_deg` degenerates (`atan2(0,0) == 0`),
        // the rotating fragment never actually moves (confirmed: distance
        // moved was exactly 0.0), and the "chirality preserved" assertion
        // was trivially true on unmoved coordinates -- not verified
        // empirically, despite this module's own doc comment (top of file)
        // claiming it was. Fixed by using a molecule where the rotated bond
        // has a real, distinct 4th atom, and by asserting the fragment
        // actually moved (a non-zero-movement guard), so this test cannot
        // silently degrade into a no-op again.
        let mol = parse("C[C@H](O)CCC").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();

        // Stereocenter is atom 1 (the @ carbon per SMILES order); its
        // neighbors are atoms 0 (CH3), 2 (O), 3 (CH2 of the propyl arm).
        let center = AtomIdx(1);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(n, _)| n).collect();
        assert!(
            neighbors.len() >= 3,
            "stereocenter must have >=3 heavy neighbors"
        );
        let signed_volume = |c: &Coords3D| -> f64 {
            let p0 = c.get(center);
            let v1 = c.get(neighbors[0]).sub(&p0);
            let v2 = c.get(neighbors[1]).sub(&p0);
            let v3 = c.get(neighbors[2]).sub(&p0);
            v1.cross(&v2).dot(&v3)
        };
        let before_sign = signed_volume(&coords).signum();

        // Rotate the distal C3-C4 bond, far from the stereocenter. Atom 5
        // (beyond atom 4) is the real, distinct 4th atom -- the rotated
        // fragment is {atom4, atom5}, a genuine 2-atom rigid-body rotation.
        assert!(
            mol.bond_between(AtomIdx(3), AtomIdx(4)).is_some(),
            "2-pentanol must have a real C3-C4 bond"
        );
        // The real "not a duplicate-index placeholder" sanity check: atom 5
        // must actually exist and be bonded to atom 4 (i.e. genuinely be the
        // further neighbor down the chain) -- `assert_ne!(AtomIdx(4),
        // AtomIdx(5))` (an earlier version of this fix) compares two
        // hardcoded literals and can never fail, so it checked nothing.
        assert!(
            mol.bond_between(AtomIdx(4), AtomIdx(5)).is_some(),
            "2-pentanol's atom 5 must be a real, distinct neighbor of atom 4, not a placeholder"
        );
        let pot = TorsionPotential {
            atoms: [AtomIdx(1), AtomIdx(3), AtomIdx(4), AtomIdx(5)],
            central_bond: (AtomIdx(3), AtomIdx(4)),
            source: TorsionKnowledgeSource::StandardExperimental,
            rule_id: "test:distal_bond".to_string(),
            terms: vec![FourierTorsionTerm::from_rdkit(1, 1, 5.0)],
            ring_size: None,
        };
        // `EmbedParameters::default()` is seeded (fixed `random_seed`), so
        // this embedding is deterministic -- convergence here is a real,
        // reproducible fact about this fixture, not a coin flip. Panicking
        // loudly on failure (rather than an earlier version's `else {
        // return; }`) matters specifically because a silent skip here is
        // the exact "test quietly does nothing" pattern this test was
        // rewritten to eliminate in the first place (see the fix note
        // above) -- an `else { return; }` guarding the one assertion this
        // test exists for would reintroduce it one line later.
        let config = TorsionOptimizationConfig::default();
        let (new_coords, report) = optimize_torsions(&mol, &coords, &[pot], &config)
            .expect("optimize_torsions must converge on this fixed, seeded embedding");

        // The whole point of using a real (non-degenerate) quadruple: the
        // rotation must have actually moved something, or this test is
        // exactly as vacuous as the duplicate-index version it replaces.
        let moved = coords.get(AtomIdx(5)).distance(&new_coords.get(AtomIdx(5)));
        assert!(
            moved > 1e-6,
            "rotation must actually move the distal atom, not silently no-op: moved={moved}, report={report:?}"
        );

        let after_sign = signed_volume(&new_coords).signum();
        assert_eq!(
            before_sign, after_sign,
            "chirality sign flipped after rotation"
        );
    }

    // -----------------------------------------------------------------
    // Negative controls (spec §13): deliberately corrupt a geometry and
    // confirm the self-check machinery this module relies on actually
    // detects it -- a check that can't fail isn't a check.
    // -----------------------------------------------------------------

    /// Negative control for "improving torsion score by breaking/tearing a
    /// ring" (spec §13, must FAIL to occur): manually stretch one bond in a
    /// copy of a real geometry (simulating what a buggy, non-rigid
    /// optimizer step might do) and confirm `max_bond_length_delta` -- the
    /// same measurement `optimize_torsions` bases its
    /// `RingIntegrityViolated` error on -- actually detects it, rather than
    /// reporting a reassuring 0.0 regardless of input.
    #[test]
    fn max_bond_length_delta_detects_an_artificially_stretched_bond() {
        let (mol, coords) = butane_coords();
        let mut corrupted = coords.clone();
        let (i, j) = {
            let (_, bond) = mol.bonds().next().unwrap();
            (bond.atom1, bond.atom2)
        };
        let pj = corrupted.get(j);
        corrupted.set(j, pj.add(&Point3::new(5.0, 0.0, 0.0))); // tear the bond
        let delta = max_bond_length_delta(&mol, &coords, &corrupted);
        assert!(delta > 1.0, "expected a large detected delta, got {delta}");
        let _ = i;
    }

    /// Same idea for ring-closure integrity: manually move one ring atom in
    /// a copy of a real cyclohexane geometry and confirm
    /// `max_ring_closure_delta` detects the resulting distance change to
    /// another ring atom, rather than silently reporting 0.0.
    #[test]
    fn max_ring_closure_delta_detects_an_artificially_broken_ring() {
        let mol = parse("C1CCCCC1").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        let rings = chematic_perception::augmented_ring_set(
            &mol,
            chematic_perception::find_sssr(&mol).rings(),
        );
        let mut corrupted = coords.clone();
        let p0 = corrupted.get(AtomIdx(0));
        corrupted.set(AtomIdx(0), p0.add(&Point3::new(3.0, 0.0, 0.0)));
        let delta = max_ring_closure_delta(rings.as_slice(), &coords, &corrupted);
        assert!(
            delta > 1.0,
            "expected a large detected ring-closure delta, got {delta}"
        );
    }

    /// Negative control for "flipping stereo and still reporting success"
    /// (spec §13, must FAIL to occur): manually mirror one substituent of a
    /// stereocenter through the center (simulating what a buggy stereo-blind
    /// transform might do) and confirm the signed-chiral-volume check this
    /// module's positive test relies on actually flips sign -- proving that
    /// check has teeth, not just that it passes on unmodified input.
    #[test]
    fn signed_volume_check_detects_an_artificial_chirality_flip() {
        let mol = parse("C[C@H](O)CC").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        let center = AtomIdx(1);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(n, _)| n).collect();
        let signed_volume = |c: &Coords3D| -> f64 {
            let p0 = c.get(center);
            let v1 = c.get(neighbors[0]).sub(&p0);
            let v2 = c.get(neighbors[1]).sub(&p0);
            let v3 = c.get(neighbors[2]).sub(&p0);
            v1.cross(&v2).dot(&v3)
        };
        let before_sign = signed_volume(&coords).signum();

        // Deliberately flip: mirror one substituent through the center.
        let mut corrupted = coords.clone();
        let center_p = corrupted.get(center);
        let sub = corrupted.get(neighbors[0]);
        let mirrored = center_p.scale(2.0).sub(&sub);
        corrupted.set(neighbors[0], mirrored);
        let after_sign = signed_volume(&corrupted).signum();

        assert_ne!(
            before_sign, after_sign,
            "the injected mirror operation must actually flip the detected sign -- otherwise this check couldn't catch a real regression"
        );
    }

    /// Negative control (spec §13, "accepting NaN/Inf energy" must FAIL to
    /// occur): `evaluate_torsion_energy` must exclude a non-finite
    /// potential's energy from the total and count it, never silently fold
    /// it in as if it were zero.
    #[test]
    fn evaluate_energy_excludes_nonfinite_terms_and_counts_them() {
        let (mol, coords) = butane_coords();
        let mut bad = butane_potential();
        bad.terms = vec![FourierTorsionTerm::new(1, 0.0, f64::NAN)];
        let good = butane_potential();
        let report = evaluate_torsion_energy(&mol, &coords, &[bad, good.clone()]).unwrap();
        assert_eq!(report.n_non_finite, 1);
        let good_only =
            evaluate_torsion_energy(&mol, &coords, std::slice::from_ref(&good)).unwrap();
        assert!(
            (report.total_energy - good_only.total_energy).abs() < 1e-9,
            "non-finite term must be excluded from total_energy, not folded in as 0"
        );
    }
}
