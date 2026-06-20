//! Distance geometry constraint satisfaction via projection.
//!
//! This module provides constraint-based coordinate correction for molecular
//! 3D geometry. It enforces bond distances and valence angles through iterative
//! geometric projection, integrated with distance geometry + force field minimization.
//!
//! **Algorithm Strategy:**
//! - Simpler than metric matrix (no eigendecomposition)
//! - O(n·k) per iteration: n atoms, k constraints
//! - Fast convergence: 5-10 iterations typical for small molecules
//! - Directly enforces user intent (bond distances, angles)
//!
//! **Pipeline:**
//! ```ignore
//! generate_coords(mol)              // rule-based 3D placement
//!   ↓
//! build_constraints(mol)            // extract bond/angle targets
//!   ↓
//! satisfy_constraints()             // iterative projection
//!   ↓
//! minimize_dreiding()               // energy minimization
//! ```

use crate::coords::{Coords3D, Point3};
use chematic_core::{AtomIdx, BondOrder, Molecule};
use core::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constraint Data Structures
// ---------------------------------------------------------------------------

/// Bond distance constraint: enforce |P_a - P_b| = target_distance ± tolerance
#[derive(Debug, Clone)]
pub struct BondConstraint {
    pub atom1: AtomIdx,
    pub atom2: AtomIdx,
    pub target_distance: f64, // Å
    pub tolerance: f64,       // Å
}

impl BondConstraint {
    /// Create a new bond constraint with default tolerance (±0.05 Å).
    pub fn new(atom1: AtomIdx, atom2: AtomIdx, target_distance: f64) -> Self {
        Self {
            atom1,
            atom2,
            target_distance,
            tolerance: 0.05,
        }
    }

    /// Check if current distance satisfies constraint.
    pub fn satisfied(&self, coords: &Coords3D) -> bool {
        let d = coords.get(self.atom1).distance(&coords.get(self.atom2));
        let lower = self.target_distance - self.tolerance;
        let upper = self.target_distance + self.tolerance;
        d >= lower && d <= upper
    }

    /// Violation magnitude; 0.0 if satisfied.
    pub fn violation(&self, coords: &Coords3D) -> f64 {
        let d = coords.get(self.atom1).distance(&coords.get(self.atom2));
        let lower = self.target_distance - self.tolerance;
        let upper = self.target_distance + self.tolerance;
        if d < lower {
            lower - d
        } else if d > upper {
            d - upper
        } else {
            0.0
        }
    }
}

/// Valence angle constraint: enforce ∠(atom1—center—atom2) ≈ target_angle ± tolerance
#[derive(Debug, Clone)]
pub struct AngleConstraint {
    pub atom1: AtomIdx,
    pub center: AtomIdx,
    pub atom2: AtomIdx,
    pub target_angle: f64, // radians
    pub tolerance: f64,    // radians
}

impl AngleConstraint {
    /// Create a new angle constraint with default tolerance (±5°).
    pub fn new(atom1: AtomIdx, center: AtomIdx, atom2: AtomIdx, target_angle: f64) -> Self {
        Self {
            atom1,
            center,
            atom2,
            target_angle,
            tolerance: 5.0_f64.to_radians(),
        }
    }

    /// Check if current angle satisfies constraint.
    pub fn satisfied(&self, coords: &Coords3D) -> bool {
        let angle = compute_angle(coords, self.atom1, self.center, self.atom2);
        let lower = self.target_angle - self.tolerance;
        let upper = self.target_angle + self.tolerance;
        angle >= lower && angle <= upper
    }

    /// Violation magnitude; 0.0 if satisfied.
    pub fn violation(&self, coords: &Coords3D) -> f64 {
        let angle = compute_angle(coords, self.atom1, self.center, self.atom2);
        let lower = self.target_angle - self.tolerance;
        let upper = self.target_angle + self.tolerance;
        if angle < lower {
            lower - angle
        } else if angle > upper {
            angle - upper
        } else {
            0.0
        }
    }
}

/// Complete set of bond and angle constraints for a molecule.
#[derive(Debug, Clone)]
pub struct ConstraintSet {
    pub bonds: Vec<BondConstraint>,
    pub angles: Vec<AngleConstraint>,
}

impl ConstraintSet {
    /// Count of constraints currently violated (outside tolerance).
    pub fn violated_count(&self, coords: &Coords3D) -> usize {
        self.bonds.iter().filter(|c| !c.satisfied(coords)).count()
            + self.angles.iter().filter(|c| !c.satisfied(coords)).count()
    }

    /// Maximum violation magnitude among all constraints.
    pub fn max_violation(&self, coords: &Coords3D) -> f64 {
        let bond_violations = self
            .bonds
            .iter()
            .map(|c| c.violation(coords))
            .fold(0.0, f64::max);
        let angle_violations = self
            .angles
            .iter()
            .map(|c| c.violation(coords))
            .fold(0.0, f64::max);
        bond_violations.max(angle_violations)
    }
}

// ---------------------------------------------------------------------------
// Ideal Parameters (Reuse from dg.rs / minimize.rs)
// ---------------------------------------------------------------------------

/// Get ideal bond length from element pair and bond order.
fn get_ideal_bond_length(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    let ea = mol.atom(a).element;
    let eb = mol.atom(b).element;

    let order = mol
        .bond_between(a, b)
        .map(|(_, bond)| bond.order)
        .unwrap_or(BondOrder::Single);

    let (lo, hi) = if ea.atomic_number() <= eb.atomic_number() {
        (ea.atomic_number(), eb.atomic_number())
    } else {
        (eb.atomic_number(), ea.atomic_number())
    };

    match (lo, hi, order) {
        // C–C
        (6, 6, BondOrder::Single) | (6, 6, BondOrder::Up) | (6, 6, BondOrder::Down) => 1.54,
        (6, 6, BondOrder::Double) => 1.34,
        (6, 6, BondOrder::Triple) => 1.20,
        (6, 6, BondOrder::Aromatic) => 1.40,
        // C–N
        (6, 7, BondOrder::Single) | (6, 7, BondOrder::Up) | (6, 7, BondOrder::Down) => 1.47,
        (6, 7, BondOrder::Double) => 1.27,
        (6, 7, BondOrder::Triple) => 1.16,
        (6, 7, BondOrder::Aromatic) => 1.34,
        // C–O
        (6, 8, BondOrder::Single) | (6, 8, BondOrder::Up) | (6, 8, BondOrder::Down) => 1.43,
        (6, 8, BondOrder::Double) => 1.22,
        (6, 8, BondOrder::Aromatic) => 1.36,
        // C–S
        (6, 16, _) => 1.82,
        // C–F
        (6, 9, _) => 1.35,
        // C–Cl
        (6, 17, _) => 1.77,
        // C–Br
        (6, 35, _) => 1.94,
        // C–I
        (6, 53, _) => 2.14,
        // C–H
        (1, 6, _) => 1.09,
        // N–H
        (1, 7, _) => 1.01,
        // O–H
        (1, 8, _) => 0.96,
        // Default
        _ => 1.54,
    }
}

/// Get ideal valence angle based on atom hybridization.
fn get_ideal_angle(mol: &Molecule, center: AtomIdx) -> f64 {
    let mut has_triple = false;
    let mut has_double_or_arom = false;

    for (_, bidx) in mol.neighbors(center) {
        match mol.bond(bidx).order {
            BondOrder::Triple => has_triple = true,
            BondOrder::Double | BondOrder::Aromatic => has_double_or_arom = true,
            _ => {}
        }
    }

    if has_triple {
        PI // 180°
    } else if has_double_or_arom {
        PI * 2.0 / 3.0 // 120°
    } else {
        109.47_f64.to_radians() // sp3
    }
}

// ---------------------------------------------------------------------------
// Constraint Assembly
// ---------------------------------------------------------------------------

/// Build constraint set from molecular topology.
///
/// Extracts:
/// - Bond constraints from all bonds (target = ideal distance)
/// - Angle constraints from all adjacent bond pairs (target = ideal angle)
pub fn build_constraints(mol: &Molecule) -> ConstraintSet {
    let mut bonds = Vec::new();
    let mut angles = Vec::new();

    // Bond constraints: one per bond
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let ideal_dist = get_ideal_bond_length(mol, a1, a2);
        bonds.push(BondConstraint::new(a1, a2, ideal_dist));
    }

    // Angle constraints: for each atom with ≥2 neighbors
    for center_idx in 0..mol.atom_count() {
        let center = AtomIdx(center_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();

        if neighbors.len() < 2 {
            continue;
        }

        let ideal_angle = get_ideal_angle(mol, center);

        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let a1 = neighbors[i];
                let a2 = neighbors[j];
                angles.push(AngleConstraint::new(a1, center, a2, ideal_angle));
            }
        }
    }

    ConstraintSet { bonds, angles }
}

// ---------------------------------------------------------------------------
// Geometric Utilities
// ---------------------------------------------------------------------------

/// Compute angle A—Center—B in radians.
pub(crate) fn compute_angle(coords: &Coords3D, a: AtomIdx, center: AtomIdx, b: AtomIdx) -> f64 {
    let pa = coords.get(a);
    let pc = coords.get(center);
    let pb = coords.get(b);

    let va = pa.sub(&pc);
    let vb = pb.sub(&pc);

    let na = va.norm();
    let nb = vb.norm();

    if na < 1e-10 || nb < 1e-10 {
        return 0.0;
    }

    let cos_angle = (va.dot(&vb) / (na * nb)).clamp(-1.0, 1.0);
    cos_angle.acos()
}

/// Return any unit vector perpendicular to `v`.
#[allow(dead_code)]
fn perpendicular_to(v: Point3) -> Point3 {
    let candidate = if v.x.abs() < 0.9 {
        Point3::new(1.0, 0.0, 0.0)
    } else {
        Point3::new(0.0, 1.0, 0.0)
    };
    let proj = v.scale(v.dot(&candidate));
    candidate.sub(&proj).normalize()
}

/// Rotate vector `v` around unit axis `axis` by angle `theta` (radians).
/// Uses Rodrigues' rotation formula.
pub(crate) fn rotate_around_axis(v: Point3, axis: Point3, theta: f64) -> Point3 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dot = axis.dot(&v);

    let term1 = v.scale(cos_t);
    let term2 = axis.cross(&v).scale(sin_t);
    let term3 = axis.scale(dot * (1.0 - cos_t));
    term1.add(&term2).add(&term3)
}

// ---------------------------------------------------------------------------
// Constraint Satisfaction via Projection
// ---------------------------------------------------------------------------

/// Project a single bond constraint: move atoms toward target distance.
///
/// Both atoms move symmetrically toward/away from their midpoint to achieve
/// the target distance. This preserves the centroid of the pair.
fn project_bond_constraint(coords: &mut Coords3D, constraint: &BondConstraint) {
    let p1 = coords.get(constraint.atom1);
    let p2 = coords.get(constraint.atom2);

    let current_dist = p1.distance(&p2);
    if current_dist < 1e-6 {
        return; // atoms coincident, can't project
    }

    let target_dist = constraint.target_distance;
    let lower = target_dist - constraint.tolerance;
    let upper = target_dist + constraint.tolerance;

    // Already satisfied?
    if current_dist >= lower && current_dist <= upper {
        return;
    }

    // Direction from p1 to p2
    let direction = p2.sub(&p1).scale(1.0 / current_dist);

    // Midpoint between atoms
    let mid = Point3::new(
        (p1.x + p2.x) / 2.0,
        (p1.y + p2.y) / 2.0,
        (p1.z + p2.z) / 2.0,
    );

    // Target distance: move to nearer boundary if out of tolerance
    let target_effective = if current_dist < lower { lower } else { upper };

    // New positions: symmetric movement from midpoint
    let offset = direction.scale(target_effective / 2.0);
    let new_p1 = mid.sub(&offset);
    let new_p2 = mid.add(&offset);

    coords.set(constraint.atom1, new_p1);
    coords.set(constraint.atom2, new_p2);
}

/// Project a single angle constraint: rotate one atom around a bond axis.
///
/// Rotates atom1 around the axis (center—atom2) to achieve the target angle.
/// atom2 is held fixed; atom1 moves.
fn project_angle_constraint(coords: &mut Coords3D, constraint: &AngleConstraint) {
    let p1 = coords.get(constraint.atom1);
    let center = coords.get(constraint.center);
    let p2 = coords.get(constraint.atom2);

    let v1 = p1.sub(&center);
    let v2 = p2.sub(&center);

    let n1 = v1.norm();
    let n2 = v2.norm();

    if n1 < 1e-10 || n2 < 1e-10 {
        return; // degenerate geometry
    }

    // Current angle
    let cos_angle = (v1.dot(&v2) / (n1 * n2)).clamp(-1.0, 1.0);
    let current_angle = cos_angle.acos();

    // Check if already satisfied
    let lower = constraint.target_angle - constraint.tolerance;
    let upper = constraint.target_angle + constraint.tolerance;

    if current_angle >= lower && current_angle <= upper {
        return;
    }

    // Axis of rotation: the bond direction center—atom2
    let axis = v2.normalize();

    // Rotation angle: delta toward target
    let delta_angle = constraint.target_angle - current_angle;

    // Rotate v1 around axis
    let v1_rotated = rotate_around_axis(v1, axis, delta_angle);
    let new_p1 = center.add(&v1_rotated);

    coords.set(constraint.atom1, new_p1);
}

/// Satisfy all constraints iteratively via geometric projection.
///
/// # Arguments
/// * `coords` - Initial coordinates
/// * `mol` - Molecule (for constraint context)
/// * `constraints` - Constraint set to satisfy
/// * `max_iterations` - Maximum projection iterations (typically 10-20)
///
/// # Returns
/// Coordinates with constraints satisfied (or as many as possible).
///
/// # Algorithm
/// For each iteration:
/// 1. Project all bond constraints (move atom pairs)
/// 2. Every other iteration: project angle constraints (avoid oscillation)
/// 3. Stop if all constraints satisfied or no progress made
pub fn satisfy_constraints(
    coords: &Coords3D,
    _mol: &Molecule,
    constraints: &ConstraintSet,
    max_iterations: usize,
) -> Coords3D {
    let mut result = coords.clone();

    for iteration in 0..max_iterations {
        let violation_before = constraints.violated_count(&result);

        // Project all bond constraints
        for constraint in &constraints.bonds {
            project_bond_constraint(&mut result, constraint);
        }

        // Project angle constraints every other iteration to avoid oscillation
        if iteration % 2 == 0 {
            for constraint in &constraints.angles {
                project_angle_constraint(&mut result, constraint);
            }
        }

        let violation_after = constraints.violated_count(&result);

        // Convergence: no constraints violated
        if violation_after == 0 {
            break;
        }

        // No progress: stop
        if iteration > 3 && (violation_before as i32 - violation_after as i32).abs() < 2 {
            break;
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_bond_constraint_creation() {
        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54);
        assert_eq!(constraint.target_distance, 1.54);
        assert_eq!(constraint.tolerance, 0.05);
    }

    #[test]
    fn test_bond_constraint_ethane_ideal_distance() {
        let mol = parse("CC").unwrap();
        let constraints = build_constraints(&mol);
        assert_eq!(constraints.bonds.len(), 1, "ethane has 1 bond");

        let bond = &constraints.bonds[0];
        assert!(
            (bond.target_distance - 1.54).abs() < 0.01,
            "C-C ideal ~1.54 Å"
        );
    }

    #[test]
    fn test_angle_constraint_creation() {
        let constraint =
            AngleConstraint::new(AtomIdx(0), AtomIdx(1), AtomIdx(2), 109.47_f64.to_radians());
        assert!((constraint.target_angle - 109.47_f64.to_radians()).abs() < 1e-6);
    }

    #[test]
    fn test_constraint_set_propane_angles() {
        let mol = parse("CCC").unwrap();
        let constraints = build_constraints(&mol);

        // Center atom (index 1) should have angle constraint with 2 neighbors
        let center_angles: Vec<_> = constraints
            .angles
            .iter()
            .filter(|a| a.center == AtomIdx(1))
            .collect();

        assert_eq!(
            center_angles.len(),
            1,
            "center of propane has 1 angle constraint"
        );
    }

    #[test]
    fn test_project_bond_constraint_too_far() {
        let _mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(5.0, 0.0, 0.0)); // too far: 5.0 Å

        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54);
        project_bond_constraint(&mut coords, &constraint);

        let d = coords.get(AtomIdx(0)).distance(&coords.get(AtomIdx(1)));
        // Should be near 1.54 (within tolerance + some projection error)
        assert!(
            d > 1.4 && d < 1.7,
            "projected distance {:.3}, expected ~1.54",
            d
        );
    }

    #[test]
    fn test_project_bond_constraint_too_close() {
        let _mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(0.5, 0.0, 0.0)); // too close: 0.5 Å

        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54);
        project_bond_constraint(&mut coords, &constraint);

        let d = coords.get(AtomIdx(0)).distance(&coords.get(AtomIdx(1)));
        assert!(
            d > 1.4 && d < 1.7,
            "projected distance {:.3}, expected ~1.54",
            d
        );
    }

    #[test]
    fn test_constraint_satisfaction_ethane() {
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(3.0, 0.0, 0.0));

        let constraints = build_constraints(&mol);
        let before_violations = constraints.violated_count(&coords);

        let satisfied = satisfy_constraints(&coords, &mol, &constraints, 10);
        let after_violations = constraints.violated_count(&satisfied);

        assert!(
            after_violations <= before_violations,
            "should reduce violations"
        );
    }

    #[test]
    fn test_constraint_set_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let constraints = build_constraints(&mol);

        assert_eq!(constraints.bonds.len(), 6, "benzene has 6 C-C bonds");
        assert_eq!(constraints.angles.len(), 6, "benzene has 6 C-C-C angles");

        // All aromatic bonds should have target ~1.40 Å
        for bond in &constraints.bonds {
            assert!((bond.target_distance - 1.40).abs() < 0.01);
        }
    }

    #[test]
    fn test_no_clashes_after_projection() {
        let mol = parse("CCCC").unwrap();
        let mut coords = Coords3D::new_zeroed(4);
        // Place in a line, too close
        for i in 0..4 {
            coords.set(AtomIdx(i as u32), Point3::new(i as f64 * 0.5, 0.0, 0.0));
        }

        let constraints = build_constraints(&mol);
        let satisfied = satisfy_constraints(&coords, &mol, &constraints, 20);

        // Check minimum distance
        let mut min_d = f64::MAX;
        for i in 0..4 {
            for j in (i + 1)..4 {
                let d = satisfied
                    .get(AtomIdx(i as u32))
                    .distance(&satisfied.get(AtomIdx(j as u32)));
                min_d = min_d.min(d);
            }
        }

        assert!(
            min_d > 0.1,
            "minimum distance {:.3}, atoms may clash",
            min_d
        );
    }

    #[test]
    fn test_compute_angle_90_degrees() {
        let coords = Coords3D::new_zeroed(3);
        let mut coords_mut = coords;
        coords_mut.set(AtomIdx(0), Point3::new(1.0, 0.0, 0.0));
        coords_mut.set(AtomIdx(1), Point3::new(0.0, 0.0, 0.0));
        coords_mut.set(AtomIdx(2), Point3::new(0.0, 1.0, 0.0));

        let angle = compute_angle(&coords_mut, AtomIdx(0), AtomIdx(1), AtomIdx(2));
        let expected = 90.0_f64.to_radians();
        assert!((angle - expected).abs() < 0.01);
    }

    #[test]
    fn test_bond_constraint_satisfied_true() {
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(1.54, 0.0, 0.0));

        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54);
        assert!(constraint.satisfied(&coords));
    }

    #[test]
    fn test_bond_constraint_satisfied_false() {
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(3.0, 0.0, 0.0));

        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54);
        assert!(!constraint.satisfied(&coords));
    }

    #[test]
    fn test_satisfy_constraints_single_bond() {
        use chematic_core::{Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        builder.add_atom(chematic_core::Atom::new(Element::C));
        builder.add_atom(chematic_core::Atom::new(Element::C));
        let _ = builder.add_bond(AtomIdx(0), AtomIdx(1), chematic_core::BondOrder::Single);
        let mol = builder.build();

        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(3.0, 0.0, 0.0)); // Far from ideal

        let target_dist = 1.54;
        let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), target_dist);
        let constraints = ConstraintSet {
            bonds: vec![constraint],
            angles: Vec::new(),
        };

        let result = satisfy_constraints(&coords, &mol, &constraints, 20);

        let dist = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(
            (dist - target_dist).abs() < 0.1,
            "Bond constraint should converge to target distance"
        );
    }

    #[test]
    fn test_satisfy_constraints_convergence() {
        use chematic_core::{Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        for _ in 0..3 {
            builder.add_atom(chematic_core::Atom::new(Element::C));
        }
        let _ = builder.add_bond(AtomIdx(0), AtomIdx(1), chematic_core::BondOrder::Single);
        let _ = builder.add_bond(AtomIdx(1), AtomIdx(2), chematic_core::BondOrder::Single);
        let mol = builder.build();

        let mut coords = Coords3D::new_zeroed(3);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(5.0, 0.0, 0.0));
        coords.set(AtomIdx(2), Point3::new(5.0, 5.0, 0.0));

        let constraints = ConstraintSet {
            bonds: vec![
                BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.54),
                BondConstraint::new(AtomIdx(1), AtomIdx(2), 1.54),
            ],
            angles: Vec::new(),
        };

        let result = satisfy_constraints(&coords, &mol, &constraints, 20);

        let dist01 = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        let dist12 = result.get(AtomIdx(1)).distance(&result.get(AtomIdx(2)));
        assert!((dist01 - 1.54).abs() < 0.15);
        assert!((dist12 - 1.54).abs() < 0.15);
    }
}
