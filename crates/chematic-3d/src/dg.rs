//! Rule-based 3D coordinate generation.
//!
//! This module implements a deterministic bond-angle-dihedral builder that
//! places heavy atoms in 3D space.  It is not a full distance-geometry solver;
//! it uses a DFS-based placement strategy with ideal bond lengths and ring
//! templates.
//!
//! Strategy:
//! 1. Find all connected components via BFS and process them independently.
//! 2. Within each component, detect rings using SSSR and place ring atoms on
//!    a regular polygon in the XY plane.
//! 3. Chain atoms are placed via DFS from a root, extending along directions
//!    chosen to approximate ideal bond angles and staggered dihedrals.
//! 4. Each component is offset along the X axis to avoid overlap.

use core::f64::consts::PI;
use std::collections::VecDeque;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Bond length lookup
// ---------------------------------------------------------------------------

/// Return the ideal bond length (angstroms) for the bond between atoms `a` and
/// `b` with the given bond order.
fn ideal_bond_len(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    let ea = mol.atom(a).element;
    let eb = mol.atom(b).element;

    // Retrieve the bond order between a and b.
    let order = mol
        .bond_between(a, b)
        .map(|(_, bond)| bond.order)
        .unwrap_or(BondOrder::Single);

    // Normalise element pair as (smaller_atomic_number, larger) for matching.
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

// ---------------------------------------------------------------------------
// Hybridisation / ideal angle
// ---------------------------------------------------------------------------

/// Rough estimate of the ideal bond angle at atom `center`, based on its
/// degree and bond orders.  Returns the angle in radians.
fn ideal_angle(mol: &Molecule, center: AtomIdx) -> f64 {
    // Check whether any neighbour bond is triple or double/aromatic.
    let mut has_triple = false;
    let mut has_double_or_arom = false;

    for (nb, bidx) in mol.neighbors(center) {
        let _ = nb;
        let order = mol.bond(bidx).order;
        match order {
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
        PI * 109.5_f64.to_radians() / PI // convert 109.5° to radians
    }
}

// ---------------------------------------------------------------------------
// Connected components
// ---------------------------------------------------------------------------

/// Return a list of connected components as atom-index lists.
fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut component: Vec<AtomIdx> = Vec::new();
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        let start_idx = AtomIdx(start as u32);
        visited[start] = true;
        queue.push_back(start_idx);

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (nb, _) in mol.neighbors(current) {
                if !visited[nb.0 as usize] {
                    visited[nb.0 as usize] = true;
                    queue.push_back(nb);
                }
            }
        }

        components.push(component);
    }

    components
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Generate 3D coordinates for all heavy atoms in `mol` using a rule-based
/// bond-angle-dihedral placement strategy.
///
/// The output coordinates are non-degenerate (no two atoms share the same
/// position) and approximate ideal bond lengths, but are not physically
/// minimised.
pub fn generate_coords(mol: &Molecule) -> Coords3D {
    let n = mol.atom_count();
    let mut coords = Coords3D::new_zeroed(n);

    if n == 0 {
        return coords;
    }

    // Single atom: place at origin.
    if n == 1 {
        coords.set(AtomIdx(0), Point3::zero());
        return coords;
    }

    // Detect rings.
    let ring_set = find_sssr(mol);

    // Process each connected component separately.
    let components = connected_components(mol);
    let mut x_offset = 0.0_f64;

    for component in &components {
        // BFS to determine component extent for the next offset.
        let placed = place_component(mol, component, &ring_set, x_offset, &mut coords);
        // Advance offset by max X extent of placed atoms + 5 Å gap.
        let max_x = placed
            .iter()
            .map(|&idx| coords.get(idx).x)
            .fold(f64::NEG_INFINITY, f64::max);
        x_offset = max_x + 5.0;
    }

    coords
}

// ---------------------------------------------------------------------------
// Component placement
// ---------------------------------------------------------------------------

/// Place all atoms in `component` starting at X = `x_offset`.
///
/// Returns the slice of atom indices that were placed (same as `component`).
fn place_component<'a>(
    mol: &Molecule,
    component: &'a [AtomIdx],
    ring_set: &chematic_perception::RingSet,
    x_offset: f64,
    coords: &mut Coords3D,
) -> &'a [AtomIdx] {
    if component.is_empty() {
        return component;
    }

    // Track which atoms have been placed.
    let total = mol.atom_count();
    let mut placed = vec![false; total];

    // First, lay out ring atoms onto polygon templates.
    place_rings(mol, component, ring_set, x_offset, coords, &mut placed);

    // Then extend non-ring atoms via DFS from the component root.
    let root = component[0];
    if !placed[root.0 as usize] {
        // No ring in this component — place root at x_offset.
        coords.set(root, Point3::new(x_offset, 0.0, 0.0));
        placed[root.0 as usize] = true;
    }

    // DFS to place remaining atoms.
    dfs_place(mol, root, &mut placed, coords);

    component
}

// ---------------------------------------------------------------------------
// Ring placement
// ---------------------------------------------------------------------------

/// Place ring atoms from SSSR onto regular polygon templates.
///
/// Each ring that contains atoms from `component` is laid out in the XY plane.
/// The first ring is centred at (`x_offset` + ring_radius, 0, 0).
/// Subsequent rings within the same component are fused to the previous ring
/// (sharing a bond edge) or offset if not fused.
fn place_rings(
    mol: &Molecule,
    component: &[AtomIdx],
    ring_set: &chematic_perception::RingSet,
    x_offset: f64,
    coords: &mut Coords3D,
    placed: &mut Vec<bool>,
) {
    let component_set: std::collections::HashSet<AtomIdx> =
        component.iter().copied().collect();

    let mut ring_cx = x_offset;
    let mut ring_cy = 0.0_f64;
    let mut first_ring = true;

    for ring in ring_set.rings() {
        // Only process rings whose atoms all belong to this component.
        if !ring.iter().all(|a| component_set.contains(a)) {
            continue;
        }

        let ring_size = ring.len();
        if ring_size == 0 {
            continue;
        }

        // Use bond length between consecutive ring atoms for the polygon side.
        let bond_len = {
            let a0 = ring[0];
            let a1 = ring[1 % ring_size];
            ideal_bond_len(mol, a0, a1)
        };

        // Circumradius of a regular polygon: r = bond_len / (2 * sin(PI / n)).
        let r = bond_len / (2.0 * (PI / ring_size as f64).sin());

        if first_ring {
            ring_cx = x_offset + r;
            ring_cy = 0.0;
            first_ring = false;
        } else {
            // Try to fuse to a previously-placed ring atom. If two atoms of
            // this ring are already placed, shift the centre to be consistent.
            let already_placed: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|a| placed[a.0 as usize])
                .collect();

            if already_placed.len() >= 2 {
                // Use the midpoint of the two most recently placed ring atoms.
                let p0 = coords.get(already_placed[0]);
                let p1 = coords.get(already_placed[1]);
                ring_cx = (p0.x + p1.x) / 2.0;
                ring_cy = (p0.y + p1.y) / 2.0 + r;
            } else if already_placed.len() == 1 {
                let p0 = coords.get(already_placed[0]);
                ring_cx = p0.x + r;
                ring_cy = p0.y;
            }
            // else keep ring_cx/ring_cy (floating ring, shouldn't happen in SSSR)
        }

        // Place each ring atom on the circle.
        for (k, &atom_idx) in ring.iter().enumerate() {
            if placed[atom_idx.0 as usize] {
                continue; // shared atom already placed by a previous ring
            }
            let angle = 2.0 * PI * k as f64 / ring_size as f64;
            let x = ring_cx + r * angle.cos();
            let y = ring_cy + r * angle.sin();
            coords.set(atom_idx, Point3::new(x, y, 0.0));
            placed[atom_idx.0 as usize] = true;
        }
    }
}

// ---------------------------------------------------------------------------
// DFS chain placement
// ---------------------------------------------------------------------------

/// DFS-based placement of atoms that have not yet been positioned.
///
/// For each unplaced neighbour of `current`, compute the ideal bond length,
/// choose a direction (bond angle from the incoming bond direction, with
/// dihedral rotated by 120° per successive neighbour to minimise clashes),
/// and recurse.
fn dfs_place(mol: &Molecule, current: AtomIdx, placed: &mut Vec<bool>, coords: &mut Coords3D) {
    let pos_current = coords.get(current);

    // Collect placed (parent-like) and unplaced neighbours.
    let mut placed_neighbors: Vec<AtomIdx> = Vec::new();
    let mut unplaced_neighbors: Vec<AtomIdx> = Vec::new();

    for (nb, _) in mol.neighbors(current) {
        if placed[nb.0 as usize] {
            placed_neighbors.push(nb);
        } else {
            unplaced_neighbors.push(nb);
        }
    }

    if unplaced_neighbors.is_empty() {
        return;
    }

    // Determine the "incoming" bond direction (from parent toward current).
    let incoming_dir: Point3 = if let Some(&parent) = placed_neighbors.first() {
        let p = coords.get(parent);
        pos_current.sub(&p).normalize()
    } else {
        // Root atom: pick +X direction.
        Point3::new(1.0, 0.0, 0.0)
    };

    // Choose a perpendicular reference vector for dihedral rotations.
    let perp = perpendicular_to(incoming_dir);

    let angle = ideal_angle(mol, current);

    for (i, &nb) in unplaced_neighbors.iter().enumerate() {
        let bond_len = ideal_bond_len(mol, current, nb);

        // Rotate incoming_dir by (PI - angle) around `perp`, then apply a
        // dihedral rotation of 120° * i around the incoming_dir axis.
        let bend_angle = PI - angle; // complement of bond angle
        let dir_bent = rotate_around_axis(incoming_dir, perp, bend_angle);

        let dihedral = (i as f64) * (2.0 * PI / 3.0); // 0°, 120°, 240°
        let dir_final = rotate_around_axis(dir_bent, incoming_dir, dihedral);

        let new_pos = pos_current.add(&dir_final.scale(bond_len));
        coords.set(nb, new_pos);
        placed[nb.0 as usize] = true;

        // Recurse.
        dfs_place(mol, nb, placed, coords);
    }
}

// ---------------------------------------------------------------------------
// Vector math helpers
// ---------------------------------------------------------------------------

/// Return any unit vector perpendicular to `v`.
fn perpendicular_to(v: Point3) -> Point3 {
    // Choose a candidate that is not parallel to v.
    let candidate = if v.x.abs() < 0.9 {
        Point3::new(1.0, 0.0, 0.0)
    } else {
        Point3::new(0.0, 1.0, 0.0)
    };
    // Gram-Schmidt: subtract projection of candidate onto v.
    let proj = v.scale(v.dot(&candidate));
    candidate.sub(&proj).normalize()
}

/// Rotate vector `v` around unit axis `axis` by angle `theta` (radians).
///
/// Uses Rodrigues' rotation formula:
///   v' = v*cos(theta) + (axis × v)*sin(theta) + axis*(axis·v)*(1 - cos(theta))
fn rotate_around_axis(v: Point3, axis: Point3, theta: f64) -> Point3 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dot = axis.dot(&v);

    // v*cos + (axis × v)*sin + axis*(dot)*(1 - cos)
    let term1 = v.scale(cos_t);
    let term2 = axis.cross(&v).scale(sin_t);
    let term3 = axis.scale(dot * (1.0 - cos_t));
    term1.add(&term2).add(&term3)
}
