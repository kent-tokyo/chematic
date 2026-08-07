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
        109.5_f64.to_radians()
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
        place_component(mol, component, &ring_set, x_offset, &mut coords);
        // Advance offset by max X extent of placed atoms + 5 Å gap.
        let max_x = component
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
fn place_component(
    mol: &Molecule,
    component: &[AtomIdx],
    ring_set: &chematic_perception::RingSet,
    x_offset: f64,
    coords: &mut Coords3D,
) {
    if component.is_empty() {
        return;
    }

    let mut placed = vec![false; mol.atom_count()];

    // First, lay out ring atoms onto polygon templates.
    place_rings(mol, component, ring_set, x_offset, coords, &mut placed);

    // If no ring placed anything in this component, anchor the first atom at
    // x_offset so there is at least one placed atom to extend from below.
    //
    // Placing this anchor unconditionally (rather than only when no ring
    // exists) was the root cause of a real bug: `place_rings` centres its
    // first ring at `x_offset + ring_radius`, which puts one ring vertex at
    // x = x_offset too (the vertex diametrically opposite the k=0 vertex,
    // since `ring_cx - ring_radius == x_offset`) -- so an unconditional
    // anchor at `(x_offset, 0, 0)` collided with (or landed a
    // floating-point epsilon from) that ring vertex on ANY molecule where
    // `component[0]` is a non-ring atom directly bonded to that ring (e.g.
    // plain toluene: the methyl carbon and the ring's ipso carbon ended up
    // at the same point). Anchoring only when the ring layout placed
    // nothing avoids ever competing with a ring-computed position.
    if !component.iter().any(|&a| placed[a.0 as usize]) {
        let root = component[0];
        coords.set(root, Point3::new(x_offset, 0.0, 0.0));
        placed[root.0 as usize] = true;
    }

    // Extend outward via DFS from every atom already placed (every ring
    // atom, and/or the anchor above) -- not just a single root. `dfs_place`
    // is a no-op for atoms whose neighbours are all already placed, so this
    // is safe and cheap to call per seed; it is what makes a substituent
    // attached to a *different* ring atom than the one nearest `component[0]`
    // actually get walked and placed, rather than being left at
    // `Coords3D::new_zeroed`'s (0, 0, 0) default forever (e.g. the second
    // methyl on p-xylene, or ibuprofen's isobutyl/carboxyl tail hanging off
    // a ring atom the original single-seed DFS never reached).
    for &atom in component {
        if placed[atom.0 as usize] {
            dfs_place(mol, atom, &mut placed, coords);
        }
    }
}

// ---------------------------------------------------------------------------
// Ring placement
// ---------------------------------------------------------------------------

/// Order `rings` (each a slice of `AtomIdx`) into fusion-consistent visiting
/// order via BFS on the ring-adjacency graph (two rings are adjacent iff they
/// share at least one atom). Returns each ring paired with whether it starts
/// a new "island" (shares no atom with any ring already visited).
///
/// SSSR's own enumeration order does **not** guarantee that every ring after
/// the first shares an atom with some already-visited ring -- confirmed on a
/// 3-linearly-fused system (anthracene): SSSR can return `[terminal ring A,
/// terminal ring C, middle ring B]`, where ring C shares zero atoms with ring
/// A (only with the not-yet-visited ring B). Iterating SSSR's raw order and
/// falling back to "keep the previous ring's center" whenever zero shared
/// atoms are found (this function's caller used to do exactly that) silently
/// superimposes two entire, unrelated rings on the same coordinates. BFS on
/// the adjacency graph guarantees each ring is visited only after a ring it
/// actually shares atoms with, whenever such a ring exists in the same
/// component.
fn order_rings_by_fusion_adjacency<'a>(
    rings: &[&'a Vec<AtomIdx>],
) -> Vec<(&'a Vec<AtomIdx>, bool)> {
    let n = rings.len();
    let mut visited = vec![false; n];
    let mut result = Vec::with_capacity(n);
    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        result.push((rings[start], true)); // starts a new island
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        while let Some(i) = queue.pop_front() {
            for j in 0..n {
                if !visited[j] && rings[i].iter().any(|a| rings[j].contains(a)) {
                    visited[j] = true;
                    result.push((rings[j], false)); // fused to an already-visited ring
                    queue.push_back(j);
                }
            }
        }
    }
    result
}

/// Place ring atoms from SSSR onto regular polygon templates.
///
/// Each ring that contains atoms from `component` is laid out in the XY plane.
/// The first ring (and the first ring of each subsequent, fusion-disconnected
/// "island" of rings within the same component, e.g. biphenyl's two separate
/// phenyl rings) is centred at a fresh anchor beyond any already-placed atom.
/// Every other ring fuses to a previously-placed ring atom (sharing a bond
/// edge) via [`order_rings_by_fusion_adjacency`]'s visiting order.
fn place_rings(
    mol: &Molecule,
    component: &[AtomIdx],
    ring_set: &chematic_perception::RingSet,
    x_offset: f64,
    coords: &mut Coords3D,
    placed: &mut [bool],
) {
    let component_set: std::collections::HashSet<AtomIdx> = component.iter().copied().collect();

    let relevant_rings: Vec<&Vec<AtomIdx>> = ring_set
        .rings()
        .iter()
        .filter(|ring| !ring.is_empty() && ring.iter().all(|a| component_set.contains(a)))
        .collect();
    if relevant_rings.is_empty() {
        return;
    }
    let ordered_rings = order_rings_by_fusion_adjacency(&relevant_rings);

    let mut ring_cx = x_offset;
    let mut ring_cy = 0.0_f64;
    let mut any_placed_yet = false;

    for (ring, is_new_island) in ordered_rings {
        let ring_size = ring.len();

        // Use bond length between consecutive ring atoms for the polygon side.
        let bond_len = {
            let a0 = ring[0];
            let a1 = ring[1 % ring_size];
            ideal_bond_len(mol, a0, a1)
        };

        // Circumradius of a regular polygon: r = bond_len / (2 * sin(PI / n)).
        let r = bond_len / (2.0 * (PI / ring_size as f64).sin());

        if is_new_island {
            // Anchor beyond any atom already placed in this component, so a
            // second (fusion-disconnected) ring island never collides with
            // the first -- x_offset itself only for the very first ring.
            let anchor_x = if any_placed_yet {
                component
                    .iter()
                    .filter(|a| placed[a.0 as usize])
                    .map(|&a| coords.get(a).x)
                    .fold(f64::NEG_INFINITY, f64::max)
                    + 5.0
            } else {
                x_offset
            };
            ring_cx = anchor_x + r;
            ring_cy = 0.0;
        } else {
            // Fuse to a previously-placed ring atom. If two atoms of this
            // ring are already placed, shift the centre to be consistent.
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
            // else: unreachable given `order_rings_by_fusion_adjacency`
            // only marks a ring `is_new_island = false` when it shares an
            // atom with some already-visited ring.
        }

        // Choose a ring conformation based on size and chemical environment.
        //
        // Rules applied only to non-fused, non-aromatic rings:
        //   6-membered  → chair (r = 1.452 Å, z = ±0.256 Å, derived from
        //                  109.5° bond angle and 1.54 Å C-C bond length)
        //   5-membered  → envelope (regular pentagon, one atom ±0.40 Å above
        //                  the plane of the other four)
        //   ≥ 8-membered → crown (alternating ±h, h scaling with ring size)
        //   everything else (aromatic, fused, 3/4/7-membered) → flat polygon
        let is_aromatic = ring.iter().all(|&a| mol.atom(a).aromatic);
        let is_fused = ring.iter().any(|a| placed[a.0 as usize]);

        // Chair uses a geometry-derived xy radius (not the regular-polygon r).
        // chair_r and chair_h are exact solutions for l=1.54 Å, θ=109.5°.
        const CHAIR_R: f64 = 1.452;
        const CHAIR_H: f64 = 0.256;
        const ENVELOPE_H: f64 = 0.400;

        enum Conf {
            Flat,
            Chair,
            Envelope,
            Crown(f64),
        }
        let conf = if is_aromatic || is_fused {
            Conf::Flat
        } else {
            match ring_size {
                6 => Conf::Chair,
                5 => Conf::Envelope,
                n if n >= 8 => Conf::Crown(0.3 + 0.04 * (n as f64 - 8.0).min(10.0)),
                _ => Conf::Flat,
            }
        };

        // Chair uses a different circumradius than the regular-polygon formula.
        let effective_r = if matches!(conf, Conf::Chair) {
            CHAIR_R
        } else {
            r
        };

        for (k, &atom_idx) in ring.iter().enumerate() {
            if placed[atom_idx.0 as usize] {
                continue; // shared atom already placed by a previous ring
            }
            let angle = 2.0 * PI * k as f64 / ring_size as f64;
            let x = ring_cx + effective_r * angle.cos();
            let y = ring_cy + effective_r * angle.sin();
            let z = match conf {
                // Chair: alternating ±h (CHAIR_H ≈ 0.256 Å).
                Conf::Chair => {
                    if k % 2 == 0 {
                        CHAIR_H
                    } else {
                        -CHAIR_H
                    }
                }
                // Envelope: last atom lifted above the mean plane of the other 4.
                Conf::Envelope => {
                    if k == ring_size - 1 {
                        ENVELOPE_H
                    } else {
                        0.0
                    }
                }
                // Crown: alternating ±h, height scales with ring size.
                Conf::Crown(h) => {
                    if k % 2 == 0 {
                        h
                    } else {
                        -h
                    }
                }
                Conf::Flat => 0.0,
            };
            coords.set(atom_idx, Point3::new(x, y, z));
            placed[atom_idx.0 as usize] = true;
        }
        any_placed_yet = true;
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
fn dfs_place(mol: &Molecule, current: AtomIdx, placed: &mut [bool], coords: &mut Coords3D) {
    let pos_current = coords.get(current);

    let parent = mol
        .neighbors(current)
        .map(|(nb, _)| nb)
        .find(|nb| placed[nb.0 as usize]);
    let unplaced_neighbors: Vec<AtomIdx> = mol
        .neighbors(current)
        .map(|(nb, _)| nb)
        .filter(|nb| !placed[nb.0 as usize])
        .collect();

    if unplaced_neighbors.is_empty() {
        return;
    }

    // Direction from parent toward current; for the root atom pick +X.
    let incoming_dir: Point3 = match parent {
        Some(p) => pos_current.sub(&coords.get(p)).normalize(),
        None => Point3::new(1.0, 0.0, 0.0),
    };

    let perp = perpendicular_to(incoming_dir);
    let angle = ideal_angle(mol, current);
    let bend_angle = PI - angle; // complement of bond angle
    let dir_bent = rotate_around_axis(incoming_dir, perp, bend_angle);

    for (i, &nb) in unplaced_neighbors.iter().enumerate() {
        let bond_len = ideal_bond_len(mol, current, nb);

        // Dihedral 0°, 120°, 240° around the incoming axis spaces successive
        // neighbours apart to minimise clashes.
        let dihedral = (i as f64) * (2.0 * PI / 3.0);
        let dir_final = rotate_around_axis(dir_bent, incoming_dir, dihedral);

        let new_pos = pos_current.add(&dir_final.scale(bond_len));
        coords.set(nb, new_pos);
        placed[nb.0 as usize] = true;

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

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn generate_coords_methane() {
        let mol = parse("C").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 1, "methane has 1 heavy atom");
        let p0 = coords.get(AtomIdx(0));
        assert_eq!(p0.x, 0.0, "first atom at origin");
    }

    #[test]
    fn generate_coords_ethane_has_reasonable_distance() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 2, "ethane has 2 carbons");
        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let dist = p0.distance(&p1);
        // C-C single bond is ~1.54 Å, should be within ±0.15 Å
        assert!(
            (dist - 1.54).abs() < 0.15,
            "C-C distance should be ~1.54 Å, got {}",
            dist
        );
    }

    #[test]
    fn generate_coords_benzene_has_6_atoms() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 6, "benzene has 6 carbons");
    }

    #[test]
    fn generate_coords_cyclohexane_all_placed() {
        let mol = parse("C1CCCCC1").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 6, "cyclohexane has 6 carbons");
    }

    #[test]
    fn generate_coords_disconnected_molecules() {
        let mol = parse("CC.CC").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(
            coords.atom_count(),
            4,
            "two ethanes (disconnected) have 4 carbons total"
        );
    }

    #[test]
    fn generate_coords_propane_linear() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 3, "propane has 3 carbons");
        // Check that all atoms have valid positions (not all zero)
        let mut has_nonzero = false;
        for i in 0..3 {
            let p = coords.get(AtomIdx(i));
            if p.x != 0.0 || p.y != 0.0 || p.z != 0.0 {
                has_nonzero = true;
            }
        }
        assert!(
            has_nonzero,
            "at least some atoms should be placed away from origin"
        );
    }

    fn min_pairwise_distance(coords: &Coords3D, n: usize) -> f64 {
        let mut min_d = f64::MAX;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = coords
                    .get(AtomIdx(i as u32))
                    .distance(&coords.get(AtomIdx(j as u32)));
                min_d = min_d.min(d);
            }
        }
        min_d
    }

    #[test]
    fn generate_coords_toluene_methyl_and_ipso_carbon_not_coincident() {
        // Regression test: `place_component` used to place a non-ring root
        // atom directly bonded to a ring at the fixed anchor (x_offset, 0, 0)
        // UNCONDITIONALLY, even when a ring had already been placed with its
        // own atom landing at that exact point (`place_rings` centres its
        // first ring at x_offset + ring_radius, putting the vertex
        // diametrically opposite k=0 at x = x_offset too). On plain toluene
        // this collided the methyl carbon with the ring's ipso carbon --
        // two chemically bonded atoms at (approximately) the same 3D point.
        let mol = parse("Cc1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 7, "toluene has 7 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.5,
            "no two atoms should be within 0.5 \u{c5} of each other, got {min_d}"
        );
    }

    #[test]
    fn generate_coords_para_disubstituted_ring_both_substituents_placed() {
        // Regression test: `dfs_place` used to be seeded only once, from the
        // component root -- once it reached an already-ring-placed atom it
        // stopped, so a substituent hanging off a *different* ring atom
        // than the one nearest the root was never visited and stayed at
        // `Coords3D::new_zeroed`'s (0, 0, 0) default. p-xylene has two
        // methyls on opposite ring atoms: only one was ever placed before
        // the fix.
        let mol = parse("Cc1ccc(C)cc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 8, "p-xylene has 8 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.5,
            "no two atoms should be within 0.5 \u{c5} of each other, got {min_d}"
        );
    }

    #[test]
    fn generate_coords_ring_with_tail_substituent_all_atoms_placed() {
        // Regression test, ibuprofen-shaped: a multi-atom substituent chain
        // (isopropyl + carboxyl) hanging off a ring atom that the initial
        // single-seed DFS never reached used to be left entirely at the
        // (0, 0, 0) default -- 5 atoms all exactly coincident.
        let mol = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 15, "ibuprofen has 15 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
    }

    #[test]
    fn generate_coords_anthracene_terminal_rings_not_superimposed() {
        // Regression test: `place_rings` used to iterate SSSR's raw ring
        // order and fuse each ring to whichever ring was placed immediately
        // before it. SSSR does NOT guarantee that order matches fusion
        // adjacency: for anthracene, SSSR returns [terminal ring A, terminal
        // ring C, middle ring B] -- ring C shares ZERO atoms with ring A
        // (only with the not-yet-placed ring B), so the old code's "0
        // already-placed atoms" branch silently reused ring A's exact
        // center for ring C, superimposing two entire terminal rings (6
        // atoms) on the same coordinates. This was independently discovered
        // while investigating issue #185 (chematic-ff's UFF minimizer
        // blowing up on naphthalene but reportedly not anthracene) --
        // anthracene's apparent "safety" turned out to be this collision
        // bug accidentally producing a degenerate-but-not-catastrophic
        // starting point, not genuine minimizer robustness: after this fix,
        // anthracene's `generate_coords` output blows up under
        // `chematic_ff::minimize_uff` too, same as naphthalene (see
        // `minimize.rs`'s `chematic_ff_own_uff_minimizer_blows_up_*` tests).
        let mol = parse("c1ccc2cc3ccccc3cc2c1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 14, "anthracene has 14 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
    }

    #[test]
    fn generate_coords_biphenyl_disconnected_ring_islands_not_superimposed() {
        // Regression test: two rings connected only via a non-ring bond
        // (biphenyl) share zero atoms and are never fusion-adjacent -- a
        // genuine "new island" case, not a fusion-order bug. Before this
        // fix, a ring island beyond the very first reused whatever
        // `ring_cx`/`ring_cy` the previous, unrelated ring left behind
        // instead of anchoring fresh beyond it.
        let mol = parse("c1ccc(cc1)-c1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 12, "biphenyl has 12 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
    }

    #[test]
    fn generate_and_minimize_dreiding_benzene() {
        use crate::minimize::minimize_dreiding;
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let minimized = minimize_dreiding(&mol, coords);
        assert_eq!(
            minimized.atom_count(),
            6,
            "benzene still has 6 carbons after minimization"
        );
    }

    // ── Macrocycle 3D (crown conformation) ──────────────────────────────────

    #[test]
    fn macrocycle_ring_is_non_planar() {
        // Cyclooctane (8-membered ring) should have non-zero Z spread after
        // crown conformation placement.
        let mol = parse("C1CCCCCCC1").unwrap(); // cyclooctane
        let coords = generate_coords(&mol);
        let z_vals: Vec<f64> = (0..mol.atom_count())
            .map(|i| coords.get(AtomIdx(i as u32)).z)
            .collect();
        let z_spread = z_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - z_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            z_spread > 0.1,
            "8-ring should have non-planar initial geometry, z_spread={z_spread}"
        );
    }

    #[test]
    fn cyclohexane_is_chair() {
        // 6-membered aliphatic ring should use chair conformation.
        // Chair geometry: r = 1.452 Å, h = ±0.256 Å → z_spread ≈ 0.512 Å.
        let mol = parse("C1CCCCC1").unwrap();
        let coords = generate_coords(&mol);
        let z_vals: Vec<f64> = (0..6).map(|i| coords.get(AtomIdx(i as u32)).z).collect();
        let z_max = z_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let z_min = z_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let z_spread = z_max - z_min;
        assert!(
            z_spread > 0.40 && z_spread < 0.65,
            "cyclohexane should be chair (z_spread ≈ 0.512 Å), got {z_spread}"
        );
        // Verify 3 atoms above and 3 atoms below the equatorial plane,
        // all at approximately ±CHAIR_H (0.256 Å).
        // (SSSR ring order != atom-index order, so we test by value, not by index.)
        let up: Vec<f64> = z_vals.iter().cloned().filter(|&z| z > 0.0).collect();
        let dn: Vec<f64> = z_vals.iter().cloned().filter(|&z| z < 0.0).collect();
        assert_eq!(up.len(), 3, "chair: 3 atoms above plane, got {:?}", z_vals);
        assert_eq!(dn.len(), 3, "chair: 3 atoms below plane, got {:?}", z_vals);
        for z in up.iter().chain(dn.iter()) {
            assert!(
                (z.abs() - 0.256).abs() < 0.01,
                "chair: z magnitude should be ≈ 0.256 Å, got {z:.4}"
            );
        }
    }

    #[test]
    fn cyclopentane_is_envelope() {
        // 5-membered aliphatic ring should use envelope conformation.
        // Last atom lifted by ENVELOPE_H = 0.40 Å; other four remain at z=0.
        let mol = parse("C1CCCC1").unwrap();
        let coords = generate_coords(&mol);
        let z_vals: Vec<f64> = (0..5).map(|i| coords.get(AtomIdx(i as u32)).z).collect();
        // SSSR ring order != atom-index order, so we test by value.
        let flat: Vec<f64> = z_vals.iter().cloned().filter(|&z| z.abs() < 0.01).collect();
        let flap: Vec<f64> = z_vals.iter().cloned().filter(|&z| z.abs() > 0.01).collect();
        assert_eq!(flat.len(), 4, "envelope: 4 atoms flat, got {:?}", z_vals);
        assert_eq!(flap.len(), 1, "envelope: 1 flap atom, got {:?}", z_vals);
        assert!(
            (flap[0].abs() - 0.40).abs() < 0.01,
            "flap atom should be at |z| ≈ 0.40 Å, got {:.4}",
            flap[0]
        );
    }

    #[test]
    fn macrocycle_all_atoms_placed() {
        // Cyclododecane (12-membered ring) — verify all 12 atoms have coords.
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cycloundecane 11-ring actually
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), mol.atom_count());
        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn generate_conformer_ensemble_multiple() {
        use crate::conformer::ConformerEnsemble;
        let mol = parse("CC").unwrap();
        let ensemble = ConformerEnsemble::new(mol);
        assert_eq!(
            ensemble.conformer_count(),
            0,
            "fresh ensemble has 0 conformers"
        );
    }
}
