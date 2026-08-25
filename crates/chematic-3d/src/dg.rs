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
//!
//! [`crate::dg_connectivity_ordered`] is a separate, independently-selectable
//! placement engine (issues #256/#255) living in its own module -- this
//! module's own `generate_coords` and its helpers below are completely
//! unaffected by it and untouched by that engine's own development.

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
///
/// `pub(crate)`: also used by [`crate::dg_connectivity_ordered`].
pub(crate) fn ideal_bond_len(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
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
///
/// `pub(crate)`: also used by [`crate::dg_connectivity_ordered`].
pub(crate) fn ideal_angle(mol: &Molecule, center: AtomIdx) -> f64 {
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
///
/// `pub(crate)`: also used by [`crate::dg_connectivity_ordered`].
pub(crate) fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
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

/// Order `rings` into a visiting order via BFS on the ring-adjacency graph,
/// where two rings are adjacent iff they share an atom (fused, e.g.
/// naphthalene) OR have a direct bond between some atom in one and some
/// atom in the other (e.g. biphenyl's two phenyls, bonded but sharing no
/// atom). Whether a *specific* ring shares an atom with whatever is
/// already placed by the time [`place_rings`] gets to it is re-checked
/// there against live placement state, not decided here -- this function's
/// only job is visiting order.
///
/// Neither adjacency alone is enough, and SSSR's own enumeration order
/// guarantees neither kind of ordering on its own:
/// - Atom-sharing only: confirmed on a 3-linearly-fused system (anthracene)
///   -- SSSR can return `[terminal ring A, terminal ring C, middle ring
///   B]`, where ring C shares zero atoms with ring A (only with the
///   not-yet-visited ring B). Falling back to "keep the previous ring's
///   center" whenever zero shared atoms are found (an earlier version of
///   this function's caller did exactly that) silently superimposes two
///   entire, unrelated rings on the same coordinates.
/// - Bond adjacency also needed: confirmed on a 3-ring direct-bond chain
///   (terphenyl, rings connected by single bonds, sharing no atoms at
///   all) -- SSSR returned `[ring1, ring3, ring2]`, i.e. the ring bonded
///   to NEITHER already-placed ring (ring3, only bonded to ring2) ahead of
///   the ring that actually connects the chain together (ring2).
///   Visiting ring3 before ring2 forces it onto an anchor with no relation
///   to ring2's real position; when ring2 is later placed and (correctly)
///   anchored via its real bond to ring3, its OTHER real bond (to ring1)
///   comes out wrong (measured: 0.66 Å, vs. an ideal ~1.5 Å single bond)
///   -- one real constraint satisfied while an earlier, arbitrarily placed
///   ring made the other unsatisfiable.
///
/// BFS over the combined adjacency guarantees each ring is visited only
/// after some ring it can be positioned relative to (by either means),
/// whenever such a ring exists in the same component.
fn order_rings_by_fusion_adjacency<'a>(
    mol: &Molecule,
    rings: &[&'a Vec<AtomIdx>],
) -> Vec<&'a Vec<AtomIdx>> {
    let n = rings.len();
    let adjacent = |i: usize, j: usize| -> bool {
        rings[i].iter().any(|a| rings[j].contains(a))
            || rings[i]
                .iter()
                .any(|&a| mol.neighbors(a).any(|(nb, _)| rings[j].contains(&nb)))
    };
    let mut visited = vec![false; n];
    let mut result = Vec::with_capacity(n);
    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        result.push(rings[start]);
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        while let Some(i) = queue.pop_front() {
            for j in 0..n {
                if !visited[j] && adjacent(i, j) {
                    visited[j] = true;
                    result.push(rings[j]);
                    queue.push_back(j);
                }
            }
        }
    }
    result
}

/// Place ring atoms from SSSR onto regular polygon templates.
///
/// Each ring that contains atoms from `component` is laid out in the XY
/// plane, visited in [`order_rings_by_fusion_adjacency`]'s order. For each
/// ring, in live placement order: if it shares an atom with something
/// already placed, fuse to that shared atom (or the midpoint of two, if
/// two are already placed). Otherwise, if it has a direct bond to an
/// already-placed atom (e.g. biphenyl's two phenyls, connected but sharing
/// no atom), anchor via that bond's real length, extending away from the
/// connected atom's own ring. Otherwise (the very first ring in the
/// component, or a ring reachable only through not-yet-placed chain atoms
/// -- see the "no direct bond" branch below), anchor at a fresh position
/// beyond any already-placed atom in the component.
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
    let ordered_rings = order_rings_by_fusion_adjacency(mol, &relevant_rings);

    let mut any_placed_yet = false;

    for ring in ordered_rings {
        let ring_size = ring.len();
        // Always assigned in every branch below before being read.
        let ring_cx: f64;
        let ring_cy: f64;

        // Use bond length between consecutive ring atoms for the polygon side.
        let bond_len = {
            let a0 = ring[0];
            let a1 = ring[1 % ring_size];
            ideal_bond_len(mol, a0, a1)
        };

        // Circumradius of a regular polygon: r = bond_len / (2 * sin(PI / n)).
        let r = bond_len / (2.0 * (PI / ring_size as f64).sin());

        // Placement strategy is decided from LIVE `placed` state here, not
        // from `order_rings_by_fusion_adjacency`'s visiting order: that
        // function only orders rings so that *some* already-placed
        // neighbour (by either adjacency kind) exists by the time a ring
        // is reached -- it does not, and should not, commit to which kind
        // applies, since that can only be known once earlier rings in the
        // same BFS traversal have actually been placed.
        let shared_atoms: Vec<AtomIdx> = ring
            .iter()
            .copied()
            .filter(|a| placed[a.0 as usize])
            .collect();

        if shared_atoms.is_empty() {
            // No shared ring atom -- may still have a real bond straight to
            // an already-placed atom -- biphenyl's two phenyl rings share
            // zero atoms but ARE directly bonded to each other. Anchoring
            // blindly at a fixed offset ignored that bond entirely and
            // left it stretched to whatever the offset was (measured:
            // exactly 5.0 Å for biphenyl, vs. an ideal ~1.4 Å aromatic C-C
            // single bond) -- both ring endpoints were already marked
            // `placed`, so `dfs_place`'s chain walk (which only visits
            // *unplaced* neighbours) never got a chance to correct it.
            // Look for such a bond first and, if found, anchor via its
            // real ideal length instead.
            let direct_bond = ring.iter().find_map(|&ring_atom| {
                mol.neighbors(ring_atom)
                    .map(|(nb, _)| nb)
                    .find(|nb| placed[nb.0 as usize])
                    .map(|anchor_atom| (ring_atom, anchor_atom))
            });

            if let Some((ring_atom, anchor_atom)) = direct_bond {
                let anchor_pos = coords.get(anchor_atom);
                // Extend away from the centroid of `anchor_atom`'s OWN
                // already-placed ring specifically, not blindly along +X
                // and not the whole component's centroid either: +X only
                // happens to work when the connecting atom sits on the
                // "outer" side of its own ring (biphenyl's own para-like
                // case), and a whole-component centroid still misleads a
                // 3+-ring chain (terphenyl) once an earlier ring has
                // already skewed the average away from the specific ring
                // being extended from. `place_rings` runs before any
                // `dfs_place` chain walk, so every already-placed atom here
                // is necessarily a ring atom -- `anchor_atom` always
                // belongs to exactly one prior entry in `relevant_rings`.
                // Measured without this: 3-phenylpyridine collapsed to
                // 0.14 Å min pairwise distance (whole-component-centroid
                // version), terphenyl's third ring still did too even with
                // it (component-wide average, not this ring's own centroid).
                let centroid_xy = relevant_rings
                    .iter()
                    .find(|other_ring| other_ring.contains(&anchor_atom))
                    .map(|other_ring| {
                        let (sx, sy, n) =
                            other_ring
                                .iter()
                                .fold((0.0, 0.0, 0.0_f64), |(sx, sy, n), &a| {
                                    let p = coords.get(a);
                                    (sx + p.x, sy + p.y, n + 1.0)
                                });
                        (sx / n, sy / n)
                    })
                    .unwrap_or((anchor_pos.x, anchor_pos.y));
                let dx = anchor_pos.x - centroid_xy.0;
                let dy = anchor_pos.y - centroid_xy.1;
                let dist = (dx * dx + dy * dy).sqrt();
                let (ux, uy) = if dist > 1e-6 {
                    (dx / dist, dy / dist)
                } else {
                    (1.0, 0.0)
                };
                let bond_len_to_ring = ideal_bond_len(mol, ring_atom, anchor_atom);
                let ring_atom_x = anchor_pos.x + ux * bond_len_to_ring;
                let ring_atom_y = anchor_pos.y + uy * bond_len_to_ring;
                let k = ring.iter().position(|&a| a == ring_atom).unwrap();
                let angle = 2.0 * PI * k as f64 / ring_size as f64;
                ring_cx = ring_atom_x - r * angle.cos();
                ring_cy = ring_atom_y - r * angle.sin();
            } else {
                // No direct bond to any already-placed atom -- this ring is
                // reachable, if at all within this component, only through
                // atoms `place_rings` hasn't placed yet (chain-bridged ring
                // islands, e.g. two phenyls linked by a -CH2CH2- bridge;
                // `place_rings` places every ring before any chain atom is
                // walked, so that bridge's own length can't be known here).
                // Known limitation, not fixed by this anchor: the fixed
                // +5 Å offset below only guarantees no collision, not a
                // correct bond length at the eventual chain-to-ring
                // junction.
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
            }
        } else {
            // Fuse to a previously-placed ring atom (shares >=1 atom with
            // something already placed). If two atoms of this ring are
            // already placed, shift the centre to be consistent.
            if shared_atoms.len() >= 2 {
                // Use the midpoint of the two most recently placed ring atoms.
                let p0 = coords.get(shared_atoms[0]);
                let p1 = coords.get(shared_atoms[1]);
                ring_cx = (p0.x + p1.x) / 2.0;
                ring_cy = (p0.y + p1.y) / 2.0 + r;
            } else {
                // shared_atoms.len() == 1, guaranteed nonempty by the
                // `else` branch of `if shared_atoms.is_empty()` above.
                let p0 = coords.get(shared_atoms[0]);
                ring_cx = p0.x + r;
                ring_cy = p0.y;
            }
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
///
/// `pub(crate)`: also used by [`crate::dg_connectivity_ordered`].
pub(crate) fn perpendicular_to(v: Point3) -> Point3 {
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
///
/// `pub(crate)`: also used by [`crate::dg_connectivity_ordered`].
pub(crate) fn rotate_around_axis(v: Point3, axis: Point3, theta: f64) -> Point3 {
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

    /// Asserts every BONDED pair in `mol` is finite and within
    /// `[min_len, max_len]` Å. Deliberately distinct from
    /// `min_pairwise_distance`, which only checks the single closest pair
    /// over ALL atoms (bonded or not) -- it can miss a specific bonded
    /// pair landing too far apart (a stretch only ever *increases*
    /// distances, so it can't move the global minimum) while some
    /// unrelated, non-bonded pair happens to be close. Panics with the
    /// offending atom indices and distance so a failure is directly
    /// actionable.
    fn assert_bonded_pairs_sane(mol: &Molecule, coords: &Coords3D, min_len: f64, max_len: f64) {
        for (_, b) in mol.bonds() {
            let p1 = coords.get(b.atom1);
            let p2 = coords.get(b.atom2);
            let d = p1.distance(&p2);
            assert!(
                d.is_finite() && d >= min_len && d <= max_len,
                "bond {}-{} has length {d:.4} \u{c5}, expected finite and within [{min_len}, {max_len}] \u{c5}",
                b.atom1.0,
                b.atom2.0
            );
        }
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
        // NOT extended with `assert_bonded_pairs_sane` like the
        // biphenyl/terphenyl/meta-linked tests below: doing so here (during
        // this PR's review) surfaced FOUR distorted bonds at this molecule's
        // ring-fusion seams -- 0.8675 Å, 1.3163 Å, 2.0365 Å, 2.2644 Å,
        // not just the single 2.2644 Å this file's own commit history
        // previously (incorrectly) described as the only imperfection.
        // Confirmed by direct comparison (temporarily disabling this PR's
        // new bond-adjacency ordering edge and re-running) that these come
        // from the PRE-EXISTING `shared_atoms.len() >= 2` fuse branch above
        // -- its `ring_cy = (p0.y + p1.y) / 2.0 + r` always extends the new
        // ring in +y from the fusion bond's midpoint, regardless of which
        // direction is actually away from the rest of the structure. This
        // predates this PR entirely (shipped in `d45f91b`, this test's own
        // prior form never checked bond lengths, only
        // `min_pairwise_distance`) and is a different code path from the
        // "new island" direct-bond anchor this PR's biphenyl/terphenyl fix
        // touches -- not fixed here, flagged for a separate decision.
    }

    #[test]
    fn generate_coords_biphenyl_disconnected_ring_islands_not_superimposed() {
        // Regression test: two rings connected only via a non-ring bond
        // (biphenyl) share zero atoms and are never fusion-adjacent -- a
        // genuine "new island" case, not a fusion-order bug. Before this
        // fix, a ring island beyond the very first reused whatever
        // `ring_cx`/`ring_cy` the previous, unrelated ring left behind
        // instead of anchoring fresh beyond it.
        //
        // The FIRST fix for this (anchor the new island at a fixed offset
        // beyond the rest of the component) traded that collision for a
        // different bug: it ignored the real single bond directly
        // connecting biphenyl's two rings entirely, stretching it to
        // exactly that fixed offset (measured: 5.0 Å, vs. an ideal ~1.5 Å
        // single bond) -- both endpoints were already `placed`, so
        // `dfs_place`'s chain walk (unplaced-neighbours only) never got a
        // chance to correct it. `min_pairwise_distance` alone could not
        // have caught this: a stretch only increases distances, so it
        // never becomes the global minimum. Fixed by detecting a direct
        // bond to an already-placed atom before falling back to the
        // fixed-offset anchor.
        let mol = parse("c1ccc(cc1)-c1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 12, "biphenyl has 12 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
        assert_bonded_pairs_sane(&mol, &coords, 1.0, 1.8);
    }

    #[test]
    fn generate_coords_terphenyl_chain_of_new_island_rings_bond_lengths_sane() {
        // Regression test: a chain of 3+ rings, each directly bonded to the
        // next with zero shared atoms (so every ring after the first is a
        // "new island" per `order_rings_by_fusion_adjacency`), needs each
        // new-island anchor to extend away from the SPECIFIC ring it is
        // bonding to, not from the whole component's running centroid.
        // Using the whole-component centroid still collapsed the third
        // ring here (measured: 0.14 \u{c5} min pairwise distance) even after
        // fixing biphenyl's simpler 2-ring case, because by the time the
        // third ring is placed the average of rings 1+2 no longer points
        // "outward" from ring 2's own extent. Fixed by finding the
        // specific already-placed ring containing the connecting atom and
        // using only its own centroid.
        let mol = parse("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 18, "terphenyl has 18 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
        assert_bonded_pairs_sane(&mol, &coords, 1.0, 1.8);
    }

    #[test]
    fn generate_coords_meta_linked_biaryl_bond_lengths_sane() {
        // Regression test: when the connecting ring atom sits on the side
        // of its ring FACING the already-placed structure (meta-linked,
        // unlike biphenyl's own para-like case), a fixed +X extension
        // direction points the new ring straight back into what's already
        // placed. Measured without the centroid-outward fix: 0.14 \u{c5}
        // min pairwise distance on this exact molecule.
        let mol = parse("c1ccc(cc1)-c1cccnc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 12, "3-phenylpyridine has 12 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
        assert_bonded_pairs_sane(&mol, &coords, 1.0, 1.8);
    }

    #[test]
    fn generate_coords_spiro_ring_adjacency_unaffected() {
        // Positive control for the two fixes above: spiro rings share
        // EXACTLY ONE atom, so `order_rings_by_fusion_adjacency` correctly
        // marks the second ring `is_new_island = false` (fused, via the
        // existing shared-atom branch this PR does not touch) -- confirms
        // the new "new island direct-bond" logic never fires for spiro
        // systems and doesn't regress them.
        let mol = parse("C1CCC2(CC1)CCCCC2").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 11, "spiro[5.5]undecane has 11 heavy atoms");
        let coords = generate_coords(&mol);
        let min_d = min_pairwise_distance(&coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
        assert_bonded_pairs_sane(&mol, &coords, 1.0, 1.8);
    }

    #[test]
    fn generate_coords_bibenzyl_chain_bridged_ring_islands_known_broken() {
        // NOT a regression test for a fix -- pins a KNOWN, still-open
        // limitation. `place_rings` places every ring in a component
        // before `dfs_place` walks any chain atom, so a ring reachable
        // from an already-placed ring only through a non-ring bridge (here
        // the -CH2CH2- of bibenzyl, PhCH2CH2Ph) has no already-placed atom
        // for the direct-bond anchor added in this PR to find -- it falls
        // through to the fixed +5 \u{c5}-offset anchor, which guarantees no
        // collision but not a correct bond length at the eventual
        // chain-to-ring junction. Fixing this needs `place_component`
        // restructured to interleave ring placement and chain DFS in true
        // graph order, not a targeted fix within `place_rings` alone --
        // out of scope here. This test exists so that future restructuring
        // has a ready-made regression fixture: it currently pins the
        // broken value so a fix is provable (this assertion should FAIL
        // once the underlying limitation is actually fixed, at which point
        // replace it with a `assert_bonded_pairs_sane` call instead of
        // deleting it).
        let mol = parse("c1ccccc1CCc1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 14, "bibenzyl has 14 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 5.0,
            "expected this known-broken chain-bridged case to still have a grossly \
             stretched bond (last measured 8.7358 \u{c5}); got worst bond {worst:.4} \u{c5} -- \
             if this now passes, the chain-bridged ring-island limitation documented above \
             was fixed and this test should be replaced with a sane-bond-length assertion"
        );
    }

    #[test]
    fn generate_coords_naphthalene_fusion_seam_known_broken() {
        // NOT a regression test for a fix -- pins a KNOWN, still-open
        // limitation (issue #255): the simplest possible 2-ring fusion
        // (exactly one shared edge) still hits `place_rings`'
        // `shared_atoms.len() >= 2` branch's fixed `+y` extension from the
        // fusion bond's midpoint, which is only correct when `+y` happens
        // to point away from the already-placed ring -- not guaranteed,
        // and not correct here. Genuinely absent from this file's test
        // suite before this addition; the file's only prior fused-ring
        // test (anthracene, 3 rings) deliberately avoided a bonded-pair
        // check for the same underlying reason (see that test's own
        // comment). This is #256 Phase 0 fixture work: pin naphthalene's
        // current broken value so the design doc's claim that a
        // connectivity-ordered rewrite should also close #255 is
        // checkable, not just asserted.
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 10, "naphthalene has 10 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 2.0,
            "expected this known-broken 2-ring fusion case to still have a distorted \
             bond (last measured 2.2644 \u{c5} vs. ~1.40 \u{c5} ideal aromatic C-C); got \
             worst bond {worst:.4} \u{c5} -- if this now passes, issue #255's fusion-seam \
             limitation was fixed and this test should be replaced with an \
             `assert_bonded_pairs_sane` call instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_quinoline_fused_heterocycle_known_broken() {
        // Same root cause and same measured value as naphthalene above
        // (identical 2-ring, 1-shared-edge fusion topology; the fused
        // pyridine nitrogen doesn't change which `place_rings` branch
        // fires) -- included as its own fixture, not merged into the
        // naphthalene test, so a future fix's coverage of at least one
        // fused *heterocycle* is checked explicitly, not inferred from an
        // all-carbon case.
        let mol = parse("c1ccc2ncccc2c1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 10, "quinoline has 10 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 2.0,
            "expected this known-broken fused-heterocycle case to still have a distorted \
             bond (last measured 2.2644 \u{c5}); got worst bond {worst:.4} \u{c5} -- if this \
             now passes, issue #255's fusion-seam limitation was fixed and this test \
             should be replaced with an `assert_bonded_pairs_sane` call instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_phenanthrene_angular_fusion_known_broken() {
        // NOT a regression test for a fix -- pins a KNOWN, still-open
        // limitation (issue #255), same fixed-`+y` fusion-seam branch as
        // naphthalene above, but ANGULAR 3-ring fusion rather than linear
        // (unlike anthracene). This specific molecule matters beyond
        // general coverage: issue #255's own history records a reverted
        // overnight fix attempt (a 2-point rigid rotation of the ring
        // template onto the real shared-edge positions) that fixed
        // anthracene exactly but broke phenanthrene with two EXACT atom
        // coincidences plus a 3.7 \u{c5} stretch -- worse than today's
        // distortion. This fixture exists so any future fix attempt is
        // checked against phenanthrene from the start, not discovered to
        // regress it after the fact.
        let mol = parse("c1ccc2c(c1)ccc1ccccc12").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 14, "phenanthrene has 14 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 3.0,
            "expected this known-broken angular-fusion case to still have a distorted \
             bond (last measured 3.3856 \u{c5}); got worst bond {worst:.4} \u{c5} -- if this \
             now passes, issue #255's fusion-seam limitation was fixed for angular fusion \
             and this test should be replaced with an `assert_bonded_pairs_sane` call \
             instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_pyrene_multi_ring_fusion_known_broken() {
        // NOT a regression test for a fix -- pins a KNOWN, still-open
        // limitation (issue #255), same fixed-`+y` fusion-seam branch,
        // 4-ring fusion. Also part of #255's reverted-fix history: the
        // same overnight attempt that broke phenanthrene also produced an
        // exact atom coincidence here. Highest-multiplicity fusion case
        // in this file's fixture set.
        let mol = parse("c1cc2ccc3cccc4ccc(c1)c2c34").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 16, "pyrene has 16 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 3.5,
            "expected this known-broken multi-ring-fusion case to still have a distorted \
             bond (last measured 3.8974 \u{c5}); got worst bond {worst:.4} \u{c5} -- if this \
             now passes, issue #255's fusion-seam limitation was fixed for multi-ring \
             fusion and this test should be replaced with an `assert_bonded_pairs_sane` \
             call instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_diphenylmethane_chain_bridge_length_1_known_broken() {
        // NOT a regression test for a fix -- pins a KNOWN, still-open
        // limitation (issue #256), same root cause as the bibenzyl test
        // below but at the SHORTEST possible ring-chain-ring bridge (a
        // single -CH2- atom, vs. bibenzyl's -CH2CH2-). Confirms #256's
        // bug isn't specific to a 2-atom bridge -- even the minimal case
        // still finds no already-placed anchor for `place_rings`' direct-
        // bond check, since the ONE bridging atom is itself unplaced at
        // `place_rings` time.
        let mol = parse("c1ccccc1Cc1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 13, "diphenylmethane has 13 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 6.0,
            "expected this known-broken chain-bridged case (bridge length 1) to still have \
             a grossly stretched bond (last measured 8.0738 \u{c5}); got worst bond \
             {worst:.4} \u{c5} -- if this now passes, issue #256's chain-bridged \
             ring-island limitation was fixed and this test should be replaced with an \
             `assert_bonded_pairs_sane` call instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_diphenylpropane_chain_bridge_length_3_known_broken() {
        // Same root cause as bibenzyl (bridge length 2) and the
        // diphenylmethane test above (bridge length 1), at bridge length
        // 3 (-CH2CH2CH2-) -- part of a bridge-length series (1, 2, 3, 4)
        // so a future fix's differential evaluation can check the fix
        // generalizes across bridge length, not just bibenzyl's single
        // data point.
        let mol = parse("c1ccccc1CCCc1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 15, "1,3-diphenylpropane has 15 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 6.0,
            "expected this known-broken chain-bridged case (bridge length 3) to still have \
             a grossly stretched bond (last measured 8.9731 \u{c5}); got worst bond \
             {worst:.4} \u{c5} -- if this now passes, issue #256's chain-bridged \
             ring-island limitation was fixed and this test should be replaced with an \
             `assert_bonded_pairs_sane` call instead of deleting it"
        );
    }

    #[test]
    fn generate_coords_diphenylbutane_chain_bridge_length_4_known_broken() {
        // Same root cause as bibenzyl (bridge length 2) and the two tests
        // above (bridge lengths 1, 3), at bridge length 4
        // (-CH2CH2CH2CH2-) -- completes the 1/2(bibenzyl)/3/4 bridge-length
        // series.
        let mol = parse("c1ccccc1CCCCc1ccccc1").unwrap();
        let n = mol.atom_count();
        assert_eq!(n, 16, "1,4-diphenylbutane has 16 heavy atoms");
        let coords = generate_coords(&mol);
        let worst = mol
            .bonds()
            .map(|(_, b)| coords.get(b.atom1).distance(&coords.get(b.atom2)))
            .fold(0.0_f64, f64::max);
        assert!(
            worst > 6.0,
            "expected this known-broken chain-bridged case (bridge length 4) to still have \
             a grossly stretched bond (last measured 8.7060 \u{c5}); got worst bond \
             {worst:.4} \u{c5} -- if this now passes, issue #256's chain-bridged \
             ring-island limitation was fixed and this test should be replaced with an \
             `assert_bonded_pairs_sane` call instead of deleting it"
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
