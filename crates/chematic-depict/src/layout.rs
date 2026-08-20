//! 2D coordinate generation for molecular depiction.
//!
//! The layout algorithm is rule-based and produces SVG pixel coordinates.
//! No physics simulation is used; atoms are placed with geometric rules.
//!
//! ## Algorithm Summary
//!
//! 1. **Ring detection**: Find SSSR (Smallest Set of Smallest Rings) via Balducci-Pearlman.
//! 2. **Ring placement**: Place each ring as a regular polygon; fused rings reflect new atoms over shared edges.
//! 3. **Chain placement**: Use DFS zigzag to place chain atoms from ring atoms or arbitrary start.
//! 4. **Fragment spacing**: Offset disconnected components horizontally to prevent overlap.
//! 5. **Collision detection**: `detect_crossings()` reports bond–bond intersections (for UI feedback).
//!
//! The algorithm prioritizes clarity (minimal crossing) over perfect physics simulation.
//! Bond angles follow tetrahedral/trigonal rules where possible.

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, BondIdx, Molecule};
use chematic_perception::find_sssr;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Bond length in SVG pixels. Scales all ring radii and chain steps.
pub const BOND_LEN: f64 = 40.0;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A 2D point in SVG coordinate space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Euclidean distance to `other`.
    pub fn dist(&self, other: &Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx * dx + dy * dy).sqrt()
    }
}

/// 2D layout result: one coordinate per atom indexed by `AtomIdx`.
pub struct Layout {
    pub coords: Vec<Point>,
}

impl Layout {
    /// Get the coordinate of atom `idx`.
    pub fn get(&self, idx: AtomIdx) -> Point {
        self.coords[idx.0 as usize]
    }

    /// Bounding box: (min_x, min_y, max_x, max_y).
    pub fn bounding_box(&self) -> (f64, f64, f64, f64) {
        if self.coords.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }
        self.coords.iter().fold(
            (f64::MAX, f64::MAX, f64::MIN, f64::MIN),
            |(min_x, min_y, max_x, max_y), p| {
                (
                    min_x.min(p.x),
                    min_y.min(p.y),
                    max_x.max(p.x),
                    max_y.max(p.y),
                )
            },
        )
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Compute a 2D layout for `mol`, returning one `Point` per atom.
///
/// Algorithm overview:
/// 1. Find ring systems (SSSR). Group rings sharing at least one atom into ring systems.
/// 2. Place each ring system: regular polygon for a single ring; fused rings are placed
///    by reflecting new atoms over the shared bond.
/// 3. Place chain atoms (not in any ring) via DFS zigzag from each ring atom or
///    from an arbitrary starting point if no rings exist.
/// 4. Offset disconnected fragments horizontally so they do not overlap.
///
/// Compute 2D layout coordinates for a molecule.
///
/// **Coordinate system:** SVG pixel space with **Y-axis pointing downward**
/// (Y increases toward the screen bottom). All returned coordinates conform to this convention.
/// This is consistent with standard SVG/canvas graphics, NOT chemical Y-up conventions.
pub fn compute_layout(mol: &Molecule) -> Layout {
    let n = mol.atom_count();
    if n == 0 {
        return Layout { coords: Vec::new() };
    }

    // Special case: single atom.
    if n == 1 {
        return Layout {
            coords: vec![Point::new(0.0, 0.0)],
        };
    }

    // Collect connected components so each can be laid out separately.
    let components = connected_components(mol);

    let mut all_coords: Vec<Option<Point>> = vec![None; n];
    let mut fragment_max_x = 0.0_f64;

    for component_atoms in &components {
        let component_set: HashSet<AtomIdx> = component_atoms.iter().copied().collect();

        // Subsets: placed will hold the coordinates for this component.
        let mut placed: Vec<Option<Point>> = vec![None; n];

        // Find SSSR for the whole molecule, then filter to this component.
        let ring_set = find_sssr(mol);
        let rings: Vec<Vec<AtomIdx>> = ring_set
            .rings()
            .iter()
            .filter(|ring| ring.iter().all(|a| component_set.contains(a)))
            .cloned()
            .collect();

        // Group rings into ring systems (connected sets of rings sharing >= 1 atom).
        let ring_systems = group_ring_systems(&rings);

        let mut atom_to_system: HashMap<AtomIdx, usize> = HashMap::new();
        for (sys_idx, system) in ring_systems.iter().enumerate() {
            for &a in system.iter().flatten() {
                atom_to_system.insert(a, sys_idx);
            }
        }
        let mut system_placed: Vec<bool> = vec![false; ring_systems.len()];

        // Seed: the ring system that anchors this component's whole
        // coordinate frame (see `seed_ring_system_index`'s doc). Every
        // other ring system gets discovered and anchored to its real
        // attachment point as the layout grows outward -- placing every
        // ring system blind at the origin (the pre-fix behavior) makes
        // unrelated ring systems of the same size collide exactly.
        if let Some(seed_idx) = seed_ring_system_index(&ring_systems) {
            place_ring_system(&ring_systems[seed_idx], None, &mut placed);
            system_placed[seed_idx] = true;
        }

        // If this component has no rings at all, seed a terminal chain atom
        // so the growth pass below has a starting point.
        seed_isolated_chain_start(mol, &component_set, &mut placed);

        // Grow everything else outward: chain atoms via DFS zigzag, newly
        // discovered ring systems anchored to their real attachment point.
        grow_layout(
            mol,
            &component_set,
            &ring_systems,
            &atom_to_system,
            &mut system_placed,
            &mut placed,
        );

        // Defensive fallbacks -- should not fire for any connected
        // component, since grow_layout's worklist reaches every atom
        // reachable from the seed via the molecule graph.
        for (sys_idx, system) in ring_systems.iter().enumerate() {
            if !system_placed[sys_idx] {
                place_ring_system(system, None, &mut placed);
            }
        }
        let mut still_unplaced: Vec<AtomIdx> = component_set
            .iter()
            .copied()
            .filter(|a| placed[a.0 as usize].is_none())
            .collect();
        still_unplaced.sort_unstable();
        let mut x = 0.0;
        for atom in still_unplaced {
            x += BOND_LEN;
            placed[atom.0 as usize] = Some(Point::new(x, 0.0));
        }

        // Offset this component to the right of the previous one.
        let x_offset = if fragment_max_x == 0.0 {
            0.0
        } else {
            fragment_max_x + 2.0 * BOND_LEN
        };

        // Find the min_x of this component so we can pack left.
        let comp_min_x = component_atoms
            .iter()
            .filter_map(|&a| placed[a.0 as usize])
            .map(|p| p.x)
            .fold(f64::MAX, f64::min);

        let shift = x_offset - comp_min_x;

        for &a in component_atoms {
            if let Some(p) = placed[a.0 as usize] {
                let shifted = Point::new(p.x + shift, p.y);
                all_coords[a.0 as usize] = Some(shifted);
                if shifted.x > fragment_max_x {
                    fragment_max_x = shifted.x;
                }
            }
        }
    }

    Layout {
        coords: all_coords
            .into_iter()
            .map(|p| p.unwrap_or(Point::new(0.0, 0.0)))
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Connected components
// ---------------------------------------------------------------------------

fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));
        visited[start] = true;

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
// Ring system grouping
// ---------------------------------------------------------------------------

/// Group rings into ring systems where rings in the same system share at least one atom.
fn group_ring_systems(rings: &[Vec<AtomIdx>]) -> Vec<Vec<Vec<AtomIdx>>> {
    if rings.is_empty() {
        return Vec::new();
    }

    // Union-Find by ring index.
    let n = rings.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut Vec<usize>, i: usize) -> usize {
        if parent[i] != i {
            parent[i] = find(parent, parent[i]);
        }
        parent[i]
    }

    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    // Two rings in the same system if they share any atom.
    for (i, ring_i) in rings.iter().enumerate() {
        let set_i: HashSet<AtomIdx> = ring_i.iter().copied().collect();
        for (j, ring_j) in rings.iter().enumerate().skip(i + 1) {
            if ring_j.iter().any(|a| set_i.contains(a)) {
                union(&mut parent, i, j);
            }
        }
    }

    // Collect into groups.
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }

    groups
        .values()
        .map(|indices| indices.iter().map(|&i| rings[i].clone()).collect())
        .collect()
}

// ---------------------------------------------------------------------------
// Ring system placement
// ---------------------------------------------------------------------------

/// Place all atoms of a ring system (a connected group of rings).
///
/// `anchor`, when `Some((entry_atom, entry_pos, dir))`, anchors the ring
/// containing `entry_atom` to `entry_pos` and extends it outward in
/// direction `dir` via [`place_first_ring_anchored`], instead of placing it
/// blind at the origin via [`place_regular_ring`]. Only the seed ring
/// system of a molecule (the one that establishes the whole layout's
/// coordinate frame) should ever pass `None` -- every other ring system
/// must be anchored to whatever already-placed atom it's attached to, or it
/// collides with unrelated geometry (see `place_regular_ring`'s doc for the
/// bug this replaces).
fn place_ring_system(
    system: &[Vec<AtomIdx>],
    anchor: Option<(AtomIdx, Point, f64)>,
    placed: &mut [Option<Point>],
) {
    if system.is_empty() {
        return;
    }

    // The anchored ring isn't necessarily system[0] -- find whichever ring
    // in this system actually contains entry_atom.
    let first_ring_idx = match anchor {
        Some((entry_atom, ..)) => system
            .iter()
            .position(|ring| ring.contains(&entry_atom))
            .unwrap_or(0),
        None => 0,
    };

    match anchor {
        Some((entry_atom, entry_pos, dir)) => {
            place_first_ring_anchored(&system[first_ring_idx], entry_atom, entry_pos, dir, placed)
        }
        None => place_regular_ring(&system[first_ring_idx], placed),
    }

    // For subsequent rings: find two atoms already placed (the shared edge),
    // then reflect unplaced atoms over that shared edge.
    let mut remaining: Vec<&Vec<AtomIdx>> = system
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != first_ring_idx)
        .map(|(_, ring)| ring)
        .collect();
    let mut iterations = 0;

    while !remaining.is_empty() && iterations < remaining.len() * 2 {
        iterations += 1;
        let mut progressed = false;

        remaining.retain(|ring| {
            // Find atoms of this ring that are already placed.
            let already_placed: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|&a| placed[a.0 as usize].is_some())
                .collect();

            if already_placed.is_empty() {
                return true; // Not ready yet.
            }

            if already_placed.len() == 1 {
                // Spiro junction: exactly one atom shared with an
                // already-placed ring. Anchor a fresh regular polygon at
                // that atom, extending away from everything placed so far
                // -- the same collision this whole function's `anchor`
                // parameter exists to avoid, reached via a different path
                // (a shared *atom* rather than a shared *edge*).
                let entry_atom = already_placed[0];
                let Some(entry_pos) = placed[entry_atom.0 as usize] else {
                    return true; // Not ready (shouldn't happen).
                };
                let dir = direction_away_from_centroid(entry_pos, placed);
                place_first_ring_anchored(ring, entry_atom, entry_pos, dir, placed);
                progressed = true;
                return false;
            }

            // Find the shared edge: two consecutive atoms in the ring that are both placed.
            let shared_edge = find_shared_edge(ring, placed);

            // Fall back: use the first two placed atoms.
            let (anchor1, anchor2) = shared_edge.unwrap_or((already_placed[0], already_placed[1]));

            // Both anchors are confirmed placed (either from find_shared_edge or already_placed).
            let (Some(p1), Some(p2)) = (placed[anchor1.0 as usize], placed[anchor2.0 as usize])
            else {
                return true; // Not ready.
            };

            // Place unplaced atoms of this ring using the regular polygon geometry,
            // anchored to the shared edge.
            place_ring_anchored(ring, anchor1, p1, anchor2, p2, placed);

            progressed = true;
            false // Remove from remaining.
        });

        if !progressed {
            break;
        }
    }

    // Any still-unplaced atoms in remaining rings: force-place them.
    // Should not fire for any ring system reachable from a shared atom or
    // edge; kept as a defensive fallback only.
    for ring in &remaining {
        place_regular_ring(ring, placed);
    }
}

/// Find the shared edge in a ring: two consecutive atoms in ring order that are both placed.
///
/// Returns `None` if no consecutive placed pair exists.
fn find_shared_edge(ring: &[AtomIdx], placed: &[Option<Point>]) -> Option<(AtomIdx, AtomIdx)> {
    let n = ring.len();
    ring.windows(2)
        .map(|w| (w[0], w[1]))
        .chain(std::iter::once((ring[n - 1], ring[0])))
        .find(|&(a, b)| placed[a.0 as usize].is_some() && placed[b.0 as usize].is_some())
}

/// Place atoms of a ring as a regular polygon centered at the origin.
fn place_regular_ring(ring: &[AtomIdx], placed: &mut [Option<Point>]) {
    let n = ring.len();
    if n == 0 {
        return;
    }

    let radius = ring_radius(n);
    // Start angle: 90 degrees (pointing up), atoms go clockwise.
    let start_angle = std::f64::consts::FRAC_PI_2;

    for (i, &atom) in ring.iter().enumerate() {
        if placed[atom.0 as usize].is_none() {
            let angle = start_angle - (2.0 * std::f64::consts::PI * i as f64) / n as f64;
            let x = radius * angle.cos();
            let y = -radius * angle.sin(); // SVG y increases downward.
            placed[atom.0 as usize] = Some(Point::new(x, y));
        }
    }
}

/// Place unplaced atoms of `ring` given that `anchor1` and `anchor2` are already placed.
///
/// The unplaced atoms are positioned so that:
/// - The ring forms a regular n-gon.
/// - The new ring extends away from the already-placed ring system.
fn place_ring_anchored(
    ring: &[AtomIdx],
    anchor1: AtomIdx,
    p1: Point,
    anchor2: AtomIdx,
    p2: Point,
    placed: &mut [Option<Point>],
) {
    let n = ring.len();
    let radius = ring_radius(n);

    // Find anchor indices in ring order.
    let idx1 = ring.iter().position(|&a| a == anchor1).unwrap_or(0);
    let idx2 = ring.iter().position(|&a| a == anchor2).unwrap_or(1);

    // Compute the new ring center:
    // it lies on the perpendicular bisector of the shared edge (p1..p2),
    // at distance = apothem (= R*cos(PI/n)) from the midpoint, on the
    // side AWAY from the existing placed atoms of this ring.

    let mid = Point::new((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let edge_len = (dx * dx + dy * dy).sqrt();
    if edge_len < 1e-10 {
        return;
    }

    // Perpendicular unit vector.
    let perp_x = -dy / edge_len;
    let perp_y = dx / edge_len;

    // Apothem: distance from midpoint of an edge to the center of a regular n-gon.
    let apothem = radius * (std::f64::consts::PI / n as f64).cos();

    // Centroid of the *entire* already-placed structure so far (not just this
    // ring's own anchor atoms). Using only this ring's placed atoms is wrong:
    // before this ring is placed, its only placed members are the two shared
    // anchor atoms, whose centroid is always exactly `mid` (the edge
    // midpoint) — equidistant from both `cand1`/`cand2` by construction, so
    // that comparison degenerates into an arbitrary tie instead of actually
    // picking the side away from the existing ring system.
    let existing_center = centroid_of_placed(placed).unwrap_or(mid);

    // Choose the candidate center farther from the existing ring centroid.
    let cand1 = Point::new(mid.x + perp_x * apothem, mid.y + perp_y * apothem);
    let cand2 = Point::new(mid.x - perp_x * apothem, mid.y - perp_y * apothem);
    let new_center = if cand1.dist(&existing_center) > cand2.dist(&existing_center) {
        cand1
    } else {
        cand2
    };

    // Angle from new_center to anchor1.
    let angle_to_a1 = (p1.y - new_center.y).atan2(p1.x - new_center.x);

    // Determine the angular step direction:
    // In the ring, going from idx1 to idx2 by +1 step in ring order should correspond to
    // going from angle_to_a1 to angle_to_a2.  We choose the sign of angle_step
    // so that the rotation from anchor1 to anchor2 (by 1 ring step) matches the geometry.
    let steps_forward = (idx2 + n - idx1) % n; // steps from idx1 to idx2 in ring order.
    let angle_to_a2 = (p2.y - new_center.y).atan2(p2.x - new_center.x);

    // Signed angle from a1 to a2 (normalized to -PI..PI).
    let mut delta = angle_to_a2 - angle_to_a1;
    while delta > std::f64::consts::PI {
        delta -= 2.0 * std::f64::consts::PI;
    }
    while delta < -std::f64::consts::PI {
        delta += 2.0 * std::f64::consts::PI;
    }

    // Expected angular step per ring step (2*PI/n in either direction).
    // If going steps_forward ring steps gives delta angle, then:
    //   angle_step = delta / steps_forward  (but might differ slightly from 2PI/n due to
    //   the fixed radius; we use the sign of delta to pick CW vs CCW).
    let angle_step = if steps_forward > 0 {
        if delta >= 0.0 {
            2.0 * std::f64::consts::PI / n as f64
        } else {
            -(2.0 * std::f64::consts::PI / n as f64)
        }
    } else {
        2.0 * std::f64::consts::PI / n as f64
    };

    // Place each unplaced atom.
    for step in 0..n {
        let ring_idx = (idx1 + step) % n;
        let atom = ring[ring_idx];

        if placed[atom.0 as usize].is_some() {
            continue;
        }

        let angle = angle_to_a1 + step as f64 * angle_step;
        let x = new_center.x + radius * angle.cos();
        let y = new_center.y + radius * angle.sin();
        placed[atom.0 as usize] = Some(Point::new(x, y));
    }
}

/// Place `ring` as a regular polygon with `entry_atom` on its circumference
/// at `entry_pos` (set there if not already placed), extending outward in
/// direction `dir` -- i.e. the ring's center sits one radius from
/// `entry_pos` in direction `dir`. Used to anchor a ring system to a real
/// attachment point (an already-placed atom for a spiro junction, or a
/// freshly-placed one bond length from a chain/exocyclic parent) instead of
/// [`place_regular_ring`]'s unconditional-origin placement.
fn place_first_ring_anchored(
    ring: &[AtomIdx],
    entry_atom: AtomIdx,
    entry_pos: Point,
    dir: f64,
    placed: &mut [Option<Point>],
) {
    let n = ring.len();
    if n == 0 {
        return;
    }
    if placed[entry_atom.0 as usize].is_none() {
        placed[entry_atom.0 as usize] = Some(entry_pos);
    }

    let radius = ring_radius(n);
    let center = Point::new(
        entry_pos.x + radius * dir.cos(),
        entry_pos.y + radius * dir.sin(),
    );
    let idx0 = ring.iter().position(|&a| a == entry_atom).unwrap_or(0);
    let angle_to_entry = (entry_pos.y - center.y).atan2(entry_pos.x - center.x);
    // Clockwise, matching place_regular_ring's convention.
    let angle_step = -2.0 * std::f64::consts::PI / n as f64;

    for step in 0..n {
        let ring_idx = (idx0 + step) % n;
        let atom = ring[ring_idx];
        if placed[atom.0 as usize].is_some() {
            continue;
        }
        let angle = angle_to_entry + step as f64 * angle_step;
        let x = center.x + radius * angle.cos();
        let y = center.y + radius * angle.sin();
        placed[atom.0 as usize] = Some(Point::new(x, y));
    }
}

/// Centroid of every currently-placed atom (across the whole molecule, not
/// just one ring/component). `None` if nothing is placed yet.
fn centroid_of_placed(placed: &[Option<Point>]) -> Option<Point> {
    let pts: Vec<Point> = placed.iter().filter_map(|p| *p).collect();
    if pts.is_empty() {
        return None;
    }
    let cx = pts.iter().map(|p| p.x).sum::<f64>() / pts.len() as f64;
    let cy = pts.iter().map(|p| p.y).sum::<f64>() / pts.len() as f64;
    Some(Point::new(cx, cy))
}

/// Direction from the centroid of everything placed so far, through
/// `entry_pos`, continued outward -- the natural "grow away from what's
/// already there" direction for a spiro-anchored ring, which (unlike a
/// chain/exocyclic attachment) has no incoming bond direction to continue.
fn direction_away_from_centroid(entry_pos: Point, placed: &[Option<Point>]) -> f64 {
    match centroid_of_placed(placed) {
        Some(c) if c.dist(&entry_pos) > 1e-9 => (entry_pos.y - c.y).atan2(entry_pos.x - c.x),
        _ => 0.0,
    }
}

/// Pick the ring system that should anchor a component's whole coordinate
/// frame: most rings, then most atoms, then lowest minimum `AtomIdx`
/// (deterministic, and prefers a fused/bridged/spiro core over a peripheral
/// single ring). Every other ring system is anchored relative to this one
/// as the layout grows outward.
fn seed_ring_system_index(ring_systems: &[Vec<Vec<AtomIdx>>]) -> Option<usize> {
    ring_systems
        .iter()
        .enumerate()
        .max_by_key(|(_, system)| {
            let n_rings = system.len();
            let mut atoms: Vec<u32> = system.iter().flatten().map(|a| a.0).collect();
            atoms.sort_unstable();
            atoms.dedup();
            let n_atoms = atoms.len();
            let min_atom = atoms.first().copied().unwrap_or(u32::MAX);
            (n_rings, n_atoms, std::cmp::Reverse(min_atom))
        })
        .map(|(i, _)| i)
}

/// Compute the circumradius for a regular n-gon with the given BOND_LEN.
fn ring_radius(n: usize) -> f64 {
    if n < 3 {
        return BOND_LEN;
    }
    BOND_LEN / (2.0 * (std::f64::consts::PI / n as f64).sin())
}

// ---------------------------------------------------------------------------
// Chain placement (DFS zigzag)
// ---------------------------------------------------------------------------

/// If nothing in `component` is placed yet (no ring system exists to seed
/// from), place a terminal atom (degree ≤1, or an arbitrary atom if none)
/// at the origin so [`grow_layout`] has a starting point.
fn seed_isolated_chain_start(
    mol: &Molecule,
    component: &HashSet<AtomIdx>,
    placed: &mut [Option<Point>],
) {
    if component.iter().any(|a| placed[a.0 as usize].is_some()) {
        return; // Already has a seed (a ring system was placed).
    }

    let mut unplaced: Vec<AtomIdx> = component.iter().copied().collect();
    unplaced.sort_unstable();
    let Some(&start) = unplaced
        .iter()
        .find(|&&a| mol.degree(a) <= 1)
        .or(unplaced.first())
    else {
        return; // Empty component (shouldn't happen).
    };

    placed[start.0 as usize] = Some(Point::new(0.0, 0.0));
}

/// Grow the layout outward from whatever's already placed (the seed ring
/// system, or the isolated chain start seeded by
/// [`seed_isolated_chain_start`]): place chain atoms via [`dfs_zigzag`], and
/// anchor any newly discovered ring system to its real attachment point via
/// [`place_ring_system`] instead of leaving it for a blind, colliding
/// placement (see `place_regular_ring`'s doc for the bug this replaces).
fn grow_layout(
    mol: &Molecule,
    component: &HashSet<AtomIdx>,
    ring_systems: &[Vec<Vec<AtomIdx>>],
    atom_to_system: &HashMap<AtomIdx, usize>,
    system_placed: &mut [bool],
    placed: &mut [Option<Point>],
) {
    let mut worklist: VecDeque<AtomIdx> = {
        let mut seeded: Vec<AtomIdx> = component
            .iter()
            .copied()
            .filter(|a| placed[a.0 as usize].is_some())
            .collect();
        seeded.sort_unstable();
        seeded.into()
    };

    while let Some(start) = worklist.pop_front() {
        if placed[start.0 as usize].is_none() {
            continue;
        }

        let mut unplaced_neighbors: Vec<AtomIdx> = mol
            .neighbors(start)
            .map(|(nb, _)| nb)
            .filter(|nb| component.contains(nb) && placed[nb.0 as usize].is_none())
            .collect();
        unplaced_neighbors.sort_unstable();

        for nb in unplaced_neighbors {
            if placed[nb.0 as usize].is_some() {
                continue; // Placed by an earlier neighbor this same pass (e.g. a shared spiro/fused atom).
            }
            // Determine outgoing direction from the already-placed atom.
            // Use a direction that avoids existing neighbors.
            let dir = best_outgoing_direction(start, placed, mol, component);
            let mut newly_ring_placed = Vec::new();
            dfs_zigzag(
                mol,
                nb,
                start,
                dir,
                placed,
                component,
                ring_systems,
                atom_to_system,
                system_placed,
                &mut newly_ring_placed,
            );
            newly_ring_placed.sort_unstable();
            worklist.extend(newly_ring_placed);
        }
    }
}

/// Compute the best outgoing direction from a ring atom to avoid collisions.
///
/// Ranks candidates by angular separation via [`ranked_candidates`] (see that
/// function's doc for why -- issue #347: this used to have its own, coarser
/// 60°-spaced candidate set with no chemistry-aware offsets, missing the
/// correct bisector for e.g. a hexagon-ring substituent by ~30°), then walks
/// them best-first and skips any candidate whose resulting position would
/// land on top of an already-placed atom elsewhere in the component --
/// angular separation from *this atom's own bonds* alone doesn't guard
/// against that (issue #347: a bridged-core molecule with a separate,
/// pre-existing placement bug has an atom sitting far from where its own
/// bonds would suggest, and the top-angular candidate can point straight at
/// it). Falls back to the top-ranked candidate if every one collides.
fn best_outgoing_direction(
    atom: AtomIdx,
    placed: &[Option<Point>],
    mol: &Molecule,
    component: &HashSet<AtomIdx>,
) -> f64 {
    let Some(origin) = placed[atom.0 as usize] else {
        return 0.0;
    };

    // Collect angles to already-placed neighbors.
    let used_angles: Vec<f64> = mol
        .neighbors(atom)
        .filter(|(nb, _)| component.contains(nb) && placed[nb.0 as usize].is_some())
        .map(|(nb, _)| {
            let p = placed[nb.0 as usize].unwrap();
            (p.y - origin.y).atan2(p.x - origin.x)
        })
        .collect();

    if used_angles.is_empty() {
        return 0.0;
    }

    let ranked = ranked_candidates(&used_angles);
    let occupies = |dir: f64| -> bool {
        let candidate = Point::new(
            origin.x + BOND_LEN * dir.cos(),
            origin.y + BOND_LEN * dir.sin(),
        );
        placed
            .iter()
            .enumerate()
            .filter(|(i, _)| component.contains(&AtomIdx(*i as u32)))
            .filter_map(|(_, p)| p.as_ref())
            .any(|p| candidate.dist(p) < BOND_LEN / 2.0)
    };

    ranked
        .iter()
        .copied()
        .find(|&dir| !occupies(dir))
        .unwrap_or(ranked[0])
}

/// Pick the direction (radians) that maximizes the minimum angular separation
/// from every angle in `used_angles`.
///
/// Thin wrapper over [`ranked_candidates`] -- see that function's doc for the
/// candidate set. Returns `0.0` if `used_angles` is empty.
fn best_direction_avoiding(used_angles: &[f64]) -> f64 {
    if used_angles.is_empty() {
        return 0.0;
    }
    ranked_candidates(used_angles)[0]
}

/// Candidate directions (radians), best-first by minimum angular separation
/// from every angle in `used_angles`.
///
/// Candidates: a 30-degree grid (12 directions covering 360°), plus, for each
/// angle already in `used_angles`: its anti-direction (α+180°), the two sp3
/// zigzag offsets (α+150°/α+210°) and the two sp2 offsets (α±120°) --
/// chemistry-aware candidates that a plain fixed grid alone can miss (issue
/// #347's ring-substituent bisector case: two ring bonds 120° apart need a
/// candidate at their exact bisector, which isn't always on a 30°-aligned
/// grid point relative to 0°, but always is one of these bond-relative
/// offsets). Shared by [`best_outgoing_direction`] and
/// [`suggest_bond_direction`] (via [`best_direction_avoiding`]) -- previously
/// duplicated with two different (one worse) candidate sets; this is the
/// single, richer implementation both now use. Returns `[0.0]` if
/// `used_angles` is empty.
fn ranked_candidates(used_angles: &[f64]) -> Vec<f64> {
    use std::f64::consts::PI;

    if used_angles.is_empty() {
        return vec![0.0];
    }

    let mut candidates: Vec<f64> = (0..12).map(|i| i as f64 * PI / 6.0).collect();
    for &a in used_angles {
        candidates.push(a + PI);
        candidates.push(a + PI - PI / 6.0);
        candidates.push(a + PI + PI / 6.0);
        candidates.push(a + 2.0 * PI / 3.0);
        candidates.push(a - 2.0 * PI / 3.0);
    }

    candidates.sort_by(|&a, &b| {
        let min_sep_a = min_angle_separation(a, used_angles);
        let min_sep_b = min_angle_separation(b, used_angles);
        min_sep_b.partial_cmp(&min_sep_a).unwrap()
    });
    candidates
}

/// Minimum angular separation between `angle` and any angle in `used`.
fn min_angle_separation(angle: f64, used: &[f64]) -> f64 {
    used.iter()
        .map(|&u| {
            let diff = (angle - u).abs();
            if diff > std::f64::consts::PI {
                2.0 * std::f64::consts::PI - diff
            } else {
                diff
            }
        })
        .fold(f64::MAX, f64::min)
}

/// Iterative zigzag placement: place `start_atom` at `BOND_LEN` from `start_parent`
/// in direction `start_dir`, then expand unplaced neighbors with alternating ±30° deflection.
///
/// The alternation is carried as a `sign` (±1.0) threaded through the DFS stack,
/// not derived from a neighbor's position in its parent's own unplaced-neighbor
/// list -- that would only ever alternate at a genuine branch point (2+ unplaced
/// neighbors), and always reapply the same deflection on an ordinary single-child
/// chain continuation (the overwhelming majority of atoms), producing a monotonic
/// per-bond rotational drift instead of a zigzag (issue #347: a plain 13-carbon
/// chain's first and last atoms landed on identical coordinates, since 12
/// consecutive -30° steps trace a full circle).
///
/// If a popped atom belongs to a not-yet-placed ring system (`atom_to_system`),
/// this anchors that whole system via [`place_ring_system`]/
/// [`place_first_ring_anchored`] instead of placing just that one atom, and
/// records every atom of the newly-placed system into `newly_ring_placed`
/// (its own neighbors, and the direction to grow from them, are the
/// caller's job -- [`grow_layout`]'s outer worklist -- not this DFS stack's,
/// since further ring growth needs `best_outgoing_direction`, not zigzag
/// deflection).
#[allow(clippy::too_many_arguments)]
fn dfs_zigzag(
    mol: &Molecule,
    start_atom: AtomIdx,
    start_parent: AtomIdx,
    start_dir: f64,
    placed: &mut [Option<Point>],
    component: &HashSet<AtomIdx>,
    ring_systems: &[Vec<Vec<AtomIdx>>],
    atom_to_system: &HashMap<AtomIdx, usize>,
    system_placed: &mut [bool],
    newly_ring_placed: &mut Vec<AtomIdx>,
) {
    let deflection = std::f64::consts::PI / 6.0;
    // 4th element: the sign to apply, and then flip, when this atom continues
    // the chain alone. -1.0 matches the pre-fix first-step behavior.
    let mut stack: Vec<(AtomIdx, AtomIdx, f64, f64)> =
        vec![(start_atom, start_parent, start_dir, -1.0)];

    while let Some((atom, parent, dir, sign)) = stack.pop() {
        if placed[atom.0 as usize].is_some() {
            continue;
        }
        let parent_pos = match placed[parent.0 as usize] {
            Some(p) => p,
            None => continue,
        };

        if let Some(&sys_idx) = atom_to_system.get(&atom) {
            let entry_pos = Point::new(
                parent_pos.x + BOND_LEN * dir.cos(),
                parent_pos.y + BOND_LEN * dir.sin(),
            );
            place_ring_system(&ring_systems[sys_idx], Some((atom, entry_pos, dir)), placed);
            system_placed[sys_idx] = true;
            newly_ring_placed.extend(ring_systems[sys_idx].iter().flatten().copied());
            continue; // Further growth from this ring's atoms is grow_layout's job.
        }

        let x = parent_pos.x + BOND_LEN * dir.cos();
        let y = parent_pos.y + BOND_LEN * dir.sin();
        placed[atom.0 as usize] = Some(Point::new(x, y));

        let unplaced: Vec<AtomIdx> = mol
            .neighbors(atom)
            .map(|(nb, _)| nb)
            .filter(|&nb| {
                nb != parent && component.contains(&nb) && placed[nb.0 as usize].is_none()
            })
            .collect();

        if unplaced.len() == 1 {
            // Ordinary chain continuation: apply this atom's sign, flip it for
            // the next step -- this is the alternation the pre-fix code missed.
            stack.push((unplaced[0], atom, dir + sign * deflection, -sign));
        } else {
            // A real branch (0, or 2+, unplaced neighbors): preserve the
            // original per-child split (alternating by position), each child
            // then continuing its own zigzag from its own starting sign.
            // Push in reverse so the first neighbor is popped first, preserving DFS order.
            for (i, nb) in unplaced.into_iter().enumerate().rev() {
                let child_sign = if i % 2 == 0 { -1.0 } else { 1.0 };
                stack.push((nb, atom, dir + child_sign * deflection, -child_sign));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Public bond-direction suggestion
// ---------------------------------------------------------------------------

/// Suggest the best direction (radians, measured from positive x-axis) for a
/// new bond leaving `atom`, given the molecule's current 2D `layout`.
///
/// Collects angles to all already-placed neighbors of `atom`, then delegates
/// candidate generation/selection to [`best_direction_avoiding`] -- see that
/// function's doc for the candidate set and selection rule.
///
/// Returns `0.0` (pointing right) when `atom` has no neighbors in `layout`.
pub fn suggest_bond_direction(mol: &Molecule, atom: AtomIdx, layout: &Layout) -> f64 {
    let origin = layout.get(atom);

    // Angles to neighbors that are already placed in the layout.
    let used_angles: Vec<f64> = mol
        .neighbors(atom)
        .filter(|(nb, _)| (nb.0 as usize) < layout.coords.len())
        .map(|(nb, _)| {
            let p = layout.get(nb);
            (p.y - origin.y).atan2(p.x - origin.x)
        })
        .collect();

    if used_angles.is_empty() {
        return 0.0;
    }

    best_direction_avoiding(&used_angles)
}

// ---------------------------------------------------------------------------
// Bond crossing detection
// ---------------------------------------------------------------------------

/// Detect which pairs of bonds have crossing 2D segments.
///
/// Returns a `Vec<(BondIdx, BondIdx)>` listing all bonds that intersect in the
/// layout. Bonds that share a common atom (adjacent bonds) are not checked.
///
/// Useful for assessing layout quality: an empty result indicates a crossing-free
/// (or at least non-crossing-bond) depiction.
pub fn detect_crossings(layout: &Layout, mol: &Molecule) -> Vec<(BondIdx, BondIdx)> {
    let bonds: Vec<(BondIdx, (Point, Point))> = mol
        .bonds()
        .map(|(bidx, bond)| {
            let p1 = layout.get(bond.atom1);
            let p2 = layout.get(bond.atom2);
            (bidx, (p1, p2))
        })
        .collect();

    let mut crossings = Vec::new();

    for i in 0..bonds.len() {
        for j in (i + 1)..bonds.len() {
            let (bidx_i, (p1_i, p2_i)) = bonds[i];
            let (bidx_j, (p1_j, p2_j)) = bonds[j];

            // Skip if bonds share an atom (adjacent bonds always "cross" at the vertex)
            let a1_i = mol.bond(bidx_i).atom1;
            let a2_i = mol.bond(bidx_i).atom2;
            let a1_j = mol.bond(bidx_j).atom1;
            let a2_j = mol.bond(bidx_j).atom2;

            if a1_i == a1_j || a1_i == a2_j || a2_i == a1_j || a2_i == a2_j {
                continue;
            }

            // Check for line segment intersection using cross product
            if segments_intersect(p1_i, p2_i, p1_j, p2_j) {
                crossings.push((bidx_i, bidx_j));
            }
        }
    }

    crossings
}

/// Check if two line segments AB and CD intersect (not including endpoints touching).
fn segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool {
    let ccw = |p1: Point, p2: Point, p3: Point| -> bool {
        (p3.y - p1.y) * (p2.x - p1.x) > (p2.y - p1.y) * (p3.x - p1.x)
    };

    ccw(a, c, d) != ccw(b, c, d) && ccw(a, b, c) != ccw(a, b, d)
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Positive control for the group_ring_systems/place_chains determinism
    /// fix: repeated calls in the SAME process on the SAME molecule must
    /// produce bit-identical coordinates. This is a stricter, different bug
    /// class than input-order dependence -- `group_ring_systems` (a plain
    /// `HashMap`) and `place_chains` (iterating a `HashSet` without sorting)
    /// used Rust's randomly-seeded default hasher, so a fresh `RandomState`
    /// derived on every `HashMap`/`HashSet` construction could reorder
    /// output even across calls with IDENTICAL input, within one run of one
    /// binary -- undetectable by any test that only varies input spelling.
    #[test]
    fn compute_layout_is_deterministic_across_repeated_calls() {
        use chematic_smiles::parse;
        // Two separate ring systems joined by a chain, so both
        // group_ring_systems (2 union-find roots) and place_chains (a
        // multi-atom unplaced chain) have real tie material to shuffle.
        let mol = parse("c1ccccc1CCCCCc1ccccc1").unwrap();

        let first = compute_layout(&mol);
        for _ in 0..100 {
            let repeat = compute_layout(&mol);
            assert_eq!(
                repeat.coords, first.coords,
                "compute_layout produced different coordinates for the same molecule \
                 across repeated calls in one process"
            );
        }
    }

    // --- Ring-system-coincidence regression tests -----------------------
    //
    // Root cause: `place_regular_ring` always centered a new ring at the
    // literal origin with no awareness of already-placed geometry, and was
    // invoked unconditionally for every ring system -- so any two ring
    // systems not fused/bridged into the same connected group (a chain- or
    // bond-mediated substituent, or a spiro junction) landed on
    // bit-for-bit identical coordinates. Fixed via `place_first_ring_anchored`
    // plus a connectivity-driven growth pass (`grow_layout`) that discovers
    // and anchors every ring system to its real attachment point.

    /// (min non-bonded pairwise distance, min bonded distance, max bonded
    /// distance) across every atom pair in `mol`'s `layout`.
    fn layout_distance_summary(mol: &Molecule, layout: &Layout) -> (f64, f64, f64) {
        let n = mol.atom_count();
        let mut min_non_bonded = f64::MAX;
        let mut min_bonded = f64::MAX;
        let mut max_bonded = f64::MIN;
        for i in 0..n {
            for j in (i + 1)..n {
                let a = AtomIdx(i as u32);
                let b = AtomIdx(j as u32);
                let d = layout.get(a).dist(&layout.get(b));
                if mol.bond_between(a, b).is_some() {
                    min_bonded = min_bonded.min(d);
                    max_bonded = max_bonded.max(d);
                } else {
                    min_non_bonded = min_non_bonded.min(d);
                }
            }
        }
        (min_non_bonded, min_bonded, max_bonded)
    }

    /// Asserts no exact/near coincidence (Tier A/B) and exact bond-length
    /// fidelity (Tier C) -- the full set, for fixtures with no pre-existing
    /// unrelated geometry bug to work around.
    fn assert_layout_clean(smiles: &str) {
        use chematic_smiles::parse;
        let mol = parse(smiles).unwrap();
        let layout = compute_layout(&mol);
        let (min_non_bonded, min_bonded, max_bonded) = layout_distance_summary(&mol, &layout);
        assert!(
            min_non_bonded > BOND_LEN / 2.0,
            "{smiles}: non-bonded atoms too close (near-collision), min_non_bonded={min_non_bonded}"
        );
        assert!(
            (min_bonded - BOND_LEN).abs() < 1e-6,
            "{smiles}: bonded distance should equal BOND_LEN, min_bonded={min_bonded}"
        );
        assert!(
            (max_bonded - BOND_LEN).abs() < 1e-6,
            "{smiles}: bonded distance should equal BOND_LEN, max_bonded={max_bonded}"
        );
    }

    #[test]
    fn ring_systems_joined_by_chain_do_not_collide() {
        // Simplest isolation of the bug: two benzene rings joined by a
        // plain chain, no shared atom, no direct bond. This SMILES is also
        // the determinism-test fixture above, which only ever checked
        // repeat-call stability -- it passed on top of the broken layout
        // for a long time.
        assert_layout_clean("c1ccccc1CCCCCc1ccccc1");
    }

    #[test]
    fn ring_systems_joined_directly_do_not_collide() {
        // Direct-bond-mediated attachment (no intervening chain atoms) --
        // exercises the outer-worklist ring-discovery path in
        // `grow_layout`, distinct from the mid-chain discovery inside
        // `dfs_zigzag` the previous test exercises.
        assert_layout_clean("c1ccc(cc1)-c1ccccc1CC");
    }

    #[test]
    fn three_substituent_rings_on_bridged_core_do_not_collide() {
        // The originally reported molecule: three separate phenyl-ring
        // substituents on a bridged bicyclic core, all landing on exactly
        // the same coordinates pre-fix (18 of 27 near-neighbor pairs at
        // distance 0.0). Tier A only (no exact/near coincidence) -- the
        // bridged core's OWN internal bond lengths have a separate,
        // pre-existing bug (`find_shared_edge` only handles a 2-atom
        // shared edge, not the 3-atom-shared case a true bridge produces --
        // e.g. this molecule's atom5-atom6 bond, a real bond, lands ~108
        // units apart instead of the expected `BOND_LEN` of 40), out of
        // scope for this fix, so Tier B/C are not asserted here.
        use chematic_smiles::parse;
        let mol = parse("C1CC2CN(CC1N2c1ccccc1)c1cccc(c1)-c1ccccc1").unwrap();
        let layout = compute_layout(&mol);
        let (min_non_bonded, _min_bonded, _max_bonded) = layout_distance_summary(&mol, &layout);
        assert!(
            min_non_bonded > 1e-6,
            "no two atoms of unrelated ring systems should land on identical coordinates: \
             min_non_bonded={min_non_bonded}"
        );
    }

    #[test]
    fn pure_chain_layout_unaffected() {
        assert_layout_clean("CCCCCCCC");
    }

    // --- Issue #347 regressions ------------------------------------------

    #[test]
    fn long_chain_does_not_wrap_onto_itself() {
        // Pre-fix, dfs_zigzag applied a constant -30°/bond drift instead of
        // alternating, so 12 consecutive bonds traced a full circle: a plain
        // 13-carbon chain's first and last atoms landed on identical
        // coordinates. 20 carbons is well past that wrap point -- a true
        // zigzag (ping-ponging between two directions) never revisits a
        // point no matter how long the chain runs, so this can't pass by
        // accident the way a shorter fixture might.
        assert_layout_clean(&"C".repeat(20));
    }

    #[test]
    fn long_chain_bond_directions_genuinely_alternate() {
        // Direct angle check, not just non-collision: asserts the actual
        // zigzag pattern (ping-ponging between exactly two directions 30°
        // apart), not merely that nothing happened to collide.
        use chematic_smiles::parse;
        let mol = parse(&"C".repeat(10)).unwrap();
        let layout = compute_layout(&mol);
        let n = mol.atom_count();
        let dirs: Vec<f64> = (0..n - 1)
            .map(|i| {
                let a = layout.get(AtomIdx(i as u32));
                let b = layout.get(AtomIdx((i + 1) as u32));
                (b.y - a.y).atan2(b.x - a.x)
            })
            .collect();
        // Consecutive bond directions must differ by exactly 30° in
        // magnitude (a real zigzag turn), never 0° (collinear) or drifting
        // to some other value.
        for w in dirs.windows(2) {
            let mut turn = (w[1] - w[0]).abs();
            if turn > std::f64::consts::PI {
                turn = 2.0 * std::f64::consts::PI - turn;
            }
            assert!(
                (turn - std::f64::consts::PI / 6.0).abs() < 1e-6,
                "expected a 30° zigzag turn between consecutive bonds, got {turn} rad \
                 (dirs={dirs:?})"
            );
        }
        // And the pattern must actually alternate (ping-pong), not drift:
        // only 2 distinct directions should appear across the whole chain.
        let mut distinct: Vec<f64> = Vec::new();
        for &d in &dirs {
            if !distinct.iter().any(|&e: &f64| (e - d).abs() < 1e-6) {
                distinct.push(d);
            }
        }
        assert_eq!(
            distinct.len(),
            2,
            "expected exactly 2 distinct bond directions (ping-pong zigzag), got {distinct:?}"
        );
    }

    #[test]
    fn ring_substituent_direction_bisects_the_open_angle() {
        // Issue #347 repro 2's exact scenario: a hexagon-ring atom with two
        // ring bonds at -30°/90° (120° apart, as any regular-hexagon vertex
        // is). The correct outward direction for an exocyclic substituent is
        // the bisector of the open 240° angle: 210°. Pre-fix,
        // best_outgoing_direction's coarse 60°-grid + broken tie-break
        // returned 180° here -- a real ~30° miss.
        let used_angles = [-std::f64::consts::PI / 6.0, std::f64::consts::PI / 2.0];
        let dir = best_direction_avoiding(&used_angles);
        let expected = 7.0 * std::f64::consts::PI / 6.0; // 210°
        let mut diff = (dir - expected).abs();
        if diff > std::f64::consts::PI {
            diff = 2.0 * std::f64::consts::PI - diff;
        }
        assert!(
            diff < 1e-6,
            "expected the 210° bisector, got {} rad ({} deg)",
            dir,
            dir.to_degrees()
        );
    }

    #[test]
    fn ring_plus_exocyclic_branch_layout_unaffected() {
        // Full end-to-end version of issue #347 repro 2.
        assert_layout_clean("C1CCCC(C(=O)CC)C1");
    }

    #[test]
    fn single_ring_layout_unaffected() {
        assert_layout_clean("c1ccccc1");
    }

    #[test]
    fn fused_and_spiro_ring_systems_unaffected() {
        // Spiro exercises `place_ring_system`'s `already_placed.len() == 1`
        // fix directly; naphthalene/decalin (fused/bridged, 2-atom shared
        // edge) confirm `place_ring_anchored`'s existing, untouched path
        // still works after `place_ring_system`'s anchor-selection refactor.
        for smiles in ["C1CCC2(CC1)CCCCC2", "c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"] {
            assert_layout_clean(smiles);
        }
    }

    #[test]
    fn test_ring_radius_hexagon() {
        // Regular hexagon: all sides == BOND_LEN, circumradius == BOND_LEN.
        let r = ring_radius(6);
        let expected = BOND_LEN; // sin(PI/6) = 0.5, so BOND_LEN / (2*0.5) = BOND_LEN
        assert!((r - expected).abs() < 1e-9, "hexagon radius = {}", r);
    }

    // --- suggest_bond_direction tests ---

    fn make_layout(coords: &[(f64, f64)]) -> Layout {
        Layout {
            coords: coords.iter().map(|&(x, y)| Point { x, y }).collect(),
        }
    }

    #[test]
    fn test_suggest_direction_no_neighbors() {
        use chematic_smiles::parse;
        // Single atom: no neighbors → default 0.0 (pointing right).
        let mol = parse("C").unwrap();
        let layout = make_layout(&[(0.0, 0.0)]);
        let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
        assert!((dir).abs() < 1e-9 || (dir - 2.0 * std::f64::consts::PI).abs() < 1e-9);
    }

    #[test]
    fn test_suggest_direction_single_bond_avoids_existing() {
        use chematic_smiles::parse;
        // Ethane C-C: atom 0 at origin, atom 1 to the right (angle = 0°).
        // Suggested direction for a third atom from atom 0 should be far from 0°.
        let mol = parse("CC").unwrap();
        let layout = make_layout(&[(0.0, 0.0), (BOND_LEN, 0.0)]);
        let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
        // Must be at least 90° away from 0° (the existing bond).
        let sep = min_angle_separation(dir, &[0.0_f64]);
        assert!(
            sep >= std::f64::consts::PI / 2.0,
            "suggested direction {dir:.3} should be ≥90° from existing bond, sep={sep:.3}"
        );
    }

    #[test]
    fn test_suggest_direction_two_bonds_finds_gap() {
        use chematic_smiles::parse;
        // Three atoms: center C bonded to left (180°) and right (0°).
        // Suggested direction should be ≈ 90° or ≈ -90° (the open gap above/below).
        let mol = parse("CCC").unwrap();
        // atom 1 (center) has neighbors at atom 0 (left) and atom 2 (right).
        let layout = make_layout(&[(-BOND_LEN, 0.0), (0.0, 0.0), (BOND_LEN, 0.0)]);
        let dir = suggest_bond_direction(&mol, AtomIdx(1), &layout);
        // The gap is ≈ 90° (top) or ≈ 270° (bottom). Both have min-sep ≈ 90° from 0° and 180°.
        let sep = {
            let used = [0.0_f64, std::f64::consts::PI];
            min_angle_separation(dir, &used)
        };
        assert!(
            sep >= std::f64::consts::PI / 2.0 - 1e-6,
            "center atom: suggested direction {dir:.3} should be ~90° from both bonds, sep={sep:.3}"
        );
    }

    #[test]
    fn test_suggest_direction_prefers_sp2_for_aromatic_ring() {
        use chematic_smiles::parse;
        // Benzene: use compute_layout to get real coordinates, then ask for
        // the exit direction from atom 0. It should be ≈120° from both ring bonds.
        let mol = parse("c1ccccc1").unwrap();
        let layout = compute_layout(&mol);
        let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
        // Atom 0 has 2 ring neighbors. The best exit is ~120° from both.
        let p0 = layout.get(AtomIdx(0));
        let used: Vec<f64> = mol
            .neighbors(AtomIdx(0))
            .map(|(nb, _)| {
                let p = layout.get(nb);
                (p.y - p0.y).atan2(p.x - p0.x)
            })
            .collect();
        let sep = min_angle_separation(dir, &used);
        // Minimum separation from both ring bonds should be ≥ 60°.
        assert!(
            sep >= std::f64::consts::PI / 3.0 - 1e-6,
            "benzene exit direction should be ≥60° from ring bonds, sep={sep:.3}"
        );
    }
}
