//! 2D coordinate generation for molecular depiction.
//!
//! The layout algorithm is rule-based and produces SVG pixel coordinates.
//! No physics simulation is used; atoms are placed with geometric rules.

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, Molecule};
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
                (min_x.min(p.x), min_y.min(p.y), max_x.max(p.x), max_y.max(p.y))
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
pub fn compute_layout(mol: &Molecule) -> Layout {
    let n = mol.atom_count();
    if n == 0 {
        return Layout { coords: Vec::new() };
    }

    // Special case: single atom.
    if n == 1 {
        return Layout { coords: vec![Point::new(0.0, 0.0)] };
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

        // Determine which atoms are in at least one ring.
        let in_ring: HashSet<AtomIdx> = rings.iter().flat_map(|r| r.iter().copied()).collect();

        // Group rings into ring systems (connected sets of rings sharing >= 1 atom).
        let ring_systems = group_ring_systems(&rings);

        // Place each ring system.
        for system in &ring_systems {
            place_ring_system(mol, system, &mut placed);
        }

        // Place chain atoms (DFS zigzag from placed ring atoms or from scratch).
        place_chains(mol, &component_set, &in_ring, &mut placed);

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
    for i in 0..n {
        let set_i: HashSet<AtomIdx> = rings[i].iter().copied().collect();
        for j in (i + 1)..n {
            if rings[j].iter().any(|a| set_i.contains(a)) {
                union(&mut parent, i, j);
            }
        }
    }

    // Collect into groups.
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
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
fn place_ring_system(_mol: &Molecule, system: &[Vec<AtomIdx>], placed: &mut Vec<Option<Point>>) {
    if system.is_empty() {
        return;
    }

    // Place the first ring as a regular polygon.
    place_regular_ring(&system[0], placed);

    // For subsequent rings: find two atoms already placed (the shared edge),
    // then reflect unplaced atoms over that shared edge.
    let mut remaining: Vec<&Vec<AtomIdx>> = system[1..].iter().collect();
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

            if already_placed.len() < 2 {
                return true; // Not ready yet.
            }

            // Find the shared edge: two consecutive atoms in the ring that are both placed.
            let shared_edge = find_shared_edge(ring, placed);

            // Fall back: use the first two placed atoms.
            let (anchor1, anchor2) = shared_edge
                .unwrap_or((already_placed[0], already_placed[1]));

            // Both anchors are confirmed placed (either from find_shared_edge or already_placed).
            let (Some(p1), Some(p2)) =
                (placed[anchor1.0 as usize], placed[anchor2.0 as usize])
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
        .find(|&(a, b)| {
            placed[a.0 as usize].is_some() && placed[b.0 as usize].is_some()
        })
}

/// Place atoms of a ring as a regular polygon centered at the origin.
fn place_regular_ring(ring: &[AtomIdx], placed: &mut Vec<Option<Point>>) {
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
    placed: &mut Vec<Option<Point>>,
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

    // Centroid of already-placed atoms in this ring (from the old ring system).
    let existing_center = {
        let pts: Vec<Point> = ring.iter().filter_map(|&a| placed[a.0 as usize]).collect();
        if pts.is_empty() {
            // No placed atoms to compare against: use the midpoint as fallback.
            mid
        } else {
            let cx = pts.iter().map(|p| p.x).sum::<f64>() / pts.len() as f64;
            let cy = pts.iter().map(|p| p.y).sum::<f64>() / pts.len() as f64;
            Point::new(cx, cy)
        }
    };

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
    while delta > std::f64::consts::PI { delta -= 2.0 * std::f64::consts::PI; }
    while delta < -std::f64::consts::PI { delta += 2.0 * std::f64::consts::PI; }

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

/// Place all atoms not yet placed (chain atoms) using DFS zigzag.
///
/// Starts from placed ring atoms that have unplaced neighbors, then
/// handles any still-unplaced atoms (isolated chains with no ring).
fn place_chains(
    mol: &Molecule,
    component: &HashSet<AtomIdx>,
    in_ring: &HashSet<AtomIdx>,
    placed: &mut Vec<Option<Point>>,
) {
    // Step 1: extend chains from ring attachment points.
    // We process ring atoms that have unplaced neighbors.
    let ring_atoms: Vec<AtomIdx> = component
        .iter()
        .copied()
        .filter(|a| in_ring.contains(a) && placed[a.0 as usize].is_some())
        .collect();

    for start in &ring_atoms {
        let unplaced_neighbors: Vec<AtomIdx> = mol
            .neighbors(*start)
            .map(|(nb, _)| nb)
            .filter(|nb| component.contains(nb) && placed[nb.0 as usize].is_none())
            .collect();

        for nb in unplaced_neighbors {
            if placed[nb.0 as usize].is_some() {
                continue;
            }
            // Determine outgoing direction from the ring atom.
            // Use a direction that avoids existing neighbors.
            let dir = best_outgoing_direction(*start, placed, mol, component);
            dfs_zigzag(mol, nb, *start, dir, placed, component);
        }
    }

    // Step 2: handle pure chain components (no ring atoms yet placed).
    // Find a terminal atom (degree 1 or degree 0) as the starting point.
    let unplaced: Vec<AtomIdx> = component
        .iter()
        .copied()
        .filter(|a| placed[a.0 as usize].is_none())
        .collect();

    if unplaced.is_empty() {
        return;
    }

    // Find a terminal atom in this component to start DFS.
    let start = unplaced
        .iter()
        .copied()
        .find(|&a| mol.degree(a) <= 1)
        .unwrap_or(unplaced[0]);

    placed[start.0 as usize] = Some(Point::new(0.0, 0.0));

    // DFS rightward at 0 degrees.
    let init_dir = 0.0_f64;
    let neighbors: Vec<AtomIdx> = mol
        .neighbors(start)
        .map(|(nb, _)| nb)
        .filter(|nb| component.contains(nb) && placed[nb.0 as usize].is_none())
        .collect();

    for (i, nb) in neighbors.into_iter().enumerate() {
        if placed[nb.0 as usize].is_some() {
            continue;
        }
        let dir = if i == 0 { init_dir } else { init_dir + std::f64::consts::PI };
        dfs_zigzag(mol, nb, start, dir, placed, component);
    }

    // Any remaining unplaced atoms: place them in a line.
    let still_unplaced: Vec<AtomIdx> = component
        .iter()
        .copied()
        .filter(|a| placed[a.0 as usize].is_none())
        .collect();

    let mut x = 0.0;
    for atom in still_unplaced {
        x += BOND_LEN;
        placed[atom.0 as usize] = Some(Point::new(x, 0.0));
    }
}

/// Compute the best outgoing direction from a ring atom to avoid collisions.
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

    // Try candidate angles at 60-degree increments and pick the one farthest from all used.
    let candidates: Vec<f64> = (0..6)
        .map(|i| i as f64 * std::f64::consts::PI / 3.0)
        .collect();

    candidates
        .into_iter()
        .max_by(|&a, &b| {
            let min_sep_a = min_angle_separation(a, &used_angles);
            let min_sep_b = min_angle_separation(b, &used_angles);
            min_sep_a.partial_cmp(&min_sep_b).unwrap()
        })
        .unwrap_or(0.0)
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

/// DFS zigzag placement: place `atom` at `BOND_LEN` from `parent` in direction `dir`,
/// then recurse on unplaced neighbors with alternating ±30° deflection.
fn dfs_zigzag(
    mol: &Molecule,
    atom: AtomIdx,
    parent: AtomIdx,
    dir: f64,
    placed: &mut Vec<Option<Point>>,
    component: &HashSet<AtomIdx>,
) {
    if placed[atom.0 as usize].is_some() {
        return;
    }

    let Some(parent_pos) = placed[parent.0 as usize] else {
        return;
    };

    let x = parent_pos.x + BOND_LEN * dir.cos();
    let y = parent_pos.y + BOND_LEN * dir.sin();
    placed[atom.0 as usize] = Some(Point::new(x, y));

    // Collect unplaced neighbors (excluding parent).
    let unplaced: Vec<AtomIdx> = mol
        .neighbors(atom)
        .map(|(nb, _)| nb)
        .filter(|&nb| nb != parent && component.contains(&nb) && placed[nb.0 as usize].is_none())
        .collect();

    // Zigzag: first neighbor turns +30°, second turns -30°, alternating.
    let deflections = [
        -std::f64::consts::PI / 6.0,
        std::f64::consts::PI / 6.0,
    ];

    for (i, nb) in unplaced.into_iter().enumerate() {
        if placed[nb.0 as usize].is_some() {
            continue;
        }
        let deflection = deflections[i % 2];
        let new_dir = dir + deflection;
        dfs_zigzag(mol, nb, atom, new_dir, placed, component);
    }
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_radius_hexagon() {
        // Regular hexagon: all sides == BOND_LEN, circumradius == BOND_LEN.
        let r = ring_radius(6);
        let expected = BOND_LEN; // sin(PI/6) = 0.5, so BOND_LEN / (2*0.5) = BOND_LEN
        assert!((r - expected).abs() < 1e-9, "hexagon radius = {}", r);
    }
}
