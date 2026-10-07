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

use std::collections::{BTreeMap, VecDeque};

// Output never depends on these maps' iteration order (the layout sorts
// where order matters; see `compute_layout_is_deterministic_across_repeated_calls`),
// so the faster non-keyed hasher gives the same coordinates.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use chematic_core::{AtomIdx, BondIdx, Molecule};

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
        // Caller-supplied layouts may intentionally omit trailing atoms;
        // `depict_data_with_coords` documents that those atoms are placed at
        // the origin. Keep every renderer on that same fail-safe contract.
        self.coords
            .get(idx.0 as usize)
            .copied()
            .unwrap_or_else(|| Point::new(0.0, 0.0))
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

    // One SSSR for the whole molecule (memoized on it), filtered per component.
    let ring_set = chematic_perception::find_sssr_shared(mol);
    for component_atoms in &components {
        let component_set: HashSet<AtomIdx> = component_atoms.iter().copied().collect();

        // The rings of this component.
        let rings: Vec<Vec<AtomIdx>> = ring_set
            .rings()
            .iter()
            .filter(|ring| ring.iter().all(|a| component_set.contains(a)))
            .cloned()
            .collect();

        // Group rings into ring systems (connected sets of rings sharing >= 1 atom).
        let ring_systems = group_ring_systems(&rings);

        let mut atom_to_system: HashMap<AtomIdx, usize> = HashMap::default();
        for (sys_idx, system) in ring_systems.iter().enumerate() {
            for &a in system.iter().flatten() {
                atom_to_system.insert(a, sys_idx);
            }
        }
        // Lays the component out with forks 60° either side of the chain
        // direction, or `forks[atom]` where that is not zero.
        let lay_forks = |forks: &[f64]| -> Vec<Option<Point>> {
            let mut placed: Vec<Option<Point>> = vec![None; n];
            let mut system_placed: Vec<bool> = vec![false; ring_systems.len()];

            // Seed: the ring system that anchors this component's whole
            // coordinate frame (see `seed_ring_system_index`'s doc). Every
            // other ring system gets discovered and anchored to its real
            // attachment point as the layout grows outward -- placing every
            // ring system blind at the origin (the pre-fix behavior) makes
            // unrelated ring systems of the same size collide exactly.
            if let Some(seed_idx) = seed_ring_system_index(&ring_systems) {
                place_ring_system(mol, &ring_systems[seed_idx], None, &mut placed);
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
                NarrowForks { angles: forks },
            );

            // Defensive fallbacks -- should not fire for any connected
            // component, since grow_layout's worklist reaches every atom
            // reachable from the seed via the molecule graph.
            for (sys_idx, system) in ring_systems.iter().enumerate() {
                if !system_placed[sys_idx] {
                    place_ring_system(mol, system, None, &mut placed);
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

            placed
        };
        let lay_with = |narrow: &[bool], angle: f64| {
            let forks: Vec<f64> = narrow
                .iter()
                .map(|&b| if b { angle } else { 0.0 })
                .collect();
            lay_forks(&forks)
        };
        let lay = |narrow: &[bool]| lay_with(narrow, std::f64::consts::PI / 6.0);
        // Forks are drawn 120° apart; where that leaves a clash or a
        // crossing the component is laid out again with the narrower ±30°
        // forks of earlier releases and the drawing with fewer defects kept
        // (a crowded peptide or sulfonamide can need the narrow ones). In a
        // component of up to 60 atoms that takes the narrow forks, only the
        // forks near the wide drawing's defects stay narrow when that does
        // as well.
        let mut placed = lay(&vec![false; n]);
        let comp_bonds = component_bonds(mol, component_atoms);
        // Clashes plus crossings of a drawing of this component.
        let defect_counts = |placed: &[Option<Point>]| {
            layout_defects(component_atoms, &comp_bonds, &Bonded(mol), placed)
        };
        let defects = |placed: &[Option<Point>]| {
            let (clashes, crossings) = defect_counts(placed);
            clashes + crossings
        };
        let wide_defects = defects(&placed);
        if wide_defects > 0 {
            let mut narrow = lay(&vec![true; n]);
            let narrow_raw = defects(&narrow);
            let mut narrow_chosen = false;
            if component_atoms.len() <= 60 {
                relieve_clashes(mol, component_atoms, &mut placed, false);
                relieve_clashes(mol, component_atoms, &mut narrow, false);
                if defects(&narrow) < defects(&placed) {
                    placed = narrow;
                    narrow_chosen = true;
                }
            } else {
                if narrow_raw < wide_defects {
                    placed = narrow;
                    narrow_chosen = true;
                }
                relieve_clashes(mol, component_atoms, &mut placed, false);
            }
            // Where the narrow forks won, keep the wide ones away from the
            // wide drawing's defects when that does as well.
            if narrow_chosen {
                let wide = lay(&vec![false; n]);
                let near = atoms_near_defects(mol, component_atoms, &wide);
                let mut local = lay(&near);
                let mut mask = vec![true; n];
                // A large component (a long peptide) keeps the local drawing
                // without relief when it needs none.
                if component_atoms.len() > 60 {
                    if defects(&local) == 0 {
                        placed = local;
                        mask = near;
                    }
                } else {
                    relieve_clashes(mol, component_atoms, &mut local, false);
                    if defects(&local) <= defects(&placed) {
                        placed = local;
                        mask = near;
                    }
                }
                // The narrow forks opened to ±45° (90° between the two
                // branches) when that does as well, after relief (large
                // components included, now that the relief stops counting a
                // move once it cannot win).
                let mut medium = lay_with(&mask, std::f64::consts::PI / 4.0);
                if defects(&medium) > 0 {
                    relieve_clashes(mol, component_atoms, &mut medium, false);
                }
                // Kept when it has no more defects and no more clashes.
                let (m, p) = (defect_counts(&medium), defect_counts(&placed));
                if m.0 <= p.0 && m.0 + m.1 <= p.0 + p.1 {
                    placed = medium;
                }
            }
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

/// Runs the layout's clash relief (see `relieve_clashes`) again on finished
/// coordinates, per connected component: for a caller that has moved atoms
/// after [`compute_layout`] (the MOL writer's E/Z reflections). No move
/// mirrors across a double bond and turns only pivot on atoms without one,
/// so every double bond's drawn geometry stays.
pub fn relieve_layout_clashes(mol: &Molecule, layout: &mut Layout) {
    if layout.coords.len() < mol.atom_count() {
        return;
    }
    let mut placed: Vec<Option<Point>> = layout.coords.iter().copied().map(Some).collect();
    for component in connected_components(mol) {
        relieve_clashes(mol, &component, &mut placed, true);
    }
    for (slot, p) in layout.coords.iter_mut().zip(placed) {
        if let Some(p) = p {
            *slot = p;
        }
    }
}

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
    mol: &Molecule,
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
    // Bounded by the starting count: `remaining` shrinks as rings are
    // placed, and a bound on its current length stopped a fused system
    // before its last ring could anchor (that ring then fell to the
    // origin-centred fallback with stretched bonds).
    let max_iterations = remaining.len() * 2;

    while !remaining.is_empty() && iterations < max_iterations {
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

            // A bridged ring can share several already-placed atoms without
            // sharing a consecutive pair in this ring's traversal order.
            // Choosing the first two then anchors the regular polygon to a
            // long diagonal and stretches one or more real bonds. Use the
            // closest placed pair as the geometric fallback; a true shared
            // bond is the shortest available pair and this remains
            // deterministic for equal distances.
            let (anchor1, anchor2) =
                shared_edge.unwrap_or_else(|| closest_placed_pair(&already_placed, placed));

            // Both anchors are confirmed placed (either from find_shared_edge or already_placed).
            let (Some(p1), Some(p2)) = (placed[anchor1.0 as usize], placed[anchor2.0 as usize])
            else {
                return true; // Not ready.
            };

            // A bridged ring whose placed atoms are one path of two or more
            // bonds: its unplaced atoms go on the regular polygon through
            // that path's ends that lands on fewer placed atoms (norbornane's
            // second ring wraps round the one-atom bridge instead of
            // retracing the first ring); otherwise the shared-edge placement.
            if !place_bridged_ring(ring, placed) {
                place_ring_anchored(ring, anchor1, p1, anchor2, p2, placed);
            }

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

    untangle_bridged_system(mol, system, anchor.map(|(a, ..)| a), placed);
}

/// A bridged ring system (two rings sharing three or more atoms), or a
/// fused or spiro system of three or more rings, whose polygon placement
/// leaves non-bonded atoms closer than half a bond (adamantane, morphinans,
/// cages fused to several rings, helicenes, a spiro naphthodioxin, a
/// chelate's fused rings) is redrawn by stress majorization on its
/// ring-graph distances, started from the polygon drawing and turned back
/// onto it, and kept only when it has fewer such pairs.
fn untangle_bridged_system(
    mol: &Molecule,
    system: &[Vec<AtomIdx>],
    anchor: Option<AtomIdx>,
    placed: &mut [Option<Point>],
) {
    let bridged = system.iter().enumerate().any(|(i, r)| {
        system[i + 1..]
            .iter()
            .any(|q| r.iter().filter(|a| q.contains(a)).count() >= 3)
    });
    if !bridged && system.len() < 3 {
        return;
    }
    let mut atoms: Vec<AtomIdx> = system.iter().flatten().copied().collect();
    atoms.sort_unstable();
    atoms.dedup();
    let n = atoms.len();
    if n > 80 || atoms.iter().any(|a| placed[a.0 as usize].is_none()) {
        return;
    }
    let pos = |a: AtomIdx| atoms.binary_search(&a).unwrap();
    let mut adj = vec![Vec::new(); n];
    for ring in system {
        for k in 0..ring.len() {
            let (i, j) = (pos(ring[k]), pos(ring[(k + 1) % ring.len()]));
            if !adj[i].contains(&j) {
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }
    // Ring systems with a template (adamantane-type cages, porphyrins) are
    // drawn from it, turned onto the polygon drawing. An anchored one keeps
    // its entry atom where the bond to the already-placed chain put it.
    // (Stress-majorized systems below keep the centred fit: shifting them
    // onto their entry atom stacked a second artemisinin unit on the first.)
    {
        let old: Vec<Point> = atoms
            .iter()
            .map(|a| placed[a.0 as usize].unwrap())
            .collect();
        if let Some(t) =
            adamantane_cage(mol, &atoms, &adj, &old).or_else(|| porphyrin_core(&adj, &old))
        {
            let mut aligned = turn_onto(&t, &old);
            if let Some(k) = anchor.and_then(|a| atoms.binary_search(&a).ok()) {
                let (sx, sy) = (old[k].x - aligned[k].x, old[k].y - aligned[k].y);
                for p in &mut aligned {
                    *p = Point::new(p.x + sx, p.y + sy);
                }
            }
            for (k, a) in atoms.iter().enumerate() {
                placed[a.0 as usize] = Some(aligned[k]);
            }
            return;
        }
    }
    // No overlapping pair (the usual case): nothing to redraw, and the
    // distance matrix below is not needed.
    let any_close = (0..n).any(|i| {
        let p = placed[atoms[i].0 as usize].unwrap();
        (i + 1..n).any(|j| {
            !adj[i].contains(&j) && p.dist(&placed[atoms[j].0 as usize].unwrap()) < 0.5 * BOND_LEN
        })
    });
    if !any_close {
        return;
    }
    // Graph distances.
    let mut dist = vec![vec![usize::MAX; n]; n];
    for (s, row) in dist.iter_mut().enumerate() {
        row[s] = 0;
        let mut queue = VecDeque::from([s]);
        while let Some(v) = queue.pop_front() {
            for &w in &adj[v] {
                if row[w] == usize::MAX {
                    row[w] = row[v] + 1;
                    queue.push_back(w);
                }
            }
        }
    }
    // How far non-bonded pairs fall inside half a bond, summed.
    let close_pairs = |x: &[Point]| {
        let mut overlap = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                if dist[i][j] > 1 {
                    overlap += (0.5 * BOND_LEN - x[i].dist(&x[j])).max(0.0);
                }
            }
        }
        overlap
    };
    let old: Vec<Point> = atoms
        .iter()
        .map(|a| placed[a.0 as usize].unwrap())
        .collect();
    let before = close_pairs(&old);
    if before == 0.0 {
        return;
    }
    let target = |i: usize, j: usize| dist[i][j] as f64 * BOND_LEN;
    // Pair targets (graph distance in bonds) and d^-2 weights, flattened;
    // a weight of 0 marks a pair with no path.
    let mut targets = vec![0.0; n * n];
    let mut base_w = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            if i != j && dist[i][j] != usize::MAX {
                let d = target(i, j);
                targets[i * n + j] = d;
                base_w[i * n + j] = 1.0 / (d * d);
            }
        }
    }
    let smacof = |mut x: Vec<Point>, bond_weight: f64, rounds: usize| -> Vec<Point> {
        // Stress majorization, one atom at a time (weights d^-2, bonds
        // `bond_weight` times that).
        // The pair weights, computed once per run (the same products the
        // loop formed for every pair in every round).
        let weights: Vec<f64> = (0..n * n)
            .map(|k| {
                let w = base_w[k];
                if dist[k / n][k % n] == 1 {
                    w * bond_weight
                } else {
                    w
                }
            })
            .collect();
        for _ in 0..rounds {
            let mut moved: f64 = 0.0;
            for i in 0..n {
                let (mut sx, mut sy, mut sw) = (0.0, 0.0, 0.0);
                let row = &weights[i * n..(i + 1) * n];
                for j in 0..n {
                    let w = row[j];
                    if w == 0.0 {
                        continue;
                    }
                    let d = targets[i * n + j];
                    let (dx, dy) = (x[i].x - x[j].x, x[i].y - x[j].y);
                    let len = (dx * dx + dy * dy).sqrt();
                    let (ux, uy) = if len > 1e-9 {
                        (dx / len, dy / len)
                    } else {
                        // Coincident atoms: separate along a fixed direction
                        // depending only on the pair's order.
                        let angle = (i * 7 + j * 13) as f64;
                        (angle.cos(), angle.sin())
                    };
                    sx += w * (x[j].x + d * ux);
                    sy += w * (x[j].y + d * uy);
                    sw += w;
                }
                if sw > 0.0 {
                    let next = Point::new(sx / sw, sy / sw);
                    moved = moved.max(next.dist(&x[i]));
                    x[i] = next;
                }
            }
            // Stop once no atom moves more than 1e-4 bond lengths.
            if moved < 1e-4 * BOND_LEN {
                break;
            }
        }
        x
    };
    let stress = |x: &[Point]| {
        let mut total = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                if dist[i][j] != usize::MAX {
                    let d = target(i, j);
                    total += (x[i].dist(&x[j]) - d).powi(2) / (d * d);
                }
            }
        }
        total
    };
    // Starts: the polygon drawing, and classical MDS of the graph distances
    // (two leading eigenvectors by power iteration), which does not inherit
    // the polygon drawing's folds.
    let mut starts = vec![old.clone()];
    {
        let mut b = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let d = if dist[i][j] == usize::MAX {
                    0.0
                } else {
                    target(i, j)
                };
                b[i][j] = -0.5 * d * d;
            }
        }
        let row: Vec<f64> = b.iter().map(|r| r.iter().sum::<f64>() / n as f64).collect();
        let all = row.iter().sum::<f64>() / n as f64;
        for i in 0..n {
            for j in 0..n {
                b[i][j] += all - row[i] - row[j];
            }
        }
        let mut vecs: Vec<(f64, Vec<f64>)> = Vec::new();
        for k in 0..2 {
            let mut v: Vec<f64> = (0..n).map(|i| 1.0 + ((i * (k + 3)) % 7) as f64).collect();
            let mut lambda = 0.0;
            for _ in 0..200 {
                let mut w: Vec<f64> = (0..n)
                    .map(|i| (0..n).map(|j| b[i][j] * v[j]).sum())
                    .collect();
                for (l, u) in &vecs {
                    let dot: f64 = (0..n).map(|i| w[i] * u[i]).sum();
                    for i in 0..n {
                        w[i] -= dot * u[i];
                    }
                    let _ = l;
                }
                let norm = w.iter().map(|t| t * t).sum::<f64>().sqrt();
                if norm < 1e-12 {
                    break;
                }
                lambda = norm;
                v = w.iter().map(|t| t / norm).collect();
            }
            vecs.push((lambda, v));
        }
        if vecs.iter().all(|(l, _)| *l > 1e-9) {
            starts.push(
                (0..n)
                    .map(|i| {
                        Point::new(
                            vecs[0].1[i] * vecs[0].0.sqrt(),
                            vecs[1].1[i] * vecs[1].0.sqrt(),
                        )
                    })
                    .collect(),
            );
        }
    }
    let mut x = old.clone();
    {
        let mut best_key = (f64::MAX, f64::MAX);
        for start in starts {
            let result = smacof(start, 1.0, 60);
            let key = (close_pairs(&result), stress(&result));
            if key.0 < best_key.0 || (key.0 == best_key.0 && key.1 < best_key.1) {
                best_key = key;
                x = result;
            }
        }
        // Then even out bond lengths, keeping the untangled arrangement when
        // that adds no overlap.
        let bond_spread = |x: &[Point]| {
            let mut worst: f64 = 0.0;
            for i in 0..n {
                for &j in &adj[i] {
                    worst = worst.max((x[i].dist(&x[j]) / BOND_LEN - 1.0).abs());
                }
            }
            worst
        };
        let refined = smacof(x.clone(), 8.0, 40);
        if close_pairs(&refined) <= best_key.0 && bond_spread(&refined) < bond_spread(&x) {
            x = refined;
        }
    }
    // Turn (no mirror) and shift back onto the polygon drawing.
    let aligned = turn_onto(&x, &old);
    if close_pairs(&aligned) < before {
        for (k, a) in atoms.iter().enumerate() {
            placed[a.0 as usize] = Some(aligned[k]);
        }
    }
}

/// The drawing of a porphyrin core: a ring system of four five-membered
/// rings, each joined to the next through one atom (the meso carbons), with
/// optionally one atom bonded to the four inner atoms (a metal). Each
/// five-membered ring is a regular pentagon with its inner atom (N) facing
/// the centre, every bond one bond length and the angles at the ring-fusion
/// carbons and the meso carbons 126°; a central atom sits at the centre.
/// The polygon placement and stress majorization drew the pentagons folded
/// inside the macrocycle. Of the drawing and its mirror image, the one
/// closer to `old` is returned. `None` for any other system.
fn porphyrin_core(adj: &[Vec<usize>], old: &[Point]) -> Option<Vec<Point>> {
    let n = adj.len();
    if n != 24 && n != 25 {
        return None;
    }
    // A central atom: bonded to exactly four atoms, each of which has two
    // other neighbours.
    let centre_atom = if n == 25 {
        let c = (0..n).find(|&i| adj[i].len() == 4)?;
        if adj[c].iter().any(|&j| adj[j].len() != 3) {
            return None;
        }
        Some(c)
    } else {
        None
    };
    // Degrees within the core (without the central atom).
    let core_nbs = |i: usize| -> Vec<usize> {
        adj[i]
            .iter()
            .copied()
            .filter(|&j| Some(j) != centre_atom)
            .collect()
    };
    let deg = |i: usize| core_nbs(i).len();
    // Pyrroles: an inner atom of core degree 2 whose two neighbours (core
    // degree 3) each have a degree-2 neighbour, those two bonded.
    let mut pyrroles: Vec<[usize; 5]> = Vec::new(); // [inner, a, b, beta_a, beta_b]
    let mut used = vec![false; n];
    for x in (0..n).filter(|&i| Some(i) != centre_atom && deg(i) == 2) {
        let nb = core_nbs(x);
        let (a, b) = (nb[0], nb[1]);
        if deg(a) != 3 || deg(b) != 3 {
            continue;
        }
        let mut found = None;
        for ba in core_nbs(a).into_iter().filter(|&t| t != x && deg(t) == 2) {
            for bb in core_nbs(b).into_iter().filter(|&t| t != x && deg(t) == 2) {
                if ba != bb && adj[ba].contains(&bb) {
                    found = Some((ba, bb));
                }
            }
        }
        let Some((ba, bb)) = found else { continue };
        let ring = [x, a, b, ba, bb];
        if ring.iter().any(|&t| used[t]) {
            return None;
        }
        for &t in &ring {
            used[t] = true;
        }
        pyrroles.push(ring);
    }
    if pyrroles.len() != 4 {
        return None;
    }
    if let Some(c) = centre_atom {
        // The central atom is bonded to the four inner atoms.
        let mut inner: Vec<usize> = pyrroles.iter().map(|r| r[0]).collect();
        let mut nbs = adj[c].clone();
        inner.sort_unstable();
        nbs.sort_unstable();
        if inner != nbs {
            return None;
        }
    }
    // The meso atom beside a fusion carbon: its remaining core neighbour.
    let meso_of = |p: &[usize; 5], side: usize| -> Option<usize> {
        let (fusion, beta) = if side == 0 {
            (p[1], p[3])
        } else {
            (p[2], p[4])
        };
        core_nbs(fusion)
            .into_iter()
            .find(|&t| t != p[0] && t != beta)
    };
    // Walk round the macrocycle: pyrrole k's side-1 fusion carbon, a meso
    // atom, then the next pyrrole's side-0 fusion carbon.
    let mut order: Vec<([usize; 5], bool)> = vec![(pyrroles[0], false)];
    let mut mesos = Vec::new();
    for _ in 0..4 {
        let (p, flipped) = *order.last().unwrap();
        let m = meso_of(&p, if flipped { 0 } else { 1 })?;
        if deg(m) != 2 || used[m] {
            return None;
        }
        let other = core_nbs(m)
            .into_iter()
            .find(|&t| t != p[if flipped { 1 } else { 2 }])?;
        let next = *pyrroles.iter().find(|q| q[1] == other || q[2] == other)?;
        mesos.push(m);
        order.push((next, next[2] == other));
    }
    if order[4].0 != order[0].0 || order[4].1 != order[0].1 {
        return None;
    }
    let mut seen_meso = mesos.clone();
    seen_meso.sort_unstable();
    seen_meso.dedup();
    if seen_meso.len() != 4 {
        return None;
    }
    // Coordinates in bond lengths.
    let r5 = 1.0 / (2.0 * 36f64.to_radians().sin());
    let (c72, s72) = (72f64.to_radians().cos(), 72f64.to_radians().sin());
    let (c36, s36) = (36f64.to_radians().cos(), 36f64.to_radians().sin());
    // Pentagon centre distance: the meso atom lies on the exterior bisector
    // of the fusion carbon and on the 45° line.
    let d = r5 * c72 + c72 + r5 * s72 + s72;
    let draw = |mirror: f64| -> Vec<Point> {
        let mut x = vec![Point::new(0.0, 0.0); n];
        for (k, &(p, flipped)) in order[..4].iter().enumerate() {
            let theta = k as f64 * std::f64::consts::FRAC_PI_2;
            let (ux, uy) = (theta.cos(), theta.sin());
            let (vx, vy) = (-uy, ux);
            let at = |a: f64, b: f64| {
                Point::new(
                    (a * ux + b * vx) * BOND_LEN,
                    mirror * (a * uy + b * vy) * BOND_LEN,
                )
            };
            // Side 1 (the walk's exit) is +v.
            let (s0, s1) = if flipped { (1.0, -1.0) } else { (-1.0, 1.0) };
            x[p[0]] = at(d - r5, 0.0);
            x[p[1]] = at(d - r5 * c72, s0 * r5 * s72);
            x[p[2]] = at(d - r5 * c72, s1 * r5 * s72);
            x[p[3]] = at(d + r5 * c36, s0 * r5 * s36);
            x[p[4]] = at(d + r5 * c36, s1 * r5 * s36);
            let m = mesos[k];
            let phi = theta + std::f64::consts::FRAC_PI_4;
            let rm = std::f64::consts::SQRT_2 * (d - r5 * c72 - c72);
            x[m] = Point::new(
                rm * phi.cos() * BOND_LEN,
                mirror * rm * phi.sin() * BOND_LEN,
            );
        }
        x
    };
    let fit = |x: &[Point]| {
        let t = turn_onto(x, old);
        t.iter()
            .zip(old)
            .map(|(a, b)| a.dist(b).powi(2))
            .sum::<f64>()
    };
    let (a, b) = (draw(1.0), draw(-1.0));
    Some(if fit(&b) < fit(&a) { b } else { a })
}

/// `x` turned (no mirror) about its centroid and shifted onto `old`'s
/// centroid, at the turn that best fits it to `old`.
fn turn_onto(x: &[Point], old: &[Point]) -> Vec<Point> {
    let n = x.len();
    let centre = |pts: &[Point]| {
        let (sx, sy) = pts.iter().fold((0.0, 0.0), |(a, b), p| (a + p.x, b + p.y));
        Point::new(sx / n as f64, sy / n as f64)
    };
    let (co, cx) = (centre(old), centre(x));
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..n {
        let (ax, ay) = (x[i].x - cx.x, x[i].y - cx.y);
        let (bx, by) = (old[i].x - co.x, old[i].y - co.y);
        num += ax * by - ay * bx;
        den += ax * bx + ay * by;
    }
    let theta = num.atan2(den);
    let (sin, cos) = theta.sin_cos();
    x.iter()
        .map(|p| {
            let (dx, dy) = (p.x - cx.x, p.y - cx.y);
            Point::new(co.x + dx * cos - dy * sin, co.y + dx * sin + dy * cos)
        })
        .collect()
}

/// RDKit's drawing of adamantane in bond lengths: the four bridgeheads, then
/// the atom joining each pair of them. The atoms joining bridgeheads 0 and 2
/// and bridgeheads 3 and 1 are inside the outline.
const ADAMANTANE_BRIDGEHEADS: [(f64, f64); 4] = [
    (0.0640, 0.8214),
    (1.6135, 0.0473),
    (0.1684, -0.9075),
    (-1.3811, -0.1334),
];
const ADAMANTANE_LINKS: [((usize, usize), (f64, f64)); 6] = [
    ((0, 2), (-0.3829, -0.0732)),
    ((0, 1), (1.0622, 0.8816)),
    ((1, 2), (1.1665, -0.8473)),
    ((2, 3), (-0.8298, -0.9678)),
    ((3, 0), (-0.9341, 0.7611)),
    ((3, 1), (-0.5467, 0.4178)),
];

/// The drawing of an adamantane-type cage (ten atoms: four bridgeheads, each
/// pair joined through one atom; adamantane, hexamine), which no polygon
/// placement or stress majorization draws without a clash: RDKit's
/// projection, with the bridgeheads assigned (of the 24 assignments and their
/// mirror images) so that no atom with a substituent lies inside the outline,
/// then closest to the current drawing `old`. `None` for any other system.
fn adamantane_cage(
    mol: &Molecule,
    atoms: &[AtomIdx],
    adj: &[Vec<usize>],
    old: &[Point],
) -> Option<Vec<Point>> {
    let n = atoms.len();
    if n != 10 {
        return None;
    }
    let heads: Vec<usize> = (0..n).filter(|&i| adj[i].len() == 3).collect();
    if heads.len() != 4 || (0..n).any(|i| adj[i].len() != 3 && adj[i].len() != 2) {
        return None;
    }
    // link[a][b]: the atom joining heads a and b.
    let mut link = [[usize::MAX; 4]; 4];
    for i in (0..n).filter(|&i| adj[i].len() == 2) {
        let a = heads.iter().position(|&h| h == adj[i][0])?;
        let b = heads.iter().position(|&h| h == adj[i][1])?;
        if a == b || link[a][b] != usize::MAX {
            return None;
        }
        link[a][b] = i;
        link[b][a] = i;
    }
    let substituted = |i: usize| mol.degree(atoms[i]) > adj[i].len();
    let centre = |pts: &[Point]| {
        let (sx, sy) = pts.iter().fold((0.0, 0.0), |(a, b), p| (a + p.x, b + p.y));
        Point::new(sx / n as f64, sy / n as f64)
    };
    let co = centre(old);
    let mut best: Option<((usize, f64), Vec<Point>)> = None;
    let mut perm = [0usize, 1, 2, 3];
    for code in 0..24 {
        // The `code`-th permutation of the four heads (factorial digits).
        let mut pool = vec![0usize, 1, 2, 3];
        let mut c = code;
        for (k, slot) in perm.iter_mut().enumerate() {
            let f = [6, 2, 1, 1][k];
            *slot = pool.remove(c / f);
            c %= f;
        }
        for mirror in [1.0, -1.0] {
            let mut x = vec![Point::new(0.0, 0.0); n];
            for (role, &(px, py)) in ADAMANTANE_BRIDGEHEADS.iter().enumerate() {
                x[heads[perm[role]]] = Point::new(mirror * px * BOND_LEN, py * BOND_LEN);
            }
            for &((a, b), (px, py)) in &ADAMANTANE_LINKS {
                x[link[perm[a]][perm[b]]] = Point::new(mirror * px * BOND_LEN, py * BOND_LEN);
            }
            let inside = [link[perm[0]][perm[2]], link[perm[3]][perm[1]]]
                .into_iter()
                .filter(|&i| substituted(i))
                .count();
            // Residual after the best turn onto `old`.
            let (mut num, mut den, mut sq) = (0.0, 0.0, 0.0);
            for i in 0..n {
                let (bx, by) = (old[i].x - co.x, old[i].y - co.y);
                num += x[i].x * by - x[i].y * bx;
                den += x[i].x * bx + x[i].y * by;
                sq += x[i].x * x[i].x + x[i].y * x[i].y + bx * bx + by * by;
            }
            let key = (inside, sq - 2.0 * num.hypot(den));
            if best
                .as_ref()
                .is_none_or(|(k, _)| key.0 < k.0 || (key.0 == k.0 && key.1 < k.1 - 1e-9))
            {
                best = Some((key, x));
            }
        }
    }
    best.map(|(_, x)| x)
}

/// Per atom, the angle each of a fork's two branches makes with the chain
/// direction when it is drawn narrower than 120° (0 for the 60° default).
#[derive(Clone, Copy)]
struct NarrowForks<'a> {
    angles: &'a [f64],
}

/// Places the unplaced atoms of `ring` when its placed atoms form one path of
/// at least three atoms in ring order and its unplaced atoms the rest: on the
/// regular polygon through the path's two ends (two centres, one on each side
/// of their chord), choosing the one whose new atoms come within half a bond
/// of fewer placed atoms, and only when it is strictly better than retracing
/// the side the path already lies on. Returns whether it placed them.
fn place_bridged_ring(ring: &[AtomIdx], placed: &mut [Option<Point>]) -> bool {
    let n = ring.len();
    let is_placed = |i: usize| placed[ring[i % n].0 as usize].is_some();
    let placed_count = (0..n).filter(|&i| is_placed(i)).count();
    if placed_count < 3 || placed_count == n {
        return false;
    }
    // The first unplaced atom after a placed one starts the unplaced run.
    let Some(start) = (0..n).find(|&i| !is_placed(i) && is_placed(i + n - 1)) else {
        return false;
    };
    let run = (0..n).take_while(|&k| !is_placed(start + k)).count();
    if run != n - placed_count {
        return false; // placed atoms are not one contiguous path
    }
    let a = ring[(start + n - 1) % n]; // placed, before the run
    let b = ring[(start + run) % n]; // placed, after the run
    let (Some(pa), Some(pb)) = (placed[a.0 as usize], placed[b.0 as usize]) else {
        return false;
    };
    let d = pa.dist(&pb);
    // a and b are `run + 1` polygon steps apart through the new atoms; the
    // polygon is scaled so they land exactly on it (tropane's six-ring on
    // the pentagon's 1,3-pair: bonds 0.93 long).
    let chord = 2.0 * ring_radius(n) * (std::f64::consts::PI * (run + 1) as f64 / n as f64).sin();
    if d < 1e-9 || chord < 1e-9 || !(0.75..=1.33).contains(&(d / chord)) {
        return false;
    }
    let radius = ring_radius(n) * d / chord;
    let mid = Point::new((pa.x + pb.x) / 2.0, (pa.y + pb.y) / 2.0);
    let h = (radius * radius - d * d / 4.0).max(0.0).sqrt();
    let (ux, uy) = ((pb.x - pa.x) / d, (pb.y - pa.y) / d);
    let (px, py) = (-uy, ux);
    let positions = |center: Point| -> Vec<Point> {
        let angle_a = (pa.y - center.y).atan2(pa.x - center.x);
        let angle_b = (pb.y - center.y).atan2(pb.x - center.x);
        // Walk from a to b the long or the short way so the run has
        // `run + 1` equal steps.
        let step = 2.0 * std::f64::consts::PI / n as f64;
        let mut ccw = angle_b - angle_a;
        while ccw < 0.0 {
            ccw += 2.0 * std::f64::consts::PI;
        }
        let want = step * (run + 1) as f64;
        let sign = if (ccw - want).abs() < (2.0 * std::f64::consts::PI - ccw - want).abs() {
            1.0
        } else {
            -1.0
        };
        (1..=run)
            .map(|k| {
                let angle = angle_a + sign * step * k as f64;
                Point::new(
                    center.x + radius * angle.cos(),
                    center.y + radius * angle.sin(),
                )
            })
            .collect()
    };
    let clashes = |pts: &[Point]| {
        pts.iter()
            .filter(|p| placed.iter().flatten().any(|q| q.dist(p) < 0.5 * BOND_LEN))
            .count()
    };
    let centers = [
        Point::new(mid.x + px * h, mid.y + py * h),
        Point::new(mid.x - px * h, mid.y - py * h),
    ];
    let options: Vec<(usize, Vec<Point>)> = centers
        .iter()
        .map(|&c| {
            let pts = positions(c);
            (clashes(&pts), pts)
        })
        .collect();
    let best = if options[0].0 <= options[1].0 { 0 } else { 1 };
    let chosen = if options[best].0 < options[1 - best].0 {
        options[best].1.clone()
    } else {
        // Both polygons retrace placed atoms (bicyclo[2.2.2]: the ends are
        // a diameter apart): bow the bridge across between the ends instead,
        // the shallowest bow on either side that clears every placed atom.
        let tie = options[best].0;
        let bowed = |bulge: f64| -> Vec<Point> {
            (1..=run)
                .map(|k| {
                    let t = k as f64 / (run + 1) as f64;
                    let off = bulge * d * (std::f64::consts::PI * t).sin();
                    Point::new(
                        pa.x + (pb.x - pa.x) * t + px * off,
                        pa.y + (pb.y - pa.y) * t + py * off,
                    )
                })
                .collect()
        };
        let Some(pts) = [0.15, -0.15, 0.25, -0.25, 0.35, -0.35]
            .iter()
            .map(|&bulge| bowed(bulge))
            .find(|pts| clashes(pts) < tie && clashes(pts) == 0)
        else {
            return false; // no gain over the ordinary placement
        };
        pts
    };
    for (k, p) in chosen.iter().enumerate() {
        placed[ring[(start + k) % n].0 as usize] = Some(*p);
    }
    true
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

fn closest_placed_pair(atoms: &[AtomIdx], placed: &[Option<Point>]) -> (AtomIdx, AtomIdx) {
    let mut best = (atoms[0], atoms[1]);
    let mut best_distance = f64::INFINITY;
    for (i, &a) in atoms.iter().enumerate() {
        for &b in &atoms[i + 1..] {
            let (Some(pa), Some(pb)) = (placed[a.0 as usize], placed[b.0 as usize]) else {
                continue;
            };
            let distance = pa.dist(&pb);
            let pair = if a <= b { (a, b) } else { (b, a) };
            let best_pair = if best.0 <= best.1 {
                (best.0, best.1)
            } else {
                (best.1, best.0)
            };
            if distance < best_distance - 1e-9
                || ((distance - best_distance).abs() <= 1e-9 && pair < best_pair)
            {
                best = (a, b);
                best_distance = distance;
            }
        }
    }
    best
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

    let cand1 = Point::new(mid.x + perp_x * apothem, mid.y + perp_y * apothem);
    let cand2 = Point::new(mid.x - perp_x * apothem, mid.y - perp_y * apothem);
    let steps_forward = (idx2 + n - idx1) % n; // steps from idx1 to idx2 in ring order.
    let candidate_geometry = |center: Point| {
        let angle_to_a1 = (p1.y - center.y).atan2(p1.x - center.x);
        let angle_to_a2 = (p2.y - center.y).atan2(p2.x - center.x);
        let mut delta = angle_to_a2 - angle_to_a1;
        while delta > std::f64::consts::PI {
            delta -= 2.0 * std::f64::consts::PI;
        }
        while delta < -std::f64::consts::PI {
            delta += 2.0 * std::f64::consts::PI;
        }
        let angle_step = if steps_forward > 0 && delta < 0.0 {
            -(2.0 * std::f64::consts::PI / n as f64)
        } else {
            2.0 * std::f64::consts::PI / n as f64
        };
        (angle_to_a1, angle_step)
    };
    let placement_error = |center: Point| {
        let (angle_to_a1, angle_step) = candidate_geometry(center);
        ring.iter()
            .enumerate()
            .filter_map(|(ring_idx, &atom)| {
                let actual = placed[atom.0 as usize]?;
                if atom == anchor1 || atom == anchor2 {
                    return None;
                }
                let steps = (ring_idx + n - idx1) % n;
                let angle = angle_to_a1 + steps as f64 * angle_step;
                let expected = Point::new(
                    center.x + radius * angle.cos(),
                    center.y + radius * angle.sin(),
                );
                Some(expected.dist(&actual))
            })
            .sum::<f64>()
    };
    let ring_has_extra_anchors = ring
        .iter()
        .any(|&atom| atom != anchor1 && atom != anchor2 && placed[atom.0 as usize].is_some());
    let new_center = if ring_has_extra_anchors {
        if placement_error(cand1) <= placement_error(cand2) {
            cand1
        } else {
            cand2
        }
    } else {
        // With only the two anchors available, choose the side whose new
        // atoms land on fewer placed atoms; when both are clear (or both
        // clash equally), the side away from the existing layout. In an
        // angular fused system the layout centroid can lie on the free side
        // of the shared edge, and the ring then lands on its neighbour.
        let clashes = |center: Point| {
            let (angle_to_a1, angle_step) = candidate_geometry(center);
            (0..n)
                .filter(|&step| placed[ring[(idx1 + step) % n].0 as usize].is_none())
                .filter(|&step| {
                    let angle = angle_to_a1 + step as f64 * angle_step;
                    let p = Point::new(
                        center.x + radius * angle.cos(),
                        center.y + radius * angle.sin(),
                    );
                    placed.iter().flatten().any(|q| q.dist(&p) < 0.5 * BOND_LEN)
                })
                .count()
        };
        let (c1, c2) = (clashes(cand1), clashes(cand2));
        let existing_center = centroid_of_placed(placed).unwrap_or(mid);
        if c1 != c2 {
            if c1 < c2 { cand1 } else { cand2 }
        } else if cand1.dist(&existing_center) > cand2.dist(&existing_center) {
            cand1
        } else {
            cand2
        }
    };
    let (angle_to_a1, angle_step) = candidate_geometry(new_center);

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
    narrow: NarrowForks<'_>,
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
            .filter(|nb| placed[nb.0 as usize].is_none())
            .collect();
        unplaced_neighbors.sort_unstable();

        // A ring atom with two ring bonds and two substituents (a
        // gem-disubstituted carbon, phenytoin's C5) puts them 30° either
        // side of the ring's outward bisector; taking the bisector for the
        // first left the second against a ring bond, crossing the
        // neighbour's substituent in five-membered rings.
        let pair_dirs: Option<[f64; 2]> = (unplaced_neighbors.len() == 2
            && atom_to_system.contains_key(&start))
        .then(|| {
            let origin = placed[start.0 as usize]?;
            let ring_nbs: Vec<Point> = mol
                .neighbors(start)
                .filter_map(|(nb, _)| placed[nb.0 as usize])
                .collect();
            if ring_nbs.len() != 2 {
                return None;
            }
            let (mx, my) = (
                (ring_nbs[0].x + ring_nbs[1].x) / 2.0,
                (ring_nbs[0].y + ring_nbs[1].y) / 2.0,
            );
            let (dx, dy) = (origin.x - mx, origin.y - my);
            if dx * dx + dy * dy < 1e-6 {
                return None;
            }
            let b = dy.atan2(dx);
            // Split the exterior angle into three equal gaps.
            let inner = {
                let (ax, ay) = (ring_nbs[0].x - origin.x, ring_nbs[0].y - origin.y);
                let (cx, cy) = (ring_nbs[1].x - origin.x, ring_nbs[1].y - origin.y);
                (ax * cy - ay * cx).abs().atan2(ax * cx + ay * cy)
            };
            let d = (2.0 * std::f64::consts::PI - inner) / 6.0;
            Some([b - d, b + d])
        })
        .flatten();
        for (k, nb) in unplaced_neighbors.into_iter().enumerate() {
            if placed[nb.0 as usize].is_some() {
                continue; // Placed by an earlier neighbor this same pass (e.g. a shared spiro/fused atom).
            }
            // Determine outgoing direction from the already-placed atom.
            // Use a direction that avoids existing neighbors.
            let dir = match pair_dirs {
                Some(dirs) if direction_is_free(start, dirs[k], placed) => dirs[k],
                _ => best_outgoing_direction(start, placed, mol),
            };
            let mut newly_ring_placed = Vec::new();
            dfs_zigzag(
                mol,
                nb,
                start,
                dir,
                placed,
                ring_systems,
                atom_to_system,
                system_placed,
                &mut newly_ring_placed,
                narrow,
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
fn best_outgoing_direction(atom: AtomIdx, placed: &[Option<Point>], mol: &Molecule) -> f64 {
    let Some(origin) = placed[atom.0 as usize] else {
        return 0.0;
    };

    // Collect angles to already-placed neighbors.
    let used_angles: Vec<f64> = mol
        .neighbors(atom)
        .filter(|(nb, _)| placed[nb.0 as usize].is_some())
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
            .filter_map(|p| p.as_ref())
            .any(|p| candidate.dist(p) < BOND_LEN / 2.0)
            || ring_centre(candidate, origin, placed)
    };

    ranked
        .iter()
        .copied()
        .find(|&dir| !occupies(dir))
        .unwrap_or(ranked[0])
}

/// Whether one bond from `atom` along `dir` is clear: no placed atom
/// within half a bond and not the middle of a placed ring.
fn direction_is_free(atom: AtomIdx, dir: f64, placed: &[Option<Point>]) -> bool {
    let Some(origin) = placed[atom.0 as usize] else {
        return false;
    };
    let candidate = Point::new(
        origin.x + BOND_LEN * dir.cos(),
        origin.y + BOND_LEN * dir.sin(),
    );
    !placed
        .iter()
        .filter_map(|p| p.as_ref())
        .any(|p| candidate.dist(p) < BOND_LEN / 2.0)
        && !ring_centre(candidate, origin, placed)
}

/// Whether `candidate` sits in the middle of a placed ring: three or more
/// placed atoms other than `origin` within 1.05 bond lengths (a hexagon's
/// centre is one bond from all six atoms, so the half-bond clash test alone
/// lets a quaternary ring atom's second substituent point into its ring).
fn ring_centre(candidate: Point, origin: Point, placed: &[Option<Point>]) -> bool {
    placed
        .iter()
        .filter_map(|p| p.as_ref())
        .filter(|p| p.dist(&origin) > 1e-6 && candidate.dist(p) < 1.05 * BOND_LEN)
        .count()
        >= 3
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

    // Separation computed once per candidate; the stable sort keeps the
    // order the per-comparison version gave.
    let mut keyed: Vec<(f64, f64)> = candidates
        .into_iter()
        .map(|c| (min_angle_separation(c, used_angles), c))
        .collect();
    keyed.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    keyed.into_iter().map(|(_, c)| c).collect()
}

/// Minimum angular separation between `angle` and any angle in `used`.
fn min_angle_separation(angle: f64, used: &[f64]) -> f64 {
    used.iter()
        .map(|&u| {
            // Reduced to [0, 2π) first: a candidate past 2π against a
            // negative used angle gave a negative separation (the exterior
            // gap of a ring-fusion atom ranked last, so a substituent went
            // into a ring).
            let diff = (angle - u).rem_euclid(2.0 * std::f64::consts::PI);
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
    ring_systems: &[Vec<Vec<AtomIdx>>],
    atom_to_system: &HashMap<AtomIdx, usize>,
    system_placed: &mut [bool],
    newly_ring_placed: &mut Vec<AtomIdx>,
    narrow: NarrowForks<'_>,
) {
    let deflection = std::f64::consts::PI / 6.0;
    // 4th element: the sign to apply, and then flip, when this atom continues
    // the chain alone. -1.0 matches the pre-fix first-step behavior.
    // 5th element: an extra turn its continuation takes once (a fork's
    // child bends back onto the chain direction, see below).
    let mut stack: Vec<(AtomIdx, AtomIdx, f64, f64, f64)> =
        vec![(start_atom, start_parent, start_dir, -1.0, 0.0)];

    while let Some((atom, parent, dir, sign, back)) = stack.pop() {
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
            // An entry atom with a second substituent (a 5,5-disubstituted
            // hydantoin, a 4,4-disubstituted piperidine) turns its ring 30°
            // so the two substituent bonds lie either side of the ring's
            // outward bisector instead of one on it; the side whose new
            // atoms land on fewer placed atoms is kept.
            let in_system = |a: AtomIdx| atom_to_system.get(&a).is_some_and(|&k| k == sys_idx);
            let (ring_nbs, outside) = mol.neighbors(atom).fold((0, 0), |(r, o), (nb, _)| {
                if in_system(nb) {
                    (r + 1, o)
                } else {
                    (r, o + 1)
                }
            });
            if ring_nbs == 2 && outside == 2 {
                let before: Vec<Option<Point>> = placed.to_vec();
                let mut best: Option<(usize, Vec<Option<Point>>)> = None;
                for turn in [std::f64::consts::PI / 6.0, -std::f64::consts::PI / 6.0] {
                    let mut trial = before.clone();
                    place_ring_system(
                        mol,
                        &ring_systems[sys_idx],
                        Some((atom, entry_pos, dir + turn)),
                        &mut trial,
                    );
                    let new_atoms: Vec<Point> = trial
                        .iter()
                        .zip(&before)
                        .filter(|(t, b)| t.is_some() && b.is_none())
                        .filter_map(|(t, _)| *t)
                        .collect();
                    let crowded = new_atoms
                        .iter()
                        .map(|q| {
                            before
                                .iter()
                                .flatten()
                                .filter(|p| p.dist(q) < BOND_LEN)
                                .count()
                        })
                        .sum::<usize>();
                    // The second substituent's spot (the parent bond mirrored
                    // across the ring's outward bisector) counts too.
                    let spot = {
                        let e = entry_pos;
                        let ring: Vec<Point> = mol
                            .neighbors(atom)
                            .filter(|&(nb, _)| in_system(nb))
                            .filter_map(|(nb, _)| trial[nb.0 as usize])
                            .collect();
                        (ring.len() == 2).then(|| {
                            let (mx, my) =
                                ((ring[0].x + ring[1].x) / 2.0, (ring[0].y + ring[1].y) / 2.0);
                            let b = (e.y - my).atan2(e.x - mx);
                            let back = dir + std::f64::consts::PI;
                            let free = 2.0 * b - back;
                            Point::new(e.x + BOND_LEN * free.cos(), e.y + BOND_LEN * free.sin())
                        })
                    };
                    let crowded = crowded
                        + spot.map_or(0, |q| {
                            trial
                                .iter()
                                .flatten()
                                .filter(|p| p.dist(&q) < BOND_LEN)
                                .count()
                        });
                    if best.as_ref().is_none_or(|(c, _)| crowded < *c) {
                        best = Some((crowded, trial));
                    }
                }
                if let Some((_, trial)) = best {
                    placed.copy_from_slice(&trial);
                }
            } else {
                place_ring_system(
                    mol,
                    &ring_systems[sys_idx],
                    Some((atom, entry_pos, dir)),
                    placed,
                );
            }
            system_placed[sys_idx] = true;
            newly_ring_placed.extend(ring_systems[sys_idx].iter().flatten().copied());
            continue; // Further growth from this ring's atoms is grow_layout's job.
        }

        // A chain atom whose zigzag position is already taken (two branches
        // meeting) moves to the nearest free direction around its parent,
        // keeping clear of the parent's other bonds; layouts without such a
        // collision are unchanged.
        let dir = free_direction(mol, parent, parent_pos, dir, placed);
        let x = parent_pos.x + BOND_LEN * dir.cos();
        let y = parent_pos.y + BOND_LEN * dir.sin();
        placed[atom.0 as usize] = Some(Point::new(x, y));

        let unplaced: Vec<AtomIdx> = mol
            .neighbors(atom)
            .map(|(nb, _)| nb)
            .filter(|&nb| nb != parent && placed[nb.0 as usize].is_none())
            .collect();

        if unplaced.len() == 1 {
            // Ordinary chain continuation: apply this atom's sign, flip it for
            // the next step -- this is the alternation the pre-fix code missed.
            stack.push((
                unplaced[0],
                atom,
                dir + back + sign * deflection,
                -sign,
                0.0,
            ));
        } else if unplaced.len() <= 2 {
            // A real branch (0 or 2 unplaced neighbors): the two children
            // 60° (30° where `narrow` marks the atom) either side of the
            // chain direction; at 60° the three
            // bonds are 120° apart (a ±30° split leaves two of them 60°
            // apart: an ester's C=O and C-C, the two aryl rings of a
            // tetrasubstituted alkene drawn on top of each other). Each
            // child's continuation turns back onto the ±30° chain direction,
            // so the chain beyond runs as it does with the ±30° split (a
            // peptide backbone does not curl).
            // Push in reverse so the first neighbor is popped first,
            // preserving DFS order.
            let fork = if narrow.angles[atom.0 as usize] != 0.0 {
                narrow.angles[atom.0 as usize]
            } else {
                2.0 * deflection
            };
            for (i, nb) in unplaced.into_iter().enumerate().rev() {
                let child_sign = if i % 2 == 0 { -1.0 } else { 1.0 };
                stack.push((
                    nb,
                    atom,
                    dir + child_sign * fork,
                    -child_sign,
                    -child_sign * (fork - deflection),
                ));
            }
        } else {
            // Three or more unplaced neighbors (a quaternary carbon, a
            // sulfonyl or phosphoryl centre): the alternating ±30° split put
            // the first and third on one point. Spread them evenly over the
            // directions away from the parent instead (a cross for three).
            let k = unplaced.len();
            let step = 2.0 * std::f64::consts::PI / (k as f64 + 1.0);
            // Three neighbours of which only one continues (a CF2 or CMe2 in
            // a chain): that one goes straight on, the terminal ones either
            // side, so a perfluoroalkyl chain does not curl back on itself.
            let mut unplaced = unplaced;
            if k == 3 {
                let going_on: Vec<usize> =
                    (0..3).filter(|&i| mol.degree(unplaced[i]) > 1).collect();
                if going_on.len() == 1 {
                    unplaced.swap(going_on[0], 1);
                }
            }
            for (i, nb) in unplaced.into_iter().enumerate().rev() {
                let offset = (i as f64 - (k as f64 - 1.0) / 2.0) * step;
                let child_sign = if offset > 0.0 { 1.0 } else { -1.0 };
                stack.push((nb, atom, dir + offset, -child_sign, 0.0));
            }
        }
    }
}

/// `dir` if the point one bond from `parent_pos` along it is free (no placed
/// atom of `component` within half a bond); otherwise the free direction in
/// 30° steps nearest to `dir` that keeps at least 30° from the parent's
/// placed bonds. `dir` itself when there is none.
fn free_direction(
    mol: &Molecule,
    parent: AtomIdx,
    parent_pos: Point,
    dir: f64,
    placed: &[Option<Point>],
) -> f64 {
    let free = |d: f64| {
        let candidate = Point::new(
            parent_pos.x + BOND_LEN * d.cos(),
            parent_pos.y + BOND_LEN * d.sin(),
        );
        !placed
            .iter()
            .filter_map(|p| p.as_ref())
            .any(|p| candidate.dist(p) < BOND_LEN / 2.0)
    };
    if free(dir) {
        return dir;
    }
    let used: Vec<f64> = mol
        .neighbors(parent)
        .filter_map(|(nb, _)| placed[nb.0 as usize])
        .map(|p| (p.y - parent_pos.y).atan2(p.x - parent_pos.x))
        .collect();
    let step = std::f64::consts::PI / 6.0;
    for k in 1..=6 {
        for sign in [1.0, -1.0] {
            let d = dir + sign * k as f64 * step;
            if min_angle_separation(d, &used) >= step - 1e-9 && free(d) {
                return d;
            }
        }
    }
    dir
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

/// Bond lookup for the clash tests (`contains(&(a, b))` with `a < b`).
struct Bonded<'a>(&'a Molecule);

impl Bonded<'_> {
    fn contains(&self, &(a, b): &(AtomIdx, AtomIdx)) -> bool {
        self.0.bond_between(a, b).is_some()
    }
}

/// Non-bonded atom pairs closer than this are a clash (the quality lane
/// counts pairs under 0.4 bond lengths).
const CLASH_DIST: f64 = 0.45 * BOND_LEN;

/// Atoms within two bonds of a clashing atom or a crossing bond's atom.
fn atoms_near_defects(mol: &Molecule, atoms: &[AtomIdx], placed: &[Option<Point>]) -> Vec<bool> {
    let at = |a: AtomIdx| placed[a.0 as usize].unwrap_or(Point::new(0.0, 0.0));
    let mut hit = vec![false; mol.atom_count()];
    // The same sweeps along x as `layout_defects`.
    let mut by_x: Vec<(f64, AtomIdx)> = atoms.iter().map(|&a| (at(a).x, a)).collect();
    by_x.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
    for (k, &(x, a)) in by_x.iter().enumerate() {
        for &(x2, b) in &by_x[k + 1..] {
            if x2 - x >= CLASH_DIST {
                break;
            }
            if at(a).dist(&at(b)) < CLASH_DIST && mol.bond_between(a, b).is_none() {
                hit[a.0 as usize] = true;
                hit[b.0 as usize] = true;
            }
        }
    }
    let mut in_component = vec![false; mol.atom_count()];
    for &a in atoms {
        in_component[a.0 as usize] = true;
    }
    let mut spans: Vec<(f64, f64, AtomIdx, AtomIdx)> = mol
        .bonds()
        .filter(|(_, e)| in_component[e.atom1.0 as usize])
        .map(|(_, e)| {
            let (xa, xb) = (at(e.atom1).x, at(e.atom2).x);
            (xa.min(xb), xa.max(xb), e.atom1, e.atom2)
        })
        .collect();
    spans.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
    for (k, &(_, hi, a, b)) in spans.iter().enumerate() {
        for &(lo2, _, c, d) in &spans[k + 1..] {
            if lo2 > hi + 1e-9 {
                break;
            }
            if a != c
                && a != d
                && b != c
                && b != d
                && boxes_meet(at(a), at(b), at(c), at(d))
                && segments_intersect(at(a), at(b), at(c), at(d))
            {
                for x in [a, b, c, d] {
                    hit[x.0 as usize] = true;
                }
            }
        }
    }
    let mut near = hit.clone();
    for _ in 0..2 {
        let frontier = near.clone();
        for (i, &f) in frontier.iter().enumerate() {
            if f {
                for (nb, _) in mol.neighbors(AtomIdx(i as u32)) {
                    near[nb.0 as usize] = true;
                }
            }
        }
    }
    near
}

/// The bonds of the component made of `atoms`.
fn component_bonds(mol: &Molecule, atoms: &[AtomIdx]) -> Vec<(AtomIdx, AtomIdx)> {
    let mut in_component = vec![false; mol.atom_count()];
    for &a in atoms {
        in_component[a.0 as usize] = true;
    }
    mol.bonds()
        .filter(|(_, e)| in_component[e.atom1.0 as usize])
        .map(|(_, e)| (e.atom1.min(e.atom2), e.atom1.max(e.atom2)))
        .collect()
}

/// Clash and crossing counts of one component's layout.
fn layout_defects(
    atoms: &[AtomIdx],
    bonds: &[(AtomIdx, AtomIdx)],
    bonded: &Bonded,
    placed: &[Option<Point>],
) -> (usize, usize) {
    let at = |a: AtomIdx| placed[a.0 as usize].unwrap_or(Point::new(0.0, 0.0));
    // Sweeps along x: a clash needs |dx| < CLASH_DIST and a crossing needs
    // overlapping x-extents, so only those pairs are tested, each with the
    // same predicate as an all-pairs count (the counts are identical).
    let mut by_x: Vec<(f64, AtomIdx)> = atoms.iter().map(|&a| (at(a).x, a)).collect();
    by_x.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
    let mut clashes = 0;
    for (k, &(x, a)) in by_x.iter().enumerate() {
        for &(x2, b) in &by_x[k + 1..] {
            if x2 - x >= CLASH_DIST {
                break;
            }
            if at(a).dist(&at(b)) < CLASH_DIST && !bonded.contains(&(a.min(b), a.max(b))) {
                clashes += 1;
            }
        }
    }
    // Margin on the x-extent so no pair the orientation test could call
    // crossing is skipped.
    const EPS: f64 = 1e-9;
    let mut spans: Vec<(f64, f64, AtomIdx, AtomIdx)> = bonds
        .iter()
        .map(|&(a, b)| {
            let (xa, xb) = (at(a).x, at(b).x);
            (xa.min(xb), xa.max(xb), a, b)
        })
        .collect();
    spans.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
    let mut crossings = 0;
    for (k, &(_, hi, a, b)) in spans.iter().enumerate() {
        for &(lo2, _, c, d) in &spans[k + 1..] {
            if lo2 > hi + EPS {
                break;
            }
            if a != c
                && a != d
                && b != c
                && b != d
                && boxes_meet(at(a), at(b), at(c), at(d))
                && segments_intersect(at(a), at(b), at(c), at(d))
            {
                crossings += 1;
            }
        }
    }
    (clashes, crossings)
}

/// Moves acyclic branches of one laid-out component to clear atom clashes
/// (two rings on one carbon, a phenyl folded back onto its neighbour): for
/// the first clashing pair, every acyclic bond on the shortest path between
/// them is a hinge, and the smaller-or-either side of it is mirrored across
/// the bond or turned 30 or 60 degrees about the hinge atom; the move that
/// lowers the clash count most without adding crossings is kept, repeated
/// while it helps. A mirror keeps every double bond's geometry and a turn is
/// only tried about an atom with no double bond, so E/Z drawings survive;
/// wedges are chosen from the coordinates afterwards.
/// With `keep_double_bonds`, a double bond is never a hinge: mirroring one
/// side across it would invert its E/Z (the layout itself sets no E/Z, so
/// the first pass may).
fn relieve_clashes(
    mol: &Molecule,
    atoms: &[AtomIdx],
    placed: &mut [Option<Point>],
    keep_double_bonds: bool,
) {
    if atoms.len() < 4 || atoms.len() > 300 {
        return;
    }
    // Moves that stack a moved atom on one that stays (within a fifth of a
    // bond) are refused: they can lower the clash count and still draw two
    // atoms as one. Some crowded chains (perfluorotributylamine) need such a
    // step that a later move clears; when the relief leaves a clash, it is
    // run again allowing them, and that drawing is kept when it has fewer
    // defects and no more stacked atoms.
    let before = placed.to_vec();
    let (clashes, crossings) = relieve_clashes_pass(mol, atoms, placed, keep_double_bonds, true);
    if clashes > 0 {
        let mut retry = before.clone();
        let (c2, x2) = relieve_clashes_pass(mol, atoms, &mut retry, keep_double_bonds, false);
        if c2 + x2 < clashes + crossings
            && c2 <= clashes
            && stacked_pairs(mol, atoms, &retry) <= stacked_pairs(mol, atoms, &before)
        {
            placed.copy_from_slice(&retry);
        }
    }
}

/// Non-bonded pairs of `atoms` drawn within a fifth of a bond.
fn stacked_pairs(mol: &Molecule, atoms: &[AtomIdx], placed: &[Option<Point>]) -> usize {
    let mut by_x: Vec<(Point, AtomIdx)> = atoms
        .iter()
        .filter_map(|&a| placed[a.0 as usize].map(|p| (p, a)))
        .collect();
    by_x.sort_unstable_by(|p, q| p.0.x.total_cmp(&q.0.x));
    let limit = 0.2 * BOND_LEN;
    let mut count = 0;
    for (k, &(p, a)) in by_x.iter().enumerate() {
        for &(q, b) in &by_x[k + 1..] {
            if q.x - p.x >= limit {
                break;
            }
            if p.dist(&q) < limit && mol.bond_between(a, b).is_none() {
                count += 1;
            }
        }
    }
    count
}

/// One clash relief run (see [`relieve_clashes`]); `no_stacking` refuses
/// moves that stack atoms. Returns the clashes and crossings left.
fn relieve_clashes_pass(
    mol: &Molecule,
    atoms: &[AtomIdx],
    placed: &mut [Option<Point>],
    keep_double_bonds: bool,
    no_stacking: bool,
) -> (usize, usize) {
    let mut atoms = atoms.to_vec();
    atoms.sort_unstable();
    let mut in_component = vec![false; mol.atom_count()];
    for &a in &atoms {
        in_component[a.0 as usize] = true;
    }
    let mut bonds: Vec<(AtomIdx, AtomIdx)> = mol
        .bonds()
        .filter(|(_, e)| in_component[e.atom1.0 as usize])
        .map(|(_, e)| (e.atom1.min(e.atom2), e.atom1.max(e.atom2)))
        .collect();
    bonds.sort_unstable();
    let bonded = Bonded(mol);
    let mut current = layout_defects(&atoms, &bonds, &bonded, placed);
    // Atoms reachable from `root` without crossing the bond to `pivot`, or
    // None when that bond is in a ring.
    let n_total = mol.atom_count();
    let mut seen_mark = vec![0u32; n_total];
    let mut seen_epoch = 0u32;
    // Atoms reachable from `root` without crossing the bond to `pivot`, or
    // None when that bond is in a ring.
    // The graph does not change between rounds, so each side is found once.
    let mut side_cache: HashMap<(AtomIdx, AtomIdx), Option<std::rc::Rc<Vec<AtomIdx>>>> =
        HashMap::default();
    let mut side_uncached = |pivot: AtomIdx, root: AtomIdx| -> Option<Vec<AtomIdx>> {
        seen_epoch += 1;
        let e = seen_epoch;
        seen_mark[root.0 as usize] = e;
        let mut stack = vec![root];
        let mut out = vec![root];
        while let Some(a) = stack.pop() {
            for (nb, _) in mol.neighbors(a) {
                if a == root && nb == pivot {
                    continue;
                }
                if nb == pivot {
                    return None;
                }
                if seen_mark[nb.0 as usize] != e {
                    seen_mark[nb.0 as usize] = e;
                    out.push(nb);
                    stack.push(nb);
                }
            }
        }
        Some(out)
    };
    let mut side = |pivot: AtomIdx, root: AtomIdx| -> Option<std::rc::Rc<Vec<AtomIdx>>> {
        side_cache
            .entry((pivot, root))
            .or_insert_with(|| side_uncached(pivot, root).map(std::rc::Rc::new))
            .clone()
    };
    // Neighbours in index order, as the breadth-first search visits them.
    let sorted_neighbours: Vec<Vec<AtomIdx>> = (0..n_total)
        .map(|i| {
            if !in_component[i] {
                return Vec::new();
            }
            let mut nbs: Vec<AtomIdx> =
                mol.neighbors(AtomIdx(i as u32)).map(|(nb, _)| nb).collect();
            nbs.sort_unstable();
            nbs
        })
        .collect();
    let mut prev = vec![u32::MAX; n_total];
    let mut shortest_path = |from: AtomIdx, to: AtomIdx| -> Vec<AtomIdx> {
        prev.iter_mut().for_each(|p| *p = u32::MAX);
        let mut queue = VecDeque::from([from]);
        prev[from.0 as usize] = from.0;
        while let Some(a) = queue.pop_front() {
            if a == to {
                break;
            }
            for &nb in &sorted_neighbours[a.0 as usize] {
                if prev[nb.0 as usize] == u32::MAX {
                    prev[nb.0 as usize] = a.0;
                    queue.push_back(nb);
                }
            }
        }
        let mut path = vec![to];
        let mut a = to;
        while a != from {
            let p = prev[a.0 as usize];
            if p == u32::MAX {
                return Vec::new();
            }
            path.push(AtomIdx(p));
            a = AtomIdx(p);
        }
        path
    };
    let multiple = |b: chematic_core::BondIdx| {
        !matches!(
            mol.bond(b).order,
            chematic_core::BondOrder::Single
                | chematic_core::BondOrder::Up
                | chematic_core::BondOrder::Down
                | chematic_core::BondOrder::Aromatic
        )
    };
    // Multiple bonds in a ring (a Kekulé-written aromatic ring, a ring
    // P=N): turning a substituent of their atoms sets no E/Z.
    let in_ring: HashSet<chematic_core::BondIdx> = mol
        .bonds()
        .filter(|&(b, e)| in_component[e.atom1.0 as usize] && multiple(b))
        .filter(|&(b, e)| {
            let mut seen = vec![false; n_total];
            let mut stack = vec![e.atom1];
            seen[e.atom1.0 as usize] = true;
            while let Some(a) = stack.pop() {
                for (nb, via) in mol.neighbors(a) {
                    if via != b && !seen[nb.0 as usize] {
                        seen[nb.0 as usize] = true;
                        stack.push(nb);
                    }
                }
            }
            seen[e.atom2.0 as usize]
        })
        .map(|(b, _)| b)
        .collect();
    // Turns pivot only on atoms without a double bond. In the layout's own
    // pass (before any E/Z is set) a ring's double bonds do not count.
    let has_double = |a: AtomIdx| {
        mol.neighbors(a)
            .any(|(_, b)| multiple(b) && (keep_double_bonds || !in_ring.contains(&b)))
    };
    for _round in 0..12 {
        if current == (0, 0) {
            break;
        }
        let snapshot: Vec<Point> = placed
            .iter()
            .map(|p| p.unwrap_or(Point::new(0.0, 0.0)))
            .collect();
        let mut in_moved = vec![false; snapshot.len()];
        let mut new_pos = snapshot.clone();
        let at = |a: AtomIdx| snapshot[a.0 as usize];
        // The first clashing atom pairs and crossing bond pairs, in atom and
        // bond order (found by sweeps along x, then sorted).
        let mut clash_pairs = Vec::new();
        let mut by_x: Vec<(f64, AtomIdx)> = atoms.iter().map(|&a| (at(a).x, a)).collect();
        by_x.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
        for (k, &(x, a)) in by_x.iter().enumerate() {
            for &(x2, b) in &by_x[k + 1..] {
                if x2 - x >= CLASH_DIST {
                    break;
                }
                let pair = (a.min(b), a.max(b));
                if at(a).dist(&at(b)) < CLASH_DIST && !bonded.contains(&pair) {
                    clash_pairs.push(pair);
                }
            }
        }
        clash_pairs.sort_unstable();
        let mut targets: Vec<(AtomIdx, AtomIdx)> = clash_pairs.iter().copied().take(4).collect();
        let all_clash_pairs = clash_pairs;
        // Crossing bond pairs: the path joins their nearer ends.
        let mut spans: Vec<(f64, f64, usize)> = bonds
            .iter()
            .enumerate()
            .map(|(k, &(a, b))| {
                let (xa, xb) = (at(a).x, at(b).x);
                (xa.min(xb), xa.max(xb), k)
            })
            .collect();
        spans.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
        let mut crossing_pairs: Vec<(usize, usize)> = Vec::new();
        for (i, &(_, hi, k)) in spans.iter().enumerate() {
            for &(lo2, _, l) in &spans[i + 1..] {
                if lo2 > hi + 1e-9 {
                    break;
                }
                let ((a, b), (c, d)) = (bonds[k], bonds[l]);
                if a != c
                    && a != d
                    && b != c
                    && b != d
                    && boxes_meet(at(a), at(b), at(c), at(d))
                    && segments_intersect(at(a), at(b), at(c), at(d))
                {
                    crossing_pairs.push((k.min(l), k.max(l)));
                }
            }
        }
        crossing_pairs.sort_unstable();
        for &(k, l) in &crossing_pairs {
            if targets.len() >= 6 {
                break;
            }
            let ((a, b), (c, d)) = (bonds[k], bonds[l]);
            targets.push((a, c));
            targets.push((b, d));
        }
        let index = RoundIndex::new(&atoms, &bonds, all_clash_pairs, crossing_pairs, at);
        let mut stamp = vec![0u32; bonds.len()];
        let mut epoch = 0u32;
        let mut best: Option<(DefectKey, Vec<(AtomIdx, Point)>)> = None;
        // A hinge reached again through another target gives the same
        // candidates, which can never beat the first evaluation.
        let mut tried: HashSet<(AtomIdx, AtomIdx)> = HashSet::default();
        let mut new: Vec<(AtomIdx, Point)> = Vec::new();
        let mut touched: Vec<usize> = Vec::new();
        for &(a, b) in &targets {
            let path = shortest_path(a, b);
            for w in path.windows(2) {
                if keep_double_bonds
                    && mol
                        .bond_between(w[0], w[1])
                        .is_some_and(|(_, e)| e.order == chematic_core::BondOrder::Double)
                {
                    continue;
                }
                for (pivot, root) in [(w[0], w[1]), (w[1], w[0])] {
                    if !tried.insert((pivot, root)) {
                        continue;
                    }
                    let Some(moved) = side(pivot, root) else {
                        continue;
                    };
                    if moved.len() * 2 > atoms.len() + 2 {
                        continue; // move the smaller side
                    }
                    let (pp, pr) = (at(pivot), at(root));
                    // The mirror across the hinge bond, then turns about the
                    // pivot (sin, cos), in this order.
                    let mut turns = [(0.0f64, 1.0f64); 10];
                    let mut n_turns = 0;
                    if !has_double(pivot) {
                        // A substituent of up to three atoms may also swing
                        // further round its hinge (a bridge atom's methyls,
                        // a tropane's N-methyl, into a free face; a carboxyl).
                        let wide: &[f64] = if moved.len() <= 3 {
                            &[
                                30.0, -30.0, 60.0, -60.0, 90.0, -90.0, 120.0, -120.0, 150.0, -150.0,
                            ]
                        } else {
                            &[30.0, -30.0, 60.0, -60.0]
                        };
                        for &deg in wide {
                            turns[n_turns] = deg.to_radians().sin_cos();
                            n_turns += 1;
                        }
                    }
                    let apply = |k: usize, q: Point| -> Point {
                        if k == 0 {
                            reflect_point(q, pp, pr)
                        } else {
                            let (sin, cos) = turns[k - 1];
                            let (dx, dy) = (q.x - pp.x, q.y - pp.y);
                            Point::new(pp.x + dx * cos - dy * sin, pp.y + dx * sin + dy * cos)
                        }
                    };
                    // The hinge's other bonds, for the turns' 30-degree check.
                    let others: Vec<f64> = if n_turns > 0 {
                        mol.neighbors(pivot)
                            .filter(|&(nb, _)| nb != root)
                            .map(|(nb, _)| {
                                let q = at(nb);
                                (q.y - pp.y).atan2(q.x - pp.x)
                            })
                            .collect()
                    } else {
                        Vec::new()
                    };
                    for k in 0..=n_turns {
                        if k > 0 {
                            // A turn must leave the hinge's bonds 30 degrees
                            // apart (`moved` starts with the root).
                            let np = apply(k, at(moved[0]));
                            let dir = (np.y - pp.y).atan2(np.x - pp.x);
                            // A three-atom group (a carboxyl) turned past 60
                            // degrees keeps them 60 degrees apart, so that it
                            // does not open a narrow fork.
                            let apart = if k > 4 && moved.len() == 3 {
                                60f64
                            } else {
                                30f64
                            };
                            if min_angle_separation(dir, &others) < apart.to_radians() - 1e-9 {
                                continue;
                            }
                        }
                        new.clear();
                        new.extend(moved.iter().map(|&m| (m, apply(k, at(m)))));
                        // Only a move that beats the current drawing and the
                        // best move so far can be kept, so the count stops
                        // once it cannot.
                        // (A move equal to the best one can still win on
                        // size; one equal to the current drawing cannot.)
                        let bound = match &best {
                            Some((k0, _)) if (k0.0, k0.1) < current => (k0.0, k0.1, true),
                            _ => (current.0, current.1, false),
                        };
                        let (clashes, crossings) = defects_after_move(
                            current,
                            bound,
                            &index,
                            &bonds,
                            &bonded,
                            &snapshot,
                            &new,
                            &mut in_moved,
                            &mut new_pos,
                            &mut stamp,
                            &mut epoch,
                            &mut touched,
                        );
                        if (clashes, crossings) < current {
                            let key = (clashes, crossings, moved.len());
                            if best.as_ref().is_none_or(|(k0, _)| key < *k0) {
                                // A move that stacks a moved atom on one that
                                // stays (within a fifth of a bond) can lower
                                // the clash count and still draw two atoms as
                                // one.
                                if no_stacking && stacks(&index, &snapshot, &new, &mut in_moved) {
                                    continue;
                                }
                                best = Some((key, new.clone()));
                            }
                        }
                    }
                }
            }
        }
        let Some(((clashes, crossings, _), new)) = best else {
            break;
        };
        for (m, q) in new {
            placed[m.0 as usize] = Some(q);
        }
        current = (clashes, crossings);
    }
    current
}

/// Whether a moved atom of `new` lands within a fifth of a bond of an atom
/// that stays.
fn stacks(
    index: &RoundIndex,
    snapshot: &[Point],
    new: &[(AtomIdx, Point)],
    in_moved: &mut [bool],
) -> bool {
    for &(m, _) in new {
        in_moved[m.0 as usize] = true;
    }
    let stacked = new.iter().any(|&(_, q)| {
        let (cx, cy) = cell(q, CLASH_DIST);
        (-1..=1).any(|dx| {
            (-1..=1).any(|dy| {
                index.atom_grid.get((cx + dx, cy + dy)).is_some_and(|list| {
                    list.iter().any(|&b| {
                        !in_moved[b.0 as usize] && q.dist(&snapshot[b.0 as usize]) < 0.2 * BOND_LEN
                    })
                })
            })
        })
    });
    for &(m, _) in new {
        in_moved[m.0 as usize] = false;
    }
    stacked
}

/// One round's lookup tables for [`defects_after_move`]: every clash and
/// crossing of the current layout, and grids of atoms (cells of
/// `CLASH_DIST`) and bond boxes (cells of a bond length).
struct RoundIndex {
    clash_pairs: Vec<(AtomIdx, AtomIdx)>,
    crossing_pairs: Vec<(usize, usize)>,
    atom_grid: DenseGrid<AtomIdx>,
    bond_grid: DenseGrid<usize>,
    /// Bond indices per atom index (empty for atoms outside the component).
    atom_bonds: Vec<Vec<usize>>,
}

/// Grid cells over a bounded range, each a list in insertion order: the
/// same lists a map from cell to list holds, without hashing (one flat
/// array, cells by offset).
struct DenseGrid<T> {
    x0: i64,
    y0: i64,
    width: i64,
    height: i64,
    /// Cell `c`'s entries are `items[starts[c]..starts[c + 1]]`.
    starts: Vec<u32>,
    items: Vec<T>,
}

impl<T: Copy> DenseGrid<T> {
    fn new(entries: &[((i64, i64), T)]) -> Self {
        let (mut x0, mut y0, mut x1, mut y1) = (i64::MAX, i64::MAX, i64::MIN, i64::MIN);
        for &((x, y), _) in entries {
            x0 = x0.min(x);
            y0 = y0.min(y);
            x1 = x1.max(x);
            y1 = y1.max(y);
        }
        if entries.is_empty() {
            return Self {
                x0: 0,
                y0: 0,
                width: 0,
                height: 0,
                starts: vec![0],
                items: Vec::new(),
            };
        }
        let (width, height) = (x1 - x0 + 1, y1 - y0 + 1);
        let index = |x: i64, y: i64| ((x - x0) * height + (y - y0)) as usize;
        let mut starts = vec![0u32; (width * height) as usize + 1];
        for &((x, y), _) in entries {
            starts[index(x, y) + 1] += 1;
        }
        for c in 1..starts.len() {
            starts[c] += starts[c - 1];
        }
        let mut fill: Vec<u32> = starts[..starts.len() - 1].to_vec();
        let mut items: Vec<T> = vec![entries[0].1; entries.len()];
        for &((x, y), v) in entries {
            let c = index(x, y);
            items[fill[c] as usize] = v;
            fill[c] += 1;
        }
        Self {
            x0,
            y0,
            width,
            height,
            starts,
            items,
        }
    }

    fn get(&self, (x, y): (i64, i64)) -> Option<&[T]> {
        let (i, j) = (x - self.x0, y - self.y0);
        if i < 0 || j < 0 || i >= self.width || j >= self.height {
            return None;
        }
        let c = (i * self.height + j) as usize;
        let list = &self.items[self.starts[c] as usize..self.starts[c + 1] as usize];
        (!list.is_empty()).then_some(list)
    }
}

const BOND_CELL: f64 = BOND_LEN;
const BOX_EPS: f64 = 1e-9;

fn cell(p: Point, size: f64) -> (i64, i64) {
    // `floor` without the library call: truncate, then step down for a
    // negative non-integer (exact for the coordinate range a layout uses).
    let floor = |v: f64| {
        let t = v as i64;
        if (t as f64) > v { t - 1 } else { t }
    };
    (floor(p.x / size), floor(p.y / size))
}

/// The grid cells a box (with a margin) overlaps.
fn box_cells(a: Point, b: Point) -> impl Iterator<Item = (i64, i64)> {
    let lo = cell(
        Point::new(a.x.min(b.x) - BOX_EPS, a.y.min(b.y) - BOX_EPS),
        BOND_CELL,
    );
    let hi = cell(
        Point::new(a.x.max(b.x) + BOX_EPS, a.y.max(b.y) + BOX_EPS),
        BOND_CELL,
    );
    (lo.0..=hi.0).flat_map(move |x| (lo.1..=hi.1).map(move |y| (x, y)))
}

fn boxes_meet(pa: Point, pb: Point, pc: Point, pd: Point) -> bool {
    pa.x.max(pb.x) + BOX_EPS >= pc.x.min(pd.x)
        && pc.x.max(pd.x) + BOX_EPS >= pa.x.min(pb.x)
        && pa.y.max(pb.y) + BOX_EPS >= pc.y.min(pd.y)
        && pc.y.max(pd.y) + BOX_EPS >= pa.y.min(pb.y)
}

impl RoundIndex {
    fn new(
        atoms: &[AtomIdx],
        bonds: &[(AtomIdx, AtomIdx)],
        clash_pairs: Vec<(AtomIdx, AtomIdx)>,
        crossing_pairs: Vec<(usize, usize)>,
        at: impl Fn(AtomIdx) -> Point,
    ) -> Self {
        let atom_entries: Vec<((i64, i64), AtomIdx)> = atoms
            .iter()
            .map(|&a| (cell(at(a), CLASH_DIST), a))
            .collect();
        let mut bond_entries: Vec<((i64, i64), usize)> = Vec::new();
        let n = atoms.iter().map(|a| a.0 as usize + 1).max().unwrap_or(0);
        let mut atom_bonds: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (k, &(a, b)) in bonds.iter().enumerate() {
            for c in box_cells(at(a), at(b)) {
                bond_entries.push((c, k));
            }
            atom_bonds[a.0 as usize].push(k);
            atom_bonds[b.0 as usize].push(k);
        }
        Self {
            clash_pairs,
            crossing_pairs,
            atom_grid: DenseGrid::new(&atom_entries),
            bond_grid: DenseGrid::new(&bond_entries),
            atom_bonds,
        }
    }
}

/// The clash and crossing counts after a rigid move of the atoms in `new`
/// (a mirror or a turn of one branch), from the counts `current` before it:
/// only atom pairs between the branch and the rest, and bond pairs with a
/// bond on or into the branch against one not inside it, change. Pairs are
/// found through `index`'s grids and tested with the same predicates as a
/// full count.
#[allow(clippy::too_many_arguments)]
fn defects_after_move(
    current: (usize, usize),
    bound: (usize, usize, bool),
    index: &RoundIndex,
    bonds: &[(AtomIdx, AtomIdx)],
    bonded: &Bonded,
    old: &[Point],
    new: &[(AtomIdx, Point)],
    in_moved: &mut [bool],
    new_pos: &mut [Point],
    stamp: &mut [u32],
    epoch: &mut u32,
    touched: &mut Vec<usize>,
) -> (usize, usize) {
    for &(m, q) in new {
        in_moved[m.0 as usize] = true;
        new_pos[m.0 as usize] = q;
    }
    let moved = |a: AtomIdx| in_moved[a.0 as usize];
    let mut clashes = current.0 as i64;
    clashes -= index
        .clash_pairs
        .iter()
        .filter(|&&(a, b)| moved(a) != moved(b))
        .count() as i64;
    // The clash count only grows from here: a move past `bound` stops early.
    let mut over = false;
    for &(m, q) in new {
        if clashes > bound.0 as i64 {
            over = true;
            break;
        }
        let (cx, cy) = cell(q, CLASH_DIST);
        for dx in -1..=1 {
            for dy in -1..=1 {
                let Some(list) = index.atom_grid.get((cx + dx, cy + dy)) else {
                    continue;
                };
                for &b in list {
                    if !moved(b)
                        && q.dist(&old[b.0 as usize]) < CLASH_DIST
                        && !bonded.contains(&(m.min(b), m.max(b)))
                    {
                        clashes += 1;
                    }
                }
            }
        }
    }
    // A move with more clashes than `bound` cannot be kept: its crossings
    // are not counted.
    let restore = |in_moved: &mut [bool], new_pos: &mut [Point]| {
        for &(m, _) in new {
            in_moved[m.0 as usize] = false;
            new_pos[m.0 as usize] = old[m.0 as usize];
        }
    };
    if over || clashes > bound.0 as i64 {
        restore(in_moved, new_pos);
        return (usize::MAX, usize::MAX);
    }
    // With as many clashes as `bound`, the count stops once the crossings
    // (which only grow from here) pass it (reach it, when a tie loses).
    let at_bound = clashes == bound.0 as i64;
    let lost = |crossings: i64| {
        at_bound && (crossings > bound.1 as i64 || (!bound.2 && crossings >= bound.1 as i64))
    };
    // Bonds on the branch: inside it, or the hinge into it.
    touched.clear();
    *epoch += 1;
    for &(m, _) in new {
        for &k in index
            .atom_bonds
            .get(m.0 as usize)
            .map_or(&[][..], Vec::as_slice)
        {
            if stamp[k] != *epoch {
                stamp[k] = *epoch;
                touched.push(k);
            }
        }
    }
    let touches = |k: usize| moved(bonds[k].0) || moved(bonds[k].1);
    let inside = |k: usize| moved(bonds[k].0) && moved(bonds[k].1);
    let changes = |k: usize, l: usize| (touches(k) || touches(l)) && !(inside(k) && inside(l));
    let mut crossings = current.1 as i64;
    crossings -= index
        .crossing_pairs
        .iter()
        .filter(|&&(k, l)| changes(k, l))
        .count() as i64;
    let share = |k: usize, l: usize| {
        let ((a, b), (c, d)) = (bonds[k], bonds[l]);
        a == c || a == d || b == c || b == d
    };
    let np = |x: AtomIdx| new_pos[x.0 as usize];
    let crosses = |k: usize, l: usize| {
        let ((a, b), (c, d)) = (bonds[k], bonds[l]);
        boxes_meet(np(a), np(b), np(c), np(d)) && segments_intersect(np(a), np(b), np(c), np(d))
    };
    for &k in touched.iter() {
        let (a, b) = bonds[k];
        // Against the bonds that stay where they are.
        *epoch += 1;
        let e = *epoch;
        for c in box_cells(np(a), np(b)) {
            let Some(list) = index.bond_grid.get(c) else {
                continue;
            };
            for &l in list {
                if stamp[l] == e || touches(l) {
                    continue;
                }
                stamp[l] = e;
                if !share(k, l) && crosses(k, l) {
                    crossings += 1;
                }
            }
        }
        if lost(crossings) {
            break;
        }
        // Against the other branch bonds, once per pair (two bonds inside
        // the branch move together).
        if !inside(k) {
            for &l in touched.iter() {
                if l != k && (inside(l) || l > k) && !share(k, l) && crosses(k, l) {
                    crossings += 1;
                }
            }
        }
        if lost(crossings) {
            break;
        }
    }
    restore(in_moved, new_pos);
    if lost(crossings) {
        return (usize::MAX, usize::MAX);
    }
    (clashes.max(0) as usize, crossings.max(0) as usize)
}

/// (clashes, crossings, atoms moved): the order in which moves are preferred.
type DefectKey = (usize, usize, usize);

/// `q` mirrored across the line through `a` and `b`.
fn reflect_point(q: Point, a: Point, b: Point) -> Point {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-12 {
        return q;
    }
    let t = ((q.x - a.x) * dx + (q.y - a.y) * dy) / len2;
    let (fx, fy) = (a.x + t * dx, a.y + t * dy);
    Point::new(2.0 * fx - q.x, 2.0 * fy - q.y)
}

/// Whether segments AB and CD cross at a point inside both (touching ends
/// and collinear segments do not count). An orientation within a relative
/// 1e-9 of zero counts as collinear: two collinear bonds of a straight chain
/// (`(…)c1cc…` drawn along one line) used to read as crossing by rounding.
fn segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool {
    let side = |p: Point, q: Point, r: Point| -> i8 {
        let (ux, uy, vx, vy) = (q.x - p.x, q.y - p.y, r.x - p.x, r.y - p.y);
        let cross = ux * vy - uy * vx;
        // |u||v| <= (|ux|+|uy|)(|vx|+|vy|): above that bound the exact test
        // below cannot call it collinear, so its square roots are skipped.
        if cross.abs() > 1e-9 * ((ux.abs() + uy.abs()) * (vx.abs() + vy.abs())) {
            return if cross > 0.0 { 1 } else { -1 };
        }
        let scale = (ux * ux + uy * uy).sqrt() * (vx * vx + vy * vy).sqrt();
        if cross.abs() <= 1e-9 * scale {
            0
        } else if cross > 0.0 {
            1
        } else {
            -1
        }
    };
    let (o1, o2) = (side(a, b, c), side(a, b, d));
    let (o3, o4) = (side(c, d, a), side(c, d, b));
    o1 * o2 < 0 && o3 * o4 < 0
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
        // distance 0.0). The bridged core shares a three-atom path between
        // rings, so placement must use all already-known ring atoms to choose
        // the correct side of the anchored regular polygon.
        use chematic_smiles::parse;
        let mol = parse("C1CC2CN(CC1N2c1ccccc1)c1cccc(c1)-c1ccccc1").unwrap();
        let layout = compute_layout(&mol);
        let (min_non_bonded, min_bonded, max_bonded) = layout_distance_summary(&mol, &layout);
        assert!(
            min_non_bonded > 1e-6,
            "no two atoms of unrelated ring systems should land on identical coordinates: \
             min_non_bonded={min_non_bonded}"
        );
        assert!(
            min_bonded > 30.0 && max_bonded < 50.0,
            "bridged-core bonds should stay near BOND_LEN: min={min_bonded}, max={max_bonded}; bonds={:?}",
            mol.bonds()
                .map(|(_, b)| (
                    b.atom1.0,
                    b.atom2.0,
                    layout.get(b.atom1).dist(&layout.get(b.atom2))
                ))
                .collect::<Vec<_>>()
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

    /// Batch 14: bridged cages, quaternary ring atoms and crowded branches
    /// are drawn with no atom within 0.4 bond lengths of a non-bonded atom,
    /// and the uncaged ones with no crossing bond (each clashed or crossed
    /// before).
    #[test]
    fn chain_forks_are_drawn_120_degrees_apart() {
        // A branch point's three bonds were drawn 150°/150°/60° (the ±30°
        // zigzag split): an ester's C=O against its C-C.
        use chematic_smiles::parse;
        for (smi, atom) in [("CCOC(=O)C", 3u32), ("CC(C)CC", 1), ("C=C(C)c1ccccc1", 1)] {
            let mol = parse(smi).unwrap();
            let layout = compute_layout(&mol);
            let c = layout.get(AtomIdx(atom));
            let mut angles: Vec<f64> = mol
                .neighbors(AtomIdx(atom))
                .map(|(nb, _)| {
                    let p = layout.get(nb);
                    (p.y - c.y)
                        .atan2(p.x - c.x)
                        .rem_euclid(std::f64::consts::TAU)
                })
                .collect();
            angles.sort_by(f64::total_cmp);
            for k in 0..3 {
                let gap = (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU);
                assert!((gap.to_degrees() - 120.0).abs() < 1e-6, "{smi}: {angles:?}");
            }
        }
    }

    #[test]
    fn angle_separation_wraps_past_two_pi() {
        // 330° against -90° is 60° apart, not negative.
        let sep = min_angle_separation(330f64.to_radians(), &[(-90f64).to_radians()]);
        assert!(
            (sep.to_degrees() - 60.0).abs() < 1e-9,
            "{}",
            sep.to_degrees()
        );
    }

    #[test]
    fn large_components_open_their_forks_too() {
        // A component of more than 60 atoms that needs narrow forks: its
        // ±45° drawing is relieved as well, and keeps no fork under 90°.
        use chematic_smiles::parse;
        let smi = "CC(C)Cc1cccc(-c2ccccc2C(C)C)c1OC1CC(CNC(=O)c2ccc(C=C3SC(=O)NC3=O)cc2)N(C(=O)c2ccccc2C(=O)c2ccc(F)cc2F)C1";
        let mol = parse(smi).unwrap();
        assert!(mol.atom_count() > 60);
        let layout = compute_layout(&mol);
        let ring_atoms: std::collections::HashSet<AtomIdx> = chematic_perception::find_sssr(&mol)
            .rings()
            .iter()
            .flatten()
            .copied()
            .collect();
        let mut narrow = 0;
        for i in 0..mol.atom_count() {
            let atom = AtomIdx(i as u32);
            if ring_atoms.contains(&atom) || mol.degree(atom) != 3 {
                continue;
            }
            let c = layout.get(atom);
            let mut angles: Vec<f64> = mol
                .neighbors(atom)
                .map(|(nb, _)| {
                    let p = layout.get(nb);
                    (p.y - c.y)
                        .atan2(p.x - c.x)
                        .rem_euclid(std::f64::consts::TAU)
                })
                .collect();
            angles.sort_by(f64::total_cmp);
            if (0..3).any(|k| {
                (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU)
                    < 90f64.to_radians() - 1e-6
            }) {
                narrow += 1;
            }
        }
        assert_eq!(narrow, 0);
        assert!(detect_crossings(&layout, &mol).is_empty());
    }

    #[test]
    fn crowded_forks_open_to_90_degrees_when_that_does_as_well() {
        // The carbonyl of a trityl phenyl ketone needs a narrower fork than
        // 120°; at ±30° its O ran 60° from a C-C, at ±45° 90°.
        use chematic_smiles::parse;
        let smi = "O=C(C1=CC=CC=C1)C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4";
        let mol = parse(smi).unwrap();
        let layout = compute_layout(&mol);
        let c = layout.get(AtomIdx(1));
        let mut angles: Vec<f64> = mol
            .neighbors(AtomIdx(1))
            .map(|(nb, _)| {
                let p = layout.get(nb);
                (p.y - c.y)
                    .atan2(p.x - c.x)
                    .rem_euclid(std::f64::consts::TAU)
            })
            .collect();
        angles.sort_by(f64::total_cmp);
        for k in 0..3 {
            let gap = (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU);
            assert!(gap.to_degrees() > 90.0 - 1e-6, "{angles:?}");
        }
        assert!(detect_crossings(&layout, &mol).is_empty());
    }

    #[test]
    fn bridged_and_crowded_layouts_have_no_clash_or_crossing() {
        use chematic_smiles::parse;
        for (smi, cage) in [
            ("CNC1CC2CCC1C2", true),        // norbornane
            ("C1C[S+]2CC[S+]1CC2", true),   // bicyclo[2.2.2]
            ("C1OCC2(CN3CCC2CC3)O1", true), // spiro quinuclidine
            // Tropane, a pinanol: the N-methyl and the bridge's methyls
            // swing into free faces instead of across a ring bond.
            ("CN1C2CCC1CC(O)C2", false),
            ("CN1C2CCC1CC(OC(=O)c1c[nH]c3ccccc13)C2", false),
            ("CC1(C)C2CCC(C)(O)C1C2", false),
            ("c1ccc(C(n2ccnc2)n2ccnc2)cc1", false), // two rings on one carbon
            ("CCC1(CC)C(=O)NC(=O)NC1=O", false),    // quaternary ring carbon
            ("OC(=O)C(F)(F)C(F)(F)C(O)=O", false),  // crowded chain
            ("c1ccc(C2(N3CCCCC3)CCCCC2)cc1", false),
            // A carboxyl stacked on a Kekulé-written ring: turned about the
            // ring atom, whose ring double bond sets no E/Z.
            (
                "OC(=O)C1=C(C=CC=C1)C(C2=CC=C(O)C(=C2)S(O)(=O)=O)=C3C=CC(=O)C(=C3)S(O)(=O)=O",
                false,
            ),
            // Aziridines on a cyclophosphazene's P (ring P=N).
            (
                "C1CCN(P2(N3CCCC3)=NP(N3CC3)(N3CC3)=NP(N3CC3)(N3CC3)=N2)C1",
                false,
            ),
            // Fused and spiro systems of three or more rings that overlap
            // as polygons (cephalotaxine core, spiro naphthodioxin).
            ("COC1=CC23CCCN2CCc2cc4c(cc2C3C1O)OCO4", false),
            ("O=C1C(O)=CC2(Oc3cccc4cccc(c34)O2)C23OC12C(O)CCC3O", false),
            // Gem-disubstituted ring entry atoms: phenytoin, a
            // 3,3-disubstituted glutarimide.
            ("O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1", false),
            ("O=C1CCC(c2ccccc2)(C2CCN(C)CC2)C(=O)N1", false),
            // Adamantane-type cages, drawn as RDKit's projection (one
            // crossing, as RDKit's), substituted on a bridgehead or a link.
            ("C1C2CC3CC1CC(C2)C3", true),
            ("C1N2CN3CN1CN(C2)C3", true),
            ("C=CC[N+]12CN3CN(CN(C3)C1)C2", true),
            ("NC12CC3CC(CC(C3)C1)C2", true),
            ("OC1C2CC3CC(C2)CC1C3", true),
            // Porphyrins (protoporphyrin IX, a zinc porphyrin): the pyrroles
            // were folded inside the macrocycle.
            (
                "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C",
                false,
            ),
            (
                "C1=CC2=CC3=CC=C4N3[Zn]35N2C1=CC1=CC=C(N13)C=C1C=CC(=N15)C=4",
                false,
            ),
            // Perfluoroalkyl chains run straight through their CF2 atoms
            // instead of curling back (perfluorotributylamine; a
            // perfluorooctyl catechol).
            (
                "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)N(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
                false,
            ),
            (
                "Oc1cccc(CCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)c1O",
                false,
            ),
            // A ring-fusion atom's substituent goes into the exterior gap
            // (it went into a ring when the gap's direction passed 2π).
            ("CC1(C)OC2C3C(COC2(COS(N)(=O)=O)O1)C3(Cl)Cl", false),
            // An aspartate carboxyl turned past 60° round its CH2, clear of
            // the neighbouring residue (three-atom groups take the wide turns).
            (
                "CCCN(NC(=O)C1CCCN1C(=O)C(NC(=O)C(NC(=O)C(CC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(NC(=O)C(CC(=O)O)NC(C)=O)C(C)O)C(C)C)C(C)C)C(=O)c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
                false,
            ),
            (
                "CC1(C)CCC2(C(=O)N3CCOCC3)CCC3(C)C(C(=O)C=C4C5(C)C=C(C#N)C(=O)C(C)(C)C5CCC43C)C2C1",
                false,
            ),
        ] {
            let mol = parse(smi).unwrap();
            let layout = compute_layout(&mol);
            let n = mol.atom_count();
            for i in 0..n {
                for j in i + 1..n {
                    let (a, b) = (AtomIdx(i as u32), AtomIdx(j as u32));
                    if mol.bond_between(a, b).is_none() {
                        let d = layout.get(a).dist(&layout.get(b));
                        assert!(d >= 0.4 * BOND_LEN, "{smi}: atoms {i} and {j} {d:.2} apart");
                    }
                }
            }
            // A drawn cage may cross itself (the tropane bridge does).
            assert!(
                cage || detect_crossings(&layout, &mol).is_empty(),
                "{smi}: crossing bonds"
            );
        }
    }
}
