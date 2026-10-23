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

#[cfg(test)]
mod issue_912_tests;
mod sp_linear;

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

mod rings;
use rings::*;
mod chains;
pub use chains::suggest_bond_direction;
use chains::*;
mod clashes;
pub use clashes::detect_crossings;
use clashes::*;
// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;

#[cfg(test)]
mod geometric_boundary_tests;
