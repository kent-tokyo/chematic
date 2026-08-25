//! Connectivity-ordered ring/chain placement engine (issues #256/#255).
//!
//! See `docs/rfcs/dg_connectivity_ordered_placement_rfc.md`. An additive,
//! parallel implementation of rule-based 3D coordinate generation --
//! [`crate::dg::generate_coords`] itself is completely untouched by this
//! module and calls none of it.
//!
//! Structural technique ported from `chematic-depict/layout.rs`'s
//! `grow_layout`/`dfs_zigzag`/`place_ring_anchored` (RFC §3): a single
//! worklist discovers and places rings and chain atoms in true connectivity
//! order (never "all rings, then all chains"), and ring-fusion center/winding
//! is derived from measurement (candidate center farther from what's already
//! placed; winding sign from the real angle between the two shared anchor
//! atoms) rather than trusting a ring template's own listed atom order.
//!
//! Per the RFC's own recorded recommendations (§4), every ring system's
//! plane is kept parallel to the global XY plane (option (b): "fixed
//! reference normal"), at whatever Z height its entry point naturally lands
//! at -- fused ring systems therefore end up coplanar even where a real 3D
//! fused system would not be. Chain placement's dihedral/bond-angle scheme is
//! [`crate::dg::generate_coords`]'s `dfs_place` unchanged (no
//! grandparent-referenced staggering, per the RFC's recorded recommendation
//! to scope this rewrite strictly to the ring-anchor defects).
//!
//! ## Status (issues #256/#255, 2026-08-25)
//!
//! Phase 1 (this engine) and Phase 2 (differential evaluation against
//! `generate_coords`, PR #387) plus the new-island ring-entry-direction fix
//! (PR #388) together cleared every go/no-go criterion recorded in
//! `ROADMAP.md`: every RFC known-broken topology recovered, a clear
//! improvement on issue #277's real-molecule population with zero new
//! breakage anywhere measured, gross clash and stereo-violation counts do
//! not regress, wall time is comparable, and quality is atom-order-stable.
//! Promoted from `pub(crate)` to `pub` on that basis -- **but `generate_coords`
//! itself, and every production caller that uses it
//! (`generate_coords_etkdg`, `embed_pipeline_v2`, ...), is deliberately
//! unchanged.** This module ships as an available, independently-selectable
//! alternative (RFC §5 Phase 3 option (c)), not a default-behavior switch or
//! topology-based routing (options (a)/(b)) -- those remain separate,
//! not-yet-made decisions.

use core::f64::consts::PI;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::find_sssr;

use crate::coords::{Coords3D, Point3};
use crate::dg::{
    connected_components, ideal_angle, ideal_bond_len, perpendicular_to, rotate_around_axis,
};

/// Group `rings` into ring systems: connected sets of rings sharing >= 1
/// atom (fused/spiro clusters). Two rings joined only by a bond, not a
/// shared atom (e.g. biphenyl's two phenyls), are separate systems -- each
/// gets its own anchored placement via the worklist below, exactly as
/// `crate::dg::generate_coords`'s own `place_rings` already distinguishes
/// "shared atom" fusion from "direct bond, new island" anchoring.
///
/// Ported from `chematic-depict/layout.rs::group_ring_systems` (same
/// algorithm, no chemistry dependency) rather than shared, since that
/// function is private to a different crate and the two engines' ring
/// representations already diverge (this one carries no SVG-specific
/// state).
fn group_ring_systems(rings: &[Vec<AtomIdx>]) -> Vec<Vec<Vec<AtomIdx>>> {
    if rings.is_empty() {
        return Vec::new();
    }
    let n = rings.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        if parent[i] != i {
            parent[i] = find(parent, parent[i]);
        }
        parent[i]
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }
    for i in 0..n {
        let set_i: HashSet<AtomIdx> = rings[i].iter().copied().collect();
        for j in (i + 1)..n {
            if rings[j].iter().any(|a| set_i.contains(a)) {
                union(&mut parent, i, j);
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }
    groups
        .values()
        .map(|idxs| idxs.iter().map(|&i| rings[i].clone()).collect())
        .collect()
}

/// Pick the ring system that anchors a component's coordinate frame: most
/// rings, then most atoms, then lowest minimum `AtomIdx` -- deterministic,
/// prefers a fused/spiro core over a peripheral single ring. Same tie-break
/// as `chematic-depict::seed_ring_system_index`.
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

/// Circumradius of a regular n-gon whose side is `bond_len`.
fn ring_circumradius(bond_len: f64, ring_size: usize) -> f64 {
    bond_len / (2.0 * (PI / ring_size as f64).sin())
}

/// Conformer Z offset at ring-position `step`, mirroring `crate::dg`'s
/// `place_rings` own rules (chair/envelope/crown for isolated, non-aromatic,
/// non-fused rings; flat otherwise). `is_fused` is passed explicitly by the
/// caller rather than inferred from `placed` state's timing (this engine's
/// placement order differs from `place_rings`', so that implicit-timing
/// trick isn't reliable here) -- true for a shared-edge/spiro continuation
/// within an already-anchored ring system, false for a system's own first
/// ring.
fn ring_conf_z(mol: &Molecule, ring: &[AtomIdx], is_fused: bool, step: usize) -> f64 {
    const CHAIR_H: f64 = 0.256;
    const ENVELOPE_H: f64 = 0.400;
    let ring_size = ring.len();
    let is_aromatic = ring.iter().all(|&a| mol.atom(a).aromatic);
    if is_aromatic || is_fused {
        return 0.0;
    }
    match ring_size {
        6 => {
            if step.is_multiple_of(2) {
                CHAIR_H
            } else {
                -CHAIR_H
            }
        }
        5 => {
            if step == ring_size - 1 {
                ENVELOPE_H
            } else {
                0.0
            }
        }
        n if n >= 8 => {
            let h = 0.3 + 0.04 * (n as f64 - 8.0).min(10.0);
            if step.is_multiple_of(2) { h } else { -h }
        }
        _ => 0.0,
    }
}

/// Chair conformation uses a geometry-derived XY radius distinct from the
/// regular-polygon circumradius (exact solution for l=1.54 Å, θ=109.5°).
fn ring_effective_radius(ring: &[AtomIdx], is_fused: bool, mol: &Molecule, r: f64) -> f64 {
    const CHAIR_R: f64 = 1.452;
    let is_aromatic = ring.iter().all(|&a| mol.atom(a).aromatic);
    if !is_aromatic && !is_fused && ring.len() == 6 {
        CHAIR_R
    } else {
        r
    }
}

/// Centroid (X, Y only) of every currently-placed atom in the component.
/// `None` if nothing is placed yet.
fn centroid_of_placed_xy(coords: &Coords3D, placed: &[bool]) -> Option<(f64, f64)> {
    let pts: Vec<(f64, f64)> = placed
        .iter()
        .enumerate()
        .filter(|&(_, &p)| p)
        .map(|(i, _)| {
            let p = coords.get(AtomIdx(i as u32));
            (p.x, p.y)
        })
        .collect();
    if pts.is_empty() {
        return None;
    }
    let n = pts.len() as f64;
    let (sx, sy) = pts
        .iter()
        .fold((0.0, 0.0), |(sx, sy), &(x, y)| (sx + x, sy + y));
    Some((sx / n, sy / n))
}

/// Direction (unit XY vector) from the centroid of everything placed so far,
/// through `entry_pos`, continued outward -- the "grow away from what's
/// already there" direction used for a spiro anchor, which (unlike a
/// chain/exocyclic attachment) has no incoming bond direction to continue.
fn direction_away_from_centroid_xy(
    entry_pos: Point3,
    coords: &Coords3D,
    placed: &[bool],
) -> (f64, f64) {
    match centroid_of_placed_xy(coords, placed) {
        Some((cx, cy)) => {
            let dx = entry_pos.x - cx;
            let dy = entry_pos.y - cy;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist > 1e-9 {
                (dx / dist, dy / dist)
            } else {
                (1.0, 0.0)
            }
        }
        None => (1.0, 0.0),
    }
}

/// Candidate roll angles (around the incoming-bond axis) tried when choosing
/// the entry direction for a newly-discovered ring system reached via a
/// direct bond (see [`choose_ring_entry_direction_xy`]).
const RING_ENTRY_CANDIDATE_COUNT: usize = 12;

/// Choose the outgoing XY direction for the first ring atom of a
/// newly-discovered ring system reached from `current` via a direct bond
/// (issue #256/#255 Phase 2's "new-island" regression: a plain chain
/// neighbor's single `i * 120°` roll has no reason to avoid the rest of the
/// already-placed structure, and post-minimisation that shows up as real
/// clashes for e.g. biphenyl/terphenyl).
///
/// Does not abandon the existing bend/dihedral scheme -- every candidate is
/// still `dir_bent` rolled around `incoming_dir`, just at
/// [`RING_ENTRY_CANDIDATE_COUNT`] angles instead of one. Ranked by: (1) XY
/// distance from the centroid of everything placed so far (farther is
/// better -- extend away from the existing structure, mirroring
/// `chematic-depict::place_ring_anchored`'s own documented reason for using
/// the whole-structure centroid rather than a single ring's own anchor,
/// which degenerates to a tie); (2) on a tie, the minimum XY distance from
/// the candidate entry position to any already-placed atom (larger is
/// better -- more room); (3) a final deterministic tie-break seeded by
/// `entry_atom`'s own `AtomIdx`, not candidate/array order -- an
/// identity-less index tie-break is wrong here (see this crate's own
/// permutation-invariance precedent).
#[allow(clippy::too_many_arguments)]
fn choose_ring_entry_direction_xy(
    entry_atom: AtomIdx,
    pos_current: Point3,
    incoming_dir: Point3,
    dir_bent: Point3,
    bond_len: f64,
    coords: &Coords3D,
    placed: &[bool],
) -> (f64, f64) {
    let centroid = centroid_of_placed_xy(coords, placed);
    let start = (entry_atom.0 as usize) % RING_ENTRY_CANDIDATE_COUNT;

    let mut best_score: Option<(f64, f64)> = None; // (centroid_dist, min_clearance)
    let mut best_dir = (1.0, 0.0);
    for offset in 0..RING_ENTRY_CANDIDATE_COUNT {
        let k = (start + offset) % RING_ENTRY_CANDIDATE_COUNT;
        let dihedral = k as f64 * (2.0 * PI / RING_ENTRY_CANDIDATE_COUNT as f64);
        let dir_final = rotate_around_axis(dir_bent, incoming_dir, dihedral);
        let xy_len = (dir_final.x * dir_final.x + dir_final.y * dir_final.y).sqrt();
        let dir_xy = if xy_len > 1e-9 {
            (dir_final.x / xy_len, dir_final.y / xy_len)
        } else {
            (1.0, 0.0)
        };
        let entry_xy = (
            pos_current.x + dir_xy.0 * bond_len,
            pos_current.y + dir_xy.1 * bond_len,
        );
        let centroid_dist = match centroid {
            Some((cx, cy)) => {
                let dx = entry_xy.0 - cx;
                let dy = entry_xy.1 - cy;
                (dx * dx + dy * dy).sqrt()
            }
            None => 0.0,
        };
        let min_clearance = placed
            .iter()
            .enumerate()
            .filter(|&(_, &p)| p)
            .map(|(idx, _)| {
                let p = coords.get(AtomIdx(idx as u32));
                let dx = entry_xy.0 - p.x;
                let dy = entry_xy.1 - p.y;
                (dx * dx + dy * dy).sqrt()
            })
            .fold(f64::INFINITY, f64::min);

        let is_better = match best_score {
            None => true,
            Some((best_centroid_dist, best_min_clearance)) => {
                if centroid_dist > best_centroid_dist + 1e-9 {
                    true
                } else if centroid_dist < best_centroid_dist - 1e-9 {
                    false
                } else {
                    min_clearance > best_min_clearance + 1e-9
                }
            }
        };
        if is_better {
            best_score = Some((centroid_dist, min_clearance));
            best_dir = dir_xy;
        }
    }
    best_dir
}

/// Two consecutive-in-ring-order atoms that are both placed (the shared
/// fusion edge), if any.
fn find_shared_edge(ring: &[AtomIdx], placed: &[bool]) -> Option<(AtomIdx, AtomIdx)> {
    let n = ring.len();
    ring.windows(2)
        .map(|w| (w[0], w[1]))
        .chain(std::iter::once((ring[n - 1], ring[0])))
        .find(|&(a, b)| placed[a.0 as usize] && placed[b.0 as usize])
}

/// Place `ring` as a regular polygon with `entry_atom` on its circumference
/// at `entry_pos`, extending outward in direction `dir_xy` (a unit vector in
/// the ring's own XY-parallel plane -- see this module's doc for why only XY
/// matters here). Ported from
/// `chematic-depict/layout.rs::place_first_ring_anchored`, with `dg.rs`'s own
/// chemistry-aware bond length and conformer puckering substituted for
/// depict's fixed `BOND_LEN`/always-flat 2D polygon.
#[allow(clippy::too_many_arguments)]
fn place_first_ring_anchored_3d(
    mol: &Molecule,
    ring: &[AtomIdx],
    entry_atom: AtomIdx,
    entry_pos: Point3,
    dir_xy: (f64, f64),
    is_fused: bool,
    coords: &mut Coords3D,
    placed: &mut [bool],
) {
    let n = ring.len();
    if n == 0 {
        return;
    }
    if !placed[entry_atom.0 as usize] {
        coords.set(entry_atom, entry_pos);
        placed[entry_atom.0 as usize] = true;
    }
    let bond_len = ideal_bond_len(mol, ring[0], ring[1 % n]);
    let r = ring_circumradius(bond_len, n);
    let effective_r = ring_effective_radius(ring, is_fused, mol, r);
    let center_x = entry_pos.x + r * dir_xy.0;
    let center_y = entry_pos.y + r * dir_xy.1;
    let base_z = entry_pos.z;
    let idx0 = ring.iter().position(|&a| a == entry_atom).unwrap_or(0);
    let angle_to_entry = (entry_pos.y - center_y).atan2(entry_pos.x - center_x);
    let angle_step = -2.0 * PI / n as f64;

    for step in 0..n {
        let ring_idx = (idx0 + step) % n;
        let atom = ring[ring_idx];
        if placed[atom.0 as usize] {
            continue;
        }
        let angle = angle_to_entry + step as f64 * angle_step;
        let x = center_x + effective_r * angle.cos();
        let y = center_y + effective_r * angle.sin();
        let z = base_z + ring_conf_z(mol, ring, is_fused, step);
        coords.set(atom, Point3::new(x, y, z));
        placed[atom.0 as usize] = true;
    }
}

/// Place unplaced atoms of `ring` given that `anchor1`/`anchor2` (a shared
/// fusion edge) are already placed. The new ring's center is whichever of
/// the shared edge's two perpendicular-bisector candidates lies farther
/// from the centroid of everything already placed, and its winding
/// direction (CW vs CCW) comes from the real, measured signed angle between
/// the two anchors around that chosen center -- never from the ring's own
/// listed atom order. This is issue #255's actual fix: `crate::dg`'s
/// `place_rings` old `shared_atoms.len() >= 2` branch used a fixed `+y`
/// extension instead, correct only when +y happened to point away from the
/// already-placed ring. Ported from
/// `chematic-depict/layout.rs::place_ring_anchored` (verified there against
/// naphthalene/phenanthrene/pyrene/anthracene via `mol.depict_data()`, see
/// the RFC §3b), with dg.rs's own chemistry-aware bond length substituted
/// for `BOND_LEN`, and a Z coordinate added (the anchor edge's shared height
/// -- both anchor atoms already sit at the same Z in every case this engine
/// produces, since ring planes are always XY-parallel).
#[allow(clippy::too_many_arguments)]
fn place_ring_anchored_3d(
    mol: &Molecule,
    ring: &[AtomIdx],
    anchor1: AtomIdx,
    p1: Point3,
    anchor2: AtomIdx,
    p2: Point3,
    coords: &mut Coords3D,
    placed: &mut [bool],
) {
    let n = ring.len();
    if n == 0 {
        return;
    }
    let bond_len = ideal_bond_len(mol, ring[0], ring[1 % n]);
    let r = ring_circumradius(bond_len, n);
    let idx1 = ring.iter().position(|&a| a == anchor1).unwrap_or(0);
    let idx2 = ring.iter().position(|&a| a == anchor2).unwrap_or(1);

    let mid = ((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let edge_len = (dx * dx + dy * dy).sqrt();
    if edge_len < 1e-10 {
        return;
    }
    let (perp_x, perp_y) = (-dy / edge_len, dx / edge_len);
    let apothem = r * (PI / n as f64).cos();

    let existing_center = centroid_of_placed_xy(coords, placed).unwrap_or(mid);
    let cand1 = (mid.0 + perp_x * apothem, mid.1 + perp_y * apothem);
    let cand2 = (mid.0 - perp_x * apothem, mid.1 - perp_y * apothem);
    let dist2 = |p: (f64, f64), q: (f64, f64)| {
        let dx = p.0 - q.0;
        let dy = p.1 - q.1;
        dx * dx + dy * dy
    };
    let new_center = if dist2(cand1, existing_center) > dist2(cand2, existing_center) {
        cand1
    } else {
        cand2
    };

    let angle_to_a1 = (p1.y - new_center.1).atan2(p1.x - new_center.0);
    let angle_to_a2 = (p2.y - new_center.1).atan2(p2.x - new_center.0);
    let mut delta = angle_to_a2 - angle_to_a1;
    while delta > PI {
        delta -= 2.0 * PI;
    }
    while delta < -PI {
        delta += 2.0 * PI;
    }
    let steps_forward = (idx2 + n - idx1) % n;
    let angle_step = if steps_forward > 0 {
        if delta >= 0.0 {
            2.0 * PI / n as f64
        } else {
            -(2.0 * PI / n as f64)
        }
    } else {
        2.0 * PI / n as f64
    };
    let base_z = (p1.z + p2.z) / 2.0;

    for step in 0..n {
        let ring_idx = (idx1 + step) % n;
        let atom = ring[ring_idx];
        if placed[atom.0 as usize] {
            continue;
        }
        let angle = angle_to_a1 + step as f64 * angle_step;
        let x = new_center.0 + r * angle.cos();
        let y = new_center.1 + r * angle.sin();
        coords.set(atom, Point3::new(x, y, base_z));
        placed[atom.0 as usize] = true;
    }
}

/// Place every ring of `system` (a connected fused/spiro cluster) not yet
/// placed, given that `system[first_ring_idx]` is already fully placed by
/// the caller (either the component's seed ring, or the first ring of a
/// chain-anchored system via [`place_first_ring_anchored_3d`]). Subsequent
/// rings are placed via a retry-until-no-progress loop over shared-edge
/// fusion or single-atom spiro, exactly like
/// `chematic-depict/layout.rs::place_ring_system`'s own loop -- this
/// sidesteps needing any SSSR-order-dependent visiting order at all (the bug
/// `crate::dg`'s `order_rings_by_fusion_adjacency` exists to work around for
/// the old engine): a ring not yet ready to place is simply retried next
/// pass.
fn grow_ring_system_3d(
    mol: &Molecule,
    system: &[Vec<AtomIdx>],
    first_ring_idx: usize,
    coords: &mut Coords3D,
    placed: &mut [bool],
) {
    let mut remaining: Vec<&Vec<AtomIdx>> = system
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != first_ring_idx)
        .map(|(_, ring)| ring)
        .collect();
    let mut iterations = 0;
    while !remaining.is_empty() && iterations < remaining.len() * 2 + 1 {
        iterations += 1;
        let mut progressed = false;
        remaining.retain(|ring| {
            let already: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|a| placed[a.0 as usize])
                .collect();
            if already.is_empty() {
                return true;
            }
            if already.len() == 1 {
                let entry_atom = already[0];
                let entry_pos = coords.get(entry_atom);
                let dir_xy = direction_away_from_centroid_xy(entry_pos, coords, placed);
                place_first_ring_anchored_3d(
                    mol, ring, entry_atom, entry_pos, dir_xy, true, coords, placed,
                );
                progressed = true;
                return false;
            }
            let (a1, a2) = find_shared_edge(ring, placed).unwrap_or((already[0], already[1]));
            let p1 = coords.get(a1);
            let p2 = coords.get(a2);
            place_ring_anchored_3d(mol, ring, a1, p1, a2, p2, coords, placed);
            progressed = true;
            false
        });
        if !progressed {
            break;
        }
    }
    // Defensive fallback -- should not fire for any ring system reachable
    // from a shared atom or edge (every ring system by construction is
    // connected via atom-sharing; `group_ring_systems` only groups rings
    // that already share an atom).
    for ring in &remaining {
        if let Some(&entry_atom) = ring.iter().find(|a| placed[a.0 as usize]) {
            let entry_pos = coords.get(entry_atom);
            place_first_ring_anchored_3d(
                mol,
                ring,
                entry_atom,
                entry_pos,
                (1.0, 0.0),
                true,
                coords,
                placed,
            );
        }
    }
}

/// Place the ring system that seeds a component's coordinate frame: its
/// first ring centered at `(x_offset + r, 0, 0)` (matching `crate::dg`'s
/// `place_rings` own convention for a component's first-placed ring), then
/// the rest of the system grown via [`grow_ring_system_3d`].
fn place_seed_ring_system_3d(
    mol: &Molecule,
    system: &[Vec<AtomIdx>],
    x_offset: f64,
    coords: &mut Coords3D,
    placed: &mut [bool],
) {
    if system.is_empty() {
        return;
    }
    let first_ring = &system[0];
    let n = first_ring.len();
    let bond_len = ideal_bond_len(mol, first_ring[0], first_ring[1 % n]);
    let r = ring_circumradius(bond_len, n);
    let effective_r = ring_effective_radius(first_ring, false, mol, r);
    let center_x = x_offset + r;
    let center_y = 0.0;
    for (k, &atom) in first_ring.iter().enumerate() {
        let angle = 2.0 * PI * k as f64 / n as f64;
        let x = center_x + effective_r * angle.cos();
        let y = center_y + effective_r * angle.sin();
        let z = ring_conf_z(mol, first_ring, false, k);
        coords.set(atom, Point3::new(x, y, z));
        placed[atom.0 as usize] = true;
    }
    grow_ring_system_3d(mol, system, 0, coords, placed);
}

/// DFS-based placement of atoms not yet positioned, ring-system-aware.
///
/// Identical bond-angle/dihedral geometry to `crate::dg`'s `dfs_place` for a
/// plain chain continuation (deliberately duplicated, not shared --
/// `dfs_place` itself must stay untouched per the RFC's "parallel
/// implementation, never an in-place rewrite" scope). The one real
/// difference: before extending toward an unplaced neighbor, check whether
/// that neighbor belongs to a not-yet-placed ring system (`atom_to_system`).
/// If so, anchor the WHOLE system there (flattening the chosen bond
/// direction's XY component only -- see this module's doc for why) and push
/// its atoms onto `worklist` for the caller to continue growing from,
/// instead of recursing directly -- a ring atom may have several of its own
/// further branches, each needing its own bond-angle/dihedral context,
/// exactly as `crate::dg`'s `place_component`'s own
/// `for atom in component { dfs_place(...) }` loop already does for every
/// ring atom today.
#[allow(clippy::too_many_arguments)]
fn dfs_place_connectivity_ordered(
    mol: &Molecule,
    current: AtomIdx,
    placed: &mut [bool],
    coords: &mut Coords3D,
    ring_systems: &[Vec<Vec<AtomIdx>>],
    atom_to_system: &HashMap<AtomIdx, usize>,
    system_placed: &mut [bool],
    worklist: &mut VecDeque<AtomIdx>,
) {
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

    let incoming_dir: Point3 = match parent {
        Some(p) => pos_current.sub(&coords.get(p)).normalize(),
        None => Point3::new(1.0, 0.0, 0.0),
    };
    let perp = perpendicular_to(incoming_dir);
    let angle = ideal_angle(mol, current);
    let bend_angle = PI - angle;
    let dir_bent = rotate_around_axis(incoming_dir, perp, bend_angle);

    for (i, &nb) in unplaced_neighbors.iter().enumerate() {
        if placed[nb.0 as usize] {
            // Already placed by an earlier iteration of this same loop (a
            // ring system spanning two of `current`'s own neighbors --
            // bridged-bicyclic territory, out of Phase 0/1 scope, but
            // guarded defensively rather than double-placed).
            continue;
        }

        if let Some(&sys_idx) = atom_to_system.get(&nb) {
            if system_placed[sys_idx] {
                continue;
            }
            let bond_len = ideal_bond_len(mol, current, nb);
            let dir_xy = choose_ring_entry_direction_xy(
                nb,
                pos_current,
                incoming_dir,
                dir_bent,
                bond_len,
                coords,
                placed,
            );
            let entry_pos = Point3::new(
                pos_current.x + dir_xy.0 * bond_len,
                pos_current.y + dir_xy.1 * bond_len,
                pos_current.z,
            );
            let system = &ring_systems[sys_idx];
            let first_ring_idx = system.iter().position(|r| r.contains(&nb)).unwrap_or(0);
            place_first_ring_anchored_3d(
                mol,
                &system[first_ring_idx],
                nb,
                entry_pos,
                dir_xy,
                false,
                coords,
                placed,
            );
            grow_ring_system_3d(mol, system, first_ring_idx, coords, placed);
            system_placed[sys_idx] = true;
            worklist.extend(system.iter().flatten().copied());
            continue;
        }

        let dihedral = (i as f64) * (2.0 * PI / 3.0);
        let dir_final = rotate_around_axis(dir_bent, incoming_dir, dihedral);
        let bond_len = ideal_bond_len(mol, current, nb);
        let new_pos = pos_current.add(&dir_final.scale(bond_len));
        coords.set(nb, new_pos);
        placed[nb.0 as usize] = true;
        dfs_place_connectivity_ordered(
            mol,
            nb,
            placed,
            coords,
            ring_systems,
            atom_to_system,
            system_placed,
            worklist,
        );
    }
}

/// Place all atoms in `component`, starting at X = `x_offset`, via the
/// connectivity-ordered engine (RFC §5, Phase 1).
fn place_component_connectivity_ordered(
    mol: &Molecule,
    component: &[AtomIdx],
    ring_set: &chematic_perception::RingSet,
    x_offset: f64,
    coords: &mut Coords3D,
) {
    if component.is_empty() {
        return;
    }
    let component_set: HashSet<AtomIdx> = component.iter().copied().collect();
    let relevant_rings: Vec<Vec<AtomIdx>> = ring_set
        .rings()
        .iter()
        .filter(|ring| !ring.is_empty() && ring.iter().all(|a| component_set.contains(a)))
        .cloned()
        .collect();
    let ring_systems = group_ring_systems(&relevant_rings);
    let mut atom_to_system: HashMap<AtomIdx, usize> = HashMap::new();
    for (i, system) in ring_systems.iter().enumerate() {
        for &a in system.iter().flatten() {
            atom_to_system.insert(a, i);
        }
    }
    let mut system_placed = vec![false; ring_systems.len()];
    let mut placed = vec![false; mol.atom_count()];
    let mut worklist: VecDeque<AtomIdx> = VecDeque::new();

    if let Some(seed_idx) = seed_ring_system_index(&ring_systems) {
        place_seed_ring_system_3d(mol, &ring_systems[seed_idx], x_offset, coords, &mut placed);
        system_placed[seed_idx] = true;
        worklist.extend(ring_systems[seed_idx].iter().flatten().copied());
    } else {
        let start = component
            .iter()
            .copied()
            .find(|&a| mol.degree(a) <= 1)
            .unwrap_or(component[0]);
        coords.set(start, Point3::new(x_offset, 0.0, 0.0));
        placed[start.0 as usize] = true;
        worklist.push_back(start);
    }

    while let Some(atom) = worklist.pop_front() {
        if !placed[atom.0 as usize] {
            continue;
        }
        dfs_place_connectivity_ordered(
            mol,
            atom,
            &mut placed,
            coords,
            &ring_systems,
            &atom_to_system,
            &mut system_placed,
            &mut worklist,
        );
    }

    // Defensive fallback -- should not fire for any connected component,
    // since the worklist above reaches every atom reachable from the seed
    // via the molecule graph (identical guarantee to
    // `chematic-depict::grow_layout`'s own analogous fallback).
    for (i, system) in ring_systems.iter().enumerate() {
        if !system_placed[i] {
            place_seed_ring_system_3d(mol, system, x_offset, coords, &mut placed);
        }
    }
    for &atom in component {
        if !placed[atom.0 as usize] {
            coords.set(atom, Point3::new(x_offset, 0.0, 0.0));
            placed[atom.0 as usize] = true;
        }
    }
}

/// Generate 3D coordinates for all heavy atoms in `mol` via the
/// connectivity-ordered placement engine.
///
/// Same contract as [`crate::dg::generate_coords`] (non-degenerate,
/// ideal-bond-length-approximating, not physically minimised) -- an
/// independently-selectable alternative implementation, not the default: no
/// existing production caller (`generate_coords_etkdg`, `embed_pipeline_v2`,
/// ...) is routed through this engine. See this module's doc for the
/// measured comparison against `generate_coords` and the go/no-go criteria
/// this engine cleared before being made public.
pub fn generate_coords_connectivity_ordered(mol: &Molecule) -> Coords3D {
    let n = mol.atom_count();
    let mut coords = Coords3D::new_zeroed(n);
    if n == 0 {
        return coords;
    }
    if n == 1 {
        coords.set(AtomIdx(0), Point3::zero());
        return coords;
    }
    let ring_set = find_sssr(mol);
    let components = connected_components(mol);
    let mut x_offset = 0.0_f64;
    for component in &components {
        place_component_connectivity_ordered(mol, component, &ring_set, x_offset, &mut coords);
        let max_x = component
            .iter()
            .map(|&idx| coords.get(idx).x)
            .fold(f64::NEG_INFINITY, f64::max);
        x_offset = max_x + 5.0;
    }
    coords
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    // Acceptance oracle for every fixture below: the RFC's own measured
    // worst-bonded-pair-distance table (docs/rfcs/
    // dg_connectivity_ordered_placement_rfc.md §2), cross-checked against
    // `crate::dg`'s own `generate_coords_*_known_broken` fixtures which pin
    // the corresponding OLD-engine distortion.

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
    /// `[min_len, max_len]` Å.
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

    /// Both checks together (`assert_bonded_pairs_sane` alone cannot catch
    /// an exact atom coincidence with otherwise-correct bond lengths -- the
    /// RFC's own history records the reverted #255 fix attempt doing
    /// exactly that on phenanthrene/pyrene, "two EXACT atom coincidences
    /// plus a 3.7 Å stretch" and "an exact atom coincidence" respectively).
    fn assert_geometry_sane(mol: &Molecule, coords: &Coords3D) {
        let n = mol.atom_count();
        let min_d = min_pairwise_distance(coords, n);
        assert!(
            min_d > 0.3,
            "no two atoms should be nearly coincident, got min pairwise distance {min_d}"
        );
        assert_bonded_pairs_sane(mol, coords, 1.0, 1.8);
    }

    #[test]
    fn connectivity_ordered_naphthalene_fusion_seam_sane() {
        // #255: `crate::dg::generate_coords` distorts this to 2.2644 Å (see
        // `generate_coords_naphthalene_fusion_seam_known_broken` in dg.rs).
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        assert_eq!(mol.atom_count(), 10, "naphthalene has 10 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_quinoline_fused_heterocycle_sane() {
        let mol = parse("c1ccc2ncccc2c1").unwrap();
        assert_eq!(mol.atom_count(), 10, "quinoline has 10 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_phenanthrene_angular_fusion_sane() {
        // The reverted #255 fix attempt produced "two EXACT atom
        // coincidences plus a 3.7 Å stretch" on this exact molecule (RFC
        // §2) -- `assert_geometry_sane`'s collision check, not just the
        // bond-length check, is what would catch a repeat of that failure.
        let mol = parse("c1ccc2c(c1)ccc1ccccc12").unwrap();
        assert_eq!(mol.atom_count(), 14, "phenanthrene has 14 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_pyrene_multi_ring_fusion_sane() {
        // The reverted #255 fix attempt produced "an exact atom
        // coincidence" here too (RFC §2) -- same reasoning as phenanthrene
        // above.
        let mol = parse("c1cc2ccc3cccc4ccc(c1)c2c34").unwrap();
        assert_eq!(mol.atom_count(), 16, "pyrene has 16 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_anthracene_linear_fusion_sane() {
        // Not in the RFC's own known-broken table (linear 3-ring fusion is
        // structurally the same kind of case as phenanthrene/pyrene above,
        // just unbranched) -- extra coverage beyond the RFC's fixture list.
        let mol = parse("c1ccc2cc3ccccc3cc2c1").unwrap();
        assert_eq!(mol.atom_count(), 14, "anthracene has 14 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_diphenylmethane_chain_bridge_length_1_sane() {
        // #256, shortest bridge: `generate_coords` distorts this to 8.0738 Å
        // (see the `..._known_broken` test in dg.rs).
        let mol = parse("c1ccccc1Cc1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 13, "diphenylmethane has 13 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_bibenzyl_chain_bridged_ring_islands_sane() {
        // #256, bridge length 2: `generate_coords` distorts this to
        // 8.7358 Å.
        let mol = parse("c1ccccc1CCc1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 14, "bibenzyl has 14 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_diphenylpropane_chain_bridge_length_3_sane() {
        let mol = parse("c1ccccc1CCCc1ccccc1").unwrap();
        assert_eq!(
            mol.atom_count(),
            15,
            "1,3-diphenylpropane has 15 heavy atoms"
        );
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_diphenylbutane_chain_bridge_length_4_sane() {
        let mol = parse("c1ccccc1CCCCc1ccccc1").unwrap();
        assert_eq!(
            mol.atom_count(),
            16,
            "1,4-diphenylbutane has 16 heavy atoms"
        );
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_spiro_ring_adjacency_sane() {
        // Positive control (RFC §4: spiro must remain sound under the new
        // traversal too, not just the old one).
        let mol = parse("C1CCC2(CC1)CCCCC2").unwrap();
        assert_eq!(
            mol.atom_count(),
            11,
            "spiro[5.5]undecane has 11 heavy atoms"
        );
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_biphenyl_new_island_sane() {
        // Positive control: two rings joined by a real bond, sharing no
        // atom -- a separate ring system per `group_ring_systems`, anchored
        // via `dfs_place_connectivity_ordered`'s direct-bond path, not
        // fusion.
        let mol = parse("c1ccc(cc1)-c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 12, "biphenyl has 12 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_terphenyl_chain_of_islands_sane() {
        let mol = parse("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 18, "terphenyl has 18 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_meta_linked_biaryl_sane() {
        // Highest-risk untested topology for this engine's ring-entry
        // direction choice: the old engine's fixed +X extension pointed
        // this ring straight back into what was already placed (measured
        // 0.14 Å min pairwise distance before its own centroid-outward
        // fix) BECAUSE the connecting ring atom sits on the side of its
        // ring FACING the already-placed structure, unlike biphenyl's
        // para-like case.
        let mol = parse("c1ccc(cc1)-c1cccnc1").unwrap();
        assert_eq!(mol.atom_count(), 12, "3-phenylpyridine has 12 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_ring_with_tail_substituent_all_atoms_placed() {
        // ibuprofen-shaped: a multi-atom substituent chain hanging off a
        // ring atom must get a real position, not the (0,0,0) default.
        let mol = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
        assert_eq!(mol.atom_count(), 15, "ibuprofen has 15 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_pure_chain_no_ring_sane() {
        // No ring anywhere in the component -- exercises the isolated
        // chain-start seeding path, not ring-system seeding. Pentane, not a
        // longer alkane: `dfs_place_connectivity_ordered`'s plain-chain
        // math is deliberately byte-identical to `crate::dg`'s `dfs_place`
        // own (see that function's doc), which already has a known,
        // pre-existing, unrelated-to-this-engine property that a long
        // enough unbranched chain's fixed bend+0°-roll-per-step placement
        // self-approaches (confirmed identical on both engines: octane
        // already measures 0.1948 Å min pairwise distance under
        // `generate_coords` too, not something this engine introduces).
        // Pentane stays at the terminal bond length (1.54 Å) on both
        // engines, so this fixture actually tests the seeding path itself
        // rather than being confounded by that separate, pre-existing
        // limitation.
        let mol = parse("CCCCC").unwrap();
        assert_eq!(mol.atom_count(), 5, "pentane has 5 heavy atoms");
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_geometry_sane(&mol, &coords);
    }

    #[test]
    fn connectivity_ordered_deterministic_across_runs() {
        // RFC §4 "atom-order reproducibility": no HashMap/HashSet iteration
        // order should leak into the output.
        let mol = parse("c1cc2ccc3cccc4ccc(c1)c2c34").unwrap(); // pyrene
        let a = generate_coords_connectivity_ordered(&mol);
        let b = generate_coords_connectivity_ordered(&mol);
        for i in 0..mol.atom_count() {
            let pa = a.get(AtomIdx(i as u32));
            let pb = b.get(AtomIdx(i as u32));
            assert_eq!(pa.x, pb.x);
            assert_eq!(pa.y, pb.y);
            assert_eq!(pa.z, pb.z);
        }
    }

    #[test]
    fn connectivity_ordered_single_atom_and_empty() {
        let mol = parse("C").unwrap();
        let coords = generate_coords_connectivity_ordered(&mol);
        assert_eq!(coords.atom_count(), 1);
        let p = coords.get(AtomIdx(0));
        assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
    }

    #[test]
    fn connectivity_ordered_all_finite_no_nan() {
        // Broad smoke check across every fixture used above, in one place.
        for smiles in [
            "c1ccc2ccccc2c1",             // naphthalene
            "c1ccc2ncccc2c1",             // quinoline
            "c1ccc2c(c1)ccc1ccccc12",     // phenanthrene
            "c1cc2ccc3cccc4ccc(c1)c2c34", // pyrene
            "c1ccc2cc3ccccc3cc2c1",       // anthracene
            "c1ccccc1Cc1ccccc1",          // diphenylmethane
            "c1ccccc1CCc1ccccc1",         // bibenzyl
            "c1ccccc1CCCc1ccccc1",        // 1,3-diphenylpropane
            "c1ccccc1CCCCc1ccccc1",       // 1,4-diphenylbutane
            "C1CCC2(CC1)CCCCC2",          // spiro[5.5]undecane
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O", // ibuprofen
        ] {
            let mol = parse(smiles).unwrap();
            let coords = generate_coords_connectivity_ordered(&mol);
            for i in 0..mol.atom_count() {
                let p = coords.get(AtomIdx(i as u32));
                assert!(
                    p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                    "non-finite coordinate for {smiles} atom {i}: {p:?}"
                );
            }
        }
    }
}
