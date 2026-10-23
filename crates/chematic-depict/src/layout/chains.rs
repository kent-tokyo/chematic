//! Acyclic growth and user-facing bond-direction suggestions.

use super::*;
// ---------------------------------------------------------------------------
// Chain placement (DFS zigzag)
// ---------------------------------------------------------------------------

/// If nothing in `component` is placed yet (no ring system exists to seed
/// from), place a terminal atom (degree ≤1, or an arbitrary atom if none)
/// at the origin so [`grow_layout`] has a starting point.
pub(super) fn seed_isolated_chain_start(
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
pub(super) fn grow_layout(
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
pub(super) fn best_outgoing_direction(
    atom: AtomIdx,
    placed: &[Option<Point>],
    mol: &Molecule,
) -> f64 {
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
pub(super) fn direction_is_free(atom: AtomIdx, dir: f64, placed: &[Option<Point>]) -> bool {
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
pub(super) fn ring_centre(candidate: Point, origin: Point, placed: &[Option<Point>]) -> bool {
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
pub(super) fn best_direction_avoiding(used_angles: &[f64]) -> f64 {
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
pub(super) fn ranked_candidates(used_angles: &[f64]) -> Vec<f64> {
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
pub(super) fn min_angle_separation(angle: f64, used: &[f64]) -> f64 {
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
pub(super) fn dfs_zigzag(
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
            // A ring hung on an atom with six or more bonds (an octahedral
            // metal with pyridines) is drawn one and a half bonds out: at
            // one bond its neighbours' rings, 60° apart, touched.
            let reach = if mol.degree(parent) >= 6 {
                1.5 * BOND_LEN
            } else {
                BOND_LEN
            };
            let entry_pos = Point::new(
                parent_pos.x + reach * dir.cos(),
                parent_pos.y + reach * dir.sin(),
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
            let next_dir =
                sp_linear::continuation_direction(mol, atom, dir, back, sign, deflection);
            stack.push((unplaced[0], atom, next_dir, -sign, 0.0));
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
pub(super) fn free_direction(
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
