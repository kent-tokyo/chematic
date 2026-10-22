//! Ring-system grouping and deterministic geometric placement.

use super::*;
// ---------------------------------------------------------------------------
// Ring system grouping
// ---------------------------------------------------------------------------

/// Group rings into ring systems where rings in the same system share at least one atom.
pub(super) fn group_ring_systems(rings: &[Vec<AtomIdx>]) -> Vec<Vec<Vec<AtomIdx>>> {
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
pub(super) fn place_ring_system(
    mol: &Molecule,
    system: &[Vec<AtomIdx>],
    anchor: Option<(AtomIdx, Point, f64)>,
    placed: &mut [Option<Point>],
) {
    if system.is_empty() {
        return;
    }

    let first_ring_idx = initialize_ring_system(system, anchor, placed);

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

fn initialize_ring_system(
    system: &[Vec<AtomIdx>],
    anchor: Option<(AtomIdx, Point, f64)>,
    placed: &mut [Option<Point>],
) -> usize {
    let first_ring_idx = anchor
        .and_then(|(entry_atom, ..)| system.iter().position(|ring| ring.contains(&entry_atom)))
        .unwrap_or(0);
    match anchor {
        Some((entry_atom, entry_pos, dir)) => {
            place_first_ring_anchored(&system[first_ring_idx], entry_atom, entry_pos, dir, placed);
        }
        None => place_regular_ring(&system[first_ring_idx], placed),
    }
    first_ring_idx
}

/// A bridged ring system (two rings sharing three or more atoms), or a
/// fused or spiro system of three or more rings, whose polygon placement
/// leaves non-bonded atoms closer than half a bond (adamantane, morphinans,
/// cages fused to several rings, helicenes, a spiro naphthodioxin, a
/// chelate's fused rings) is redrawn by stress majorization on its
/// ring-graph distances, started from the polygon drawing and turned back
/// onto it, and kept only when it has fewer such pairs.
pub(super) fn untangle_bridged_system(
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
pub(super) fn porphyrin_core(adj: &[Vec<usize>], old: &[Point]) -> Option<Vec<Point>> {
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
pub(super) fn turn_onto(x: &[Point], old: &[Point]) -> Vec<Point> {
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
pub(super) fn adamantane_cage(
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
pub(super) struct NarrowForks<'a> {
    pub(super) angles: &'a [f64],
}

/// Places the unplaced atoms of `ring` when its placed atoms form one path of
/// at least three atoms in ring order and its unplaced atoms the rest: on the
/// regular polygon through the path's two ends (two centres, one on each side
/// of their chord), choosing the one whose new atoms come within half a bond
/// of fewer placed atoms, and only when it is strictly better than retracing
/// the side the path already lies on. Returns whether it placed them.
pub(super) fn place_bridged_ring(ring: &[AtomIdx], placed: &mut [Option<Point>]) -> bool {
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
pub(super) fn find_shared_edge(
    ring: &[AtomIdx],
    placed: &[Option<Point>],
) -> Option<(AtomIdx, AtomIdx)> {
    let n = ring.len();
    ring.windows(2)
        .map(|w| (w[0], w[1]))
        .chain(std::iter::once((ring[n - 1], ring[0])))
        .find(|&(a, b)| placed[a.0 as usize].is_some() && placed[b.0 as usize].is_some())
}

pub(super) fn closest_placed_pair(
    atoms: &[AtomIdx],
    placed: &[Option<Point>],
) -> (AtomIdx, AtomIdx) {
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
pub(super) fn place_regular_ring(ring: &[AtomIdx], placed: &mut [Option<Point>]) {
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
///
/// The candidate centre lies on the shared edge's perpendicular bisector at
/// the regular polygon's apothem. Existing placed atoms decide which side of
/// that edge preserves the established drawing.
pub(super) fn place_ring_anchored(
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

    place_unplaced_ring_atoms(
        ring,
        idx1,
        new_center,
        angle_to_a1,
        angle_step,
        radius,
        placed,
    );
}

fn place_unplaced_ring_atoms(
    ring: &[AtomIdx],
    idx1: usize,
    center: Point,
    angle_to_a1: f64,
    angle_step: f64,
    radius: f64,
    placed: &mut [Option<Point>],
) {
    let n = ring.len();
    for step in 0..n {
        let ring_idx = (idx1 + step) % n;
        let atom = ring[ring_idx];

        if placed[atom.0 as usize].is_some() {
            continue;
        }

        let angle = angle_to_a1 + step as f64 * angle_step;
        let x = center.x + radius * angle.cos();
        let y = center.y + radius * angle.sin();
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
pub(super) fn place_first_ring_anchored(
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
pub(super) fn centroid_of_placed(placed: &[Option<Point>]) -> Option<Point> {
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
pub(super) fn direction_away_from_centroid(entry_pos: Point, placed: &[Option<Point>]) -> f64 {
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
pub(super) fn seed_ring_system_index(ring_systems: &[Vec<Vec<AtomIdx>>]) -> Option<usize> {
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
pub(super) fn ring_radius(n: usize) -> f64 {
    if n < 3 {
        return BOND_LEN;
    }
    BOND_LEN / (2.0 * (std::f64::consts::PI / n as f64).sin())
}
