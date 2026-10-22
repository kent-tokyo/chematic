//! Crossing detection and deterministic local clash relief.

use super::*;
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
pub(super) struct Bonded<'a>(pub(super) &'a Molecule);

impl Bonded<'_> {
    pub(super) fn contains(&self, &(a, b): &(AtomIdx, AtomIdx)) -> bool {
        self.0.bond_between(a, b).is_some()
    }
}

/// Non-bonded atom pairs closer than this are a clash (the quality lane
/// counts pairs under 0.4 bond lengths).
const CLASH_DIST: f64 = 0.45 * BOND_LEN;

/// Atoms within two bonds of a clashing atom or a crossing bond's atom.
pub(super) fn atoms_near_defects(
    mol: &Molecule,
    atoms: &[AtomIdx],
    placed: &[Option<Point>],
) -> Vec<bool> {
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
pub(super) fn component_bonds(mol: &Molecule, atoms: &[AtomIdx]) -> Vec<(AtomIdx, AtomIdx)> {
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
pub(super) fn layout_defects(
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
pub(super) fn relieve_clashes(
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
pub(super) fn stacked_pairs(mol: &Molecule, atoms: &[AtomIdx], placed: &[Option<Point>]) -> usize {
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
pub(super) fn relieve_clashes_pass(
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
pub(super) fn stacks(
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
pub(super) struct RoundIndex {
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
pub(super) struct DenseGrid<T> {
    x0: i64,
    y0: i64,
    width: i64,
    height: i64,
    /// Cell `c`'s entries are `items[starts[c]..starts[c + 1]]`.
    starts: Vec<u32>,
    items: Vec<T>,
}

impl<T: Copy> DenseGrid<T> {
    pub(super) fn new(entries: &[((i64, i64), T)]) -> Self {
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

    pub(super) fn get(&self, (x, y): (i64, i64)) -> Option<&[T]> {
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

pub(super) fn cell(p: Point, size: f64) -> (i64, i64) {
    // `floor` without the library call: truncate, then step down for a
    // negative non-integer (exact for the coordinate range a layout uses).
    let floor = |v: f64| {
        let t = v as i64;
        if (t as f64) > v { t - 1 } else { t }
    };
    (floor(p.x / size), floor(p.y / size))
}

/// The grid cells a box (with a margin) overlaps.
pub(super) fn box_cells(a: Point, b: Point) -> impl Iterator<Item = (i64, i64)> {
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

pub(super) fn boxes_meet(pa: Point, pb: Point, pc: Point, pd: Point) -> bool {
    pa.x.max(pb.x) + BOX_EPS >= pc.x.min(pd.x)
        && pc.x.max(pd.x) + BOX_EPS >= pa.x.min(pb.x)
        && pa.y.max(pb.y) + BOX_EPS >= pc.y.min(pd.y)
        && pc.y.max(pd.y) + BOX_EPS >= pa.y.min(pb.y)
}

impl RoundIndex {
    pub(super) fn new(
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
pub(super) fn defects_after_move(
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
pub(super) fn reflect_point(q: Point, a: Point, b: Point) -> Point {
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
pub(super) fn segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool {
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
