//! 2D coordinates and wedge bonds that carry a molecule's stereo into a MOL
//! block.
//!
//! MDL files express tetrahedral stereo as a wedge or hash bond on 2D
//! coordinates and E/Z stereo by the coordinates alone. A molecule parsed
//! from SMILES has neither, and writing it with zero coordinates dropped
//! every tetrahedral centre (and wrote SMILES `/`/`\` markers as wedge codes
//! that only chematic read back). [`stereo_depiction`] lays the molecule out,
//! reflects substituents so each stereo double bond has its E/Z geometry,
//! and picks one wedge or hash per tetrahedral centre, each checked by
//! re-reading it the way the MOL reader does.

use std::collections::HashMap;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder};

/// A wedge (`Up`) or hash (`Down`) bond drawn from `start` (the stereocentre).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wedge {
    pub start: AtomIdx,
    pub order: BondOrder,
}

/// Coordinates (Å, y up) and wedges for writing `mol` with its stereo.
#[derive(Debug, Clone, Default)]
pub struct StereoDepiction {
    pub coords: Vec<(f64, f64)>,
    pub wedges: HashMap<BondIdx, Wedge>,
    /// Tetrahedral centres no wedge on these coordinates could express.
    pub unexpressed_centres: Vec<AtomIdx>,
    /// Stereo double bonds whose geometry could not be set (ring bonds).
    pub unexpressed_double_bonds: Vec<BondIdx>,
    /// Double bonds that could carry E/Z but have none declared. Written as
    /// "either" so a reader does not take the drawn geometry as stereo.
    pub unspecified_double_bonds: Vec<BondIdx>,
}

/// `mol` without stereo markers on bonds: `Up`/`Down` become single bonds
/// and no bond directions are kept; atoms and stereo neighbour orders are
/// unchanged. Bond indices are preserved.
fn without_bond_marks(mol: &Molecule, wedges: &HashMap<BondIdx, Wedge>) -> Molecule {
    let mut b = MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        b.add_atom(atom.clone());
    }
    for (idx, bond) in mol.bonds() {
        let (a1, a2, order) = match wedges.get(&idx) {
            Some(w) => {
                let other = if bond.atom1 == w.start {
                    bond.atom2
                } else {
                    bond.atom1
                };
                (w.start, other, w.order)
            }
            None => match bond.order {
                BondOrder::Up | BondOrder::Down => (bond.atom1, bond.atom2, BondOrder::Single),
                order => (bond.atom1, bond.atom2, order),
            },
        };
        let _ = b.add_bond(a1, a2, order);
    }
    b.copy_stereo_from(mol);
    b.copy_stereo_groups_from(mol);
    b.build()
}

/// E/Z labels of `mol`'s stereo double bonds, keyed by bond.
fn ez_labels(mol: &Molecule) -> HashMap<BondIdx, chematic_core::atom::CipCode> {
    chematic_chem::assign_ez_bonds(mol).into_iter().collect()
}

/// The substituent of `end` (double-bonded to `other`) that carries a
/// `/`/`\\` marker in `mol`, and whether the marker puts it above `end`.
/// A stored direction (read from its anchor, atom1 when none) wins over a
/// literal `Up`/`Down` order (read from atom1), as in the SMILES writer.
fn marked_substituent(mol: &Molecule, end: AtomIdx, other: AtomIdx) -> Option<(AtomIdx, bool)> {
    mol.neighbors(end).find_map(|(sub, bi)| {
        if sub == other {
            return None;
        }
        let bond = mol.bond(bi);
        let (dir, anchor) = match (mol.bond_direction(bi), bond.order) {
            (Some(dir), _) => (dir, mol.bond_direction_anchor(bi).unwrap_or(bond.atom1)),
            (None, BondOrder::Up | BondOrder::Down) => (bond.order, bond.atom1),
            _ => return None,
        };
        // `Up` read from `anchor`: the far atom is above the anchor.
        let far_above = dir == BondOrder::Up;
        Some((sub, if anchor == end { far_above } else { !far_above }))
    })
}

/// Whether `coords` draw the stereo double bond `bond` with the geometry
/// its markers in `mol` declare; `None` when there are no markers on both
/// ends or a marked substituent lies on the double-bond axis.
fn ez_drawn_as_declared(mol: &Molecule, coords: &[(f64, f64)], bond: BondIdx) -> Option<bool> {
    let (p, q) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    let (s1, up1) = marked_substituent(mol, p, q)?;
    let (s2, up2) = marked_substituent(mol, q, p)?;
    let (pp, qq) = (coords[p.0 as usize], coords[q.0 as usize]);
    let axis = (qq.0 - pp.0, qq.1 - pp.1);
    let side = |end: AtomIdx, sub: AtomIdx| {
        let (e, s) = (coords[end.0 as usize], coords[sub.0 as usize]);
        let cross = axis.0 * (s.1 - e.1) - axis.1 * (s.0 - e.0);
        (cross.abs() > 1e-6).then_some(cross > 0.0)
    };
    let (g1, g2) = (side(p, s1)?, side(q, s2)?);
    Some((g1 == g2) == (up1 == up2))
}

/// Double bonds outside small rings whose ends each have one or two
/// distinct substituents (an implicit H or a lone pair making up the
/// rest), with no E/Z in `declared`.
fn unspecified_double_bonds(
    mol: &Molecule,
    declared: &HashMap<BondIdx, chematic_core::atom::CipCode>,
) -> Vec<BondIdx> {
    let classes = chematic_smiles::topological_equivalence_classes(mol);
    let end_ok = |end: AtomIdx, other: AtomIdx, db: BondIdx| {
        let subs: Vec<AtomIdx> = mol
            .neighbors(end)
            .filter(|&(nb, _)| nb != other)
            .map(|(nb, _)| nb)
            .collect();
        let cumulated = mol
            .neighbors(end)
            .any(|(_, bi)| bi != db && mol.bond(bi).order == BondOrder::Double);
        let h = chematic_core::implicit_hcount(mol, end) as usize;
        !cumulated
            && !subs.is_empty()
            && subs.len() + h <= 2
            && (subs.len() < 2 || classes[subs[0].0 as usize] != classes[subs[1].0 as usize])
    };
    let mut out: Vec<BondIdx> = mol
        .bonds()
        .filter(|&(bi, b)| {
            b.order == BondOrder::Double
                && !declared.contains_key(&bi)
                && (!mol.atom(b.atom1).aromatic || !mol.atom(b.atom2).aromatic)
                && !in_ring_smaller_than(mol, bi, 8)
                && end_ok(b.atom1, b.atom2, bi)
                && end_ok(b.atom2, b.atom1, bi)
        })
        .map(|(bi, _)| bi)
        .collect();
    out.sort_by_key(|b| b.0);
    out
}

/// Whether `bond` lies on a ring of fewer than `size` atoms.
fn in_ring_smaller_than(mol: &Molecule, bond: BondIdx, size: usize) -> bool {
    let (start, goal) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    let mut seen = vec![false; mol.atom_count()];
    seen[start.0 as usize] = true;
    let mut frontier = vec![start];
    for _ in 0..size.saturating_sub(2) {
        let mut next = Vec::new();
        for &a in &frontier {
            for (nb, bi) in mol.neighbors(a) {
                if bi == bond {
                    continue;
                }
                if nb == goal {
                    return true;
                }
                if !seen[nb.0 as usize] {
                    seen[nb.0 as usize] = true;
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    false
}

/// Atoms on `end`'s side of the double bond `end`=`other`, or `None` when
/// the bond is in a ring (both sides connected).
fn side(mol: &Molecule, end: AtomIdx, other: AtomIdx) -> Option<Vec<AtomIdx>> {
    let mut seen = vec![false; mol.atom_count()];
    seen[end.0 as usize] = true;
    let mut stack = vec![end];
    let mut out = vec![end];
    while let Some(a) = stack.pop() {
        for (nb, _) in mol.neighbors(a) {
            if a == end && nb == other {
                continue;
            }
            if nb == other {
                return None;
            }
            if !seen[nb.0 as usize] {
                seen[nb.0 as usize] = true;
                out.push(nb);
                stack.push(nb);
            }
        }
    }
    Some(out)
}

fn reflect(coords: &mut [(f64, f64)], atoms: &[AtomIdx], p: (f64, f64), q: (f64, f64)) {
    let (dx, dy) = (q.0 - p.0, q.1 - p.1);
    let len2 = dx * dx + dy * dy;
    if len2 == 0.0 {
        return;
    }
    for a in atoms {
        let (x, y) = coords[a.0 as usize];
        let t = ((x - p.0) * dx + (y - p.1) * dy) / len2;
        let (fx, fy) = (p.0 + t * dx, p.1 + t * dy);
        coords[a.0 as usize] = (2.0 * fx - x, 2.0 * fy - y);
    }
}

/// Whether the wedge set gives `centre` the configuration `mol` declares.
fn wedge_matches(
    mol: &Molecule,
    coords: &[(f64, f64)],
    wedges: &HashMap<BondIdx, Wedge>,
    centre: AtomIdx,
) -> Option<bool> {
    // Neighbours drawn in (nearly) the same direction leave the parity to
    // each reader's tie-breaking; do not rely on such a drawing.
    let c = coords[centre.0 as usize];
    let angles: Vec<f64> = mol
        .neighbors(centre)
        .map(|(nb, _)| {
            let p = coords[nb.0 as usize];
            (p.1 - c.1).atan2(p.0 - c.0)
        })
        .collect();
    for (i, a) in angles.iter().enumerate() {
        for b in &angles[i + 1..] {
            let d = (a - b).rem_euclid(std::f64::consts::TAU);
            if d.min(std::f64::consts::TAU - d) < 15f64.to_radians() {
                return None;
            }
        }
    }
    // Four bonds with the three plain ones inside a half-plane and the
    // wedge between two of them is read differently by different toolkits.
    let wedged = wedges.iter().find_map(|(&b, w)| {
        (w.start == centre).then(|| {
            let bond = mol.bond(b);
            if bond.atom1 == centre {
                bond.atom2
            } else {
                bond.atom1
            }
        })
    })?;
    if angles.len() == 4 {
        let angle_of = |a: AtomIdx| {
            let p = coords[a.0 as usize];
            (p.1 - c.1)
                .atan2(p.0 - c.0)
                .rem_euclid(std::f64::consts::TAU)
        };
        let w = angle_of(wedged);
        let mut plain: Vec<f64> = mol
            .neighbors(centre)
            .filter(|&(nb, _)| nb != wedged)
            .map(|(nb, _)| angle_of(nb))
            .collect();
        plain.sort_by(f64::total_cmp);
        for i in 0..3 {
            let (lo, hi) = if i < 2 {
                (plain[i], plain[i + 1])
            } else {
                (plain[2], plain[0] + std::f64::consts::TAU)
            };
            let in_gap = (lo..hi).contains(&w) || (lo..hi).contains(&(w + std::f64::consts::TAU));
            if hi - lo > std::f64::consts::PI && !in_gap {
                return None;
            }
        }
    }
    let drawn = without_bond_marks(mol, wedges);
    let (chirality, order) =
        chematic_perception::stereo2d_local::local_parity_from_wedges(&drawn, coords, centre)?;
    let declared = mol.stereo_neighbor_order(centre)?;
    let want = mol.atom(centre).chirality;
    if !want.is_tetrahedral() || !chirality.is_tetrahedral() || declared.len() != order.len() {
        return None;
    }
    let odd = match (
        <[u32; 4]>::try_from(declared),
        <[u32; 4]>::try_from(order.as_slice()),
    ) {
        (Ok(a), Ok(b)) => chematic_core::remap_tetrahedral_parity(a, b).ok()?,
        _ => permutation_is_odd(declared, &order)?,
    };
    Some((chirality == want) != odd)
}

fn permutation_is_odd(a: &[u32], b: &[u32]) -> Option<bool> {
    let mut perm: Vec<usize> = b
        .iter()
        .map(|x| a.iter().position(|y| y == x))
        .collect::<Option<_>>()?;
    let mut swaps = 0;
    for i in 0..perm.len() {
        while perm[i] != i {
            let j = perm[i];
            if j >= perm.len() {
                return None;
            }
            perm.swap(i, j);
            swaps += 1;
        }
    }
    Some(swaps % 2 == 1)
}

/// Lay `mol` out and choose wedges for its stereo; see the module docs.
/// `coords`, when non-empty and not all zero, are used instead of a fresh
/// layout (and are not modified).
pub fn stereo_depiction(mol: &Molecule, coords: &[(f64, f64)]) -> StereoDepiction {
    let n = mol.atom_count();
    let given = coords.len() >= n && coords.iter().any(|&(x, y)| x != 0.0 || y != 0.0);
    let layout = if given {
        coords[..n].to_vec()
    } else {
        // Depiction layouts are in pixels with y down; MOL wants Å with y up
        // (RDKit's 1.5 Å bond length).
        let scale = 1.5 / chematic_depict::layout::BOND_LEN;
        chematic_depict::compute_layout(mol)
            .coords
            .iter()
            .map(|p| (p.x * scale, -p.y * scale))
            .collect()
    };
    let mut out = StereoDepiction {
        coords: layout,
        ..Default::default()
    };

    // E/Z: reflect one side of each double bond whose geometry disagrees.
    // A reflection across bond B's axis mirrors every double bond on the
    // moved side as a whole (an atom bonded to B lies on the mirror line),
    // so it never changes a bond already set.
    let target = ez_labels(mol);
    if !target.is_empty() {
        let mut bonds: Vec<BondIdx> = target.keys().copied().collect();
        bonds.sort_by_key(|b| b.0);
        for &bond in &bonds {
            if given || ez_drawn_as_declared(mol, &out.coords, bond) != Some(false) {
                continue;
            }
            let e = mol.bond(bond);
            let sides = (side(mol, e.atom1, e.atom2), side(mol, e.atom2, e.atom1));
            if let (Some(a), Some(b)) = sides {
                let atoms = if a.len() <= b.len() { a } else { b };
                let (p, q) = (
                    out.coords[e.atom1.0 as usize],
                    out.coords[e.atom2.0 as usize],
                );
                reflect(&mut out.coords, &atoms, p, q);
            }
        }
        out.unexpressed_double_bonds = bonds
            .into_iter()
            .filter(|&b| ez_drawn_as_declared(mol, &out.coords, b) != Some(true))
            .collect();
    }

    out.unspecified_double_bonds = unspecified_double_bonds(mol, &target);

    // Tetrahedral centres: one wedge or hash each, checked by re-reading.
    let centres: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(i, a)| a.chirality.is_tetrahedral() && mol.stereo_neighbor_order(*i).is_some())
        .map(|(i, _)| i)
        .collect();
    if !given {
        place_terminal_neighbours(mol, &mut out.coords, &centres);
    }
    (out.wedges, out.unexpressed_centres) = place_wedges(mol, &out.coords, &centres);

    // The layout can put bridge atoms on top of others, where no wedge
    // gives a parity. Pull coincident atoms apart and try again.
    if !given && !out.unexpressed_centres.is_empty() {
        let mut nudged = out.coords.clone();
        if separate_coincident_atoms(mol, &mut nudged) {
            let (wedges, unexpressed) = place_wedges(mol, &nudged, &centres);
            let ez_still_ok = target
                .keys()
                .filter(|b| !out.unexpressed_double_bonds.contains(b))
                .all(|&b| ez_drawn_as_declared(mol, &nudged, b) == Some(true));
            if unexpressed.len() < out.unexpressed_centres.len() && ez_still_ok {
                out.coords = nudged;
                out.wedges = wedges;
                out.unexpressed_centres = unexpressed;
            }
        }
    }
    out
}

/// Put each one-bond neighbour of a stereocentre (a methyl, OH, halogen,
/// explicit H) in the middle of the widest gap between the centre's other
/// bonds, so a wedge to it reads unambiguously.
fn place_terminal_neighbours(mol: &Molecule, coords: &mut [(f64, f64)], centres: &[AtomIdx]) {
    use std::f64::consts::TAU;
    for &centre in centres {
        let c = coords[centre.0 as usize];
        let terminal: Vec<AtomIdx> = mol
            .neighbors(centre)
            .map(|(nb, _)| nb)
            .filter(|&nb| mol.degree(nb) == 1)
            .collect();
        for t in terminal {
            let mut angles: Vec<f64> = mol
                .neighbors(centre)
                .filter(|&(nb, _)| nb != t)
                .map(|(nb, _)| {
                    let p = coords[nb.0 as usize];
                    (p.1 - c.1).atan2(p.0 - c.0).rem_euclid(TAU)
                })
                .collect();
            if angles.is_empty() {
                continue;
            }
            angles.sort_by(f64::total_cmp);
            let (mut best, mut gap) = (angles[0] + TAU / 2.0, 0.0);
            for (i, &a) in angles.iter().enumerate() {
                let next = if i + 1 < angles.len() {
                    angles[i + 1]
                } else {
                    angles[0] + TAU
                };
                if next - a > gap {
                    gap = next - a;
                    best = a + gap / 2.0;
                }
            }
            coords[t.0 as usize] = (c.0 + 1.5 * best.cos(), c.1 + 1.5 * best.sin());
        }
    }
}

/// Move every atom that sits on an earlier atom (closer than 0.05 Å)
/// half a bond length towards the centroid of its neighbours, or off to
/// one side when that centroid is the same point. Returns whether any atom
/// moved.
fn separate_coincident_atoms(mol: &Molecule, coords: &mut [(f64, f64)]) -> bool {
    let mut moved = false;
    for i in 0..coords.len() {
        let clash = (0..i).any(|j| {
            let (dx, dy) = (coords[i].0 - coords[j].0, coords[i].1 - coords[j].1);
            dx * dx + dy * dy < 0.0025
        });
        if !clash {
            continue;
        }
        let nbs: Vec<(f64, f64)> = mol
            .neighbors(AtomIdx(i as u32))
            .map(|(nb, _)| coords[nb.0 as usize])
            .collect();
        let (cx, cy) = if nbs.is_empty() {
            coords[i]
        } else {
            let n = nbs.len() as f64;
            (
                nbs.iter().map(|p| p.0).sum::<f64>() / n,
                nbs.iter().map(|p| p.1).sum::<f64>() / n,
            )
        };
        let (dx, dy) = (cx - coords[i].0, cy - coords[i].1);
        let len = (dx * dx + dy * dy).sqrt();
        let (ux, uy) = if len > 1e-6 {
            (dx / len, dy / len)
        } else {
            (0.6, 0.8)
        };
        coords[i] = (coords[i].0 + 0.75 * ux, coords[i].1 + 0.75 * uy);
        moved = true;
    }
    moved
}

/// One wedge or hash per centre, each verified by re-reading; the centres
/// no single wedge could express.
fn place_wedges(
    mol: &Molecule,
    coords: &[(f64, f64)],
    centres: &[AtomIdx],
) -> (HashMap<BondIdx, Wedge>, Vec<AtomIdx>) {
    let mut wedges: HashMap<BondIdx, Wedge> = HashMap::new();
    let mut unexpressed = Vec::new();
    let is_centre = |a: AtomIdx| centres.contains(&a);
    let in_ring = ring_bonds(mol);
    for &centre in centres {
        let mut candidates: Vec<(u32, AtomIdx, BondIdx)> = mol
            .neighbors(centre)
            .filter(|&(_, b)| {
                matches!(
                    mol.bond(b).order,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down
                ) && !wedges.contains_key(&b)
            })
            .map(|(nb, b)| {
                // Prefer a terminal, non-ring, non-stereocentre neighbour
                // that is not on a double bond (a wedge there would sit next
                // to an E/Z marker).
                let mut rank = 0;
                if mol
                    .neighbors(nb)
                    .any(|(_, nbb)| mol.bond(nbb).order == BondOrder::Double)
                {
                    rank += 8;
                }
                if is_centre(nb) {
                    rank += 4;
                }
                if in_ring[b.0 as usize] {
                    rank += 2;
                }
                if mol.degree(nb) > 1 {
                    rank += 1;
                }
                (rank, nb, b)
            })
            .collect();
        candidates.sort_by_key(|&(rank, nb, _)| (rank, nb.0));
        let mut placed = false;
        for (_, _, bond) in candidates {
            for order in [BondOrder::Up, BondOrder::Down] {
                wedges.insert(
                    bond,
                    Wedge {
                        start: centre,
                        order,
                    },
                );
                if wedge_matches(mol, coords, &wedges, centre) == Some(true) {
                    placed = true;
                    break;
                }
                wedges.remove(&bond);
            }
            if placed {
                break;
            }
        }
        if !placed {
            unexpressed.push(centre);
        }
    }
    (wedges, unexpressed)
}

/// Per bond, whether it lies in a ring (is not a bridge), by an iterative
/// Tarjan low-link search.
fn ring_bonds(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut order = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut in_ring = vec![true; mol.bond_count()];
    let mut counter = 0;
    for root in 0..n {
        if order[root] != usize::MAX {
            continue;
        }
        // (atom, parent bond, neighbour cursor)
        let mut stack: Vec<(usize, Option<BondIdx>, usize)> = vec![(root, None, 0)];
        order[root] = counter;
        low[root] = counter;
        counter += 1;
        while let Some(&(v, parent, cursor)) = stack.last() {
            let nbs: Vec<(AtomIdx, BondIdx)> = mol.neighbors(AtomIdx(v as u32)).collect();
            if cursor < nbs.len() {
                if let Some(top) = stack.last_mut() {
                    top.2 += 1;
                }
                let (w, b) = nbs[cursor];
                if Some(b) == parent {
                    continue;
                }
                let w = w.0 as usize;
                if order[w] == usize::MAX {
                    order[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push((w, Some(b), 0));
                } else {
                    low[v] = low[v].min(order[w]);
                }
            } else {
                stack.pop();
                if let Some(&(u, _, _)) = stack.last() {
                    low[u] = low[u].min(low[v]);
                    if low[v] > order[u]
                        && let Some(b) = parent
                    {
                        in_ring[b.0 as usize] = false;
                    }
                }
            }
        }
    }
    in_ring
}

/// Whether writing `mol` needs [`stereo_depiction`]: it has a tetrahedral
/// centre, or SMILES-style `/`/`\\` markers next to a double bond that only
/// coordinates can express. A wedge/hash on a molecule with neither (e.g. a
/// MOL record whose wedges assign no centre) is written as read.
pub(crate) fn needs_stereo_depiction(mol: &Molecule) -> bool {
    let next_to_double = |a: AtomIdx| {
        mol.neighbors(a)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double)
    };
    mol.atoms().any(|(_, a)| a.chirality.is_tetrahedral())
        || mol.bonds().any(|(i, b)| {
            mol.bond_direction(i).is_some()
                || (matches!(b.order, BondOrder::Up | BondOrder::Down)
                    && (next_to_double(b.atom1) || next_to_double(b.atom2)))
        })
}

#[cfg(test)]
mod tests {
    use crate::mol2000::{MolMetadata, read_mol_with_diagnostics, write_mol};
    use crate::mol3000::{read_mol_v3000_with_diagnostics, write_mol_v3000};
    use chematic_smiles::{canonical_smiles, parse};

    const CASES: &[&str] = &[
        "N[C@@H](C)C(=O)O",
        "N[C@H](C)C(=O)O",
        "F/C=C/F",
        "F/C=C\\F",
        "C/C=C/[C@H](O)CC",
        "C[C@H]1CC[C@@H](O)CC1",
        "C[C@@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
        "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "C/C(F)=C(\\Cl)Br",
        "CC(=O)O[C@@H]1C[C@H](N)C1",
    ];

    #[test]
    fn mol_blocks_without_coords_keep_tetrahedral_and_ez_stereo() {
        for smi in CASES {
            let mol = parse(smi).unwrap();
            let want = canonical_smiles(&mol);
            let v2 = write_mol(&mol, &MolMetadata::default());
            let got2 = canonical_smiles(&read_mol_with_diagnostics(&v2).unwrap().mol);
            assert_eq!(got2, want, "V2000 {smi}\n{v2}");
            let v3 = write_mol_v3000(&mol, &MolMetadata::default(), &[]);
            let got3 = canonical_smiles(&read_mol_v3000_with_diagnostics(&v3).unwrap().mol);
            assert_eq!(got3, want, "V3000 {smi}\n{v3}");
        }
    }

    #[test]
    fn conjugated_and_branched_ez_round_trip() {
        // Methyl-branched polyene (retinoid) and an amidine whose ends are
        // conjugated: the reader used to reject both.
        for smi in [
            "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C(=O)O)C(C)(C)CCC1",
            "CN/C(=N/C1=CC(C)(C)Oc2ccc(C#N)cc21)NC#N",
            "C/C=C/C=C\\C=C\\C",
        ] {
            let mol = parse(smi).unwrap();
            let v2 = write_mol(&mol, &MolMetadata::default());
            let back = read_mol_with_diagnostics(&v2).unwrap().mol;
            assert_eq!(
                canonical_smiles(&back),
                canonical_smiles(&mol),
                "{smi}\n{v2}"
            );
        }
    }

    #[test]
    fn ring_double_bond_stereo_is_written_as_either_not_guessed() {
        // A macrocyclic E/Z bond cannot be set by reflecting one side; it is
        // written "either" (bond stereo 3) instead of with a wrong geometry.
        let mol = parse("C1CCCC/C=C/CCCCC1").unwrap();
        let d = super::stereo_depiction(&mol, &[]);
        let v2 = write_mol(&mol, &MolMetadata::default());
        for &b in &d.unexpressed_double_bonds {
            let line = v2.lines().nth(4 + mol.atom_count() + b.0 as usize).unwrap();
            assert_eq!(line[9..12].trim(), "3", "{v2}");
        }
        let back = canonical_smiles(&read_mol_with_diagnostics(&v2).unwrap().mol);
        assert!(
            back == canonical_smiles(&mol) || !back.contains('/'),
            "{back}"
        );
    }

    #[test]
    fn unspecified_stereo_double_bond_is_written_as_either() {
        // No E/Z declared: the drawn geometry must not become stereo.
        let mol = parse("CC=CC[C@H](N)C(=O)O").unwrap();
        let v2 = write_mol(&mol, &MolMetadata::default());
        let back = read_mol_with_diagnostics(&v2).unwrap().mol;
        assert_eq!(canonical_smiles(&back), canonical_smiles(&mol), "{v2}");
    }

    #[test]
    fn one_wedge_per_centre_and_no_wedge_codes_for_ez() {
        let mol = parse("F/C=C/F").unwrap();
        let v2 = write_mol(&mol, &MolMetadata::default());
        let bond_lines: Vec<&str> = v2
            .lines()
            .skip(4 + mol.atom_count())
            .take(mol.bond_count())
            .collect();
        assert!(bond_lines.iter().all(|l| l[9..12].trim() == "0"), "{v2}");

        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let v2 = write_mol(&mol, &MolMetadata::default());
        let bond_lines: Vec<&str> = v2
            .lines()
            .skip(4 + mol.atom_count())
            .take(mol.bond_count())
            .collect();
        let wedged: Vec<&&str> = bond_lines
            .iter()
            .filter(|l| l[9..12].trim() != "0")
            .collect();
        assert_eq!(wedged.len(), 1, "{v2}");
        // drawn from the stereocentre (atom 2, 1-based)
        assert_eq!(wedged[0][0..3].trim(), "2", "{v2}");
    }
}
