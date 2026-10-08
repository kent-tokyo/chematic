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
    if wedges.is_empty()
        && let Some(mut copy) = mol.atoms_bonds_and_stereo_copy()
    {
        for (idx, bond) in mol.bonds() {
            if matches!(bond.order, BondOrder::Up | BondOrder::Down) {
                copy.set_bond_order(idx, BondOrder::Single);
            }
        }
        return copy;
    }
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
    // The substituents of `end` when it can carry E/Z, before the check
    // that two substituents differ (classes are computed only if needed).
    let end_subs = |end: AtomIdx, other: AtomIdx, db: BondIdx| {
        let subs: Vec<AtomIdx> = mol
            .neighbors(end)
            .filter(|&(nb, _)| nb != other)
            .map(|(nb, _)| nb)
            .collect();
        let cumulated = mol
            .neighbors(end)
            .any(|(_, bi)| bi != db && mol.bond(bi).order == BondOrder::Double);
        let h = chematic_core::implicit_hcount(mol, end) as usize;
        (!cumulated && !subs.is_empty() && subs.len() + h <= 2).then_some(subs)
    };
    let classes = std::cell::OnceCell::new();
    let distinct = |subs: &[AtomIdx]| {
        subs.len() < 2 || {
            let classes =
                classes.get_or_init(|| chematic_smiles::topological_equivalence_classes(mol));
            classes[subs[0].0 as usize] != classes[subs[1].0 as usize]
        }
    };
    let mut out: Vec<BondIdx> = mol
        .bonds()
        .filter(|&(bi, b)| {
            if b.order != BondOrder::Double
                || declared.contains_key(&bi)
                || (mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic)
            {
                return false;
            }
            let (Some(s1), Some(s2)) = (
                end_subs(b.atom1, b.atom2, bi),
                end_subs(b.atom2, b.atom1, bi),
            ) else {
                return false;
            };
            !in_ring_smaller_than(mol, bi, 8) && distinct(&s1) && distinct(&s2)
        })
        .map(|(bi, _)| bi)
        .collect();
    out.sort_by_key(|b| b.0);
    out
}

/// Double bonds of a molecule without declared stereo that a reader would
/// take as E/Z from drawn coordinates. Writers mark them "either" (V2000
/// stereo 3 / V3000 `CFG=2`) whenever they write real coordinates, so a
/// layout never invents a configuration the input did not have.
pub(crate) fn undeclared_stereo_double_bonds(mol: &Molecule) -> Vec<BondIdx> {
    unspecified_double_bonds(mol, &HashMap::new())
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

/// Ring double bonds (macrocycles: rings below eight atoms carry no E/Z)
/// cannot be set by reflecting a side, and a ring drawn as a regular polygon
/// makes every ring double bond cis. Move one end of a bond that needs to be
/// trans, with its own substituents, across the line from its other ring
/// neighbour to the far end of the double bond (both bond lengths kept),
/// when that sets the bond, no bond already drawn as declared changes, and
/// nothing lands on another atom.
fn flip_ring_double_bonds(mol: &Molecule, coords: &mut [(f64, f64)], bonds: &[BondIdx]) {
    let n = mol.atom_count();
    for &bond in bonds {
        if ez_drawn_as_declared(mol, coords, bond) != Some(false) {
            continue;
        }
        let ok_before: Vec<BondIdx> = bonds
            .iter()
            .copied()
            .filter(|&b| ez_drawn_as_declared(mol, coords, b) == Some(true))
            .collect();
        let e = mol.bond(bond);
        for (end, far) in [(e.atom1, e.atom2), (e.atom2, e.atom1)] {
            // `end`'s ring neighbour: the one still connected to `far` when
            // `end` is removed.
            let ring_nbs: Vec<AtomIdx> = mol
                .neighbors(end)
                .map(|(nb, _)| nb)
                .filter(|&nb| nb != far && connected_without(mol, nb, far, end))
                .collect();
            let [pivot] = ring_nbs[..] else {
                continue;
            };
            // `end` and its own substituents.
            let mut group = vec![end];
            let mut seen = vec![false; n];
            for a in [end, pivot, far] {
                seen[a.0 as usize] = true;
            }
            let mut stack = vec![end];
            while let Some(a) = stack.pop() {
                for (nb, _) in mol.neighbors(a) {
                    if !seen[nb.0 as usize] {
                        seen[nb.0 as usize] = true;
                        group.push(nb);
                        stack.push(nb);
                    }
                }
            }
            if group.len() * 2 > n {
                continue;
            }
            let mut moved = coords.to_vec();
            let (p, q) = (coords[pivot.0 as usize], coords[far.0 as usize]);
            reflect(&mut moved, &group, p, q);
            let in_group = |i: usize| group.iter().any(|g| g.0 as usize == i);
            let clear = group.iter().all(|g| {
                let (gx, gy) = moved[g.0 as usize];
                (0..n).all(|i| {
                    let (dx, dy) = (moved[i].0 - gx, moved[i].1 - gy);
                    in_group(i) || dx * dx + dy * dy > 0.25
                })
            });
            if clear
                && ez_drawn_as_declared(mol, &moved, bond) == Some(true)
                && ok_before
                    .iter()
                    .all(|&b| ez_drawn_as_declared(mol, &moved, b) == Some(true))
            {
                coords.copy_from_slice(&moved);
                break;
            }
        }
    }
}

/// Whether `from` reaches `to` without passing through `blocked`.
fn connected_without(mol: &Molecule, from: AtomIdx, to: AtomIdx, blocked: AtomIdx) -> bool {
    let mut seen = vec![false; mol.atom_count()];
    seen[blocked.0 as usize] = true;
    seen[from.0 as usize] = true;
    let mut stack = vec![from];
    while let Some(a) = stack.pop() {
        if a == to {
            return true;
        }
        for (nb, _) in mol.neighbors(a) {
            if !seen[nb.0 as usize] {
                seen[nb.0 as usize] = true;
                stack.push(nb);
            }
        }
    }
    false
}

/// Non-bonded pairs closer than 0.45 bonds (0.675 Å) between `moved` and
/// the other atoms.
fn moved_clashes(mol: &Molecule, coords: &[(f64, f64)], moved: &[AtomIdx]) -> usize {
    let mut in_moved = vec![false; coords.len()];
    for a in moved {
        in_moved[a.0 as usize] = true;
    }
    let limit = 0.675 * 0.675;
    let mut count = 0;
    for &a in moved {
        let (x, y) = coords[a.0 as usize];
        for (j, &(u, v)) in coords.iter().enumerate() {
            if in_moved[j] {
                continue;
            }
            let d2 = (x - u) * (x - u) + (y - v) * (y - v);
            if d2 < limit && mol.bond_between(a, AtomIdx(j as u32)).is_none() {
                count += 1;
            }
        }
    }
    count
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

/// Whether the wedges in `drawn` (a copy of `mol` without its own bond
/// marks, see [`place_wedges`]) give `centre` the configuration `mol`
/// declares; `wedged` is the neighbour at the end of `centre`'s wedge.
fn wedge_matches(
    mol: &Molecule,
    drawn: &Molecule,
    coords: &[(f64, f64)],
    centre: AtomIdx,
    wedged: AtomIdx,
) -> Option<bool> {
    wedge_matches_with(mol, drawn, coords, centre, wedged, NEIGHBOUR_SEPARATION_DEG)
}

/// Neighbours of a wedged centre drawn closer than this (degrees) leave the
/// parity to each reader's tie-breaking.
const NEIGHBOUR_SEPARATION_DEG: f64 = 15.0;
/// The separation accepted when no wedge passes [`NEIGHBOUR_SEPARATION_DEG`]
/// (a bridgehead squeezed by the layout): RDKit reads such a drawing as
/// declared, while without a wedge every reader loses the centre.
const FALLBACK_NEIGHBOUR_SEPARATION_DEG: f64 = 5.0;

fn wedge_matches_with(
    mol: &Molecule,
    drawn: &Molecule,
    coords: &[(f64, f64)],
    centre: AtomIdx,
    wedged: AtomIdx,
    min_separation_deg: f64,
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
            if d.min(std::f64::consts::TAU - d) < min_separation_deg.to_radians() {
                return None;
            }
        }
    }
    // Four bonds with the three plain ones inside a half-plane and the
    // wedge between two of them is read differently by different toolkits.
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
    let (chirality, order) =
        chematic_perception::stereo2d_local::local_parity_from_wedges(drawn, coords, centre)?;
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
    if given {
        return depict_on(mol, coords[..n].to_vec(), true);
    }
    let first = depict_on(mol, layout_angstrom(mol, false), false);
    let lost = |d: &StereoDepiction| d.unexpressed_centres.len() + d.unexpressed_double_bonds.len();
    if first.unexpressed_centres.is_empty() && first.unexpressed_double_bonds.is_empty() {
        return first;
    }
    // The layout of fused and bridged systems depends on atom order, and a
    // stretched ring closure or two bonds in one direction can leave a
    // centre undrawable, or a ring double bond that cannot be flipped to its
    // declared geometry. Lay the atoms out in reverse order as well and
    // keep whichever drawing loses less.
    let second = depict_on(mol, layout_angstrom(mol, true), false);
    if lost(&second) < lost(&first) {
        second
    } else {
        first
    }
}

/// `compute_layout` in Å with y up (depiction layouts are in pixels with y
/// down; RDKit's 1.5 Å bond length), optionally laid out with the atoms in
/// reverse order and mapped back.
fn layout_angstrom(mol: &Molecule, reversed: bool) -> Vec<(f64, f64)> {
    let scale = 1.5 / chematic_depict::layout::BOND_LEN;
    let to_angstrom = |p: &chematic_depict::layout::Point| (p.x * scale, -p.y * scale);
    if !reversed {
        return chematic_depict::compute_layout(mol)
            .coords
            .iter()
            .map(to_angstrom)
            .collect();
    }
    let n = mol.atom_count();
    let mut b = MoleculeBuilder::new();
    for i in (0..n).rev() {
        b.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }
    let flip = |a: AtomIdx| AtomIdx((n - 1 - a.0 as usize) as u32);
    let mut bonds: Vec<_> = mol.bonds().map(|(_, bond)| bond.clone()).collect();
    bonds.reverse();
    for bond in bonds {
        let order = match bond.order {
            BondOrder::Up | BondOrder::Down => BondOrder::Single,
            order => order,
        };
        let _ = b.add_bond(flip(bond.atom1), flip(bond.atom2), order);
    }
    let layout = chematic_depict::compute_layout(&b.build());
    (0..n)
        .map(|i| to_angstrom(&layout.coords[n - 1 - i]))
        .collect()
}

/// [`chematic_depict::relieve_layout_clashes`] on Å coordinates (y up).
fn relieve_after_reflection(mol: &Molecule, coords: &mut [(f64, f64)]) {
    let scale = 1.5 / chematic_depict::layout::BOND_LEN;
    let mut layout = chematic_depict::Layout {
        coords: coords
            .iter()
            .map(|&(x, y)| chematic_depict::Point::new(x / scale, -y / scale))
            .collect(),
    };
    chematic_depict::relieve_layout_clashes(mol, &mut layout);
    for (c, p) in coords.iter_mut().zip(&layout.coords) {
        *c = (p.x * scale, -p.y * scale);
    }
}

/// The E/Z reflection, terminal-neighbour placement and wedge choice on
/// `layout`; `given` coordinates are kept as they are.
fn depict_on(mol: &Molecule, layout: Vec<(f64, f64)>, given: bool) -> StereoDepiction {
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
                let (p, q) = (
                    out.coords[e.atom1.0 as usize],
                    out.coords[e.atom2.0 as usize],
                );
                // Either side gives the declared geometry; reflect the one
                // that lands on fewer atoms (two aryl rings of a
                // tetrasubstituted alkene end up side by side otherwise),
                // the smaller one on a tie.
                let mut with_a = out.coords.clone();
                reflect(&mut with_a, &a, p, q);
                let mut with_b = out.coords.clone();
                reflect(&mut with_b, &b, p, q);
                let (ca, cb) = (
                    moved_clashes(mol, &with_a, &a),
                    moved_clashes(mol, &with_b, &b),
                );
                out.coords = if ca < cb || (ca == cb && a.len() <= b.len()) {
                    with_a
                } else {
                    with_b
                };
            }
        }
        if !given {
            flip_ring_double_bonds(mol, &mut out.coords, &bonds);
            // A reflected side can land on the rest of the drawing; clear
            // such clashes again (moves that keep every E/Z geometry).
            relieve_after_reflection(mol, &mut out.coords);
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
    let mut ez_bonds: Vec<BondIdx> = target.keys().copied().collect();
    ez_bonds.sort_by_key(|b| b.0);
    (out.wedges, out.unexpressed_centres) = place_wedges(mol, &out.coords, &centres, &ez_bonds);

    // The layout can put bridge atoms on top of others, where no wedge
    // gives a parity. Pull coincident atoms apart and try again.
    if !given && !out.unexpressed_centres.is_empty() {
        let mut nudged = out.coords.clone();
        if separate_coincident_atoms(mol, &mut nudged) {
            let (wedges, unexpressed) = place_wedges(mol, &nudged, &centres, &ez_bonds);
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

    // A bridge drawn as a regular polygon can leave two bonds of a centre
    // in nearly one direction. Move a two-bond neighbour (a bridge atom)
    // into the widest gap between the centre's other bonds, stretching its
    // far bond, when that lets more centres be drawn.
    let mut attempts = 0;
    while !given && !out.unexpressed_centres.is_empty() && attempts < 8 {
        let mut improved = false;
        for &centre in &out.unexpressed_centres.clone() {
            let movable: Vec<AtomIdx> = mol
                .neighbors(centre)
                .map(|(nb, _)| nb)
                .filter(|&nb| mol.degree(nb) == 2 && !centres.contains(&nb))
                .collect();
            for nb in movable {
                attempts += 1;
                let Some(dir) = widest_gap_direction(mol, &out.coords, centre, nb) else {
                    continue;
                };
                let mut moved = out.coords.clone();
                let c = moved[centre.0 as usize];
                moved[nb.0 as usize] = (c.0 + 1.5 * dir.cos(), c.1 + 1.5 * dir.sin());
                let clear = (0..moved.len()).all(|i| {
                    let (dx, dy) = (
                        moved[i].0 - moved[nb.0 as usize].0,
                        moved[i].1 - moved[nb.0 as usize].1,
                    );
                    i == nb.0 as usize || dx * dx + dy * dy > 0.25
                });
                let ez_still_ok = target
                    .keys()
                    .filter(|b| !out.unexpressed_double_bonds.contains(b))
                    .all(|&b| ez_drawn_as_declared(mol, &moved, b) == Some(true));
                if !clear || !ez_still_ok {
                    continue;
                }
                let (wedges, unexpressed) = place_wedges(mol, &moved, &centres, &ez_bonds);
                if unexpressed.len() < out.unexpressed_centres.len() {
                    out.coords = moved;
                    out.wedges = wedges;
                    out.unexpressed_centres = unexpressed;
                    improved = true;
                    break;
                }
            }
            if improved {
                break;
            }
        }
        if !improved {
            break;
        }
    }
    out
}

/// The direction (radians) through the middle of the widest gap between
/// `centre`'s bonds other than the one to `skip`.
fn widest_gap_direction(
    mol: &Molecule,
    coords: &[(f64, f64)],
    centre: AtomIdx,
    skip: AtomIdx,
) -> Option<f64> {
    gap_directions(mol, coords, centre, skip).into_iter().next()
}

/// The middles of the gaps between `centre`'s bonds other than the one to
/// `skip`, widest gap first.
fn gap_directions(
    mol: &Molecule,
    coords: &[(f64, f64)],
    centre: AtomIdx,
    skip: AtomIdx,
) -> Vec<f64> {
    use std::f64::consts::TAU;
    let c = coords[centre.0 as usize];
    let mut angles: Vec<f64> = mol
        .neighbors(centre)
        .filter(|&(nb, _)| nb != skip)
        .map(|(nb, _)| {
            let p = coords[nb.0 as usize];
            (p.1 - c.1).atan2(p.0 - c.0).rem_euclid(TAU)
        })
        .collect();
    if angles.is_empty() {
        return Vec::new();
    }
    angles.sort_by(f64::total_cmp);
    let mut gaps: Vec<(f64, f64)> = angles
        .iter()
        .enumerate()
        .map(|(i, &a)| {
            let next = angles.get(i + 1).copied().unwrap_or(angles[0] + TAU);
            (next - a, a + (next - a) / 2.0)
        })
        .collect();
    // Stable: equal gaps keep their angular order, as the widest-gap scan did.
    gaps.sort_by(|a, b| b.0.total_cmp(&a.0));
    gaps.into_iter().map(|(_, mid)| mid).collect()
}

/// Whether segments `a`-`b` and `c`-`d` cross properly.
fn segments_cross(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    let orient = |p: (f64, f64), q: (f64, f64), r: (f64, f64)| {
        (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)
    };
    let (d1, d2) = (orient(c, d, a), orient(c, d, b));
    let (d3, d4) = (orient(a, b, c), orient(a, b, d));
    d1 * d2 < -1e-9 && d3 * d4 < -1e-9
}

/// Clashes (non-bonded atoms closer than 0.45 bonds) and bond crossings the
/// terminal atom `t` at `coords[t]`, bonded to `centre`, makes.
fn terminal_defects(mol: &Molecule, coords: &[(f64, f64)], centre: AtomIdx, t: AtomIdx) -> usize {
    let (c, p) = (coords[centre.0 as usize], coords[t.0 as usize]);
    let clashes = moved_clashes(mol, coords, &[t]);
    let crossings = mol
        .bonds()
        .filter(|(_, b)| ![b.atom1, b.atom2].iter().any(|&a| a == centre || a == t))
        .filter(|(_, b)| {
            segments_cross(c, p, coords[b.atom1.0 as usize], coords[b.atom2.0 as usize])
        })
        .count();
    clashes + crossings
}

/// Put each one-bond neighbour of a stereocentre (a methyl, OH, halogen,
/// explicit H) in the middle of the widest gap between the centre's other
/// bonds, so a wedge to it reads unambiguously. When that spot lands on
/// another atom or across a bond (a crowded fused ring), the next widest gap
/// is tried, and the atom stays where the layout put it when no gap is
/// clear and its own spot is.
fn place_terminal_neighbours(mol: &Molecule, coords: &mut [(f64, f64)], centres: &[AtomIdx]) {
    for &centre in centres {
        let terminal: Vec<AtomIdx> = mol
            .neighbors(centre)
            .map(|(nb, _)| nb)
            .filter(|&nb| mol.degree(nb) == 1)
            .collect();
        for t in terminal {
            let c = coords[centre.0 as usize];
            let original = coords[t.0 as usize];
            let original_defects = terminal_defects(mol, coords, centre, t);
            let mut best: Option<((f64, f64), usize)> = None;
            for dir in gap_directions(mol, coords, centre, t) {
                let spot = (c.0 + 1.5 * dir.cos(), c.1 + 1.5 * dir.sin());
                coords[t.0 as usize] = spot;
                let defects = terminal_defects(mol, coords, centre, t);
                if best.is_none_or(|(_, d)| defects < d) {
                    best = Some((spot, defects));
                }
                if defects == 0 {
                    break;
                }
            }
            coords[t.0 as usize] = match best {
                Some((spot, defects)) if defects <= original_defects => spot,
                _ => original,
            };
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
    ez_bonds: &[BondIdx],
) -> (HashMap<BondIdx, Wedge>, Vec<AtomIdx>) {
    let mut wedges: HashMap<BondIdx, Wedge> = HashMap::new();
    let mut unexpressed = Vec::new();
    // One unmarked copy with the accepted wedges set in place. The MOL
    // reader takes a wedge's direction from its first atom, so a wedge from
    // `centre` on a bond stored the other way round is the opposite order
    // there (it reads the same at both ends).
    let mut drawn = without_bond_marks(mol, &wedges);
    // Centre -> the neighbour its wedge points at.
    let mut wedge_ends: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let stored_order = |bond: BondIdx, w: Wedge| {
        if mol.bond(bond).atom1 == w.start {
            w.order
        } else if w.order == BondOrder::Up {
            BondOrder::Down
        } else {
            BondOrder::Up
        }
    };
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
                // that is not on a stereo double bond (a wedge there would
                // sit next to an E/Z marker).
                let mut rank = 0;
                if mol.neighbors(nb).any(|(_, nbb)| ez_bonds.contains(&nbb)) {
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
        for separation in [NEIGHBOUR_SEPARATION_DEG, FALLBACK_NEIGHBOUR_SEPARATION_DEG] {
            for &(_, nb, bond) in &candidates {
                for order in [BondOrder::Up, BondOrder::Down] {
                    let wedge = Wedge {
                        start: centre,
                        order,
                    };
                    drawn.set_bond_order(bond, stored_order(bond, wedge));
                    // A wedge also shows at its wide end, so a centre there
                    // that already has its wedge must still read as declared.
                    let other_ok = || {
                        !is_centre(nb)
                            || wedge_ends.get(&nb).is_none_or(|&end| {
                                wedge_matches(mol, &drawn, coords, nb, end) == Some(true)
                            })
                    };
                    if wedge_matches_with(mol, &drawn, coords, centre, nb, separation) == Some(true)
                        && other_ok()
                    {
                        wedge_ends.insert(centre, nb);
                        wedges.insert(bond, wedge);
                        placed = true;
                        break;
                    }
                    drawn.set_bond_order(bond, BondOrder::Single);
                }
                if placed {
                    break;
                }
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
            if let Some((w, b)) = mol.neighbors(AtomIdx(v as u32)).nth(cursor) {
                if let Some(top) = stack.last_mut() {
                    top.2 += 1;
                }
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
/// What to draw for `mol`: a copy with the MOL writer's wedges (SMILES `/`
/// `\` marks become plain bonds) and its stereo layout in depiction units
/// (E/Z double bonds drawn with their declared geometry), or `mol` itself
/// with the ordinary layout when it has no stereo.
pub fn depiction_with_stereo(mol: &Molecule) -> (Option<Molecule>, chematic_depict::Layout) {
    if !needs_stereo_depiction(mol) {
        return (None, chematic_depict::compute_layout(mol));
    }
    let depiction = stereo_depiction(mol, &[]);
    let scale = 1.5 / chematic_depict::layout::BOND_LEN;
    let layout = chematic_depict::Layout {
        coords: depiction
            .coords
            .iter()
            .map(|&(x, y)| chematic_depict::Point::new(x / scale, -y / scale))
            .collect(),
    };
    (Some(without_bond_marks(mol, &depiction.wedges)), layout)
}

/// Runs `draw` on what [`depiction_with_stereo`] gives for `mol`.
pub fn with_stereo_depiction<R>(
    mol: &Molecule,
    draw: impl FnOnce(&Molecule, &chematic_depict::Layout) -> R,
) -> R {
    let (copy, layout) = depiction_with_stereo(mol);
    draw(copy.as_ref().unwrap_or(mol), &layout)
}

/// 2D coordinates (Å, y up) for writing `mol` to a MOL block: its stereo
/// depiction when it has stereo, the ordinary layout otherwise.
pub fn mol_block_coords(mol: &Molecule) -> Vec<(f64, f64)> {
    if needs_stereo_depiction(mol) {
        stereo_depiction(mol, &[]).coords
    } else {
        layout_angstrom(mol, false)
    }
}

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
        // Macrocycle trans double bonds (a regular polygon draws them cis).
        "C1CCCCC/C=C/CCCC1",
        "O=C1CCCCC/C=C/C=C\\CCCCO1",
        // Bridged systems whose layout puts two bonds of a centre in nearly
        // one direction; neighbouring centres whose wedges meet.
        "CC12N[C@H](CC3=C1CCCC3)c1ccccc12",
        "O=C(O)[C@@H]1[C@H](C(=O)O)[C@H]2CC[C@@H]1O2",
        // Pentacyclic: the last fused ring used to fall to the fallback
        // placement with two stretched bonds.
        "CCOC(=O)[C@]12CCC(C)(C)CC1C1C(=O)C=C3[C@@]4(C)C=C(C#N)C(=O)C(C)(C)[C@@H]4CC[C@@]3(C)[C@]1(C)CC2",
        // A bridgehead and a cage spiro centre whose layout leaves two
        // neighbours 9 and 12 degrees apart: wedged with the narrower
        // fallback separation (RDKit 2026.03.6 reads both as declared).
        "O=C1OC(=O)[C@@H]2[C@H]1[C@@H]1O[C@H]2C[C@@H]1COC(=O)[C@@H]1C[C@H]1c1ccccc1",
        "CC1CC(=O)OC[C@]23C[C@@H](O)[C@H](C)C[C@H]2O[C@@H]2C[C@@H](OC(=O)/C=C\\C=C\\C(C(C)O)O[C@@H](O)C1)[C@@]3(C)C21CO1",
        // Bryostatin (ChEMBL 5k row 3533): the macrocycle's trans ring double
        // bond cannot be flipped on the first layout, only on the layout of
        // the reversed atom order.
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2CC([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
    ];

    #[test]
    fn stereo_report_lists_what_a_block_cannot_carry() {
        let meta = MolMetadata::default();
        for smi in CASES {
            let mol = parse(smi).unwrap();
            let (block, loss) = crate::mol2000::write_mol_with_stereo_report(&mol, &meta, &[]);
            assert!(loss.is_empty(), "{smi}: {loss}");
            assert_eq!(block, write_mol(&mol, &meta));
            let (block3, loss3) =
                crate::mol3000::write_mol_v3000_with_stereo_report(&mol, &meta, &[]);
            assert!(loss3.is_empty(), "{smi}: {loss3}");
            assert_eq!(block3, write_mol_v3000(&mol, &meta, &[]));
        }
        let sp = parse("F[Pt@SP1](Cl)(Br)I").unwrap();
        let loss = crate::mol2000::write_mol_with_stereo_report(&sp, &meta, &[]).1;
        assert_eq!(
            loss.non_tetrahedral_centres,
            vec![chematic_core::AtomIdx(1)]
        );
        assert!(!loss.is_empty());
        assert!(loss.to_string().contains("square-planar"), "{loss}");

        // Coordinates that put every bond of the centre on one line: no
        // wedge can express it, and the report says so.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let line: Vec<(f64, f64)> = (0..mol.atom_count())
            .map(|i| (i as f64 * 1.5, 0.0))
            .collect();
        let loss = crate::mol2000::write_mol_with_stereo_report(&mol, &meta, &line).1;
        assert_eq!(loss.centres, vec![chematic_core::AtomIdx(1)], "{loss}");

        // V2000 has no field for enhanced stereo groups; V3000 does.
        let v3 = "\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 5 4 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1.3 0.75 0 0\nM  V30 3 F 2.6 0 0 0\nM  V30 4 Cl 1.3 2.25 0 0\nM  V30 5 Br 0 1.5 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 2 4 CFG=1\nM  V30 4 1 2 5\nM  V30 END BOND\nM  V30 BEGIN COLLECTION\nM  V30 MDLV30/STERAC1 ATOMS=(1 2)\nM  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n";
        let grouped = read_mol_v3000_with_diagnostics(v3).unwrap().mol;
        assert!(!grouped.stereo_groups().is_empty());
        let v2_loss = crate::mol2000::write_mol_with_stereo_report(&grouped, &meta, &[]).1;
        assert!(v2_loss.stereo_groups_dropped, "{v2_loss}");
        let v3_loss = crate::mol3000::write_mol_v3000_with_stereo_report(&grouped, &meta, &[]).1;
        assert!(!v3_loss.stereo_groups_dropped, "{v3_loss}");
    }

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
    fn laid_out_molecule_without_stereo_marks_undeclared_double_bonds_either() {
        // A molecule with no stereo at all takes the plain layout, which
        // still draws each double bond cis or trans. Readers (RDKit among
        // them) take that drawn geometry as E/Z unless the bond is "either":
        // 533 of the exposed 10k rows read back with an invented E/Z.
        use crate::mol2000::{write_laid_out_mol, write_mol_with_coords};
        for smi in [
            "CC=CC",
            "CN(C)N=Nc1ccccc1",
            "O=C(Nc1ccccc1)C(=NNc1ccccc1)N=Nc1ccccc1",
        ] {
            let mol = parse(smi).unwrap();
            assert!(!super::needs_stereo_depiction(&mol));
            let (v2, loss) = write_laid_out_mol(&mol, &MolMetadata::default());
            assert!(loss.is_empty(), "{smi}");
            let either = super::undeclared_stereo_double_bonds(&mol);
            assert!(!either.is_empty(), "{smi}");
            for &b in &either {
                let line = v2.lines().nth(4 + mol.atom_count() + b.0 as usize).unwrap();
                assert_eq!(line[9..12].trim(), "3", "{v2}");
            }
            let back = read_mol_with_diagnostics(&v2).unwrap().mol;
            assert_eq!(canonical_smiles(&back), canonical_smiles(&mol), "{v2}");
            // Caller-supplied coordinates are data: that contract is unchanged.
            let coords = super::mol_block_coords(&mol);
            let plain = write_mol_with_coords(&mol, &MolMetadata::default(), &coords);
            for &b in &either {
                let line = plain
                    .lines()
                    .nth(4 + mol.atom_count() + b.0 as usize)
                    .unwrap();
                assert_eq!(line[9..12].trim(), "0", "{plain}");
            }
        }
        // Symmetric ends and small rings are not stereo bonds.
        for smi in ["CC(C)=C(C)C", "C1=CCCCC1", "C=CC"] {
            let mol = parse(smi).unwrap();
            assert!(
                super::undeclared_stereo_double_bonds(&mol).is_empty(),
                "{smi}"
            );
        }
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

    /// The depiction (SVG, depict data) shows the declared E/Z geometry and
    /// a wedge from the centre, and draws SMILES `/` `\` as plain bonds.
    #[test]
    fn stereo_depiction_does_not_stack_atoms() {
        use chematic_core::AtomIdx;
        // Tamoxifen's cis aryl rings (drawn on each other when the alkene's
        // substituents were 60° apart) and a terminal OH of a fused-ring
        // centre (moved into the widest gap, onto a ring atom).
        for smi in [
            "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
            "Cc1cc(=O)n(C(=O)OC(C)(C)C)c2c3c(ccc12)OC(C)(C)[C@H](O)[C@@H]3O",
        ] {
            let mol = parse(smi).unwrap();
            let d = super::stereo_depiction(&mol, &[]);
            assert!(d.unexpressed_centres.is_empty() && d.unexpressed_double_bonds.is_empty());
            let n = mol.atom_count();
            for i in 0..n {
                for j in i + 1..n {
                    let (a, b) = (AtomIdx(i as u32), AtomIdx(j as u32));
                    let (p, q) = (d.coords[i], d.coords[j]);
                    let dist = ((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2)).sqrt();
                    assert!(
                        mol.bond_between(a, b).is_some() || dist >= 0.6,
                        "{smi}: atoms {i} and {j} {dist:.2} Å apart"
                    );
                }
            }
        }
    }

    #[test]
    fn depiction_shows_declared_stereo() {
        use chematic_core::{AtomIdx, BondOrder};
        let side =
            |q: chematic_depict::Point, a: chematic_depict::Point, b: chematic_depict::Point| {
                (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x)
            };
        for (smi, cis) in [("F/C=C\\F", true), ("F/C=C/F", false)] {
            let mol = parse(smi).unwrap();
            let (copy, layout) = super::depiction_with_stereo(&mol);
            let copy = copy.unwrap();
            let p = |i: u32| layout.coords[i as usize];
            assert_eq!(
                side(p(0), p(1), p(2)) * side(p(3), p(1), p(2)) > 0.0,
                cis,
                "{smi}"
            );
            assert!(
                copy.bonds()
                    .all(|(_, b)| !matches!(b.order, BondOrder::Up | BondOrder::Down)),
                "{smi}"
            );
        }
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let (copy, _) = super::depiction_with_stereo(&mol);
        let wedges: Vec<_> = copy
            .unwrap()
            .bonds()
            .filter(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
            .map(|(_, b)| b.atom1)
            .collect();
        assert_eq!(wedges, vec![AtomIdx(1)]);
        assert!(
            super::depiction_with_stereo(&parse("CCO").unwrap())
                .0
                .is_none()
        );
    }
}
