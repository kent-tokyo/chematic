//! Consistent SMILES `/`/`\` markers for double-bond (E/Z) stereo.
//!
//! A double bond's geometry is the relative side of one substituent at each
//! end. SMILES spells it with direction markers on substituent bonds, read
//! "up" or "down" from the end. Only the relative direction of a double
//! bond's two markers matters, so all of one double bond's markers may flip
//! together; but a single bond between two double bonds (conjugation) is a
//! substituent of both, and its one marker must suit both. When a marker is
//! moved (an H atom removed, a template writing new markers next to carried
//! ones) an alkene end can end up with both substituents marked on the same
//! side, which readers treat as a conflict and drop, or read for the wrong
//! double bond. This module rewrites the markers of a set of geometry facts
//! with a two-colouring of those shared bonds, as RDKit does for conjugated
//! polyenes.

use std::collections::HashMap;

use crate::{AtomIdx, BondIdx, BondOrder, Molecule};

/// One end of an [`EzFact`]: the substituent bond used as reference and its
/// side, "up" as seen from the end atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EzFactEnd {
    pub atom: AtomIdx,
    pub reference: BondIdx,
    pub up: bool,
}

/// The geometry of one double bond: its two reference substituents and
/// their sides. Equal sides mean the references are cis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EzFact {
    pub double: BondIdx,
    pub ends: [EzFactEnd; 2],
}

/// The raw direction of `bidx` and the atom it is read from. A stored
/// direction ([`Molecule::bond_direction`]) wins over a literal `Up`/`Down`
/// order, as in the SMILES writers.
pub fn raw_direction(mol: &Molecule, bidx: BondIdx) -> Option<(BondOrder, AtomIdx)> {
    if let Some(direction) = mol.bond_direction(bidx) {
        let anchor = mol
            .bond_direction_anchor(bidx)
            .unwrap_or(mol.bond(bidx).atom1);
        return Some((direction, anchor));
    }
    let bond = mol.bond(bidx);
    matches!(bond.order, BondOrder::Up | BondOrder::Down).then_some((bond.order, bond.atom1))
}

/// `true` when `direction`, read from `anchor`, puts the bond "up" as seen
/// from `end`.
pub fn direction_is_up(direction: BondOrder, anchor: AtomIdx, end: AtomIdx) -> bool {
    match direction {
        BondOrder::Up => anchor == end,
        BondOrder::Down => anchor != end,
        _ => false,
    }
}

/// The side of `bidx` seen from `end`, if the bond carries a marker.
pub fn marker_side(mol: &Molecule, bidx: BondIdx, end: AtomIdx) -> Option<bool> {
    raw_direction(mol, bidx).map(|(d, a)| direction_is_up(d, a, end))
}

/// Substituent (non-double) bonds of `end`, other than `double`.
fn substituents(mol: &Molecule, end: AtomIdx, double: BondIdx) -> Vec<BondIdx> {
    mol.neighbors(end)
        .filter(|&(_, b)| b != double && mol.bond(b).order != BondOrder::Double)
        .map(|(_, b)| b)
        .collect()
}

/// Rewrite the markers of every double bond in `facts` so each keeps its
/// geometry and no bond is asked for two directions. Marks each fact's two
/// reference bonds; markers on other substituents of those ends are
/// removed. Returns `false`, leaving `mol` unchanged, when the facts cannot
/// be satisfied together or a bond that would change is not a plain single
/// bond carrying a literal marker (stored directions on aromatic ring bonds
/// are left alone).
pub fn write_ez_facts(mol: &mut Molecule, facts: &[EzFact]) -> bool {
    if facts.is_empty() {
        return true;
    }
    // A reference that cannot carry a marker (a zero-order bond, which
    // RDKit keeps as a stereo atom's bond) hands its side to the end's
    // other substituent.
    let markable = |mol: &Molecule, b: BondIdx| {
        matches!(
            mol.bond(b).order,
            BondOrder::Single | BondOrder::Up | BondOrder::Down
        )
    };
    let mut swapped: Vec<EzFact> = facts.to_vec();
    for fact in &mut swapped {
        for end in &mut fact.ends {
            if markable(mol, end.reference) {
                continue;
            }
            if let Some(sibling) = substituents(mol, end.atom, fact.double)
                .into_iter()
                .find(|&b| b != end.reference && markable(mol, b))
            {
                end.reference = sibling;
                end.up = !end.up;
            }
        }
    }
    let facts = &swapped[..];
    let marked: std::collections::HashSet<BondIdx> = facts
        .iter()
        .flat_map(|f| f.ends.iter().map(|e| e.reference))
        .collect();
    // Side wanted for each marked bond, per fact, read from its atom1.
    let mut wanted: HashMap<BondIdx, Vec<(usize, bool)>> = HashMap::new();
    let mut cleared: Vec<BondIdx> = Vec::new();
    for (i, fact) in facts.iter().enumerate() {
        for end in &fact.ends {
            for b in substituents(mol, end.atom, fact.double) {
                if mol.bond_direction(b).is_some() {
                    return false;
                }
                if marked.contains(&b) {
                    let up = if b == end.reference { end.up } else { !end.up };
                    let from_atom1 = if mol.bond(b).atom1 == end.atom {
                        up
                    } else {
                        !up
                    };
                    wanted.entry(b).or_default().push((i, from_atom1));
                } else if matches!(mol.bond(b).order, BondOrder::Up | BondOrder::Down) {
                    cleared.push(b);
                }
            }
        }
    }
    for &b in &marked {
        if !matches!(
            mol.bond(b).order,
            BondOrder::Single | BondOrder::Up | BondOrder::Down
        ) {
            return false;
        }
    }
    // Parity union-find: flip[i] mirrors all of fact i's markers.
    let mut parent: Vec<(usize, bool)> = (0..facts.len()).map(|i| (i, false)).collect();
    fn find(parent: &mut [(usize, bool)], x: usize) -> (usize, bool) {
        let (p, parity) = parent[x];
        if p == x {
            return (x, false);
        }
        let (root, rp) = find(parent, p);
        parent[x] = (root, parity ^ rp);
        (root, parity ^ rp)
    }
    let mut keys: Vec<BondIdx> = wanted.keys().copied().collect();
    keys.sort_unstable();
    for b in &keys {
        let list = &wanted[b];
        for pair in list.windows(2) {
            let ((a, ua), (c, uc)) = (pair[0], pair[1]);
            let (ra, pa) = find(&mut parent, a);
            let (rc, pc) = find(&mut parent, c);
            let need = ua ^ uc;
            if ra == rc {
                if pa ^ pc != need {
                    return false;
                }
            } else {
                parent[ra] = (rc, pa ^ pc ^ need);
            }
        }
    }
    for b in cleared {
        mol.set_bond_order(b, BondOrder::Single);
    }
    for b in keys {
        let (i, from_atom1) = wanted[&b][0];
        let (_, flip) = find(&mut parent, i);
        let direction = if from_atom1 ^ flip {
            BondOrder::Up
        } else {
            BondOrder::Down
        };
        if mol.bond(b).order != direction {
            mol.set_bond_order(b, direction);
        }
    }
    true
}

/// The geometry facts the markers of `mol` spell, one per stereo double
/// bond (both ends marked). An end whose two substituents are marked on the
/// same side holds markers of two double bonds (its own, and a conjugated
/// neighbour's on the shared single bond); its own is taken to be the one
/// on the bond no other double bond reads. Returns `None` when such a
/// conflict cannot be attributed (both or neither substituent shared), and
/// otherwise the facts and whether any end had a conflict.
pub fn read_ez_facts(mol: &Molecule) -> Option<(Vec<EzFact>, bool)> {
    let n = mol.atom_count();
    let mut double_ends = vec![false; n];
    for (_, bond) in mol.bonds() {
        if bond.order == BondOrder::Double {
            double_ends[bond.atom1.0 as usize] = true;
            double_ends[bond.atom2.0 as usize] = true;
        }
    }
    let mut conflict = false;
    let mut facts: Vec<EzFact> = Vec::new();
    for (bidx, bond) in mol.bonds() {
        if bond.order != BondOrder::Double
            || mol.atom(bond.atom1).aromatic && mol.atom(bond.atom2).aromatic
        {
            continue;
        }
        let mut ends: Vec<EzFactEnd> = Vec::with_capacity(2);
        for end in [bond.atom1, bond.atom2] {
            let subs = substituents(mol, end, bidx);
            if subs.is_empty() || subs.len() > 2 {
                break;
            }
            let marks: Vec<(BondIdx, bool, bool)> = subs
                .iter()
                .filter_map(|&b| {
                    let side = marker_side(mol, b, end)?;
                    let far = mol.bond(b).other(end)?;
                    Some((b, side, double_ends[far.0 as usize]))
                })
                .collect();
            let chosen = match marks.as_slice() {
                [] => None,
                [one] => Some(*one),
                [a, b] if a.1 != b.1 => Some(if a.2 && !b.2 { *b } else { *a }),
                [a, b] => {
                    conflict = true;
                    match (a.2, b.2) {
                        (false, true) => Some(*a),
                        (true, false) => Some(*b),
                        _ => return None,
                    }
                }
                _ => None,
            };
            let Some((reference, up, _)) = chosen else {
                break;
            };
            ends.push(EzFactEnd {
                atom: end,
                reference,
                up,
            });
        }
        if let [a, b] = ends.as_slice() {
            facts.push(EzFact {
                double: bidx,
                ends: [*a, *b],
            });
        }
    }
    Some((facts, conflict))
}

/// Repair markers that ask an alkene end for the same side twice (see
/// [`read_ez_facts`] for how they are attributed), rewriting every stereo
/// double bond with [`write_ez_facts`]. Returns `true` when markers
/// changed; a molecule without such a conflict, or one whose conflict cannot
/// be attributed, is left unchanged.
pub fn reconcile_ez_markers(mol: &mut Molecule) -> bool {
    let directional = |order: BondOrder| matches!(order, BondOrder::Up | BondOrder::Down);
    if !mol.has_bond_directions() && !mol.bonds().any(|(_, b)| directional(b.order)) {
        return false;
    }
    let Some((facts, true)) = read_ez_facts(mol) else {
        return false;
    };
    let mut trial = mol.clone();
    if write_ez_facts(&mut trial, &facts) {
        *mol = trial;
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Element, MoleculeBuilder};

    /// F-C1=C2-C3=C4-C, with markers given per bond as (from, to, order).
    fn diene(markers: &[(usize, usize, BondOrder)]) -> Molecule {
        let mut b = MoleculeBuilder::new();
        let ids: Vec<AtomIdx> = [
            Element::F,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
        ]
        .into_iter()
        .map(|e| b.add_atom(Atom::new(e)))
        .collect();
        let order_of = |x: usize, y: usize| {
            markers
                .iter()
                .find(|&&(a, c, _)| (a, c) == (x, y))
                .map(|&(_, _, o)| o)
        };
        for (x, y, base) in [
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 5, BondOrder::Single),
        ] {
            b.add_bond(ids[x], ids[y], order_of(x, y).unwrap_or(base))
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn no_conflict_is_left_alone() {
        let mut m = diene(&[
            (0, 1, BondOrder::Up),
            (2, 3, BondOrder::Up),
            (4, 5, BondOrder::Up),
        ]);
        let before = m.clone();
        assert!(!reconcile_ez_markers(&mut m));
        for (b, bond) in before.bonds() {
            assert_eq!(m.bond(b).order, bond.order);
        }
    }

    #[test]
    fn write_facts_ties_conjugated_double_bonds() {
        // Facts: F trans C3 across C1=C2; C2 trans C5 across C3=C4. The
        // shared C2-C3 bond is the reference of the second fact and a
        // substituent of the first.
        let mut m = diene(&[]);
        let b01 = m.bond_between(AtomIdx(0), AtomIdx(1)).unwrap().0;
        let b23 = m.bond_between(AtomIdx(2), AtomIdx(3)).unwrap().0;
        let b45 = m.bond_between(AtomIdx(4), AtomIdx(5)).unwrap().0;
        let d12 = m.bond_between(AtomIdx(1), AtomIdx(2)).unwrap().0;
        let d34 = m.bond_between(AtomIdx(3), AtomIdx(4)).unwrap().0;
        let facts = [
            EzFact {
                double: d12,
                ends: [
                    EzFactEnd {
                        atom: AtomIdx(1),
                        reference: b01,
                        up: true,
                    },
                    EzFactEnd {
                        atom: AtomIdx(2),
                        reference: b23,
                        up: false,
                    },
                ],
            },
            EzFact {
                double: d34,
                ends: [
                    EzFactEnd {
                        atom: AtomIdx(3),
                        reference: b23,
                        up: true,
                    },
                    EzFactEnd {
                        atom: AtomIdx(4),
                        reference: b45,
                        up: false,
                    },
                ],
            },
        ];
        assert!(write_ez_facts(&mut m, &facts));
        // Both facts hold on the written markers.
        for f in &facts {
            let s0 = marker_side(&m, f.ends[0].reference, f.ends[0].atom).unwrap();
            let s1 = marker_side(&m, f.ends[1].reference, f.ends[1].atom).unwrap();
            assert_eq!(s0 == s1, f.ends[0].up == f.ends[1].up);
        }
    }
}
