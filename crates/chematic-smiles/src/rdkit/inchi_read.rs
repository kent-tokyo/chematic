//! RDKit's `InchiToMol` post-processing (RDKit 2026.03.1
//! `External/INCHI-API/inchi.cpp`): the molecule `Chem.MolFromInchi` builds
//! from the IUPAC library's `GetStructFromINCHI` output, before `removeHs`,
//! sanitization and stereo perception. The library call itself lives in
//! chematic-inchi (native feature); this module only sees its output as
//! plain data, so it needs no C library.
//!
//! The result is returned as a chematic [`Molecule`] encoding RDKit's
//! pre-sanitization state (bracket H counts, RDKit's bond numbering, chiral
//! tags as neighbour orders, `/` `\` single bonds), so that
//! [`crate::rdkit_canonical_smiles`] on it runs the rest of RDKit's
//! pipeline.

use std::collections::{BTreeSet, HashSet, VecDeque};

use chematic_core::{
    Atom as CAtom, AtomIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder,
    STEREO_H_SENTINEL,
};

use super::RdkitSmilesError;
use super::mol::{Atom, Bond, BondDir, BondType, ChiralTag, Mol};
use super::periodic;
use super::stereo::assign_atom_cip_ranks;

/// `ISOTOPIC_SHIFT_FLAG` (inchi_api.h).
const ISOTOPIC_SHIFT_FLAG: i32 = 10000;
const INCHI_BOND_TYPE_TRIPLE: i8 = 3;
const INCHI_BOND_TYPE_ALTERN: i8 = 4;
const INCHI_PARITY_NONE: i8 = 0;
const INCHI_PARITY_ODD: i8 = 1;
const INCHI_PARITY_EVEN: i8 = 2;
const INCHI_PARITY_UNDEFINED: i8 = 4;
const INCHI_STEREO_TYPE_DOUBLE_BOND: i8 = 1;
const INCHI_STEREO_TYPE_TETRAHEDRAL: i8 = 2;

/// One `inchi_Atom` of `GetStructFromINCHI`'s output.
#[derive(Debug, Clone, Default)]
pub struct InchiOutputAtom {
    /// `elname`.
    pub element: String,
    /// `(neighbor, bond_type, bond_stereo)` for each of `num_bonds`.
    pub bonds: Vec<(usize, i8, i8)>,
    /// `num_iso_H` (index 0: non-isotopic H).
    pub num_iso_h: [i8; 4],
    /// `isotopic_mass`.
    pub isotopic_mass: i16,
    /// `radical`.
    pub radical: i8,
    /// `charge`.
    pub charge: i8,
}

/// One `inchi_Stereo0D` of `GetStructFromINCHI`'s output.
#[derive(Debug, Clone, Default)]
pub struct InchiOutputStereo0D {
    pub neighbor: [i16; 4],
    pub central_atom: i16,
    pub stereo_type: i8,
    pub parity: i8,
}

fn unsupported(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(what.into())
}

/// `assignBondDirs`: bond directions satisfying "same direction" (`z`) and
/// "different direction" (`e`) pairs; false when they conflict.
fn assign_bond_dirs(mol: &mut Mol, z: &[(usize, usize)], e: &[(usize, usize)]) -> bool {
    let mut pending: BTreeSet<usize> = BTreeSet::new();
    for &(a, b) in z.iter().chain(e) {
        pending.insert(a);
        pending.insert(b);
    }
    let mut queue: VecDeque<(usize, BondDir)> = VecDeque::new();
    while !pending.is_empty() || !queue.is_empty() {
        let Some((cur, dir)) = queue.pop_front() else {
            let first = *pending.iter().next().expect("pending");
            queue.push_back((first, BondDir::EndUpRight));
            continue;
        };
        if mol.bonds[cur].dir != BondDir::None {
            if mol.bonds[cur].dir != dir {
                return false;
            }
            continue;
        }
        mol.bonds[cur].dir = dir;
        pending.remove(&cur);
        for (rules, d) in [(z, dir), (e, dir.flipped())] {
            for &(a, b) in rules {
                let other = if a == cur {
                    Some(b)
                } else if b == cur {
                    Some(a)
                } else {
                    None
                };
                if let Some(o) = other
                    && o != cur
                {
                    if mol.bonds[o].dir != BondDir::None {
                        if mol.bonds[o].dir != d {
                            return false;
                        }
                    } else {
                        queue.push_back((o, d));
                    }
                }
            }
        }
    }
    true
}

/// `InchiToMol` up to (not including) its `cleanUp`, `removeHs`,
/// sanitization and stereo assignment, returned as a chematic molecule
/// (see the module documentation).
pub fn rdkit_molecule_from_inchi_output(
    atoms: &[InchiOutputAtom],
    stereo: &[InchiOutputStereo0D],
) -> Result<Molecule, RdkitSmilesError> {
    let mut m = Mol::default();
    let mut isotopes: Vec<(u32, usize, i8)> = Vec::new();
    for (i, ia) in atoms.iter().enumerate() {
        let anum = periodic::atomic_number(&ia.element)
            .ok_or_else(|| unsupported(format!("element {}", ia.element)))?;
        let mut atom = Atom::new(anum);
        let ref_weight = (periodic::atomic_weight(anum) + 0.5) as i32;
        if ia.isotopic_mass != 0 {
            let iso = i32::from(ia.isotopic_mass) - ISOTOPIC_SHIFT_FLAG;
            if iso != 0 {
                atom.isotope = (iso + ref_weight).max(0) as u32;
            }
        }
        atom.charge = i32::from(ia.charge);
        if ia.radical == 2 || ia.radical == 3 {
            atom.radicals = (ia.radical - 1) as u32;
        }
        atom.num_explicit_hs = ia.num_iso_h[0].max(0) as u32;
        if ia.num_iso_h[1] != 0 {
            isotopes.push((1, i, ia.num_iso_h[1]));
        } else if ia.num_iso_h[2] != 0 {
            isotopes.push((2, i, ia.num_iso_h[2]));
        } else if ia.num_iso_h[3] != 0 {
            isotopes.push((3, i, ia.num_iso_h[3]));
        }
        atom.no_implicit = true;
        m.add_atom(atom);
    }
    let mut register: HashSet<(usize, usize)> = HashSet::new();
    for (i, ia) in atoms.iter().enumerate() {
        for &(nbr, bt, bst) in &ia.bonds {
            if register.contains(&(i, nbr)) || register.contains(&(nbr, i)) {
                continue;
            }
            register.insert((i, nbr));
            let (bt, aromatic) = match bt {
                1 => (BondType::Single, false),
                2 => (BondType::Double, false),
                INCHI_BOND_TYPE_TRIPLE => (BondType::Triple, false),
                INCHI_BOND_TYPE_ALTERN => (BondType::Aromatic, true),
                other => return Err(unsupported(format!("InChI bond type {other}"))),
            };
            if bst != 0 {
                return Err(unsupported("InChI 2D bond stereo"));
            }
            let mut b = Bond::new(i, nbr, bt);
            b.aromatic = aromatic;
            m.add_bond(b);
        }
    }
    for (iso, aid, repeat) in isotopes {
        for _ in 0..repeat.max(0) {
            let mut h = Atom::new(1);
            h.isotope = iso;
            let j = m.add_atom(h);
            m.add_bond(Bond::new(j, aid, BondType::Single));
        }
    }
    m.update_property_cache(false)?;

    let mut z_pairs: Vec<(usize, usize)> = Vec::new();
    let mut e_pairs: Vec<(usize, usize)> = Vec::new();
    if !stereo.is_empty() {
        let mut ranks = Vec::new();
        assign_atom_cip_ranks(&m, &mut ranks);
        for s in stereo {
            if s.parity == INCHI_PARITY_NONE || s.parity == INCHI_PARITY_UNDEFINED {
                continue;
            }
            match s.stereo_type {
                INCHI_STEREO_TYPE_DOUBLE_BOND => {
                    let left = s.neighbor[1] as usize;
                    let right = s.neighbor[2] as usize;
                    let orig_left = s.neighbor[0] as usize;
                    let orig_right = s.neighbor[3] as usize;
                    if m.bond_between(left, right).is_none() {
                        continue; // cumulene stereo is ignored
                    }
                    let find_nbrs = |r: usize| -> (Option<usize>, Option<usize>) {
                        let (mut nbr, mut extra): (Option<usize>, Option<usize>) = (None, None);
                        let mut cip: i64 = -1;
                        for &b in &m.atom_bonds[r] {
                            if !matches!(m.bonds[b].bt, BondType::Single | BondType::Aromatic) {
                                continue;
                            }
                            let a = m.bonds[b].other(r);
                            let c = i64::from(ranks[a]);
                            if c > cip {
                                if nbr.is_some() {
                                    extra = nbr;
                                }
                                nbr = Some(a);
                                cip = c;
                            } else {
                                extra = Some(a);
                            }
                        }
                        (nbr, extra)
                    };
                    let (Some(left_nbr), extra_left) = find_nbrs(left) else {
                        continue;
                    };
                    let (Some(right_nbr), extra_right) = find_nbrs(right) else {
                        continue;
                    };
                    let switch_ez = (orig_left == left_nbr) != (orig_right == right_nbr);
                    let mut parity = s.parity;
                    if switch_ez {
                        if parity == INCHI_PARITY_ODD {
                            parity = INCHI_PARITY_EVEN;
                        } else if parity == INCHI_PARITY_EVEN {
                            parity = INCHI_PARITY_ODD;
                        }
                    }
                    let mut pair_bonds = |r: usize, nbr: usize, extra: Option<usize>| -> usize {
                        let bond = m.bond_between(r, nbr).expect("bonded");
                        if let Some(x) = extra {
                            let mut modifier = -1;
                            if m.bonds[bond].begin != r {
                                modifier = -modifier;
                            }
                            let extra_bond = m.bond_between(r, x).expect("bonded");
                            if m.bonds[extra_bond].begin != r {
                                modifier = -modifier;
                            }
                            if modifier == 1 {
                                z_pairs.push((bond, extra_bond));
                            } else {
                                e_pairs.push((bond, extra_bond));
                            }
                        }
                        bond
                    };
                    let left_bond = pair_bonds(left, left_nbr, extra_left);
                    let right_bond = pair_bonds(right, right_nbr, extra_right);
                    let mut modifier = -1;
                    if m.bonds[left_bond].begin != left {
                        modifier = -modifier;
                    }
                    if m.bonds[right_bond].begin != right {
                        modifier = -modifier;
                    }
                    if parity == INCHI_PARITY_ODD {
                        if modifier == 1 {
                            e_pairs.push((left_bond, right_bond));
                        } else {
                            z_pairs.push((left_bond, right_bond));
                        }
                    } else if parity == INCHI_PARITY_EVEN {
                        if modifier == 1 {
                            z_pairs.push((left_bond, right_bond));
                        } else {
                            e_pairs.push((left_bond, right_bond));
                        }
                    }
                }
                INCHI_STEREO_TYPE_TETRAHEDRAL => {
                    let c = s.central_atom as usize;
                    let mut odd = false;
                    let mut nid = 0;
                    if s.neighbor[0] == s.central_atom {
                        nid = 1;
                        if m.degree(c) == 3 {
                            odd = true;
                        }
                    }
                    let mut nbr_bonds = Vec::new();
                    for k in nid..4 {
                        let end = s.neighbor[k] as usize;
                        let b = m
                            .bond_between(c, end)
                            .ok_or_else(|| unsupported("stereo neighbour not bonded"))?;
                        nbr_bonds.push(b);
                    }
                    if nbr_bonds.len() != m.degree(c) {
                        return Err(unsupported("tetrahedral stereo with a non-bonded ligand"));
                    }
                    if m.perturbation_is_odd(c, &nbr_bonds) {
                        odd = !odd;
                    }
                    let mut tag = if s.parity == INCHI_PARITY_ODD {
                        ChiralTag::Ccw
                    } else {
                        ChiralTag::Cw
                    };
                    if odd {
                        tag = if tag == ChiralTag::Ccw {
                            ChiralTag::Cw
                        } else {
                            ChiralTag::Ccw
                        };
                    }
                    m.atoms[c].chiral = tag;
                }
                _ => {}
            }
        }
        // A failure only logs a warning in RDKit.
        let _ = assign_bond_dirs(&mut m, &z_pairs, &e_pairs);
    }
    inchi_clean_up(&mut m);
    to_chematic(&m)
}

/// The chematic molecule `parse::from_chematic` turns back into `m`.
fn to_chematic(m: &Mol) -> Result<Molecule, RdkitSmilesError> {
    let mut b = MoleculeBuilder::new();
    let mut orders: Vec<(AtomIdx, Vec<u32>)> = Vec::new();
    for (i, a) in m.atoms.iter().enumerate() {
        let element = Element::from_symbol(periodic::symbol(a.anum))
            .ok_or_else(|| unsupported(format!("element {}", a.anum)))?;
        let mut atom = CAtom::new(element);
        atom.isotope = (a.isotope != 0).then_some(a.isotope as u16);
        atom.charge = a.charge as i8;
        atom.hydrogen_count = Some(a.num_explicit_hs as u8);
        if a.chiral != ChiralTag::Unspecified {
            // The tag is relative to the bond order; `from_chematic` reads a
            // neighbour order equal to the bond order as the written tag,
            // then applies the parser's degree-3 inversion for a SMILES
            // start atom with one H, which is undone here. The implicit H
            // goes last, as in RDKit's convention.
            let mut order: Vec<u32> = m.nbrs(i).map(|x| x as u32).collect();
            let is_start = !m.nbrs(i).any(|x| x < i);
            let mut tag = a.chiral;
            if m.degree(i) == 3 && is_start && a.num_explicit_hs == 1 {
                tag = if tag == ChiralTag::Cw {
                    ChiralTag::Ccw
                } else {
                    ChiralTag::Cw
                };
            }
            if a.num_explicit_hs > 0 {
                order.push(STEREO_H_SENTINEL);
            }
            atom.chirality = if tag == ChiralTag::Ccw {
                Chirality::CounterClockwise
            } else {
                Chirality::Clockwise
            };
            orders.push((AtomIdx(i as u32), order));
        }
        b.add_atom(atom);
    }
    for bond in &m.bonds {
        let order = match (bond.bt, bond.dir) {
            (BondType::Single, BondDir::EndUpRight) => BondOrder::Up,
            (BondType::Single, BondDir::EndDownRight) => BondOrder::Down,
            (BondType::Single, BondDir::None) => BondOrder::Single,
            (BondType::Double, _) => BondOrder::Double,
            (BondType::Triple, _) => BondOrder::Triple,
            (BondType::Aromatic, _) => BondOrder::Aromatic,
            (bt, _) => return Err(unsupported(format!("bond type {bt:?}"))),
        };
        b.add_bond(AtomIdx(bond.begin as u32), AtomIdx(bond.end as u32), order)
            .map_err(|e| unsupported(e.to_string()))?;
    }
    for (idx, order) in orders {
        b.set_stereo_neighbor_order(idx, order);
    }
    Ok(b.build())
}

/// `findAlternatingBonds`: the shortest path of alternating bonds from
/// `start` to an atom of `anum` with `charge`, as a stack of bonds whose top
/// (last element) is the bond at `start`.
#[allow(clippy::too_many_arguments)]
fn find_alternating_bonds(
    m: &Mol,
    current: usize,
    anum: u32,
    charge: i32,
    next_bt: Option<BondType>,
    ending_bt: BondType,
    cur_len: usize,
    max_len: usize,
    last_bond: Option<usize>,
    path: &mut Vec<usize>,
    visited: &mut HashSet<usize>,
) -> Option<usize> {
    if last_bond.is_none() {
        visited.clear();
        path.clear();
    }
    visited.insert(current);
    if let Some(lb) = last_bond
        && m.atoms[current].anum == anum
        && m.bonds[lb].bt == ending_bt
        && m.atoms[current].charge == charge
    {
        if path.is_empty() || path.len() > cur_len {
            path.clear();
            path.push(lb);
            return Some(current);
        }
        return None;
    }
    if max_len <= cur_len {
        return None;
    }
    let mut target = None;
    for &b in &m.atom_bonds[current] {
        let nbr = m.bonds[b].other(current);
        if visited.contains(&nbr) {
            continue;
        }
        let bt = m.bonds[b].bt;
        if Some(bt) == next_bt {
            let next = if next_bt == Some(BondType::Single) {
                BondType::Double
            } else {
                BondType::Single
            };
            if let Some(t) = find_alternating_bonds(
                m,
                nbr,
                anum,
                charge,
                Some(next),
                ending_bt,
                cur_len + 1,
                max_len,
                Some(b),
                path,
                visited,
            ) {
                target = Some(t);
            }
        } else if ending_bt != BondType::Single
            && ending_bt != BondType::Double
            && bt == ending_bt
            && let Some(t) = find_alternating_bonds(
                m,
                nbr,
                anum,
                charge,
                None,
                ending_bt,
                cur_len + 1,
                0,
                Some(b),
                path,
                visited,
            )
        {
            target = Some(t);
        }
    }
    if target.is_some()
        && let Some(lb) = last_bond
    {
        path.push(lb);
    }
    target
}

fn find_alt(
    m: &Mol,
    start: usize,
    anum: u32,
    charge: i32,
    next_bt: BondType,
    ending_bt: BondType,
    max_len: usize,
) -> Option<(usize, Vec<usize>)> {
    let mut path = Vec::new();
    let mut visited = HashSet::new();
    find_alternating_bonds(
        m,
        start,
        anum,
        charge,
        Some(next_bt),
        ending_bt,
        0,
        max_len,
        None,
        &mut path,
        &mut visited,
    )
    .map(|t| (t, path))
}

/// Swap single and double bonds along `path`.
fn flip_path(m: &mut Mol, path: &[usize]) {
    for &b in path.iter().rev() {
        m.bonds[b].bt = if m.bonds[b].bt == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
    }
}

/// `SubstructMatch(mol, query)` for a small plain query (atoms by atomic
/// number, bonds by type, `None`: any type), unique by atom set.
fn substruct_matches(
    m: &Mol,
    q_atoms: &[u32],
    q_bonds: &[(usize, usize, Option<BondType>)],
) -> Vec<Vec<usize>> {
    fn extend(
        m: &Mol,
        q_atoms: &[u32],
        q_bonds: &[(usize, usize, Option<BondType>)],
        map: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
        seen: &mut HashSet<Vec<usize>>,
    ) {
        let k = map.len();
        if k == q_atoms.len() {
            let mut key = map.clone();
            key.sort_unstable();
            if seen.insert(key) {
                out.push(map.clone());
            }
            return;
        }
        for t in 0..m.atoms.len() {
            if m.atoms[t].anum != q_atoms[k] || map.contains(&t) {
                continue;
            }
            let ok = q_bonds.iter().all(|&(a, b, bt)| {
                let (x, y) = if a == k && b < k {
                    (t, map[b])
                } else if b == k && a < k {
                    (t, map[a])
                } else {
                    return true;
                };
                match m.bond_between(x, y) {
                    Some(bi) => bt.is_none_or(|bt| m.bonds[bi].bt == bt),
                    None => false,
                }
            });
            if ok {
                map.push(t);
                extend(m, q_atoms, q_bonds, map, out, seen);
                map.pop();
            }
        }
    }
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    extend(m, q_atoms, q_bonds, &mut Vec::new(), &mut out, &mut seen);
    out
}

fn set_bt(m: &mut Mol, a: usize, b: usize, bt: BondType) {
    if let Some(x) = m.bond_between(a, b) {
        m.bonds[x].bt = bt;
    }
}

fn ev(m: &mut Mol, a: usize) -> i32 {
    m.calc_explicit_valence(a, false).unwrap_or(-1)
}

use BondType::{Double as D, Single as S, Triple as T};

/// `_Valence4NCleanUp1`: C1=NN=[N-]=N1.
fn v4n_1(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 7 || m.atoms[a].charge != -1 || ev(m, a) != 4 {
        return false;
    }
    m.atoms[a].anum = 50;
    m.atoms[a].charge = 0;
    let mt = substruct_matches(
        m,
        &[6, 7, 50, 7, 7],
        &[
            (0, 1, Some(S)),
            (1, 2, Some(D)),
            (2, 3, Some(D)),
            (3, 4, Some(S)),
            (4, 0, Some(D)),
        ],
    );
    if mt.len() != 1 {
        m.atoms[a].anum = 7;
        m.atoms[a].charge = -1;
        return false;
    }
    let p = &mt[0];
    set_bt(m, p[0], p[1], D);
    set_bt(m, p[1], p[2], S);
    set_bt(m, p[2], p[3], S);
    set_bt(m, p[3], p[4], D);
    set_bt(m, p[4], p[0], S);
    m.atoms[a].anum = 7;
    m.atoms[a].charge = -1;
    true
}

fn v4n_2(m: &mut Mol, a: usize) -> bool {
    let Some((t, path)) = find_alt(m, a, 7, 0, D, D, 1) else {
        return false;
    };
    let top = *path.last().expect("path");
    m.bonds[top].bt = S;
    m.atoms[a].charge = 0;
    m.atoms[t].charge = -1;
    true
}

fn v5n_1(m: &mut Mol, a: usize) -> bool {
    let Some((t, path)) = find_alt(m, a, 7, 1, D, D, 5) else {
        return false;
    };
    m.atoms[t].charge = 0;
    ev(m, t);
    flip_path(m, &path);
    m.atoms[a].charge = 1;
    true
}

fn v5n_2(m: &mut Mol, a: usize) -> bool {
    let Some((t, mut path)) = find_alt(m, a, 7, -1, T, S, 2) else {
        return false;
    };
    let b = path.pop().expect("path");
    m.bonds[b].bt = S;
    let other = m.bonds[b].other(a);
    m.atoms[other].charge = -1;
    let b2 = *path.last().expect("path");
    m.bonds[b2].bt = D;
    m.atoms[t].charge = 0;
    ev(m, t);
    ev(m, a);
    true
}

fn v5n_3(m: &mut Mol, a: usize) -> bool {
    let Some((t, path)) = find_alt(m, a, 7, 0, D, D, 1) else {
        return false;
    };
    if find_alt(m, a, 8, 0, D, D, 1).is_none() {
        m.atoms[t].charge = -1;
        ev(m, t);
        let top = *path.last().expect("path");
        m.bonds[top].bt = S;
        m.atoms[a].charge = 1;
        ev(m, a);
    }
    true
}

fn v5n_4(m: &mut Mol, a: usize) -> bool {
    let mut found: Vec<(usize, usize)> = Vec::new();
    for &b in &m.atom_bonds[a] {
        let nbr = m.bonds[b].other(a);
        if m.atoms[nbr].anum == 14 && m.atoms[nbr].charge == -1 && m.bonds[b].bt == D {
            if found.len() >= 2 {
                return false;
            }
            found.push((nbr, b));
        }
    }
    if found.len() != 2 {
        return false;
    }
    for (nbr, b) in found {
        m.atoms[nbr].charge = 0;
        m.bonds[b].bt = S;
    }
    true
}

fn v5n_5(m: &mut Mol, a: usize, anum: u32) -> bool {
    // RDKit shares one visited set between the two searches.
    let mut visited = HashSet::new();
    let mut path_u = Vec::new();
    let mut path_c = Vec::new();
    let unch = find_alternating_bonds(
        m,
        a,
        anum,
        0,
        Some(D),
        D,
        0,
        7,
        None,
        &mut path_u,
        &mut visited,
    );
    let chg = find_alternating_bonds(
        m,
        a,
        anum,
        1,
        Some(D),
        D,
        0,
        7,
        None,
        &mut path_c,
        &mut visited,
    );
    if unch.is_none() && chg.is_none() {
        return false;
    }
    let path = if unch.is_none() { path_c } else { path_u };
    if let (Some(_), Some(c)) = (unch, chg) {
        m.atoms[c].charge = 0;
        m.atoms[c].num_explicit_hs = 0;
    }
    m.atoms[a].charge = 1;
    flip_path(m, &path);
    match (unch, chg) {
        (Some(u), Some(_)) => m.atoms[u].num_explicit_hs = 1,
        (Some(u), None) => m.atoms[u].charge = -1,
        (None, Some(c)) => m.atoms[c].charge = 0,
        (None, None) => {}
    }
    if let Some(c) = chg {
        ev(m, c);
    }
    if let Some(u) = unch {
        ev(m, u);
    }
    true
}

fn v5n_6(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 7 || m.atoms[a].charge != 0 || ev(m, a) != 5 {
        return false;
    }
    m.atoms[a].anum = 50;
    let mt = substruct_matches(
        m,
        &[6, 6, 50, 6, 6, 7, 6],
        &[
            (0, 1, Some(S)),
            (1, 2, Some(D)),
            (2, 3, Some(D)),
            (3, 4, None),
            (4, 5, Some(S)),
            (5, 0, Some(D)),
            (2, 6, Some(S)),
        ],
    );
    if mt.len() != 1 {
        m.atoms[a].anum = 7;
        return false;
    }
    let p = &mt[0];
    set_bt(m, p[0], p[1], D);
    set_bt(m, p[1], p[2], S);
    set_bt(m, p[4], p[5], D);
    set_bt(m, p[5], p[0], S);
    m.atoms[a].anum = 7;
    m.atoms[a].charge = 1;
    true
}

fn v5n_7(m: &mut Mol, a: usize) -> bool {
    let Some((t, path)) = find_alt(m, a, 8, 0, D, D, 5) else {
        return false;
    };
    if m.atoms[a].anum != 7 || m.atoms[a].charge != 0 || ev(m, a) != 5 {
        return false;
    }
    m.atoms[a].anum = 50;
    let mt = substruct_matches(
        m,
        &[6, 6, 50, 7, 6, 8, 6],
        &[
            (0, 1, None),
            (1, 2, Some(D)),
            (2, 3, Some(D)),
            (3, 4, Some(S)),
            (4, 5, Some(S)),
            (5, 0, Some(S)),
            (2, 6, Some(S)),
        ],
    );
    if mt.len() != 1 {
        m.atoms[a].anum = 7;
        return false;
    }
    let p = &mt[0];
    set_bt(m, p[1], p[2], S);
    flip_path(m, &path);
    m.atoms[t].charge = -1;
    m.atoms[a].anum = 7;
    true
}

fn v5n_8(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 7 || m.atoms[a].charge != 0 || ev(m, a) != 5 {
        return false;
    }
    m.atoms[a].anum = 50;
    let mt = substruct_matches(
        m,
        &[6, 7, 6, 7, 7, 50],
        &[
            (0, 1, Some(S)),
            (1, 2, Some(D)),
            (2, 3, Some(S)),
            (3, 4, Some(D)),
            (4, 0, Some(S)),
            (5, 0, Some(D)),
        ],
    );
    if mt.len() != 1 {
        m.atoms[a].anum = 7;
        return false;
    }
    let p = &mt[0];
    set_bt(m, p[1], p[2], S);
    set_bt(m, p[2], p[3], D);
    set_bt(m, p[3], p[4], S);
    set_bt(m, p[4], p[0], D);
    set_bt(m, p[5], p[0], S);
    m.atoms[p[1]].charge = -1;
    m.atoms[a].anum = 7;
    m.atoms[a].charge = 1;
    true
}

fn v5n_9(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 7 || m.atoms[a].charge != 0 || ev(m, a) != 5 {
        return false;
    }
    m.atoms[a].anum = 50;
    let mt = substruct_matches(
        m,
        &[6, 7, 7, 6, 6, 50],
        &[
            (0, 1, Some(S)),
            (1, 2, Some(D)),
            (2, 3, Some(S)),
            (3, 4, Some(D)),
            (4, 0, Some(S)),
            (5, 0, Some(D)),
        ],
    );
    if mt.len() != 1 {
        m.atoms[a].anum = 7;
        return false;
    }
    let p = &mt[0];
    set_bt(m, p[0], p[1], D);
    set_bt(m, p[1], p[2], S);
    set_bt(m, p[5], p[0], S);
    m.atoms[p[2]].charge = -1;
    m.atoms[a].anum = 7;
    m.atoms[a].charge = 1;
    true
}

fn v5n_a(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 7 || m.atoms[a].charge != 0 || ev(m, a) != 5 {
        return false;
    }
    let matches = substruct_matches(m, &[7, 7], &[(0, 1, Some(D))]);
    if matches.is_empty() {
        return false;
    }
    let mut best: Option<Vec<usize>> = None;
    for mt in matches {
        if mt[0] == a || mt[1] == a {
            continue;
        }
        m.atoms[mt[0]].anum = 50;
        m.atoms[mt[1]].anum = 50;
        if let Some((_, path)) = find_alt(m, a, 50, 0, D, D, 9)
            && best.as_ref().is_none_or(|b| path.len() < b.len())
        {
            best = Some(path);
        }
        m.atoms[mt[0]].anum = 7;
        m.atoms[mt[1]].anum = 7;
    }
    match best {
        Some(path) if !path.is_empty() => {
            flip_path(m, &path);
            m.atoms[a].charge = 1;
            ev(m, a);
            true
        }
        _ => false,
    }
}

fn v5n_b(m: &mut Mol, a: usize) -> bool {
    let Some((t, path)) = find_alt(m, a, 6, 0, D, D, 1) else {
        return false;
    };
    m.atoms[t].charge = -1;
    ev(m, t);
    let top = *path.last().expect("path");
    m.bonds[top].bt = S;
    m.atoms[a].charge = 1;
    ev(m, a);
    true
}

fn v7s_1(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 16 || m.atoms[a].charge != -1 || ev(m, a) != 7 {
        return false;
    }
    let (mut n_c, mut n_o) = (0, 0);
    let mut last_o: Option<usize> = None;
    for &b in &m.atom_bonds[a] {
        let o = m.bonds[b].other(a);
        match m.atoms[o].anum {
            8 => {
                if m.bonds[b].bt != D {
                    n_o = 100;
                    break;
                }
                last_o = Some(o);
                n_o += 1;
            }
            6 => {
                if m.bonds[b].bt != S {
                    n_c = 100;
                    break;
                }
                n_c += 1;
            }
            _ => {
                n_c = 100;
                break;
            }
        }
    }
    match last_o {
        Some(o) if n_c == 1 || n_o == 3 => {
            set_bt(m, o, a, S);
            m.atoms[o].charge = -1;
            m.atoms[a].charge = 0;
            ev(m, o);
            ev(m, a);
            true
        }
        _ => false,
    }
}

fn v7s_2(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 16 || m.atoms[a].charge != -1 || ev(m, a) != 7 {
        return false;
    }
    let Some((_, path)) = find_alt(m, a, 7, 0, D, T, 3) else {
        return false;
    };
    for &b in path.iter().rev() {
        m.bonds[b].bt = match m.bonds[b].bt {
            S => D,
            D => S,
            T => D,
            other => other,
        };
    }
    m.atoms[a].charge = 0;
    ev(m, a);
    true
}

fn v7s_3(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 16 || m.atoms[a].charge != -1 || ev(m, a) != 7 {
        return false;
    }
    let Some((t, path)) = find_alt(m, a, 7, 0, D, D, 1) else {
        return false;
    };
    let top = *path.last().expect("path");
    m.bonds[top].bt = S;
    m.atoms[t].charge = -1;
    m.atoms[a].charge = 0;
    ev(m, a);
    true
}

fn v8s_1(m: &mut Mol, a: usize) -> bool {
    if m.atoms[a].anum != 16 || m.atoms[a].charge != -1 || ev(m, a) != 7 {
        return false;
    }
    let Some((t, path)) = find_alt(m, a, 7, 0, D, D, 9) else {
        return false;
    };
    flip_path(m, &path);
    m.atoms[t].charge = -1;
    ev(m, t);
    m.atoms[t].num_explicit_hs = 0;
    m.atoms[a].charge = 0;
    ev(m, a);
    true
}

fn v8cl_1(m: &mut Mol, a: usize) -> bool {
    if ev(m, a) != 8 || m.atoms[a].charge != -1 {
        return false;
    }
    if !m.nbrs(a).all(|o| m.atoms[o].anum == 8) {
        return false;
    }
    m.atoms[a].charge = 3;
    for b in m.atom_bonds[a].clone() {
        if m.bonds[b].bt == D {
            m.bonds[b].bt = S;
            let o = m.bonds[b].other(a);
            m.atoms[o].charge = -1;
            ev(m, o);
        }
    }
    ev(m, a);
    true
}

fn v5cl_1(m: &mut Mol, a: usize) -> bool {
    if ev(m, a) != 6 || m.atoms[a].charge != 1 {
        return false;
    }
    let Some((t, path)) = find_alt(m, a, 8, -1, S, S, 1) else {
        return false;
    };
    let top = *path.last().expect("path");
    m.bonds[top].bt = D;
    m.atoms[a].charge = 0;
    m.atoms[t].charge = 0;
    ev(m, a);
    true
}

fn v3cl_1(m: &mut Mol, a: usize) -> bool {
    if ev(m, a) != 3 || m.atoms[a].charge != 0 {
        return false;
    }
    let Some((_, path)) = find_alt(m, a, 16, 0, T, T, 1) else {
        return false;
    };
    let top = *path.last().expect("path");
    m.bonds[top].bt = S;
    ev(m, a);
    true
}

/// `cleanUp` of RDKit's InChI reader (not `MolOps::cleanUp`).
fn inchi_clean_up(m: &mut Mol) {
    for a in 0..m.atoms.len() {
        match m.atoms[a].anum {
            7 => {
                if ev(m, a) == 4 {
                    if v4n_1(m, a) {
                        continue;
                    }
                    if m.atoms[a].charge == -1 {
                        v4n_2(m, a);
                    }
                    continue;
                }
                if m.atoms[a].charge != 0 {
                    continue;
                }
                let arom = m.atoms[a].aromatic;
                m.atoms[a].aromatic = false;
                if ev(m, a) == 5 {
                    let _ = v5n_6(m, a)
                        || v5n_7(m, a)
                        || v5n_8(m, a)
                        || v5n_9(m, a)
                        || v5n_a(m, a)
                        || v5n_1(m, a)
                        || v5n_2(m, a)
                        || v5n_3(m, a)
                        || v5n_4(m, a)
                        || v5n_5(m, a, 8)
                        || v5n_5(m, a, 16)
                        || v5n_5(m, a, 9)
                        || v5n_5(m, a, 17)
                        || v5n_b(m, a);
                }
                if arom {
                    m.atoms[a].aromatic = true;
                }
            }
            17 => {
                if ev(m, a) == 8 && v8cl_1(m, a) {
                    continue;
                }
                if ev(m, a) == 5 && v5cl_1(m, a) {
                    continue;
                }
                if ev(m, a) == 3 {
                    v3cl_1(m, a);
                }
            }
            16 => {
                let e = ev(m, a);
                if e == 7 {
                    if v7s_1(m, a) || v7s_2(m, a) || v7s_3(m, a) {
                        continue;
                    }
                    v8s_1(m, a);
                } else if e == 8 {
                    v8s_1(m, a);
                }
            }
            35 => {
                let e = ev(m, a);
                if e == 3 && m.atoms[a].charge == 0 && m.degree(a) == 1 {
                    let b = m.atom_bonds[a][0];
                    let o = m.bonds[b].other(a);
                    if m.atoms[o].anum == 34 {
                        m.bonds[b].bt = S;
                    }
                }
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(atoms: &[(u32, i32)], bonds: &[(usize, usize, BondType)]) -> Mol {
        let mut mol = Mol::default();
        for &(anum, charge) in atoms {
            let mut atom = Atom::new(anum);
            atom.charge = charge;
            atom.no_implicit = true;
            mol.add_atom(atom);
        }
        for &(a, b, bt) in bonds {
            mol.add_bond(Bond::new(a, b, bt));
        }
        mol
    }

    #[test]
    fn bond_direction_constraints_detect_conflicts() {
        let mut mol = graph(&[(6, 0); 4], &[(0, 1, S), (1, 2, S), (2, 3, S)]);
        assert!(assign_bond_dirs(&mut mol, &[(0, 1)], &[(1, 2)]));
        assert_eq!(mol.bonds[0].dir, mol.bonds[1].dir);
        assert_ne!(mol.bonds[1].dir, mol.bonds[2].dir);

        let mut conflict = graph(&[(6, 0); 3], &[(0, 1, S), (1, 2, S)]);
        assert!(!assign_bond_dirs(&mut conflict, &[(0, 1)], &[(0, 1)]));
    }

    #[test]
    fn alternating_path_and_substructure_helpers_are_deterministic() {
        let mut mol = graph(
            &[(7, 0), (6, 0), (7, 1), (8, 0)],
            &[(0, 1, D), (1, 2, S), (1, 3, S)],
        );
        let (target, path) = find_alt(&mol, 0, 7, 1, D, S, 3).unwrap();
        assert_eq!(target, 2);
        assert_eq!(path.len(), 2);
        flip_path(&mut mol, &path);
        assert_eq!(mol.bonds[0].bt, S);
        assert_eq!(mol.bonds[1].bt, D);

        let matches = substruct_matches(&mol, &[7, 6, 8], &[(0, 1, Some(S)), (1, 2, None)]);
        assert_eq!(matches, vec![vec![0, 1, 3]]);
        assert!(find_alt(&mol, 3, 17, 0, S, S, 1).is_none());
    }

    #[test]
    fn inchi_cleanup_normalizes_common_hypervalent_groups() {
        let mut perchlorate = graph(
            &[(17, -1), (8, 0), (8, 0), (8, 0), (8, 0)],
            &[(0, 1, D), (0, 2, D), (0, 3, D), (0, 4, D)],
        );
        inchi_clean_up(&mut perchlorate);
        assert_eq!(perchlorate.atoms[0].charge, 3);
        assert!(perchlorate.atoms[1..].iter().all(|a| a.charge == -1));
        assert!(perchlorate.bonds.iter().all(|b| b.bt == S));

        let mut sulfonate = graph(
            &[(16, -1), (6, 0), (8, 0), (8, 0), (8, 0)],
            &[(0, 1, S), (0, 2, D), (0, 3, D), (0, 4, D)],
        );
        inchi_clean_up(&mut sulfonate);
        assert_eq!(sulfonate.atoms[0].charge, 0);
        assert_eq!(sulfonate.atoms.iter().filter(|a| a.charge == -1).count(), 1);

        let mut azide_like = graph(
            &[(7, -1), (7, 0), (6, 0), (6, 0)],
            &[(0, 1, D), (0, 2, S), (0, 3, S)],
        );
        inchi_clean_up(&mut azide_like);
        assert_eq!(azide_like.atoms[0].charge, 0);
        assert_eq!(azide_like.atoms[1].charge, -1);
        assert_eq!(azide_like.bonds[0].bt, S);

        let mut selenium_bromide = graph(&[(35, 0), (34, 0)], &[(0, 1, T)]);
        inchi_clean_up(&mut selenium_bromide);
        assert_eq!(selenium_bromide.bonds[0].bt, S);
    }

    #[test]
    fn inchi_output_preserves_isotopes_radicals_and_stereo() {
        let alkene = [
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(1, 1, 0)],
                num_iso_h: [3, 0, 0, 0],
                ..Default::default()
            },
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(0, 1, 0), (2, 2, 0)],
                num_iso_h: [1, 0, 0, 0],
                ..Default::default()
            },
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(1, 2, 0), (3, 1, 0)],
                num_iso_h: [1, 0, 0, 0],
                ..Default::default()
            },
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(2, 1, 0)],
                num_iso_h: [3, 0, 0, 0],
                ..Default::default()
            },
        ];
        let stereo = [InchiOutputStereo0D {
            neighbor: [0, 1, 2, 3],
            central_atom: -1,
            stereo_type: INCHI_STEREO_TYPE_DOUBLE_BOND,
            parity: INCHI_PARITY_ODD,
        }];
        let mol = rdkit_molecule_from_inchi_output(&alkene, &stereo).unwrap();
        let written = crate::canonical_smiles(&mol);
        assert!(written.contains('/') || written.contains('\\'));

        let isotopic = [InchiOutputAtom {
            element: "C".into(),
            isotopic_mass: (ISOTOPIC_SHIFT_FLAG + 1) as i16,
            radical: 3,
            num_iso_h: [1, 0, 1, 0],
            ..Default::default()
        }];
        let mol = rdkit_molecule_from_inchi_output(&isotopic, &[]).unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert!(mol.atoms().any(|(_, a)| a.isotope == Some(2)));
    }

    #[test]
    fn inchi_output_rejects_bond_stereo_and_bad_tetrahedral_neighbors() {
        let bad_bond_stereo = [
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(1, 1, 1)],
                ..Default::default()
            },
            InchiOutputAtom {
                element: "C".into(),
                ..Default::default()
            },
        ];
        assert!(rdkit_molecule_from_inchi_output(&bad_bond_stereo, &[]).is_err());

        let tetra = [
            InchiOutputAtom {
                element: "C".into(),
                bonds: vec![(1, 1, 0)],
                ..Default::default()
            },
            InchiOutputAtom {
                element: "F".into(),
                bonds: vec![(0, 1, 0)],
                ..Default::default()
            },
        ];
        let stereo = [InchiOutputStereo0D {
            neighbor: [0, 1, 2, 3],
            central_atom: 0,
            stereo_type: INCHI_STEREO_TYPE_TETRAHEDRAL,
            parity: INCHI_PARITY_EVEN,
        }];
        assert!(rdkit_molecule_from_inchi_output(&tetra, &stereo).is_err());
    }
}
