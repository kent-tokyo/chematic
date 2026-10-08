//! `Chirality::findPotentialStereo(mol, cleanIt=false, flagPossible=true)`
//! (RDKit 2026.03.1 `FindStereo.cpp`), as `Chem.FindPotentialStereo(mol)`
//! calls it: the new stereo perception's iterative symmetry-class
//! refinement, with every unspecified potential stereo atom/bond made
//! distinguishable.

use chematic_perception::{RdkitRankAtom, RdkitRankBond, rdkit_rank_fragment_atoms_with_symbols};

use super::mol::{BondStereo, BondType, ChiralTag, Mol, count_swaps};
use super::periodic;
use super::stereo::{atom_nonzero_degree, is_atom_potential_tetrahedral_center};

/// `Atom::NOATOM`.
pub(crate) const NOATOM: usize = usize::MAX;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StereoType {
    AtomTetrahedral,
    AtomTrigonalBipyramidal,
    AtomOctahedral,
    BondDouble,
}

/// `Chirality::StereoSpecified`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub(crate) enum Specified {
    Unspecified,
    Specified,
    Unknown,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Descriptor {
    None,
    TetCw,
    TetCcw,
    BondCis,
    BondTrans,
}

#[derive(Clone, Debug)]
pub(crate) struct StereoInfo {
    pub kind: StereoType,
    pub centered_on: usize,
    pub controlling: Vec<usize>,
    pub specified: Specified,
    pub descriptor: Descriptor,
}

/// `isAtomPotentialNontetrahedralCenter` (no non-tetrahedral tags exist in
/// the model).
fn is_atom_potential_nontetrahedral_center(mol: &Mol, a: usize) -> bool {
    let tnz = atom_nonzero_degree(mol, a) + mol.total_num_hs(a) as usize;
    let anum = mol.atoms[a].anum;
    if !(2..=6).contains(&tnz) || (anum < 12 && anum != 4) {
        return false;
    }
    mol.atoms[a].chiral == ChiralTag::Unspecified && tnz >= 4
}

fn is_atom_potential_stereo_atom(mol: &Mol, a: usize) -> bool {
    is_atom_potential_tetrahedral_center(mol, a) || is_atom_potential_nontetrahedral_center(mol, a)
}

/// `getTotalNumHs(includeNeighbors=true)`.
fn total_num_hs_with_nbrs(mol: &Mol, a: usize) -> usize {
    mol.total_num_hs(a) as usize + mol.nbrs(a).filter(|&nb| mol.atoms[nb].anum == 1).count()
}

/// `isBondPotentialStereoBond`.
fn is_bond_potential_stereo_bond(mol: &Mol, b: usize) -> bool {
    let bond = &mol.bonds[b];
    if bond.bt != BondType::Double {
        return false;
    }
    let (beg, end) = (bond.begin, bond.end);
    let (bd, ed) = (mol.total_degree(beg), mol.total_degree(end));
    if bd > 1
        && bd < 4
        && ed > 1
        && ed < 4
        && total_num_hs_with_nbrs(mol, beg) < 2
        && total_num_hs_with_nbrs(mol, end) < 2
    {
        let ri = mol.ring_info();
        !ri.bond_rings
            .iter()
            .any(|bring| bring.len() < 8 && bring.contains(&b))
    } else {
        false
    }
}

/// `getStereoInfo(const Atom *)`.
fn atom_stereo_info(mol: &Mol, a: usize) -> StereoInfo {
    let orig: Vec<usize> = mol.nbrs(a).collect();
    let mut controlling = orig.clone();
    controlling.sort_unstable();
    let mut info = StereoInfo {
        kind: StereoType::AtomTetrahedral,
        centered_on: a,
        controlling,
        specified: Specified::Unspecified,
        descriptor: Descriptor::None,
    };
    match mol.atoms[a].chiral {
        ChiralTag::Cw | ChiralTag::Ccw => {
            info.specified = Specified::Specified;
            let mut cw = mol.atoms[a].chiral == ChiralTag::Cw;
            if count_swaps(&orig, &info.controlling) % 2 == 1 {
                cw = !cw;
            }
            info.descriptor = if cw {
                Descriptor::TetCw
            } else {
                Descriptor::TetCcw
            };
        }
        ChiralTag::Unspecified => {
            if is_atom_potential_nontetrahedral_center(mol, a) {
                info.kind = match mol.total_degree(a) {
                    5 => StereoType::AtomTrigonalBipyramidal,
                    6 => StereoType::AtomOctahedral,
                    _ => StereoType::AtomTetrahedral,
                };
            }
        }
    }
    info
}

/// `getStereoInfo(const Bond *)` for a double bond.
fn bond_stereo_info(mol: &Mol, b: usize) -> StereoInfo {
    let bond = &mol.bonds[b];
    let mut controlling = Vec::with_capacity(4);
    for atom in [bond.begin, bond.end] {
        for &nb in &mol.atom_bonds[atom] {
            if nb != b {
                controlling.push(mol.bonds[nb].other(atom));
            }
        }
        let pad = 3usize.saturating_sub(mol.degree(atom));
        controlling.extend(std::iter::repeat_n(NOATOM, pad));
    }
    let mut info = StereoInfo {
        kind: StereoType::BondDouble,
        centered_on: b,
        controlling,
        specified: Specified::Unspecified,
        descriptor: Descriptor::None,
    };
    match bond.stereo {
        BondStereo::Any => info.specified = Specified::Unknown,
        BondStereo::E | BondStereo::Z if bond.stereo_atoms.len() == 2 => {
            info.specified = Specified::Specified;
            let mut cis = bond.stereo == BondStereo::Z;
            let first_at_begin = bond.stereo_atoms[0] == info.controlling[0];
            let first_at_end = bond.stereo_atoms[1] == info.controlling[2];
            if first_at_begin ^ first_at_end {
                cis = !cis;
            }
            info.descriptor = if cis {
                Descriptor::BondCis
            } else {
                Descriptor::BondTrans
            };
        }
        _ => {}
    }
    info
}

/// `getBondSymbol`.
fn bond_symbol(mol: &Mol, b: usize) -> String {
    let bond = &mol.bonds[b];
    if bond.aromatic {
        return ":".into();
    }
    match bond.bt {
        BondType::Single => "-",
        BondType::Double => "=",
        BondType::Triple => "#",
        BondType::Aromatic => ":",
        _ => "?",
    }
    .into()
}

/// `getAtomCompareSymbol`.
fn atom_compare_symbol(mol: &Mol, a: usize) -> String {
    let atom = &mol.atoms[a];
    format!(
        "{}{}{}",
        atom.isotope,
        periodic::symbol(atom.anum),
        atom.charge
    )
}

struct State {
    known_atoms: Vec<bool>,
    possible_atoms: Vec<bool>,
    fixed_atoms: Vec<bool>,
    atom_symbols: Vec<String>,
    known_bonds: Vec<bool>,
    possible_bonds: Vec<bool>,
    fixed_bonds: Vec<bool>,
    bond_symbols: Vec<String>,
    ring_atoms: Vec<u32>,
    ring_bonds: Vec<u32>,
}

/// `flagRingStereo` (with the possible sets, `cleanIt=false`).
fn flag_ring_stereo(mol: &Mol, st: &mut State) {
    let ri = mol.ring_info();
    let n = mol.atoms.len();
    let mut in_ring = vec![false; n];
    for ridx in 0..ri.atom_rings.len() {
        let aring = &ri.atom_rings[ridx];
        let bring = &ri.bond_rings[ridx];
        let mut n_here = 0u32;
        let sz = aring.len();
        let half = sz / 2 + sz % 2;
        in_ring.iter_mut().for_each(|x| *x = false);
        for ai in 0..sz {
            let aidx = aring[ai];
            if !st.known_atoms[aidx] && !st.possible_atoms[aidx] {
                continue;
            }
            // (RDKit's `continue` here only moves on to the next divisor;
            // the fused-ring exploration below always runs.)
            for div in [2usize, 3] {
                if !sz.is_multiple_of(div) {
                    continue;
                }
                let inc = sz / div;
                let mut by_bond = 0usize;
                let mut by_atom = 0usize;
                let mut k = inc;
                while k < sz {
                    let other = aring[(ai + k) % sz];
                    for &bidx in &mol.atom_bonds[other] {
                        if (st.known_bonds[bidx] || st.possible_bonds[bidx])
                            && !bring.contains(&bidx)
                        {
                            by_bond += 1;
                            break;
                        }
                    }
                    if by_bond == 0 && (st.known_atoms[other] || st.possible_atoms[other]) {
                        by_atom += 1;
                    }
                    k += inc;
                }
                if by_bond == div - 1 || by_atom == div - 1 {
                    n_here += 1 + by_bond as u32;
                    let mut k = 0;
                    while k < sz {
                        in_ring[aring[(ai + k) % sz]] = true;
                        k += inc;
                    }
                }
            }
            if ri.num_atom_rings(aidx) > 1 {
                let mut prev = aidx;
                for step in 1..=half {
                    let other = aring[(ai + step) % sz];
                    let bnd = mol.bond_between(prev, other).expect("ring bond");
                    if ri.num_bond_rings(bnd) < 2 {
                        break;
                    }
                    if st.known_atoms[other] || st.possible_atoms[other] {
                        n_here += 2;
                        in_ring[aidx] = true;
                        in_ring[other] = true;
                        break;
                    }
                    prev = other;
                }
            }
        }
        if n_here > 1 {
            for &a in aring {
                if in_ring[a] {
                    st.ring_atoms[a] += 1;
                }
            }
            for &b in bring {
                st.ring_bonds[b] += 1;
            }
        }
    }
}

/// `areStereobondControllingAtomsDupes`.
fn controlling_atoms_dupes(
    mol: &Mol,
    b: usize,
    c1: usize,
    c2: usize,
    ranks: &[u32],
    st: &State,
) -> bool {
    if ranks[c1] != ranks[c2] {
        return false;
    }
    let ri = mol.ring_info();
    let (m1, m2) = (&ri.atom_members[c1], &ri.atom_members[c2]);
    let (mut i1, mut i2) = (0, 0);
    while i1 < m1.len() && i2 < m2.len() {
        if m1[i1] < m2[i2] {
            i1 += 1;
            continue;
        } else if m1[i1] > m2[i2] {
            i2 += 1;
            continue;
        }
        let ring = &ri.atom_rings[m1[i1]];
        i1 += 1;
        i2 += 1;
        if ring.len() % 2 == 1 {
            continue;
        }
        for bond_end in [mol.bonds[b].begin, mol.bonds[b].end] {
            if let Some(pos) = ring.iter().position(|&x| x == bond_end) {
                let opposite = ring[(pos + ring.len() / 2) % ring.len()];
                if st.possible_atoms[opposite] || st.known_atoms[opposite] {
                    return false;
                }
                if mol.degree(opposite) == 3 {
                    for &nb in &mol.atom_bonds[opposite] {
                        let out = mol.bonds[nb].other(opposite);
                        if !ring.contains(&out) && (st.possible_bonds[nb] || st.known_bonds[nb]) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

fn update_atoms(mol: &Mol, ranks: &[u32], st: &mut State, res: &mut Vec<StereoInfo>) -> bool {
    let mut need = false;
    let ri = mol.ring_info();
    for aidx in 0..mol.atoms.len() {
        if !(st.known_atoms[aidx] || st.possible_atoms[aidx]) {
            continue;
        }
        let mut info = atom_stereo_info(mol, aidx);
        if st.fixed_atoms[aidx] {
            res.push(info);
            continue;
        }
        let mut nbrs: Vec<u32> = Vec::with_capacity(info.controlling.len());
        let mut dupe = false;
        if info.kind == StereoType::AtomTetrahedral {
            for &nbr in &info.controlling {
                let rnk = ranks[nbr];
                if nbrs.contains(&rnk) {
                    if st.ring_atoms[aidx] != 0 {
                        match mol.bond_between(aidx, nbr) {
                            Some(b) if st.ring_bonds[b] != 0 => {}
                            _ => {
                                dupe = true;
                                break;
                            }
                        }
                    } else {
                        dupe = true;
                        break;
                    }
                } else {
                    nbrs.push(rnk);
                }
            }
        }
        if !dupe {
            let mut acs = st.atom_symbols[aidx].clone();
            if !st.possible_atoms[aidx] {
                let mut sorted = nbrs.clone();
                sorted.sort_unstable();
                if info.kind == StereoType::AtomTetrahedral {
                    let swaps = count_swaps_u32(&nbrs, &sorted);
                    if swaps % 2 == 1 {
                        info.descriptor = match info.descriptor {
                            Descriptor::TetCcw => Descriptor::TetCw,
                            Descriptor::TetCw => Descriptor::TetCcw,
                            d => d,
                        };
                    }
                    match info.descriptor {
                        Descriptor::TetCw => acs = atom_compare_symbol(mol, aidx) + "_CW",
                        Descriptor::TetCcw => acs = atom_compare_symbol(mol, aidx) + "_CCW",
                        _ => {}
                    }
                }
                st.fixed_atoms[aidx] = true;
            }
            if st.atom_symbols[aidx] != acs {
                st.atom_symbols[aidx] = acs;
                need = true;
            }
            res.push(info);
        } else {
            need |= st.possible_atoms[aidx];
            st.possible_atoms[aidx] = false;
            st.atom_symbols[aidx] = atom_compare_symbol(mol, aidx);
            if st.ring_atoms[aidx] != 0 {
                st.ring_atoms[aidx] = 0;
                need = true;
                for ridx in 0..ri.atom_rings.len() {
                    let aring = &ri.atom_rings[ridx];
                    let mut n_here = 0;
                    for &ra in aring {
                        st.fixed_atoms[ra] = false;
                        n_here += (st.ring_atoms[ra] > 0) as u32;
                    }
                    if n_here <= 1 {
                        if n_here == 1 {
                            for &ra in aring {
                                if st.ring_atoms[ra] != 0 {
                                    st.ring_atoms[ra] -= 1;
                                    break;
                                }
                            }
                        }
                        for &rb in &ri.bond_rings[ridx] {
                            if st.ring_bonds[rb] != 0 {
                                st.ring_bonds[rb] -= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    need
}

fn count_swaps_u32(reference: &[u32], probe: &[u32]) -> usize {
    let r: Vec<usize> = reference.iter().map(|&x| x as usize).collect();
    let p: Vec<usize> = probe.iter().map(|&x| x as usize).collect();
    count_swaps(&r, &p)
}

fn update_bonds(mol: &Mol, ranks: &[u32], st: &mut State, res: &mut Vec<StereoInfo>) -> bool {
    let mut need = false;
    for bidx in 0..mol.bonds.len() {
        if !(st.known_bonds[bidx] || st.possible_bonds[bidx]) {
            continue;
        }
        if mol.bonds[bidx].bt != BondType::Double {
            continue;
        }
        let mut info = bond_stereo_info(mol, bidx);
        let c = &info.controlling;
        if (c[0] == NOATOM && c[1] == NOATOM) || (c[2] == NOATOM && c[3] == NOATOM) {
            st.fixed_bonds[bidx] = true;
        }
        if st.fixed_bonds[bidx] {
            res.push(info);
            continue;
        }
        let mut dupe = false;
        let mut needs_swap = false;
        for (x, y) in [(0usize, 1usize), (2, 3)] {
            let (cx, cy) = (info.controlling[x], info.controlling[y]);
            if cx != NOATOM && cy != NOATOM {
                if controlling_atoms_dupes(mol, bidx, cx, cy, ranks, st) {
                    dupe = true;
                } else if ranks[cx] < ranks[cy] {
                    info.controlling.swap(x, y);
                    needs_swap = !needs_swap;
                }
            }
        }
        if !dupe {
            if needs_swap {
                info.descriptor = match info.descriptor {
                    Descriptor::BondCis => Descriptor::BondTrans,
                    Descriptor::BondTrans => Descriptor::BondCis,
                    d => d,
                };
            }
            let mut gbs = st.bond_symbols[bidx].clone();
            match info.specified {
                Specified::Specified => match info.descriptor {
                    Descriptor::BondCis => gbs.push_str("_cis"),
                    Descriptor::BondTrans => gbs.push_str("_trans"),
                    _ => {}
                },
                Specified::Unknown => gbs.push_str("_unk"),
                Specified::Unspecified => {}
            }
            if st.bond_symbols[bidx] != gbs {
                st.bond_symbols[bidx] = gbs;
                need = true;
            }
            if !st.possible_bonds[bidx] {
                st.fixed_bonds[bidx] = true;
            }
            res.push(info);
        } else if st.possible_bonds[bidx] {
            st.possible_bonds[bidx] = false;
            st.bond_symbols[bidx] = bond_symbol(mol, bidx);
            need = true;
        }
    }
    need
}

fn ranks(mol: &Mol, st: &State) -> Vec<u32> {
    let atoms: Vec<RdkitRankAtom> = (0..mol.atoms.len())
        .map(|i| RdkitRankAtom {
            atomic_num: mol.atoms[i].anum,
            formal_charge: mol.atoms[i].charge,
            total_num_hs: mol.total_num_hs(i),
            num_rings: mol.num_atom_rings(i) as u32,
            ..RdkitRankAtom::default()
        })
        .collect();
    let bonds: Vec<RdkitRankBond> = mol
        .bonds
        .iter()
        .map(|b| RdkitRankBond {
            begin: b.begin as u32,
            end: b.end as u32,
            bond_type: if b.aromatic {
                BondType::Aromatic as u32
            } else {
                b.bt as u32
            },
            ..RdkitRankBond::default()
        })
        .collect();
    rdkit_rank_fragment_atoms_with_symbols(&atoms, &bonds, &st.atom_symbols, &st.bond_symbols)
}

fn refine(mol: &Mol, st: &mut State) -> (Vec<StereoInfo>, Vec<u32>) {
    let mut res = Vec::new();
    let mut aranks = Vec::new();
    let mut need = true;
    while need {
        res.clear();
        aranks = ranks(mol, st);
        need = update_atoms(mol, &aranks, st, &mut res);
        need |= update_bonds(mol, &aranks, st, &mut res);
    }
    (res, aranks)
}

/// `Chem.FindPotentialStereo(mol)` (`cleanIt=False, flagPossible=True`):
/// atoms (index order) then double bonds (index order). Bond entries carry
/// their controlling atoms ordered by RDKit's symmetry ranks.
pub(crate) fn find_potential_stereo(mol: &Mol) -> Vec<StereoInfo> {
    let n = mol.atoms.len();
    let nb = mol.bonds.len();
    let mut st = State {
        known_atoms: vec![false; n],
        possible_atoms: vec![false; n],
        fixed_atoms: vec![false; n],
        atom_symbols: vec![String::new(); n],
        known_bonds: vec![false; nb],
        possible_bonds: vec![false; nb],
        fixed_bonds: vec![false; nb],
        bond_symbols: vec![String::new(); nb],
        ring_atoms: vec![0; n],
        ring_bonds: vec![0; nb],
    };
    // initAtomInfo
    for a in 0..n {
        st.atom_symbols[a] = atom_compare_symbol(mol, a);
        if is_atom_potential_stereo_atom(mol, a) {
            let info = atom_stereo_info(mol, a);
            match info.specified {
                Specified::Unknown => {
                    st.known_atoms[a] = true;
                    st.atom_symbols[a] += &a.to_string();
                }
                Specified::Specified => {
                    st.known_atoms[a] = true;
                    st.atom_symbols[a] += match info.descriptor {
                        Descriptor::TetCcw => "_CCW",
                        Descriptor::TetCw => "_CW",
                        _ => "_STEREO",
                    };
                }
                Specified::Unspecified => {
                    st.possible_atoms[a] = true;
                    st.atom_symbols[a] += &format!("_{a}");
                }
            }
        }
    }
    // initBondInfo
    for b in 0..nb {
        st.bond_symbols[b] = bond_symbol(mol, b);
        if is_bond_potential_stereo_bond(mol, b) {
            let info = bond_stereo_info(mol, b);
            match info.specified {
                Specified::Unknown => {
                    st.known_bonds[b] = true;
                    st.bond_symbols[b] += &format!("_{b}");
                }
                Specified::Specified => {
                    st.known_bonds[b] = true;
                    st.bond_symbols[b] += match info.descriptor {
                        Descriptor::BondCis => "_cis",
                        Descriptor::BondTrans => "_trans",
                        _ => "_STEREO",
                    };
                }
                Specified::Unspecified => {
                    st.possible_bonds[b] = true;
                    st.bond_symbols[b] += &format!("_{b}");
                }
            }
        }
    }
    let orig_possible_atoms = st.possible_atoms.clone();
    let orig_possible_bonds = st.possible_bonds.clone();
    flag_ring_stereo(mol, &mut st);
    let (mut res, _) = refine(mol, &mut st);
    if st.possible_atoms != orig_possible_atoms || st.possible_bonds != orig_possible_bonds {
        st.possible_atoms = orig_possible_atoms;
        for i in 0..n {
            if !st.fixed_atoms[i] && st.known_atoms[i] {
                st.possible_atoms[i] = true;
                st.known_atoms[i] = false;
            }
            if st.possible_atoms[i] {
                st.atom_symbols[i] += &format!("_{i}");
            }
        }
        st.possible_bonds = orig_possible_bonds;
        for i in 0..nb {
            if !st.fixed_bonds[i] && st.known_bonds[i] {
                st.possible_bonds[i] = true;
                st.known_bonds[i] = false;
            }
            if st.possible_bonds[i] {
                st.bond_symbols[i] += &format!("_{i}");
            }
        }
        flag_ring_stereo(mol, &mut st);
        res = refine(mol, &mut st).0;
    }
    res
}
