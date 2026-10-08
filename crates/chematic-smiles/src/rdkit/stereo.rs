//! RDKit's legacy stereo perception, `MolOps::assignStereochemistry` with
//! `Chirality::getUseLegacyStereoPerception() == true` (RDKit 2026.03.1
//! `Chirality.cpp`: `legacyStereoPerception`, `assignAtomCIPRanks`,
//! `assignAtomChiralCodes`, `assignBondStereoCodes`, `rerankAtoms`,
//! `findChiralAtomSpecialCases`; `FindStereo.cpp`:
//! `isAtomPotentialTetrahedralCenter`).

use std::collections::{BTreeSet, VecDeque};

use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Hybridization, Mol};
use super::periodic;
use super::sanitize::atom_has_conjugated_bond;

const MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO: usize = 8;

/// `shouldDetectDoubleBondStereo`.
pub(crate) fn should_detect_double_bond_stereo(mol: &Mol, b: usize) -> bool {
    let ri = mol.ring_info();
    ri.num_bond_rings(b) == 0 || ri.min_bond_ring_size(b) >= MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
}

/// `Chirality::detail::bondAffectsAtomChirality`.
fn bond_affects_atom_chirality(mol: &Mol, b: usize, a: usize) -> bool {
    let bond = &mol.bonds[b];
    !(bond.bt == BondType::Dative && bond.begin == a)
}

/// `Chirality::detail::getAtomNonzeroDegree`.
pub(crate) fn atom_nonzero_degree(mol: &Mol, a: usize) -> usize {
    mol.atom_bonds[a]
        .iter()
        .filter(|&&b| bond_affects_atom_chirality(mol, b, a))
        .count()
}

/// `Chirality::detail::has_protium_neighbor`.
fn has_protium_neighbor(mol: &Mol, a: usize) -> bool {
    mol.nbrs(a)
        .any(|nb| mol.atoms[nb].anum == 1 && mol.atoms[nb].isotope == 0)
}

/// `queryIsAtomBridgehead`.
pub(crate) fn is_atom_bridgehead(mol: &Mol, a: usize) -> bool {
    if mol.degree(a) < 3 {
        return false;
    }
    let Some(ri) = mol.rings.as_ref() else {
        return false;
    };
    let mut atom_ring_bonds = vec![false; mol.bonds.len()];
    let mut count = 0;
    for &b in &mol.atom_bonds[a] {
        if ri.num_bond_rings(b) != 0 {
            atom_ring_bonds[b] = true;
            count += 1;
        }
    }
    if count < 3 {
        return false;
    }
    let nr = ri.bond_rings.len();
    let mut rings_overlap = vec![false; nr];
    for i in 0..nr {
        let mut in_ring_i = vec![false; mol.bonds.len()];
        let mut atom_in_ring_i = false;
        for &b in &ri.bond_rings[i] {
            in_ring_i[b] = true;
            if atom_ring_bonds[b] {
                atom_in_ring_i = true;
            }
        }
        if !atom_in_ring_i {
            continue;
        }
        for j in i + 1..nr {
            let mut overlap = 0;
            let mut atom_in_ring_j = false;
            for &b in &ri.bond_rings[j] {
                if atom_ring_bonds[b] {
                    atom_in_ring_j = true;
                }
                if in_ring_i[b] {
                    overlap += 1;
                }
                if overlap >= 2 && atom_in_ring_j {
                    rings_overlap[i] = true;
                    rings_overlap[j] = true;
                    break;
                }
            }
        }
        if !rings_overlap[i] {
            return false;
        }
    }
    true
}

/// `Chirality::detail::isAtomPotentialTetrahedralCenter`.
pub(crate) fn is_atom_potential_tetrahedral_center(mol: &Mol, a: usize) -> bool {
    let nz = atom_nonzero_degree(mol, a);
    let atom = &mol.atoms[a];
    let tnz = nz + mol.total_num_hs(a) as usize;
    if tnz > 4 {
        return false;
    }
    if nz == 4 {
        true
    } else if nz <= 1 {
        false
    } else if nz < 3 && atom.anum != 15 && atom.anum != 33 {
        false
    } else if atom.anum == 15 || atom.anum == 33 {
        true
    } else if nz == 3 {
        if mol.total_num_hs(a) == 1 {
            !has_protium_neighbor(mol, a)
        } else if (atom.anum == 16 || atom.anum == 34)
            && (atom.explicit_valence == 4 || (atom.explicit_valence == 3 && atom.charge == 1))
        {
            true
        } else if atom.anum == 7 {
            atom.hybrid == Hybridization::Sp3
                && !atom_has_conjugated_bond(mol, a)
                && (mol.ring_info().is_atom_in_ring_of_size(a, 3) || is_atom_bridgehead(mol, a))
        } else {
            false
        }
    } else {
        false
    }
}

/// `isAtomPotentialChiralCenter`: (legal center, has duplicate neighbour
/// ranks); `nbrs` gets (rank, bond index) for every bond when ranks are
/// given.
fn is_atom_potential_chiral_center(
    mol: &Mol,
    a: usize,
    ranks: &[u32],
    nbrs: &mut Vec<(u32, usize)>,
) -> (bool, bool) {
    let atom = &mol.atoms[a];
    let mut legal = true;
    let mut dupes = false;
    let nz = atom_nonzero_degree(mol, a);
    let tnz = nz + mol.total_num_hs(a) as usize;
    if tnz > 4 {
        legal = false;
    } else {
        if tnz < 3 {
            legal = false;
        } else if nz < 3 && atom.anum != 15 && atom.anum != 33 {
            legal = false;
        } else if nz == 3 {
            if mol.total_num_hs(a) == 1 {
                if has_protium_neighbor(mol, a) {
                    legal = false;
                }
            } else {
                legal = false;
                if atom.anum == 7 {
                    if atom.hybrid == Hybridization::Sp3
                        && !atom_has_conjugated_bond(mol, a)
                        && (mol.ring_info().is_atom_in_ring_of_size(a, 3)
                            || is_atom_bridgehead(mol, a))
                    {
                        legal = true;
                    }
                } else if atom.anum == 15 || atom.anum == 33 {
                    legal = true;
                } else if (atom.anum == 16 || atom.anum == 34)
                    && (atom.explicit_valence == 4
                        || (atom.explicit_valence == 3 && atom.charge == 1))
                {
                    legal = true;
                }
            }
        }
        if legal && !ranks.is_empty() {
            let mut seen = vec![false; mol.atoms.len()];
            for &b in &mol.atom_bonds[a] {
                let other = mol.bonds[b].other(a);
                nbrs.push((ranks[other], b));
                if !bond_affects_atom_chirality(mol, b, a) {
                    continue;
                }
                let r = ranks[other] as usize;
                if seen[r] {
                    dupes = true;
                    break;
                }
                seen[r] = true;
            }
        }
    }
    (legal, dupes)
}

/// `buildCIPInvariants`.
fn build_cip_invariants(mol: &Mol) -> Vec<i64> {
    const N_MASS_BITS: i64 = 10;
    const MAX_MASS: i64 = 1 << N_MASS_BITS;
    mol.atoms
        .iter()
        .map(|atom| {
            let num = i64::from(atom.anum % 128);
            let mut mass: i64 = 0;
            if atom.isotope != 0 {
                mass =
                    i64::from(atom.isotope) - i64::from(periodic::most_common_isotope(atom.anum));
                if mass >= 0 {
                    mass += 1;
                }
            }
            mass += MAX_MASS / 2;
            if mass < 0 {
                mass = 0;
            } else {
                mass %= MAX_MASS;
            }
            let mut invariant = num;
            invariant = (invariant << N_MASS_BITS) | mass;
            let mapnum = atom.map.map_or(-1i64, i64::from);
            let mapnum = (mapnum + 1) % 1024;
            invariant = (invariant << 10) | mapnum;
            invariant
        })
        .collect()
}

/// `findSegmentsToResort`: assigns ranks along `sorted`, returns tied
/// sections and the number of distinct ranks.
fn find_segments_to_resort(
    sorted: &[usize],
    cip: &[Vec<i32>],
    curr_rank: &mut [u32],
) -> (Vec<(usize, usize)>, usize) {
    let mut res: Vec<(usize, usize)> = Vec::new();
    let mut num_independent = sorted.len();
    let mut current = 0usize;
    let mut running = 0u32;
    curr_rank[sorted[0]] = 0;
    let mut in_equal = false;
    for i in 1..sorted.len() {
        if cip[sorted[current]] == cip[sorted[i]] {
            curr_rank[sorted[i]] = running;
            num_independent -= 1;
            if !in_equal {
                in_equal = true;
                res.push((i - 1, 0));
            }
        } else {
            running += 1;
            curr_rank[sorted[i]] = running;
            current = i;
            if in_equal {
                res.last_mut().expect("open section").1 = i;
                in_equal = false;
            }
        }
    }
    if in_equal {
        res.last_mut().expect("open section").1 = sorted.len() - 1;
    }
    (res, num_independent)
}

/// `iterateCIPRanks`.
fn iterate_cip_ranks(mol: &Mol, invars: &[i64], ranks: &mut Vec<u32>, seed_with_invars: bool) {
    let n = mol.atoms.len();
    ranks.resize(n, 0);
    if n == 0 {
        return;
    }
    // Entry buffers are reused across calls (RDKit reserves 16 per entry).
    let mut pool: Vec<Vec<i32>> =
        CIP_ENTRY_POOL.with(|pool| std::mem::take(&mut *pool.borrow_mut()));
    while pool.len() < n {
        pool.push(Vec::with_capacity(16));
    }
    let cip = &mut pool[..n];
    for (e, &v) in cip.iter_mut().zip(invars) {
        e.clear();
        e.push(v as i32);
    }
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by(|&x, &y| cip[x].cmp(&cip[y]));
    let mut curr_rank = vec![0u32; n];
    let (mut needs_sorting, mut num_ranks) = find_segments_to_resort(&sorted, cip, &mut curr_rank);
    ranks.copy_from_slice(&curr_rank);
    for i in 0..n {
        if seed_with_invars {
            cip[i][0] = invars[i] as i32;
        } else {
            cip[i][0] = mol.atoms[i].anum as i32;
            cip[i].push(ranks[i] as i32);
        }
    }
    let cip_rank_index = if seed_with_invars { 1 } else { 2 };
    let max_its = n / 2 + 1;
    let mut num_its = 0;
    let mut last_num_ranks: i64 = -1;
    // Bond features: per atom, (count, neighbour) in bond order, stored
    // back to back (atom `a`'s are `features[feature_start[a]..feature_start[a + 1]]`).
    let mut feature_start: Vec<usize> = Vec::with_capacity(n + 1);
    let mut features: Vec<(u32, usize)> = Vec::with_capacity(2 * mol.bonds.len());
    feature_start.push(0);
    for a in 0..n {
        for &b in &mol.atom_bonds[a] {
            let bond = &mol.bonds[b];
            let nbr = bond.other(a);
            let special = bond.bt == BondType::Double
                && mol.atoms[nbr].anum == 15
                && matches!(mol.degree(nbr), 3 | 4);
            let count = if special { 1 } else { bond.bt.twice() };
            features.push((count, nbr));
        }
        feature_start.push(features.len());
    }
    let total_hs: Vec<usize> = (0..n).map(|a| mol.total_num_hs(a) as usize).collect();
    let mut vals: Vec<i32> = Vec::with_capacity(16);
    while !needs_sorting.is_empty()
        && num_its < max_its
        && (last_num_ranks < 0 || (last_num_ranks as usize) < num_ranks)
    {
        for index in 0..n {
            vals.clear();
            for &(count, nbr) in &features[feature_start[index]..feature_start[index + 1]] {
                for _ in 0..count {
                    vals.push(ranks[nbr] as i32 + 1);
                }
            }
            vals.sort_unstable_by(|x, y| y.cmp(x));
            let entry = &mut cip[index];
            entry.reserve(vals.len() + total_hs[index]);
            entry.extend_from_slice(&vals);
            entry.resize(entry.len() + total_hs[index], 0);
        }
        last_num_ranks = num_ranks as i64;
        for &(first, last) in &needs_sorting {
            sorted[first..=last].sort_by(|&x, &y| cip[x].cmp(&cip[y]));
        }
        let (ns, nr) = find_segments_to_resort(&sorted, cip, &mut curr_rank);
        needs_sorting = ns;
        num_ranks = nr;
        ranks.copy_from_slice(&curr_rank);
        if last_num_ranks as usize != num_ranks {
            for i in 0..n {
                cip[i].resize(cip_rank_index + 1, 0);
                cip[i][cip_rank_index] = ranks[i] as i32;
            }
        }
        num_its += 1;
    }
    CIP_ENTRY_POOL.with(|p| *p.borrow_mut() = pool);
}

thread_local! {
    static CIP_ENTRY_POOL: std::cell::RefCell<Vec<Vec<i32>>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// `assignAtomCIPRanks`.
pub(crate) fn assign_atom_cip_ranks(mol: &Mol, ranks: &mut Vec<u32>) {
    let invars = build_cip_invariants(mol);
    iterate_cip_ranks(mol, &invars, ranks, false);
}

/// `findAtomNeighborDirHelper` (pairs of neighbour, direction relative to
/// the double bond atom; empty when no direction is set).
fn find_atom_neighbor_dir_helper(
    mol: &Mol,
    a: usize,
    ref_bond: usize,
    ranks: &[u32],
) -> Vec<(usize, BondDir)> {
    let mut seen_dir = false;
    let mut neighbors = Vec::new();
    for &b in &mol.atom_bonds[a] {
        let bond = &mol.bonds[b];
        let mut dir = bond.dir;
        if b != ref_bond {
            if dir.is_set() {
                seen_dir = true;
                if a != bond.begin {
                    dir = dir.flipped();
                }
            }
            neighbors.push((bond.other(a), dir));
        }
    }
    if !seen_dir {
        neighbors.clear();
    } else if neighbors.len() == 2 && ranks[neighbors[0].0] == ranks[neighbors[1].0] {
        neighbors.clear();
    } else if !neighbors[0].1.is_set() {
        neighbors[0].1 = neighbors[1].1.flipped();
    } else if neighbors.len() > 1 && !neighbors[1].1.is_set() {
        neighbors[1].1 = neighbors[0].1.flipped();
    }
    neighbors
}

/// Neighbour directions for a double bond with a requested cis/trans
/// configuration, as `SetDoubleBondNeighborDirections` would set them and
/// [`find_atom_neighbor_dir_helper`] would then read them: the requested
/// stereo atoms get one direction (the same one for cis), the other
/// neighbour on each end the opposite.
type NeighborDirections = Vec<(usize, BondDir)>;
type DoubleBondNeighborDirections = (NeighborDirections, NeighborDirections);

fn requested_neighbor_dirs(
    mol: &Mol,
    b: usize,
    (stereo_begin, stereo_end, trans): (usize, usize, bool),
    ranks: &[u32],
) -> DoubleBondNeighborDirections {
    let side = |a: usize, stereo_atom: usize, dir: BondDir| {
        // Only single and aromatic bonds take a direction.
        let mut out: Vec<(usize, BondDir)> = mol.atom_bonds[a]
            .iter()
            .filter(|&&nb| {
                nb != b && matches!(mol.bonds[nb].bt, BondType::Single | BondType::Aromatic)
            })
            .map(|&nb| {
                let o = mol.bonds[nb].other(a);
                (o, if o == stereo_atom { dir } else { dir.flipped() })
            })
            .collect();
        if out.len() == 2 && ranks[out[0].0] == ranks[out[1].0] {
            out.clear();
        }
        out
    };
    let (beg, end) = (mol.bonds[b].begin, mol.bonds[b].end);
    let up = BondDir::EndUpRight;
    (
        side(beg, stereo_begin, up),
        side(end, stereo_end, if trans { up.flipped() } else { up }),
    )
}

/// `assignAtomChiralCodes`: (unassigned atoms remain, any assigned).
fn assign_atom_chiral_codes(
    mol: &mut Mol,
    ranks: &mut Vec<u32>,
    flag_possible: bool,
) -> (bool, bool) {
    let mut atom_changed = false;
    let mut unassigned = 0;
    for a in 0..mol.atoms.len() {
        let tag = mol.atoms[a].chiral;
        if !(flag_possible || tag != ChiralTag::Unspecified) {
            continue;
        }
        if mol.atoms[a].cip_code.is_some() {
            continue;
        }
        if ranks.is_empty() {
            assign_atom_cip_ranks(mol, ranks);
        }
        let mut nbrs = Vec::new();
        let (legal, dupes) = is_atom_potential_chiral_center(mol, a, ranks, &mut nbrs);
        if legal {
            unassigned += 1;
        }
        if legal && !dupes && flag_possible {
            mol.atoms[a].chirality_possible = true;
        }
        if legal && !dupes && tag != ChiralTag::Unspecified {
            atom_changed = true;
            unassigned -= 1;
            // std::sort with Rankers::pairLess (ranks are distinct here,
            // except for bonds that do not affect chirality).
            nbrs.sort_by_key(|p| p.0);
            let indices: Vec<usize> = nbrs.iter().map(|p| p.1).collect();
            let mut odd = mol.perturbation_is_odd(a, &indices);
            if indices.len() == 3 && mol.total_num_hs(a) == 1 {
                odd = !odd;
            }
            let mut t = tag;
            if odd {
                t = if t == ChiralTag::Ccw {
                    ChiralTag::Cw
                } else {
                    ChiralTag::Ccw
                };
            }
            mol.atoms[a].cip_code = Some(if t == ChiralTag::Ccw { b'S' } else { b'R' });
        }
    }
    (unassigned > 0, atom_changed)
}

/// `assignBondStereoCodes`: (unassigned bonds remain, any assigned).
fn assign_bond_stereo_codes(mol: &mut Mol, ranks: &mut Vec<u32>) -> (bool, bool) {
    let mut assigned = false;
    let mut unassigned = 0;
    let mut to_clear = vec![false; mol.bonds.len()];
    for b in 0..mol.bonds.len() {
        if mol.bonds[b].bt != BondType::Double || mol.bonds[b].stereo != BondStereo::None {
            continue;
        }
        if ranks.is_empty() {
            assign_atom_cip_ranks(mol, ranks);
        }
        mol.bonds[b].stereo_atoms.clear();
        if !should_detect_double_bond_stereo(mol, b) {
            continue;
        }
        let (beg, end) = (mol.bonds[b].begin, mol.bonds[b].end);
        if !matches!(mol.degree(beg), 2 | 3) || !matches!(mol.degree(end), 2 | 3) {
            continue;
        }
        unassigned += 1;
        let (beg_n, end_n) = match mol.bonds[b].requested {
            Some(req) => requested_neighbor_dirs(mol, b, req, ranks),
            None => (
                find_atom_neighbor_dir_helper(mol, beg, b, ranks),
                find_atom_neighbor_dir_helper(mol, end, b, ranks),
            ),
        };
        if beg_n.is_empty() || end_n.is_empty() {
            continue;
        }
        let pick = |n: &Vec<(usize, BondDir)>| {
            if n.len() == 1 || ranks[n[0].0] > ranks[n[1].0] {
                n[0]
            } else {
                n[1]
            }
        };
        let (beg_nbr, beg_dir) = pick(&beg_n);
        let (end_nbr, end_dir) = pick(&end_n);
        let conflicting_begin = beg_n.len() == 2 && beg_n[0].1 == beg_n[1].1;
        let conflicting_end = end_n.len() == 2 && end_n[0].1 == end_n[1].1;
        if conflicting_begin || conflicting_end {
            mol.bonds[b].stereo = BondStereo::None;
            assigned = true;
            if conflicting_begin {
                for k in 0..2 {
                    if let Some(x) = mol.bond_between(beg_n[k].0, beg) {
                        to_clear[x] = true;
                    }
                }
            }
            if conflicting_end {
                for k in 0..2 {
                    if let Some(x) = mol.bond_between(end_n[k].0, end) {
                        to_clear[x] = true;
                    }
                }
            }
        } else {
            mol.bonds[b].stereo_atoms = vec![beg_nbr, end_nbr];
            mol.bonds[b].stereo = if beg_dir == end_dir {
                BondStereo::Z
            } else {
                BondStereo::E
            };
            assigned = true;
        }
        unassigned -= 1;
    }
    for (b, clear) in to_clear.into_iter().enumerate() {
        if clear {
            mol.bonds[b].dir = BondDir::None;
        }
    }
    (unassigned > 0, assigned)
}

/// `rerankAtoms`.
fn rerank_atoms(mol: &Mol, ranks: &mut Vec<u32>) {
    let n = mol.atoms.len();
    let mut factor: i64 = 100;
    while factor < n as i64 {
        factor *= 10;
    }
    let mut invars = vec![0i64; n];
    for i in 0..n {
        invars[i] = i64::from(ranks[i]) * factor;
        match mol.atoms[i].cip_code {
            Some(b'S') => invars[i] += 10,
            Some(b'R') => invars[i] += 20,
            _ => {}
        }
        for &b in &mol.atom_bonds[i] {
            let bond = &mol.bonds[b];
            if bond.bt == BondType::Double {
                match bond.stereo {
                    BondStereo::E => invars[i] += 1,
                    BondStereo::Z => invars[i] += 2,
                    _ => {}
                }
            }
        }
    }
    iterate_cip_ranks(mol, &invars, ranks, true);
}

/// `atomIsCandidateForRingStereochem`.
fn atom_is_candidate_for_ring_stereochem(mol: &mut Mol, a: usize, ranks: &[u32]) -> bool {
    if let Some(v) = mol.atoms[a].ring_stereochem_cand {
        return v;
    }
    let mut res = false;
    let ri = mol.ring_info();
    if ri.num_atom_rings(a) != 0 {
        if mol.atoms[a].anum == 7
            && mol.total_degree(a) == 3
            && !ri.is_atom_in_ring_of_size(a, 3)
            && !is_atom_bridgehead(mol, a)
        {
            // Not cached (RDKit returns before setting the property).
            return false;
        }
        let mut non_ring: Vec<usize> = Vec::new();
        let mut ring_nbrs: Vec<usize> = Vec::new();
        let mut ring_nbr_ranks: BTreeSet<u32> = BTreeSet::new();
        for &b in &mol.atom_bonds[a] {
            let nbr = mol.bonds[b].other(a);
            if ri.num_bond_rings(b) == 0 {
                non_ring.push(nbr);
            } else {
                ring_nbrs.push(nbr);
                ring_nbr_ranks.insert(ranks[nbr]);
            }
        }
        res = match non_ring.len() {
            2 => {
                ranks[non_ring[0]] != ranks[non_ring[1]] && ring_nbrs.len() != ring_nbr_ranks.len()
            }
            1 => ring_nbrs.len() > ring_nbr_ranks.len(),
            0 => {
                (ring_nbrs.len() == 4 && ring_nbr_ranks.len() == 3)
                    || (ring_nbrs.len() == 3 && ring_nbr_ranks.len() == 2)
            }
            _ => false,
        };
    }
    mol.atoms[a].ring_stereochem_cand = Some(res);
    res
}

/// `findChiralAtomSpecialCases`.
fn find_chiral_atom_special_cases(mol: &mut Mol, ranks: &[u32]) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut possible = vec![false; n];
    let mut atoms_seen = vec![false; n];
    let mut atoms_used = vec![false; n];
    let mut bonds_seen = vec![false; mol.bonds.len()];
    for a in 0..n {
        if atoms_seen[a] {
            continue;
        }
        if mol.atoms[a].chiral == ChiralTag::Unspecified
            || mol.atoms[a].cip_code.is_some()
            || mol.num_atom_rings(a) == 0
            || !atom_is_candidate_for_ring_stereochem(mol, a, ranks)
        {
            continue;
        }
        let mut next_atoms: VecDeque<usize> = VecDeque::new();
        for &b in &mol.atom_bonds[a] {
            if !bonds_seen[b] {
                bonds_seen[b] = true;
                if mol.num_bond_rings(b) != 0 {
                    let o = mol.bonds[b].other(a);
                    if !atoms_seen[o] {
                        next_atoms.push_back(o);
                        atoms_used[o] = true;
                    }
                }
            }
        }
        let mut ring_stereo_atoms: Vec<i32> = Vec::new();
        if !next_atoms.is_empty()
            && let Some(v) = &mol.atoms[a].ring_stereo_atoms
        {
            ring_stereo_atoms = v.clone();
        }
        while let Some(r) = next_atoms.pop_front() {
            atoms_seen[r] = true;
            if mol.atoms[r].chiral != ChiralTag::Unspecified
                && mol.atoms[r].cip_code.is_none()
                && atom_is_candidate_for_ring_stereochem(mol, r, ranks)
            {
                let same = if mol.atoms[r].chiral == mol.atoms[a].chiral {
                    1
                } else {
                    -1
                };
                ring_stereo_atoms.push(same * (r as i32 + 1));
                let mut o = mol.atoms[r].ring_stereo_atoms.clone().unwrap_or_default();
                o.push(same * (a as i32 + 1));
                mol.atoms[r].ring_stereo_atoms = Some(o);
                possible[r] = true;
                possible[a] = true;
            }
            for &b in &mol.atom_bonds[r] {
                if !bonds_seen[b] {
                    bonds_seen[b] = true;
                    if mol.num_bond_rings(b) != 0 {
                        let o = mol.bonds[b].other(r);
                        if !atoms_seen[o] && !atoms_used[o] {
                            next_atoms.push_back(o);
                            atoms_used[o] = true;
                        }
                    }
                }
            }
        }
        if !ring_stereo_atoms.is_empty() {
            mol.atoms[a].ring_stereo_atoms = Some(ring_stereo_atoms.clone());
            let idx_of = |e: i32| -> usize { (e.unsigned_abs() - 1) as usize };
            for (k, &entry) in ring_stereo_atoms.iter().enumerate() {
                let ridx = idx_of(entry);
                let mut lring = mol.atoms[ridx]
                    .ring_stereo_atoms
                    .clone()
                    .unwrap_or_default();
                for &oentry in &ring_stereo_atoms[k + 1..] {
                    let oidx = idx_of(oentry);
                    let different = (entry < 0) ^ (oentry < 0);
                    lring.push(if different {
                        -(oidx as i32 + 1)
                    } else {
                        oidx as i32 + 1
                    });
                    let mut olring = mol.atoms[oidx]
                        .ring_stereo_atoms
                        .clone()
                        .unwrap_or_default();
                    olring.push(if different {
                        -(ridx as i32 + 1)
                    } else {
                        ridx as i32 + 1
                    });
                    mol.atoms[oidx].ring_stereo_atoms = Some(olring);
                }
                mol.atoms[ridx].ring_stereo_atoms = Some(lring);
            }
        } else {
            possible[a] = false;
        }
        atoms_seen[a] = true;
    }
    possible
}

/// `legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters)`.
///
/// Returns the atoms' final `_CIPRank` values (empty where RDKit sets no
/// `_CIPRank`: no CIP ranking was needed).
pub(crate) fn legacy_stereo_perception(
    mol: &mut Mol,
    clean_it: bool,
    flag_possible: bool,
) -> Vec<u32> {
    legacy_stereo_perception_impl(mol, clean_it, flag_possible, false)
}

/// `legacy_stereo_perception(mol, true, true)` for callers that read neither
/// `chirality_possible` nor the returned ranks. On a molecule without
/// stereo atoms or bonds the flagging pass only sets `chirality_possible`
/// (no CIP codes, bond stereo or ring-stereo data survive it), so it is
/// skipped there, and with it the CIP ranking.
pub(crate) fn legacy_stereo_perception_unflagged(mol: &mut Mol) {
    legacy_stereo_perception_impl(mol, true, true, true);
}

fn legacy_stereo_perception_impl(
    mol: &mut Mol,
    clean_it: bool,
    flag_possible: bool,
    skip_flag_only_pass: bool,
) -> Vec<u32> {
    let mut has_stereo_atoms = false;
    let mut has_potential_stereo_atoms = false;
    for a in 0..mol.atoms.len() {
        if clean_it {
            mol.atoms[a].cip_code = None;
            mol.atoms[a].chirality_possible = false;
        }
        if !has_stereo_atoms && mol.atoms[a].chiral != ChiralTag::Unspecified {
            has_stereo_atoms = true;
        } else if !has_potential_stereo_atoms {
            let mut nbrs = Vec::new();
            has_potential_stereo_atoms = is_atom_potential_chiral_center(mol, a, &[], &mut nbrs).0;
        }
    }
    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    for b in 0..mol.bonds.len() {
        if clean_it {
            let bt = mol.bonds[b].bt;
            if (bt == BondType::Double || bt == BondType::Aromatic)
                && !should_detect_double_bond_stereo(mol, b)
            {
                if mol.bonds[b].stereo != BondStereo::None {
                    mol.bonds[b].stereo = BondStereo::None;
                    mol.bonds[b].stereo_atoms.clear();
                }
                continue;
            } else if bt == BondType::Double && mol.bonds[b].stereo != BondStereo::Any {
                mol.bonds[b].stereo = BondStereo::None;
                mol.bonds[b].stereo_atoms.clear();
            }
        }
        if !has_stereo_bonds && mol.bonds[b].bt == BondType::Double {
            let mut is_specified = false;
            if mol.bonds[b].requested.is_some() {
                has_stereo_bonds = true;
                is_specified = true;
            }
            for end in [mol.bonds[b].begin, mol.bonds[b].end] {
                if mol.atom_bonds[end]
                    .iter()
                    .any(|&nb| mol.bonds[nb].dir.is_set())
                {
                    has_stereo_bonds = true;
                    is_specified = true;
                    break;
                }
            }
            if !has_potential_stereo_bonds
                && !is_specified
                && should_detect_double_bond_stereo(mol, b)
            {
                has_potential_stereo_bonds = true;
            }
        }
        if !clean_it && has_stereo_bonds && has_potential_stereo_bonds {
            break;
        }
    }
    let mut ranks: Vec<u32> = Vec::new();
    let mut keep_going = has_stereo_atoms | has_stereo_bonds;
    if !keep_going && !skip_flag_only_pass {
        keep_going = flag_possible && (has_potential_stereo_atoms || has_potential_stereo_bonds);
    }
    while keep_going {
        let changed_atoms;
        let changed_bonds;
        if has_stereo_atoms || has_potential_stereo_atoms {
            let (h, c) = assign_atom_chiral_codes(mol, &mut ranks, flag_possible);
            has_stereo_atoms = h;
            changed_atoms = c;
        } else {
            changed_atoms = false;
        }
        if has_stereo_bonds || has_potential_stereo_bonds {
            let (h, c) = assign_bond_stereo_codes(mol, &mut ranks);
            has_stereo_bonds = h;
            changed_bonds = c;
        } else {
            changed_bonds = false;
        }
        keep_going = (has_stereo_atoms || has_stereo_bonds) && (changed_atoms || changed_bonds);
        if keep_going {
            rerank_atoms(mol, &mut ranks);
        }
    }
    if !clean_it {
        return ranks;
    }
    for atom in &mut mol.atoms {
        atom.ring_stereochem_cand = None;
        atom.ring_stereo_atoms = None;
    }
    let possible =
        if ranks.is_empty() && !mol.atoms.iter().any(|a| a.chiral != ChiralTag::Unspecified) {
            vec![false; mol.atoms.len()]
        } else {
            if ranks.is_empty() {
                assign_atom_cip_ranks(mol, &mut ranks);
            }
            find_chiral_atom_special_cases(mol, &ranks)
        };
    for a in 0..mol.atoms.len() {
        let atom = &mol.atoms[a];
        if atom.chiral != ChiralTag::Unspecified
            && atom.chiral.nontet().is_none()
            && atom.cip_code.is_none()
            && (!possible[a] || atom.ring_stereo_atoms.is_none())
        {
            mol.atoms[a].chiral = ChiralTag::Unspecified;
            let atom = &mol.atoms[a];
            if atom.num_explicit_hs == 1 && atom.charge == 0 && !atom.aromatic {
                mol.atoms[a].num_explicit_hs = 0;
                mol.atoms[a].no_implicit = false;
                // calcExplicitValence(false) / calcImplicitValence(false)
                if let Ok(v) = mol.calc_explicit_valence_value(a, false) {
                    mol.atoms[a].explicit_valence = v;
                }
                if let Ok(v) = mol.calc_implicit_valence_value(a, false) {
                    mol.atoms[a].implicit_valence = v;
                }
            }
        }
    }
    for b in 0..mol.bonds.len() {
        let (beg, end) = (mol.bonds[b].begin, mol.bonds[b].end);
        if mol.bonds[b].bt == BondType::Double
            && mol.bonds[b].stereo == BondStereo::Any
            && (mol.degree(beg) == 1 || mol.degree(end) == 1)
        {
            mol.bonds[b].stereo = BondStereo::None;
        }
        if mol.bonds[b].bt == BondType::Double
            && matches!(mol.bonds[b].stereo, BondStereo::Any | BondStereo::None)
        {
            for batom in [beg, end] {
                for nb_i in mol.atom_bonds[batom].clone() {
                    if nb_i == b {
                        continue;
                    }
                    let nbond = &mol.bonds[nb_i];
                    if nbond.dir.is_set()
                        && matches!(nbond.bt, BondType::Single | BondType::Aromatic)
                    {
                        let far = nbond.other(batom);
                        let ok_to_clear = !mol.atom_bonds[far].iter().any(|&nb_j| {
                            mol.bonds[nb_j].bt == BondType::Double
                                && !matches!(
                                    mol.bonds[nb_j].stereo,
                                    BondStereo::Any | BondStereo::None
                                )
                        });
                        if ok_to_clear {
                            mol.bonds[nb_i].dir = BondDir::None;
                        }
                    }
                }
            }
        }
    }
    ranks
}
