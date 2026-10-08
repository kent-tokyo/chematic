//! `MolOps::setAromaticity(mol, AROMATICITY_RDKIT)` (RDKit 2026.03.1
//! `Aromaticity.cpp`, `aromaticityHelper(mol, srings, 0, 0, true)`).

use std::collections::{BTreeMap, HashSet};

use super::kekulize::{make_ring_neighbor_map, pick_fused_rings};
use super::mol::{BondType, Mol};
use super::periodic;

/// Rings larger than this are not fused neighbours (`maxFusedAromaticRingSize`).
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Edon {
    Vacant,
    One,
    Two,
    Any,
    No,
}

/// `MolOps::countAtomElec`.
pub(crate) fn count_atom_elec(mol: &Mol, a: usize) -> i32 {
    let atom = &mol.atoms[a];
    let dv = periodic::default_valence(atom.anum);
    if dv <= 1 {
        return -1;
    }
    let mut degree = mol.degree(a) as i32 + mol.total_num_hs(a) as i32;
    for &b in &mol.atom_bonds[a] {
        if mol.bonds[b].valence_contrib(a) == 0.0 {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let nlp = (periodic::n_outer_elecs(atom.anum) - dv - atom.charge).max(0);
    let mut res = (dv - degree) + nlp - atom.radicals as i32;
    if res > 1 {
        let n_unsat = atom.explicit_valence - mol.degree(a) as i32;
        if n_unsat > 1 {
            res = 1;
        }
    }
    res
}

fn incident_non_cyclic_multiple_bond(mol: &Mol, a: usize) -> Option<usize> {
    for &b in &mol.atom_bonds[a] {
        if mol.num_bond_rings(b) == 0 && mol.bonds[b].valence_contrib(a) >= 2.0 {
            return Some(mol.bonds[b].other(a));
        }
    }
    None
}

fn incident_cyclic_multiple_bond(mol: &Mol, a: usize) -> bool {
    mol.atom_bonds[a]
        .iter()
        .any(|&b| mol.num_bond_rings(b) != 0 && mol.bonds[b].valence_contrib(a) >= 2.0)
}

fn incident_multiple_bond(mol: &Mol, a: usize) -> bool {
    let mut deg = mol.degree(a) as i32 + mol.atoms[a].num_explicit_hs as i32;
    for &b in &mol.atom_bonds[a] {
        if mol.bonds[b].valence_contrib(a).round() as i32 == 0 {
            deg -= 1;
        }
    }
    mol.atoms[a].explicit_valence != deg
}

fn get_atom_donor_type_arom(mol: &Mol, a: usize) -> Edon {
    let atom = &mol.atoms[a];
    if atom.anum == 0 {
        return if incident_cyclic_multiple_bond(mol, a) {
            Edon::One
        } else {
            Edon::Any
        };
    }
    let mut nelec = count_atom_elec(mol, a);
    if nelec < 0 {
        Edon::No
    } else if nelec == 0 {
        if incident_non_cyclic_multiple_bond(mol, a).is_some() {
            Edon::Vacant
        } else if incident_cyclic_multiple_bond(mol, a) {
            Edon::One
        } else {
            Edon::No
        }
    } else if nelec == 1 {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, a) {
            if periodic::more_electro_negative(mol.atoms[who].anum, atom.anum) {
                Edon::Vacant
            } else {
                Edon::One
            }
        } else if incident_multiple_bond(mol, a) {
            Edon::One
        } else if atom.charge == 1 {
            Edon::Vacant
        } else {
            Edon::No
        }
    } else {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, a)
            && periodic::more_electro_negative(mol.atoms[who].anum, atom.anum)
        {
            nelec -= 1;
        }
        if nelec % 2 == 1 { Edon::One } else { Edon::Two }
    }
}

fn is_atom_cand_for_arom(mol: &Mol, a: usize, edon: Edon) -> bool {
    let atom = &mol.atoms[a];
    if atom.anum > 18 && atom.anum != 34 && atom.anum != 52 {
        return false;
    }
    if edon == Edon::No {
        return false;
    }
    let def_val = periodic::default_valence(atom.anum);
    if def_val > 0
        && mol.total_valence(a)
            > periodic::default_valence((atom.anum as i32 - atom.charge).max(0) as u32)
    {
        return false;
    }
    if atom.radicals != 0 && (atom.anum != 6 || atom.charge != 0) {
        return false;
    }
    let n_unsat = atom.explicit_valence - mol.degree(a) as i32;
    if n_unsat > 1 {
        let mut n_mult = 0;
        for &b in &mol.atom_bonds[a] {
            match mol.bonds[b].bt {
                BondType::Double | BondType::Triple => n_mult += 1,
                _ => {}
            }
            if n_mult > 1 {
                break;
            }
        }
        if n_mult > 1 {
            return false;
        }
    }
    true
}

fn min_max_elecs(t: Edon) -> (i32, i32) {
    match t {
        Edon::Any => (1, 2),
        Edon::One => (1, 1),
        Edon::Two => (2, 2),
        Edon::No | Edon::Vacant => (0, 0),
    }
}

fn apply_huckel(ring: &[usize], edon: &[Edon]) -> bool {
    let (mut rlw, mut rup) = (0, 0);
    let mut n_any = 0;
    for &idx in ring {
        if edon[idx] == Edon::Any {
            n_any += 1;
            if n_any > 1 {
                return false;
            }
        }
        let (lw, up) = min_max_elecs(edon[idx]);
        rlw += lw;
        rup += up;
    }
    if rup >= 6 {
        (rlw..=rup).any(|rie| (rie - 2) % 4 == 0)
    } else {
        rup == 2
    }
}

/// `RDKit::nextCombination`.
fn next_combination(comb: &mut [usize], tot: usize) -> Option<usize> {
    let nelem = comb.len();
    let mut celem = nelem as isize - 1;
    while comb[celem as usize] == tot - nelem + celem as usize {
        celem -= 1;
        if celem < 0 {
            return None;
        }
    }
    let c = celem as usize;
    comb[c] += 1;
    for i in c + 1..nelem {
        comb[i] = comb[i - 1] + 1;
    }
    Some(c)
}

/// `RingUtils::checkFused`.
fn check_fused(rids: &[usize], neigh: &[Vec<usize>]) -> bool {
    let mut done = vec![true; neigh.len()];
    for &r in rids {
        done[r] = false;
    }
    let mut fused = Vec::new();
    pick_fused_rings(rids[0], neigh, &mut fused, &mut done);
    fused.len() == rids.len()
}

#[allow(clippy::too_many_arguments)]
fn apply_huckel_to_fused(
    mol: &mut Mol,
    srings: &[Vec<usize>],
    brings: &[Vec<usize>],
    fused: &[usize],
    edon: &[Edon],
    neigh: &[Vec<usize>],
    max_num_fused: usize,
) {
    let nrings = fused.len();
    let mut fused_bonds = vec![false; mol.bonds.len()];
    for &r in fused {
        for &b in &brings[r] {
            fused_bonds[b] = true;
        }
    }
    let n_ring_bonds = fused_bonds.iter().filter(|&&x| x).count();
    let mut done_bonds: HashSet<usize> = HashSet::new();
    let mut cur_size = 0usize;
    let mut comb: Vec<usize> = Vec::new();
    let mut pos: Option<usize> = None;
    loop {
        if pos.is_none() {
            if cur_size == 2 && nrings > 300 {
                break;
            }
            cur_size += 1;
            if cur_size > nrings.min(max_num_fused) || done_bonds.len() >= n_ring_bonds {
                break;
            }
            comb = (0..cur_size).collect();
            pos = Some(0);
        } else {
            pos = next_combination(&mut comb, nrings);
        }
        if pos.is_none() {
            continue;
        }
        let cur_rs: Vec<usize> = comb.iter().map(|&i| fused[i]).collect();
        if !neigh.is_empty() && !check_fused(&cur_rs, neigh) {
            continue;
        }
        let mut ats_in_system = vec![0u32; mol.atoms.len()];
        for &r in &cur_rs {
            for &a in &srings[r] {
                ats_in_system[a] += 1;
            }
        }
        let unon: Vec<usize> = (0..mol.atoms.len())
            .filter(|&i| ats_in_system[i] == 1 || ats_in_system[i] == 2)
            .collect();
        if apply_huckel(&unon, edon) {
            // markAtomsBondsArom
            let mut bnd_cntr: BTreeMap<usize, u32> = BTreeMap::new();
            for &r in &cur_rs {
                for &b in &brings[r] {
                    *bnd_cntr.entry(b).or_insert(0) += 1;
                }
            }
            for (&b, &cnt) in &bnd_cntr {
                if cnt == 1 {
                    let bond = &mut mol.bonds[b];
                    bond.aromatic = true;
                    if matches!(bond.bt, BondType::Single | BondType::Double) {
                        bond.bt = BondType::Aromatic;
                        let (x, y) = (bond.begin, bond.end);
                        mol.atoms[x].aromatic = true;
                        mol.atoms[y].aromatic = true;
                    }
                    done_bonds.insert(b);
                }
            }
        }
    }
}

pub(crate) fn set_aromaticity(mol: &mut Mol) {
    let srings = mol.ring_info().atom_rings.clone();
    let n = mol.atoms.len();
    let mut acands = vec![false; n];
    let mut aseen = vec![false; n];
    let mut edon = vec![Edon::No; n];
    let mut c_rings: Vec<Vec<usize>> = Vec::new();
    for sring in &srings {
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &a in sring {
            if all_dummy && mol.atoms[a].anum != 0 {
                all_dummy = false;
            }
            if aseen[a] {
                if !acands[a] {
                    all_aromatic = false;
                }
                continue;
            }
            aseen[a] = true;
            edon[a] = get_atom_donor_type_arom(mol, a);
            acands[a] = is_atom_cand_for_arom(mol, a, edon[a]);
            if !acands[a] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            c_rings.push(sring.clone());
        }
    }
    let brings: Vec<Vec<usize>> = c_rings
        .iter()
        .map(|ring| {
            (0..ring.len())
                .map(|k| {
                    mol.bond_between(ring[k], ring[(k + 1) % ring.len()])
                        .expect("ring atoms are bonded")
                })
                .collect()
        })
        .collect();
    let neigh = make_ring_neighbor_map(&brings, MAX_FUSED_AROMATIC_RING_SIZE, 1);
    let cnrs = c_rings.len();
    let mut fus_done = vec![false; cnrs];
    let mut curr = 0;
    while curr < cnrs {
        let mut fused = Vec::new();
        pick_fused_rings(curr, &neigh, &mut fused, &mut fus_done);
        apply_huckel_to_fused(mol, &c_rings, &brings, &fused, &edon, &neigh, 6);
        match (0..cnrs).find(|&r| !fus_done[r]) {
            Some(r) => curr = r,
            None => break,
        }
    }
}
