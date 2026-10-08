//! `MolOps::Kekulize(mol, markAtomsBonds=true, canonical)` (RDKit
//! 2026.03.1 `Kekulize.cpp`): [`kekulize`] is the `canonical=false` call
//! `sanitizeMol` makes (the atom ranks are the atom indices),
//! [`kekulize_ranked`] takes the atom ranks of a `canonical=true` call.

use std::collections::{HashMap, VecDeque};

use super::RdkitSmilesError;
use super::mol::{BondType, Mol};
use super::periodic;

const MAX_BACKTRACKS: usize = 100;

/// `RingUtils::pickFusedRings` (recursive, depth first).
pub(crate) fn pick_fused_rings(
    curr: usize,
    neigh: &[Vec<usize>],
    res: &mut Vec<usize>,
    done: &mut [bool],
) {
    done[curr] = true;
    res.push(curr);
    for &n in &neigh[curr] {
        if !done[n] {
            pick_fused_rings(n, neigh, res, done);
        }
    }
}

/// `RingUtils::makeRingNeighborMap(brings, neighMap, maxSize, maxOverlapSize)`.
pub(crate) fn make_ring_neighbor_map(
    brings: &[Vec<usize>],
    max_size: usize,
    max_overlap: usize,
) -> Vec<Vec<usize>> {
    let n = brings.len();
    let mut neigh = vec![Vec::new(); n];
    for i in 0..n {
        if max_size != 0 && brings[i].len() > max_size {
            continue;
        }
        for j in i + 1..n {
            if max_size != 0 && brings[j].len() > max_size {
                continue;
            }
            let inter = brings[i].iter().filter(|b| brings[j].contains(b)).count();
            if inter > 0 && (max_overlap == 0 || inter <= max_overlap) {
                neigh[i].push(j);
                neigh[j].push(i);
            }
        }
    }
    neigh
}

pub(crate) fn kekulize(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    kekulize_ranked(mol, None)
}

/// `Kekulize` with the traversal ordered by `ranks` (`None`: atom indices).
/// No bond carries a wedge here, so the wedge preferences never apply.
pub(crate) fn kekulize_ranked(
    mol: &mut Mol,
    ranks: Option<&[u32]>,
) -> Result<(), RdkitSmilesError> {
    let n_atoms = mol.atoms.len();
    let mut found_aromatic = mol.bonds.iter().any(|b| b.aromatic);
    let mut valences = vec![0i32; n_atoms];
    let mut dummy = vec![false; n_atoms];
    for a in 0..n_atoms {
        mol.calc_implicit_valence(a, false)?;
        valences[a] = mol.total_valence(a);
        if mol.is_aromatic_atom(a) {
            found_aromatic = true;
        }
        if mol.atoms[a].anum == 0 {
            dummy[a] = true;
        }
    }
    if !found_aromatic {
        return Ok(());
    }
    if !mol.bonds.is_empty() {
        let allrings: Vec<Vec<usize>> = mol.ring_info().atom_rings.clone();
        let arings: Vec<Vec<usize>> = allrings
            .into_iter()
            .filter(|ring| ring.iter().any(|&a| !dummy[a]))
            .collect();
        let brings: Vec<Vec<usize>> = arings
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
        let neigh = make_ring_neighbor_map(&brings, 0, 0);
        let cnrs = arings.len();
        let mut fus_done = vec![false; cnrs];
        let mut curr = 0;
        while curr < cnrs {
            let mut fused = Vec::new();
            pick_fused_rings(curr, &neigh, &mut fused, &mut fus_done);
            let frings: Vec<&Vec<usize>> = fused.iter().map(|&r| &arings[r]).collect();
            kekulize_fused(mol, &frings, ranks)?;
            match (0..cnrs).find(|&r| !fus_done[r]) {
                Some(r) => curr = r,
                None => break,
            }
        }
    }
    // markAtomsBonds
    for b in &mut mol.bonds {
        b.aromatic = false;
    }
    for a in 0..n_atoms {
        if mol.atoms[a].aromatic {
            if mol.num_atom_rings(a) == 0 {
                return Err(RdkitSmilesError::Sanitization(format!(
                    "non-ring atom {a} marked aromatic"
                )));
            }
            mol.atoms[a].aromatic = false;
            let atom = &mol.atoms[a];
            if (atom.anum == 7 || atom.anum == 15) && atom.charge == 0 && atom.num_explicit_hs == 1
            {
                mol.atoms[a].no_implicit = false;
                mol.atoms[a].num_explicit_hs = 0;
                mol.update_atom_property_cache(a, false)?;
            }
        }
    }
    for a in 0..n_atoms {
        if mol.total_valence(a) != valences[a] {
            return Err(RdkitSmilesError::Sanitization(format!(
                "Kekulization somehow screwed up valence on {a}"
            )));
        }
    }
    Ok(())
}

fn kekulize_fused(
    mol: &mut Mol,
    arings: &[&Vec<usize>],
    ranks: Option<&[u32]>,
) -> Result<(), RdkitSmilesError> {
    let mut all_atms: Vec<usize> = Vec::new();
    for ring in arings {
        for &a in ring.iter() {
            if !all_atms.contains(&a) {
                all_atms.push(a);
            }
        }
    }
    let n = mol.atoms.len();
    let mut d_bnd_cands = vec![false; n];
    let mut questions = Vec::new();
    let mut done = Vec::new();
    mark_dbond_cands(mol, &all_atms, &mut d_bnd_cands, &mut questions, &mut done);
    let d_bnd_adds = vec![false; mol.bonds.len()];
    let mut kekulized =
        kekulize_worker(mol, &all_atms, d_bnd_cands.clone(), d_bnd_adds, done, ranks);
    if !kekulized && !questions.is_empty() {
        kekulized = permute_dummies_and_kekulize(mol, &all_atms, &d_bnd_cands, &questions, ranks);
    }
    if !kekulized {
        return Err(RdkitSmilesError::Sanitization("Can't kekulize mol.".into()));
    }
    Ok(())
}

fn mark_dbond_cands(
    mol: &mut Mol,
    all_atms: &[usize],
    d_bnd_cands: &mut [bool],
    questions: &mut Vec<usize>,
    done: &mut Vec<usize>,
) {
    let has_arom_or_dummy = all_atms
        .iter()
        .any(|&a| mol.atoms[a].anum == 0 || mol.is_aromatic_atom(a));
    if !has_arom_or_dummy {
        return;
    }
    let rings = mol.ring_info();
    let mut is_ring_not_cand = vec![false; rings.atom_rings.len()];
    for (ri, aring) in rings.atom_rings.iter().enumerate() {
        is_ring_not_cand[ri] = true;
        for &ai in aring {
            if mol.is_aromatic_atom(ai) && rings.num_atom_rings(ai) == 1 {
                is_ring_not_cand[ri] = false;
                break;
            }
        }
    }
    let mut make_single: Vec<usize> = Vec::new();
    let mut in_all = vec![false; mol.atoms.len()];
    for &a in all_atms {
        in_all[a] = true;
        if mol.atoms[a].anum != 0 && !mol.is_aromatic_atom(a) {
            done.push(a);
            continue;
        }
        let mut sbo: i32 = 0;
        let mut n_to_ignore: i32 = 0;
        let mut non_ar_non_dummy_nbr = 0usize;
        for &b in &mol.atom_bonds[a] {
            let bond = &mol.bonds[b];
            let other = bond.other(a);
            if mol.atoms[other].anum != 0 && !mol.atoms[other].aromatic && in_all[other] {
                non_ar_non_dummy_nbr += 1;
            }
            if bond.aromatic
                && matches!(
                    bond.bt,
                    BondType::Single | BondType::Double | BondType::Aromatic
                )
            {
                sbo += 1;
                make_single.push(b);
            } else {
                let contrib = bond.valence_contrib(a).round() as i32;
                sbo += contrib;
                if contrib == 0 {
                    n_to_ignore += 1;
                }
            }
        }
        let rinfo = mol.ring_info();
        let num_atom_rings = rinfo.num_atom_rings(a);
        let num_non_cand = rinfo
            .atom_members(a)
            .iter()
            .filter(|&&ri| is_ring_not_cand[ri])
            .count();
        if mol.atoms[a].anum == 0
            && non_ar_non_dummy_nbr < num_atom_rings
            && num_non_cand < num_atom_rings
        {
            d_bnd_cands[a] = true;
            questions.push(a);
        } else {
            sbo += mol.total_num_hs(a) as i32;
            let atom = &mol.atoms[a];
            let mut dv = periodic::default_valence(atom.anum);
            let mut chrg = atom.charge;
            if periodic::is_early_atom(atom.anum) {
                chrg = -chrg;
            }
            if atom.anum == 6 && chrg > 0 {
                chrg = -chrg;
            }
            dv += chrg;
            let tbo = mol.total_valence(a);
            let n_radicals = atom.radicals as i32;
            let total_degree = mol.degree(a) as i32 + mol.implicit_valence(a) - n_to_ignore;
            let val_list = periodic::valence_list(atom.anum);
            let mut vi = 1;
            while tbo > dv && vi < val_list.len() && val_list[vi] > 0 {
                dv = i32::from(val_list[vi]) + chrg;
                vi += 1;
            }
            if tbo == 5
                && sbo == 4
                && dv == 3
                && total_degree == 3
                && n_radicals == 0
                && chrg == 0
                && mol.total_num_hs(a) == 0
                && matches!(atom.anum, 7 | 15 | 33)
            {
                dv = 5;
            }
            if total_degree + n_radicals >= dv {
                continue;
            }
            if dv == sbo + 1 + n_radicals || (n_radicals == 0 && atom.no_implicit && dv == sbo + 2)
            {
                d_bnd_cands[a] = true;
            }
        }
    }
    for b in make_single {
        mol.bonds[b].bt = BondType::Single;
    }
}

/// `done` with `done_count[a]`, the number of times `a` is in it.
struct Done {
    list: Vec<usize>,
    count: Vec<u32>,
}

impl Done {
    fn new(list: Vec<usize>, n: usize) -> Self {
        let mut count = vec![0u32; n];
        for &a in &list {
            count[a] += 1;
        }
        Self { list, count }
    }

    fn contains(&self, a: usize) -> bool {
        self.count[a] > 0
    }

    fn push(&mut self, a: usize) {
        self.list.push(a);
        self.count[a] += 1;
    }

    fn truncate(&mut self, len: usize) {
        for &a in &self.list[len..] {
            self.count[a] -= 1;
        }
        self.list.truncate(len);
    }
}

fn back_track(
    mol: &mut Mol,
    last_opt: usize,
    done: &mut Done,
    aqueue: &mut VecDeque<usize>,
    d_bnd_cands: &mut [bool],
    d_bnd_adds: &mut [bool],
) {
    let first = done
        .list
        .iter()
        .position(|&x| x == last_opt)
        .unwrap_or(done.list.len());
    let last = done
        .list
        .iter()
        .rposition(|&x| x == last_opt)
        .expect("last option was visited");
    for &x in done.list[last..].iter().rev() {
        aqueue.push_front(x);
    }
    // `done` becomes its first `first` atoms.
    done.truncate(first);
    for bi in 0..mol.bonds.len() {
        if d_bnd_adds[bi] {
            let (a1, a2) = (mol.bonds[bi].begin, mol.bonds[bi].end);
            if !done.contains(a1) && !done.contains(a2) {
                d_bnd_adds[bi] = false;
                mol.bonds[bi].bt = BondType::Single;
                d_bnd_cands[a1] = true;
                d_bnd_cands[a2] = true;
            }
        }
    }
}

fn kekulize_worker(
    mol: &mut Mol,
    all_atms: &[usize],
    mut d_bnd_cands: Vec<bool>,
    mut d_bnd_adds: Vec<bool>,
    done: Vec<usize>,
    ranks: Option<&[u32]>,
) -> bool {
    let n = mol.atoms.len();
    let mut done = Done::new(done, n);
    let mut astack: VecDeque<usize> = VecDeque::new();
    let mut options: HashMap<usize, VecDeque<usize>> = HashMap::new();
    let mut last_opt: Option<usize> = None;
    let mut local_added = vec![false; mol.bonds.len()];
    let mut in_all = vec![false; n];
    for &a in all_atms {
        in_all[a] = true;
    }
    // `lessByRank`: by rank, then index (canonical=false: ranks are the
    // indices); no wedged bonds.
    let key = |a: usize| (ranks.map_or(a as u32, |r| r[a]), a);
    let mut sorted = all_atms.to_vec();
    sorted.sort_unstable_by_key(|&a| key(a));
    let mut btmoves: VecDeque<usize> = VecDeque::new();
    let mut num_bt = 0usize;
    let mut nbrs: Vec<usize> = Vec::new();
    let mut lstack: Vec<usize> = Vec::new();
    while done.list.len() < sorted.len() || !astack.is_empty() {
        let curr = if let Some(c) = astack.pop_front() {
            c
        } else {
            match sorted.iter().copied().find(|&a| !done.contains(a)) {
                Some(c) => c,
                None => return false,
            }
        };
        done.push(curr);
        let c_cand = d_bnd_cands[curr];
        let mut opts: VecDeque<usize>;
        if let Some(o) = options.get(&curr) {
            opts = o.clone();
        } else {
            opts = VecDeque::new();
            nbrs.clear();
            nbrs.extend(
                mol.nbrs(curr)
                    .filter(|&nb| in_all[nb] && !done.contains(nb)),
            );
            nbrs.sort_unstable_by_key(|&a| key(a));
            lstack.clear();
            for &nb in &nbrs {
                let nb_bond = mol.bond_between(curr, nb).expect("bonded");
                if !astack.contains(&nb) {
                    lstack.push(nb);
                }
                if c_cand
                    && d_bnd_cands[nb]
                    && (mol.bonds[nb_bond].aromatic
                        || mol.atoms[curr].anum == 0
                        || mol.atoms[nb].anum == 0)
                {
                    opts.push_back(nb);
                }
            }
            astack.extend(lstack.iter().copied());
        }
        if c_cand {
            if let Some(ncnd) = opts.pop_front() {
                let bnd = mol.bond_between(curr, ncnd).expect("bonded");
                mol.bonds[bnd].bt = BondType::Double;
                mol.bonds[bnd].dir = super::mol::BondDir::None;
                d_bnd_cands[curr] = false;
                d_bnd_cands[ncnd] = false;
                d_bnd_adds[bnd] = true;
                local_added[bnd] = true;
                if let Some(stored) = options.get_mut(&curr) {
                    if opts.is_empty() {
                        options.remove(&curr);
                        btmoves.pop_back();
                        last_opt = btmoves.back().copied();
                    } else {
                        *stored = opts;
                    }
                } else if !opts.is_empty() {
                    last_opt = Some(curr);
                    btmoves.push_back(curr);
                    options.insert(curr, opts);
                }
            } else if mol.atoms[curr].anum != 0 {
                match last_opt {
                    Some(lo) if num_bt < MAX_BACKTRACKS => {
                        back_track(
                            mol,
                            lo,
                            &mut done,
                            &mut astack,
                            &mut d_bnd_cands,
                            &mut d_bnd_adds,
                        );
                        num_bt += 1;
                    }
                    _ => {
                        for bi in 0..mol.bonds.len() {
                            if local_added[bi] {
                                mol.bonds[bi].bt = BondType::Single;
                            }
                        }
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn permute_dummies_and_kekulize(
    mol: &mut Mol,
    all_atms: &[usize],
    d_bnd_cands: &[bool],
    questions: &[usize],
    ranks: Option<&[u32]>,
) -> bool {
    let mut in_play = vec![false; mol.atoms.len()];
    for &a in all_atms {
        in_play[a] = true;
    }
    let mut pos: u32 = 1;
    let mut kekulized = false;
    while !kekulized && !questions.is_empty() {
        for b in &mut mol.bonds {
            if b.aromatic && b.bt != BondType::Single && in_play[b.begin] && in_play[b.end] {
                b.bt = BondType::Single;
            }
        }
        if questions.len() >= 32 || pos >= (1u32 << questions.len()) {
            break;
        }
        let mut t_cands = d_bnd_cands.to_vec();
        for (i, &q) in questions.iter().enumerate() {
            if pos & (1u32 << i) != 0 {
                t_cands[q] = false;
            }
        }
        pos += 1;
        kekulized = kekulize_worker(
            mol,
            all_atms,
            t_cands,
            vec![false; mol.bonds.len()],
            Vec::new(),
            ranks,
        );
    }
    kekulized
}
