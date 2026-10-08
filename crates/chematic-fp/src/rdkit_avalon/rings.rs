// Ported from the Avalon Cheminformatics Toolkit, Copyright 2001-2011 Novartis
// Pharma AG, BSD-3-Clause license; see THIRD_PARTY_NOTICES.md.
//! Ring perception of the Avalon toolkit, ported operation for operation.
//!
//! Sources (Avalon Cheminformatics Toolkit, `AvalonToolkit_2.0.5-pre.3`,
//! BSD-3-Clause, see `THIRD_PARTY_NOTICES.md`): `set.c` (bit sets),
//! `graph.c` (`RingList`, `SortRings`, `CombineRings`, `ProperRingPairs`),
//! `utilities.c` (`RingState`) and `perceive.c` (`SetRingSizeFlags`,
//! `PerceiveAromaticBonds`, `PerceiveDYAromaticity`).
//!
//! The ring lists are singly linked, prepend-built lists in C; here they are
//! `Vec`s holding the list in head-first order. `CombineRings` breaks ties
//! with glibc's `rand()` after `srand(1)`, which [`GlibcRand`] reproduces.

use super::{AROMATIC, AvalonMolecule, DOUBLE, SINGLE, TRIPLE};

/// `bit_set_t`: a set of the integers `0..=max_member`.
#[derive(Clone, Debug)]
pub(crate) struct BitSet {
    max_member: usize,
    words: Vec<u64>,
}

impl BitSet {
    pub(crate) fn new(max_member: usize) -> Self {
        Self {
            max_member,
            words: vec![0; max_member / 64 + 1],
        }
    }

    pub(crate) fn is_member(&self, member: usize) -> bool {
        member <= self.max_member && self.words[member / 64] & (1u64 << (member % 64)) != 0
    }

    fn put(&mut self, member: usize) {
        if member > self.max_member {
            // `PutMember` grows the set; members never exceed the bond count
            // here, but mirror the C behaviour anyway.
            self.max_member = member;
            self.words.resize(member / 64 + 1, 0);
        }
        self.words[member / 64] |= 1u64 << (member % 64);
    }

    pub(crate) fn cardinality(&self) -> i32 {
        self.words.iter().map(|w| w.count_ones() as i32).sum()
    }

    fn xor_with(&mut self, other: &BitSet) {
        for (a, b) in self.words.iter_mut().zip(other.words.iter()) {
            *a ^= *b;
        }
    }

    fn and_with(&mut self, other: &BitSet) {
        for (i, a) in self.words.iter_mut().enumerate() {
            *a &= other.words.get(i).copied().unwrap_or(0);
        }
    }

    fn intersection_is_empty(&self, other: &BitSet) -> bool {
        self.words
            .iter()
            .zip(other.words.iter())
            .all(|(a, b)| a & b == 0)
    }

    /// Members in increasing order (`NextMember` iteration).
    pub(crate) fn members(&self) -> impl Iterator<Item = usize> + '_ {
        (0..=self.max_member).filter(move |&m| self.is_member(m))
    }
}

/// `bond_set_node`: one ring of a ring list.
#[derive(Clone, Debug)]
pub(crate) struct RingNode {
    pub(crate) bond_set: BitSet,
    pub(crate) cardinality: i32,
}

/// glibc's `rand()` (the `TYPE_3` additive feedback generator) after
/// `srand(1)`.
pub(crate) struct GlibcRand {
    r: Vec<u32>,
}

impl GlibcRand {
    pub(crate) fn new(seed: u32) -> Self {
        let mut r: Vec<u32> = Vec::with_capacity(400);
        let seed = if seed == 0 { 1 } else { seed };
        r.push(seed);
        let mut word = seed as i64;
        for _ in 1..31 {
            let hi = word / 127_773;
            let lo = word % 127_773;
            word = 16_807 * lo - 2_836 * hi;
            if word < 0 {
                word += 2_147_483_647;
            }
            r.push(word as u32);
        }
        for i in 31..34 {
            let v = r[i - 31];
            r.push(v);
        }
        let mut g = Self { r };
        // glibc discards the first 310 outputs.
        for _ in 0..310 {
            g.step();
        }
        g
    }

    fn step(&mut self) -> u32 {
        let n = self.r.len();
        let v = self.r[n - 31].wrapping_add(self.r[n - 3]);
        self.r.push(v);
        if self.r.len() > 64 {
            self.r.drain(0..self.r.len() - 34);
        }
        v
    }

    pub(crate) fn rand(&mut self) -> i32 {
        (self.step() >> 1) as i32
    }
}

/// `RingList`: basis rings of the graph `bonds` (pairs of 1-based atom
/// numbers; pairs containing atom 0 are ignored), head first.
pub(crate) fn ring_list(bonds: &[[u32; 2]]) -> Vec<RingNode> {
    const UNLINKED: i32 = -1;
    let nbonds = bonds.len();
    let mut natoms = 0u32;
    for b in bonds {
        natoms = natoms.max(b[0]).max(b[1]);
    }
    let natoms = natoms as usize;
    let mut tree_color = vec![0i32; natoms + 1];
    let mut tree_link = vec![UNLINKED; natoms + 1];
    let other =
        |b: usize, atom: usize| -> usize { (bonds[b][0] as usize + bonds[b][1] as usize) - atom };
    // Built in discovery order; the C list is prepend-built.
    let mut found: Vec<RingNode> = Vec::new();
    let mut color = 0i32;
    for b in 0..nbonds {
        let mut at1 = bonds[b][0] as usize;
        let mut at2 = bonds[b][1] as usize;
        if at1 == 0 || at2 == 0 {
            continue;
        }
        if tree_color[at1] == 0 && tree_color[at2] == 0 {
            color += 1;
            tree_link[at2] = b as i32;
            tree_color[at1] = color;
            tree_color[at2] = color;
        } else if tree_color[at1] == 0 {
            tree_color[at1] = tree_color[at2];
            tree_link[at1] = b as i32;
        } else if tree_color[at2] == 0 {
            tree_color[at2] = tree_color[at1];
            tree_link[at2] = b as i32;
        } else if tree_color[at1] != tree_color[at2] {
            let new_color = tree_color[at1];
            let old_color = tree_color[at2];
            for c in tree_color.iter_mut() {
                if *c == old_color {
                    *c = new_color;
                }
            }
            let mut tmp_link = tree_link[at2];
            tree_link[at2] = b as i32;
            while tmp_link != UNLINKED {
                at2 = other(tmp_link as usize, at2);
                std::mem::swap(&mut tmp_link, &mut tree_link[at2]);
            }
        } else {
            let mut level1 = 0;
            let mut trace = at1;
            while tree_link[trace] != UNLINKED {
                trace = other(tree_link[trace] as usize, trace);
                level1 += 1;
            }
            let mut level2 = 0;
            trace = at2;
            while tree_link[trace] != UNLINKED {
                trace = other(tree_link[trace] as usize, trace);
                level2 += 1;
            }
            if level1 > level2 {
                std::mem::swap(&mut level1, &mut level2);
                std::mem::swap(&mut at1, &mut at2);
            }
            let mut set = BitSet::new(nbonds);
            let mut cardinality = 1;
            set.put(b);
            for _ in 0..(level2 - level1) {
                let l = tree_link[at2] as usize;
                set.put(l);
                cardinality += 1;
                at2 = other(l, at2);
            }
            while at1 != at2 {
                let l1 = tree_link[at1] as usize;
                set.put(l1);
                cardinality += 1;
                at1 = other(l1, at1);
                let l2 = tree_link[at2] as usize;
                set.put(l2);
                cardinality += 1;
                at2 = other(l2, at2);
            }
            found.push(RingNode {
                bond_set: set,
                cardinality,
            });
        }
    }
    found.reverse();
    found
}

/// `SortRings`: in-place exchange sort by increasing cardinality.
fn sort_rings(list: &mut [RingNode]) {
    for i in 0..list.len() {
        for j in i + 1..list.len() {
            if list[i].cardinality > list[j].cardinality {
                list.swap(i, j);
            }
        }
    }
}

/// `CombineRings`: XOR pairs of rings until no ring gets smaller.
pub(crate) fn combine_rings(mut list: Vec<RingNode>) -> Vec<RingNode> {
    if list.is_empty() {
        return list;
    }
    let mut rng = GlibcRand::new(1);
    let mut ntoggle = 0;
    loop {
        let mut changed = false;
        sort_rings(&mut list);
        for i in 0..list.len() {
            for j in i + 1..list.len() {
                if list[i].bond_set.intersection_is_empty(&list[j].bond_set) {
                    continue;
                }
                let mut set = list[i].bond_set.clone();
                set.xor_with(&list[j].bond_set);
                let size = set.cardinality();
                if size > 0 && (size <= list[i].cardinality || size <= list[j].cardinality) {
                    let target = if list[i].cardinality > list[j].cardinality {
                        i
                    } else {
                        j
                    };
                    if list[target].cardinality > size || (rng.rand() / 10) % 2 != 0 {
                        if list[target].cardinality == size {
                            ntoggle += 1;
                        } else {
                            ntoggle = 0;
                        }
                        changed = true;
                        list[target].cardinality = size;
                        list[target].bond_set = set;
                    }
                }
            }
        }
        if ntoggle > 4 {
            changed = false;
        }
        if !changed {
            break;
        }
    }
    list
}

/// `ProperRingPairs`: XORs of base-ring pairs that share exactly one path,
/// head first (the C list is prepend-built).
pub(crate) fn proper_ring_pairs(
    base_rings: &[RingNode],
    maxnode: usize,
    bonds: &[[u32; 2]],
) -> Vec<RingNode> {
    let mut found = Vec::new();
    if base_rings.is_empty() {
        return found;
    }
    let mut atom_touched = vec![0i32; maxnode + 1];
    for i in 0..base_rings.len() {
        for j in i + 1..base_rings.len() {
            let (p1, p2) = (&base_rings[i], &base_rings[j]);
            if p1.bond_set.intersection_is_empty(&p2.bond_set) {
                continue;
            }
            let mut set = p1.bond_set.clone();
            set.and_with(&p2.bond_set);
            atom_touched.iter_mut().for_each(|t| *t = 0);
            for b in 0..set.max_member {
                if set.is_member(b) {
                    atom_touched[bonds[b][0] as usize] += 1;
                    atom_touched[bonds[b][1] as usize] += 1;
                }
            }
            let nterminal = atom_touched.iter().filter(|&&t| t == 1).count();
            if nterminal == 2 {
                let mut x = p1.bond_set.clone();
                x.xor_with(&p2.bond_set);
                let cardinality = x.cardinality();
                found.push(RingNode {
                    bond_set: x,
                    cardinality,
                });
            }
        }
    }
    found.reverse();
    found
}

fn bond_graph(mol: &AvalonMolecule) -> Vec<[u32; 2]> {
    mol.bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect()
}

/// `RingState`: per atom the number of attached ring bonds, per bond the
/// number of basis rings containing it.
pub(crate) fn ring_state(mol: &AvalonMolecule) -> (Vec<i32>, Vec<i32>) {
    let mut atom_status = vec![0i32; mol.atoms.len()];
    let mut bond_status = vec![0i32; mol.bonds.len()];
    if mol.bonds.is_empty() {
        return (atom_status, bond_status);
    }
    let rings = combine_rings(ring_list(&bond_graph(mol)));
    for r in &rings {
        for (i, s) in bond_status.iter_mut().enumerate() {
            if r.bond_set.is_member(i) {
                *s += 1;
            }
        }
    }
    for (i, b) in mol.bonds.iter().enumerate() {
        if bond_status[i] > 0 {
            atom_status[(b.atoms[0] - 1) as usize] += 1;
            atom_status[(b.atoms[1] - 1) as usize] += 1;
        }
    }
    (atom_status, bond_status)
}

/// Neighbourhood of one atom (`neighbourhood_t`): 0-based atom and bond
/// indices in bond order.
#[derive(Clone, Debug, Default)]
pub(crate) struct Neighbourhood {
    pub(crate) atoms: Vec<usize>,
    pub(crate) bonds: Vec<usize>,
}

/// `SetupNeighbourhood`.
pub(crate) fn setup_neighbourhood(mol: &AvalonMolecule) -> Vec<Neighbourhood> {
    let mut nbp = vec![Neighbourhood::default(); mol.atoms.len()];
    for (i, b) in mol.bonds.iter().enumerate() {
        let a0 = (b.atoms[0] - 1) as usize;
        let a1 = (b.atoms[1] - 1) as usize;
        nbp[a0].atoms.push(a1);
        nbp[a0].bonds.push(i);
        nbp[a1].atoms.push(a0);
        nbp[a1].bonds.push(i);
    }
    nbp
}

#[allow(clippy::too_many_arguments)]
fn mark_recursive(
    atom_rsize: &mut [u32],
    bond_rsize: &mut [u32],
    touched_atoms: &mut [bool],
    touched_bonds: &mut [bool],
    start_index: usize,
    path_length: u32,
    current_index: usize,
    max_size: u32,
    nbp: &[Neighbourhood],
) {
    for i in 0..nbp[current_index].atoms.len() {
        let ai = nbp[current_index].atoms[i];
        if ai < start_index {
            continue;
        }
        if ai == start_index {
            if path_length < 3 {
                continue;
            }
            for (j, &t) in touched_atoms.iter().enumerate() {
                if t {
                    atom_rsize[j] |= 1 << path_length;
                }
            }
            for (j, &t) in touched_bonds.iter().enumerate() {
                if t {
                    bond_rsize[j] |= 1 << path_length;
                }
            }
            continue;
        }
        if touched_atoms[ai] {
            continue;
        }
        if path_length + 1 > max_size {
            continue;
        }
        if atom_rsize[ai] == 0 {
            continue;
        }
        let bi = nbp[current_index].bonds[i];
        if bond_rsize[bi] == 0 {
            continue;
        }
        touched_atoms[ai] = true;
        touched_bonds[bi] = true;
        mark_recursive(
            atom_rsize,
            bond_rsize,
            touched_atoms,
            touched_bonds,
            start_index,
            path_length + 1,
            ai,
            max_size,
            nbp,
        );
        touched_atoms[ai] = false;
        touched_bonds[bi] = false;
    }
}

/// `SetRingSizeFlags`: returns `(atom rsize_flags, bond rsize_flags)`.
pub(crate) fn set_ring_size_flags(
    mol: &AvalonMolecule,
    max_size: u32,
    nbp: &[Neighbourhood],
) -> (Vec<u32>, Vec<u32>) {
    let (atom_status, bond_status) = ring_state(mol);
    let mut atom_rsize: Vec<u32> = atom_status.iter().map(|&s| u32::from(s > 0)).collect();
    let mut bond_rsize: Vec<u32> = bond_status.iter().map(|&s| u32::from(s > 0)).collect();
    let mut touched_atoms = vec![false; mol.atoms.len()];
    let mut touched_bonds = vec![false; mol.bonds.len()];
    for i in 0..mol.atoms.len() {
        if atom_rsize[i] == 0 {
            continue;
        }
        touched_atoms[i] = true;
        mark_recursive(
            &mut atom_rsize,
            &mut bond_rsize,
            &mut touched_atoms,
            &mut touched_bonds,
            i,
            1,
            i,
            max_size,
            nbp,
        );
        touched_atoms[i] = false;
    }
    (atom_rsize, bond_rsize)
}

/// `PerceiveAromaticBonds`: sets `bond_types` to `AROMATIC` for bonds in
/// Hückel-sized rings of sp2 atoms.
pub(crate) fn perceive_aromatic_bonds(mol: &AvalonMolecule, bond_types: &mut [i32]) {
    let graph = bond_graph(mol);
    let n_atoms = mol.atoms.len();
    let mut bond_is_in_ring = vec![false; mol.bonds.len()];
    let base = combine_rings(ring_list(&graph));
    for r in &base {
        for m in r.bond_set.members() {
            bond_is_in_ring[m] = true;
        }
    }
    let pairs = proper_ring_pairs(&base, n_atoms, &graph);
    // Prepending the pair list node by node reverses it.
    let mut rings: Vec<RingNode> = pairs.into_iter().rev().collect();
    rings.extend(base);

    let mut sp_count = vec![0i32; n_atoms + 1];
    loop {
        let mut changed = false;
        for r in &rings {
            sp_count.iter_mut().for_each(|c| *c = 0);
            let mut ndouble = 0;
            let mut nsingle = 0;
            for i in r.bond_set.members() {
                if bond_types[i] == SINGLE {
                    nsingle += 1;
                }
            }
            let mut is_cumulene = false;
            for i in r.bond_set.members() {
                let b = &mol.bonds[i];
                if bond_types[i] == DOUBLE {
                    ndouble += 1;
                    sp_count[b.atoms[0] as usize] += 1;
                    if sp_count[b.atoms[0] as usize] > 1 {
                        is_cumulene = true;
                    }
                    sp_count[b.atoms[1] as usize] += 1;
                    if sp_count[b.atoms[1] as usize] > 1 {
                        is_cumulene = true;
                    }
                } else if bond_types[i] == TRIPLE {
                    is_cumulene = true;
                }
            }
            for i in r.bond_set.members() {
                let b = &mol.bonds[i];
                if bond_types[i] == AROMATIC {
                    if sp_count[b.atoms[0] as usize] == 0 {
                        sp_count[b.atoms[0] as usize] += 1;
                    }
                    if sp_count[b.atoms[1] as usize] == 0 {
                        sp_count[b.atoms[1] as usize] += 1;
                    }
                }
            }
            let mut ring_is_aromatic = !is_cumulene && (r.cardinality - 2) % 4 == 0;
            for i in r.bond_set.members() {
                let b = &mol.bonds[i];
                if sp_count[b.atoms[0] as usize] != 1 || sp_count[b.atoms[1] as usize] != 1 {
                    ring_is_aromatic = false;
                }
            }
            if ring_is_aromatic && (ndouble > 0 || nsingle > 0) {
                for i in r.bond_set.members() {
                    if bond_is_in_ring[i] && bond_types[i] != AROMATIC {
                        changed = true;
                        bond_types[i] = AROMATIC;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
}

fn symbol_in(symbol: &str, list: &[&str]) -> bool {
    list.contains(&symbol)
}

/// `PerceiveDYAromaticity`: sets `bond_types` to `AROMATIC` for the rings
/// Daylight perceives as aromatic.
pub(crate) fn perceive_dy_aromaticity(
    mol: &AvalonMolecule,
    nbp: &[Neighbourhood],
    bond_types: &mut [i32],
) {
    let n_atoms = mol.atoms.len();
    let n_bonds = mol.bonds.len();
    let mut graph = bond_graph(mol);
    let first = ring_list(&graph);
    if first.is_empty() {
        return;
    }
    let mut bond_is_in_ring = vec![false; n_bonds];
    for r in &first {
        for m in r.bond_set.members() {
            bond_is_in_ring[m] = true;
        }
    }
    let mut atom_is_in_ring = vec![false; n_atoms];
    let mut candidate: Vec<bool> = mol.atoms.iter().map(|a| a.symbol != "C").collect();
    for (i, b) in mol.bonds.iter().enumerate() {
        if bond_types[i] > SINGLE && bond_types[i] != TRIPLE {
            candidate[(b.atoms[0] - 1) as usize] = true;
            candidate[(b.atoms[1] - 1) as usize] = true;
        }
    }
    for i in 0..n_bonds {
        let b = &mol.bonds[i];
        graph[i] = [b.atoms[0] as u32, b.atoms[1] as u32];
        if !bond_is_in_ring[i] {
            graph[i] = [0, 0];
        } else {
            atom_is_in_ring[(graph[i][0] - 1) as usize] = true;
            atom_is_in_ring[(graph[i][1] - 1) as usize] = true;
            if !candidate[(graph[i][0] - 1) as usize] {
                graph[i][0] = 0;
            }
            if !candidate[(graph[i][1] - 1) as usize] {
                graph[i][1] = 0;
            }
        }
    }
    let rl = ring_list(&graph);
    if rl.is_empty() {
        return;
    }
    let base = combine_rings(rl);
    let pairs = proper_ring_pairs(&base, n_atoms, &graph);
    let mut ring_list_all: Vec<RingNode> = pairs.into_iter().rev().collect();
    ring_list_all.extend(base);

    // Fused pairs that could also be aromatic, prepended to the candidates.
    let mut fused: Vec<RingNode> = Vec::new();
    for i in 0..ring_list_all.len() {
        for j in i + 1..ring_list_all.len() {
            let (p1, p2) = (&ring_list_all[i], &ring_list_all[j]);
            if p1.bond_set.intersection_is_empty(&p2.bond_set) {
                continue;
            }
            let mut set = p1.bond_set.clone();
            set.xor_with(&p2.bond_set);
            let c = set.cardinality();
            if c == p1.cardinality + p2.cardinality - 2 {
                fused.push(RingNode {
                    bond_set: set,
                    cardinality: c,
                });
            }
        }
    }
    let mut candidates: Vec<RingNode> = fused.into_iter().rev().collect();
    candidates.extend(ring_list_all);

    loop {
        let mut changed = false;
        for r in &candidates {
            let mut npi = 0;
            let mut conjugated = true;
            for i in 0..n_atoms {
                if !atom_is_in_ring[i] {
                    continue;
                }
                let ap = &mol.atoms[i];
                let mut in_ring_double = 0;
                let mut in_ring_aromatic = 0;
                let mut exo_pull = false;
                let mut is_in_ring = false;
                let mut local_pi = 0;
                for j in 0..nbp[i].atoms.len() {
                    let alp = &mol.atoms[nbp[i].atoms[j]];
                    let bj = nbp[i].bonds[j];
                    if bond_is_in_ring[bj] && r.bond_set.is_member(bj) {
                        is_in_ring = true;
                        if bond_types[bj] == AROMATIC {
                            in_ring_aromatic += 1;
                        }
                        if bond_types[bj] == DOUBLE {
                            in_ring_double += 1;
                        }
                    } else if bond_types[bj] == DOUBLE
                        && ap.symbol == "C"
                        && !bond_is_in_ring[bj]
                        && symbol_in(&alp.symbol, &["O", "S", "P", "N", "L"])
                    {
                        exo_pull = true;
                    }
                }
                if !is_in_ring {
                    continue;
                }
                if (in_ring_aromatic >= 1 || in_ring_double == 1)
                    && (symbol_in(&ap.symbol, &["C", "N", "A", "*"]) || ap.symbol == "L")
                {
                    local_pi = 1;
                } else if in_ring_aromatic == 0
                    && in_ring_double == 0
                    && ap.charge == 0
                    && symbol_in(&ap.symbol, &["N", "S", "O"])
                {
                    // Symbol-list ("L") atoms never come from an RDKit
                    // molecule, so the list branch of the C condition is
                    // not reachable here.
                    local_pi = 2;
                } else if in_ring_aromatic == 0
                    && in_ring_double == 0
                    && ap.charge == 0
                    && exo_pull
                    && ap.symbol == "C"
                {
                    local_pi = 0;
                } else {
                    conjugated = false;
                }
                if ap.charge < 0 && symbol_in(&ap.symbol, &["C", "N"]) {
                    conjugated = false;
                }
                npi += local_pi;
            }
            if !conjugated {
                continue;
            }
            if npi % 4 != 2 {
                continue;
            }
            for i in 0..n_bonds {
                if bond_is_in_ring[i] && r.bond_set.is_member(i) && bond_types[i] != AROMATIC {
                    bond_types[i] = AROMATIC;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    let mut arom_atom = vec![false; n_atoms];
    for (i, b) in mol.bonds.iter().enumerate() {
        if bond_types[i] == AROMATIC {
            arom_atom[(b.atoms[0] - 1) as usize] = true;
            arom_atom[(b.atoms[1] - 1) as usize] = true;
        }
    }
    for (i, b) in mol.bonds.iter().enumerate() {
        if arom_atom[(b.atoms[0] - 1) as usize]
            && arom_atom[(b.atoms[1] - 1) as usize]
            && bond_is_in_ring[i]
            && bond_types[i] == SINGLE
        {
            bond_types[i] = AROMATIC;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GlibcRand;

    #[test]
    fn glibc_rand_matches_srand_1() {
        // First outputs of glibc `rand()` after `srand(1)`.
        let mut g = GlibcRand::new(1);
        let got: Vec<i32> = (0..5).map(|_| g.rand()).collect();
        assert_eq!(
            got,
            vec![1804289383, 846930886, 1681692777, 1714636915, 1957747793]
        );
    }
}
