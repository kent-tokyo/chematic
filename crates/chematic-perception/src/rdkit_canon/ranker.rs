const STEREONONE: u32 = 0;
const STEREOANY: u32 = 1;
const STEREOZ: u32 = 2;
const STEREOE: u32 = 3;
const STEREOCIS: u32 = 4;
const STEREOTRANS: u32 = 5;
/// RDKit `Bond::AROMATIC`.
const BT_AROMATIC: u32 = 12;

#[derive(Clone)]
struct BondHolder {
    bond_type: u32,
    bond_stereo: u32,
    nbr_sym_class: u32,
    nbr_idx: u32,
    stype: u32,
    ctrl: [Option<u32>; 4],
    /// Bond index (for `Input::bond_symbols`).
    bidx: u32,
}

#[derive(Default)]
struct CanonAtom {
    index: u32,
    /// `Ranker::epoch` when `bonds` was last refreshed by `update_nbr_index`.
    nbr_epoch: u64,
    degree: u32,
    total_hs: u32,
    has_ring_nbr: bool,
    is_ring_stereo: bool,
    // Inline for the usual degree <= 4 (no per-atom allocation).
    nbr_ids: SmallVec<[u32; 4]>,
    bonds: Vec<BondHolder>,
    neighbor_num: Vec<i32>,
    revisited: Vec<i32>,
}

/// A bond as the ranking sees it: its atoms, RDKit bond type, RDKit bond
/// stereo and, for cis/trans stereo, its stereo atoms.
type InputBond = (u32, u32, u32, u32, Option<(u32, u32)>);

/// Per-atom and per-bond invariants the ranking reads, in RDKit's terms.
struct Input {
    anum: Vec<u32>,
    isotope: Vec<u32>,
    charge: Vec<i32>,
    map: Vec<i32>,
    /// 0 none, 1 `CHI_TETRAHEDRAL_CW`, 2 `CHI_TETRAHEDRAL_CCW` (relative to
    /// `nbrs` order), 3 any other tag.
    chiral: Vec<u8>,
    nrings: Vec<u32>,
    /// RDKit's `_ringStereoAtoms` on a `CW`/`CCW` atom (ring stereo).
    ring_stereo: Vec<bool>,
    /// Adjacency in bond insertion order.
    nbrs: Vec<SmallVec<[u32; 4]>>,
    bonds: Vec<InputBond>,
    /// `rankFragmentAtoms`' `atomSymbols` / `bondSymbols`: when given, atoms
    /// compare by (class, degree, symbol) only and bonds by symbol first.
    atom_symbols: Option<Vec<String>>,
    bond_symbols: Option<Vec<String>>,
    /// `breakTies`.
    break_ties: bool,
    /// Run the special-chirality refinement on ties (`useChirality &&
    /// includeRingStereo`).
    special_chirality: bool,
    /// `rankFragmentAtoms`' `atomsInPlay` / `bondsInPlay`.
    in_play: Option<(Vec<bool>, Vec<bool>)>,
}

fn sign(v: i64) -> i32 {
    v.signum() as i32
}

fn cmp_u(a: u32, b: u32) -> i32 {
    sign(i64::from(a) - i64::from(b))
}

struct Ranker<'a> {
    inp: &'a Input,
    atoms: Vec<CanonAtom>,
    use_nbrs: bool,
    /// Classes touched since the last activation (see `activate_touched`).
    touched_list: Vec<usize>,
    /// Reused partition/scratch buffer for `hanoi`.
    scratch: Vec<u32>,
    /// Bumped whenever an atom's `index` changes: a neighbour list refreshed
    /// at the current epoch is already up to date.
    epoch: u64,
}

#[derive(Clone, Copy)]
enum Functor {
    Atom,
    SpecialChirality,
    SpecialSymmetry,
}

/// RDKit `countSwapsToInterconvert(ref, probe)`.
fn count_swaps(reference: &[u32], probe: &[u32]) -> usize {
    let mut probe: SmallVec<[u32; 4]> = SmallVec::from_slice(probe);
    let mut n = 0;
    for i in 0..reference.len().min(probe.len()) {
        if probe[i] != reference[i]
            && let Some(j) = (i..probe.len()).find(|&j| probe[j] == reference[i])
        {
            probe.swap(i, j);
            n += 1;
        }
    }
    n
}

impl Ranker<'_> {
    fn new(inp: &Input) -> Ranker<'_> {
        let n = inp.anum.len();
        let mut atoms: Vec<CanonAtom> = (0..n)
            .map(|i| CanonAtom {
                index: i as u32,
                ..CanonAtom::default()
            })
            .collect();
        // `isRingStereoAtom` / `hasRingNbr` (`advancedInitCanonAtom`).
        // Atoms outside `atomsInPlay` are left as `initFragmentCanonAtoms`
        // leaves them (no `advancedInitCanonAtom`).
        for (i, atom) in atoms.iter_mut().enumerate() {
            if inp.in_play.as_ref().is_some_and(|(a, _)| !a[i]) {
                continue;
            }
            atom.bonds.reserve_exact(inp.nbrs[i].len());
            atom.is_ring_stereo = inp.ring_stereo[i];
            atom.has_ring_nbr = inp.nbrs[i].iter().any(|&nb| inp.ring_stereo[nb as usize]);
        }
        let mut ranker = Ranker {
            inp,
            atoms,
            use_nbrs: false,
            touched_list: Vec::new(),
            scratch: Vec::new(),
            epoch: 1,
        };
        for bidx in 0..inp.bonds.len() {
            let (a, b, ..) = inp.bonds[bidx];
            if let Some((atoms_in, bonds_in)) = &inp.in_play
                && !(bonds_in[bidx] && atoms_in[a as usize] && atoms_in[b as usize])
            {
                continue;
            }
            for (x, y) in [(a, b), (b, a)] {
                let h = ranker.make_holder(bidx, y);
                let at = &mut ranker.atoms[x as usize];
                at.nbr_ids.push(y);
                at.degree += 1;
                at.bonds.push(h);
            }
        }
        for i in 0..n {
            ranker.atoms[i].total_hs = 0;
            let mut bonds = std::mem::take(&mut ranker.atoms[i].bonds);
            // std::sort with bondholder::greater; equal holders are
            // interchangeable for every later comparison.
            bonds.sort_by(|x, y| ranker.bh_compare(y, x).cmp(&0));
            ranker.atoms[i].bonds = bonds;
        }
        ranker
    }

    fn make_holder(&self, bidx: usize, other: u32) -> BondHolder {
        let (a, b, bt, st, satoms) = self.inp.bonds[bidx];
        let mut h = BondHolder {
            bond_type: bt,
            bond_stereo: st,
            nbr_sym_class: 0,
            nbr_idx: other,
            stype: st,
            ctrl: [None; 4],
            bidx: bidx as u32,
        };
        if (st == STEREOCIS || st == STEREOTRANS)
            && let Some((s0, s1)) = satoms
        {
            h.ctrl[0] = Some(s0);
            h.ctrl[2] = Some(s1);
            let nb = &self.inp.nbrs;
            if nb[a as usize].len() > 2 {
                for &x in &nb[a as usize] {
                    if x != b && x != s0 {
                        h.ctrl[1] = Some(x);
                    }
                }
            }
            if nb[b as usize].len() > 2 {
                for &x in &nb[b as usize] {
                    if x != a && x != s1 {
                        h.ctrl[3] = Some(x);
                    }
                }
            }
        }
        h
    }

    fn flip_if_needed(&self, st: u32, ctrl: &[Option<u32>; 4]) -> u32 {
        let idx = |c: Option<u32>| c.map(|a| self.atoms[a as usize].index);
        let mut flip = false;
        if let (Some(c1), Some(c0)) = (idx(ctrl[1]), idx(ctrl[0]))
            && c1 > c0
        {
            flip = !flip;
        }
        if let (Some(c3), Some(c2)) = (idx(ctrl[3]), idx(ctrl[2]))
            && c3 > c2
        {
            flip = !flip;
        }
        match (flip, st) {
            (true, STEREOCIS) => STEREOTRANS,
            (true, STEREOTRANS) => STEREOCIS,
            _ => st,
        }
    }

    fn compare_stereo(&self, x: &BondHolder, y: &BondHolder) -> i32 {
        let (st1, st2) = (x.stype, y.stype);
        if st1 == STEREONONE {
            return if st2 == STEREONONE { 0 } else { -1 };
        }
        if st2 == STEREONONE {
            return 1;
        }
        if st1 == STEREOANY {
            return if st2 == STEREOANY { 0 } else { -1 };
        }
        if st2 == STEREOANY {
            return 1;
        }
        if matches!(st1, STEREOE | STEREOZ) && matches!(st2, STEREOE | STEREOZ) {
            return cmp_u(st1, st2);
        }
        cmp_u(
            self.flip_if_needed(st1, &x.ctrl),
            self.flip_if_needed(st2, &y.ctrl),
        )
    }

    fn bh_compare(&self, x: &BondHolder, y: &BondHolder) -> i32 {
        if let Some(sym) = &self.inp.bond_symbols {
            match sym[x.bidx as usize].cmp(&sym[y.bidx as usize]) {
                std::cmp::Ordering::Less => return -1,
                std::cmp::Ordering::Greater => return 1,
                std::cmp::Ordering::Equal => {}
            }
        }
        let c = cmp_u(x.bond_type, y.bond_type);
        if c != 0 {
            return c;
        }
        let c = cmp_u(x.bond_stereo, y.bond_stereo);
        if c != 0 {
            return c;
        }
        let c = cmp_u(x.nbr_sym_class, y.nbr_sym_class);
        if c != 0 {
            return c;
        }
        if x.bond_stereo != 0 && y.bond_stereo != 0 {
            return self.compare_stereo(x, y);
        }
        0
    }

    fn set_index(&mut self, i: usize, value: u32) {
        if self.atoms[i].index != value {
            self.atoms[i].index = value;
            self.epoch += 1;
        }
    }

    /// `updateAtomNeighborIndex`: refresh neighbour classes, re-sort
    /// descending. The result depends only on the atoms' current `index`
    /// values, so a list refreshed since the last change is kept.
    fn update_nbr_index(&mut self, i: usize) {
        if self.atoms[i].nbr_epoch == self.epoch {
            return;
        }
        self.atoms[i].nbr_epoch = self.epoch;
        let mut nbrs = std::mem::take(&mut self.atoms[i].bonds);
        for nb in &mut nbrs {
            nb.nbr_sym_class = self.atoms[nb.nbr_idx as usize].index;
        }
        for k in 1..nbrs.len() {
            if self.bh_compare(&nbrs[k], &nbrs[k - 1]) <= 0 {
                continue;
            }
            let v = nbrs[k].clone();
            let mut j = k;
            loop {
                nbrs[j] = nbrs[j - 1].clone();
                j -= 1;
                if !(j > 0 && self.bh_compare(&v, &nbrs[j - 1]) > 0) {
                    break;
                }
            }
            nbrs[j] = v;
        }
        self.atoms[i].bonds = nbrs;
    }

    fn chiral_rank(&self, i: usize) -> u32 {
        let mut perm: SmallVec<[u32; 4]> = SmallVec::new();
        for &x in &self.inp.nbrs[i] {
            let r = self.atoms[x as usize].index;
            if perm.contains(&r) {
                break;
            }
            perm.push(r);
        }
        if perm.len() != self.inp.nbrs[i].len() {
            return 0;
        }
        let ct = self.inp.chiral[i];
        if ct != 1 && ct != 2 {
            return 0;
        }
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        let swaps = count_swaps(&perm, &sorted);
        let mut res = if ct == 1 { 2 } else { 1 };
        if swaps % 2 == 1 {
            res = if res == 2 { 1 } else { 2 };
        }
        res
    }

    fn ring_nbr_code(&self, i: usize) -> u32 {
        if !self.atoms[i].has_ring_nbr {
            return 0;
        }
        let mut code: u32 = 0;
        for &j in &self.atoms[i].nbr_ids {
            if self.atoms[j as usize].is_ring_stereo {
                code = code.wrapping_add(self.atoms[j as usize].index.wrapping_mul(10000) + 1);
            }
        }
        code
    }

    fn basecomp(&self, i: usize, j: usize) -> i32 {
        let (a, b) = (&self.atoms[i], &self.atoms[j]);
        let inp = self.inp;
        if let Some(sym) = &inp.atom_symbols {
            for c in [
                cmp_u(a.index, b.index),
                sign(i64::from(inp.map[i]) - i64::from(inp.map[j])),
                cmp_u(a.degree, b.degree),
            ] {
                if c != 0 {
                    return c;
                }
            }
            return match sym[i].cmp(&sym[j]) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Greater => 1,
                std::cmp::Ordering::Equal => 0,
            };
        }
        let checks = [
            cmp_u(a.index, b.index),
            sign(i64::from(inp.map[i]) - i64::from(inp.map[j])),
            cmp_u(a.degree, b.degree),
            cmp_u(inp.anum[i], inp.anum[j]),
            cmp_u(inp.isotope[i], inp.isotope[j]),
            cmp_u(a.total_hs, b.total_hs),
            // RDKit compares the formal charge as `unsigned int`.
            cmp_u(inp.charge[i] as u32, inp.charge[j] as u32),
        ];
        for c in checks {
            if c != 0 {
                return c;
            }
        }
        let (ci, cj) = (inp.chiral[i] != 0, inp.chiral[j] != 0);
        if ci != cj {
            return if ci { 1 } else { -1 };
        }
        if ci && cj {
            let c = cmp_u(self.chiral_rank(i), self.chiral_rank(j));
            if c != 0 {
                return c;
            }
        }
        cmp_u(self.ring_nbr_code(i), self.ring_nbr_code(j))
    }

    fn compare_bond_lists(&self, i: usize, j: usize, size_too: bool) -> i32 {
        let (a, b) = (&self.atoms[i].bonds, &self.atoms[j].bonds);
        for (x, y) in a.iter().zip(b.iter()) {
            let c = self.bh_compare(x, y);
            if c != 0 {
                return c;
            }
        }
        if size_too {
            return sign(a.len() as i64 - b.len() as i64);
        }
        0
    }

    fn nbr_num_swaps(&self, atom_idx: usize) -> Vec<(u32, u32)> {
        let inp = self.inp;
        let is_ring = inp.nrings[atom_idx] > 0;
        let mut res = Vec::new();
        for nb in &self.atoms[atom_idx].bonds {
            let n = nb.nbr_idx as usize;
            if is_ring && inp.chiral[n] != 0 {
                let mut seen: Vec<u32> = Vec::new();
                let mut too_many = false;
                let mut reference = Vec::new();
                for &nn in &self.atoms[n].nbr_ids {
                    reference.push(nn);
                    if nn as usize != atom_idx {
                        let cls = self.atoms[nn as usize].index;
                        if seen.contains(&cls) {
                            too_many = true;
                        } else {
                            seen.push(cls);
                        }
                    }
                }
                let mut probe = vec![atom_idx as u32];
                probe.extend(
                    self.atoms[n]
                        .bonds
                        .iter()
                        .map(|b| b.nbr_idx)
                        .filter(|&x| x as usize != atom_idx),
                );
                if too_many {
                    res.push((nb.nbr_sym_class, 0));
                } else {
                    let odd = count_swaps(&reference, &probe) % 2 == 1;
                    match inp.chiral[n] {
                        1 => res.push((nb.nbr_sym_class, if odd { 2 } else { 1 })),
                        2 => res.push((nb.nbr_sym_class, if odd { 1 } else { 2 })),
                        _ => {}
                    }
                }
            } else {
                res.push((nb.nbr_sym_class, 0));
            }
        }
        res.sort_unstable();
        res
    }

    fn out_of_play(&self, i: usize) -> bool {
        self.inp.in_play.as_ref().is_some_and(|(a, _)| !a[i])
    }

    fn compare(&mut self, f: Functor, i: usize, j: usize) -> i32 {
        match f {
            _ if self.out_of_play(i) && self.out_of_play(j) => 0,
            Functor::Atom => {
                let v = self.basecomp(i, j);
                if v != 0 {
                    return v;
                }
                if self.use_nbrs {
                    self.update_nbr_index(i);
                    self.update_nbr_index(j);
                    return self.compare_bond_lists(i, j, true);
                }
                0
            }
            Functor::SpecialChirality => {
                self.update_nbr_index(i);
                self.update_nbr_index(j);
                let v = self.compare_bond_lists(i, j, false);
                if v != 0 {
                    return v;
                }
                let si = self.nbr_num_swaps(i);
                let sj = self.nbr_num_swaps(j);
                for (x, y) in si.iter().zip(sj.iter()) {
                    let c = cmp_u(x.1, y.1);
                    if c != 0 {
                        return c;
                    }
                }
                0
            }
            Functor::SpecialSymmetry => {
                let (a, b) = (&self.atoms[i], &self.atoms[j]);
                match a.neighbor_num.cmp(&b.neighbor_num) {
                    std::cmp::Ordering::Less => return -1,
                    std::cmp::Ordering::Greater => return 1,
                    _ => {}
                }
                match a.revisited.cmp(&b.revisited) {
                    std::cmp::Ordering::Less => return -1,
                    std::cmp::Ordering::Greater => return 1,
                    _ => {}
                }
                self.update_nbr_index(i);
                self.update_nbr_index(j);
                self.compare_bond_lists(i, j, true)
            }
        }
    }

    /// `compareRingAtomsConcerningNumNeighbors`.
    fn ring_atoms_num_neighbors(&mut self) {
        let n = self.atoms.len();
        let nrings = &self.inp.nrings;
        let mut cur = vec![false; n];
        let mut revis = vec![0i32; n];
        for idx in 0..n {
            if nrings[idx] < 1 {
                continue;
            }
            let mut visited = vec![false; n];
            let mut last = vec![false; n];
            let mut neighbors: VecDeque<usize> = VecDeque::from([idx]);
            let mut neighbor_num = Vec::new();
            let mut rv = Vec::new();
            while !neighbors.is_empty() {
                let mut num_level = 0i32;
                let mut next_level: Vec<usize> = Vec::new();
                while let Some(nidx) = neighbors.pop_front() {
                    if nrings[nidx] < 1 {
                        continue;
                    }
                    last[nidx] = true;
                    visited[nidx] = true;
                    for &iidx in &self.atoms[nidx].nbr_ids {
                        let iidx = iidx as usize;
                        if !visited[iidx] {
                            cur[iidx] = true;
                            num_level += 1;
                            visited[iidx] = true;
                            next_level.push(iidx);
                        }
                    }
                }
                let mut modified: Vec<usize> = Vec::new();
                for &i2 in &next_level {
                    for &jidx in &self.atoms[i2].nbr_ids {
                        let jidx = jidx as usize;
                        if cur[jidx] || last[jidx] {
                            if revis[jidx] == 0 {
                                modified.push(jidx);
                            }
                            revis[jidx] += 1;
                        }
                    }
                }
                last.iter_mut().for_each(|x| *x = false);
                for &i2 in &next_level {
                    last[i2] = true;
                    cur[i2] = false;
                }
                let mut tmp: Vec<i32> = modified.iter().map(|&i2| revis[i2]).collect();
                tmp.sort_unstable();
                tmp.push(-1);
                rv.extend(tmp);
                for &i2 in &modified {
                    revis[i2] = 0;
                }
                neighbor_num.push(num_level);
                neighbor_num.push(-1);
                neighbors.extend(next_level);
            }
            self.atoms[idx].neighbor_num = neighbor_num;
            self.atoms[idx].revisited = rv;
        }
    }

    /// RDKit `detail::hanoi` over `buf`, whose `[0, len)` is the partition
    /// and `[len, 2 len)` the scratch area. Returns true when the result
    /// sits in the scratch area.
    #[allow(clippy::too_many_arguments)]
    fn hanoi(
        &mut self,
        f: Functor,
        buf: &mut [u32],
        base: usize,
        nel: usize,
        temp: usize,
        count: &mut [i32],
        changed: &mut [bool],
    ) -> bool {
        if nel == 1 {
            count[buf[base] as usize] = 1;
            return false;
        }
        if nel == 2 {
            let (n1, n2) = (buf[base] as usize, buf[base + 1] as usize);
            let stat = if changed[n1] || changed[n2] {
                self.compare(f, n1, n2)
            } else {
                0
            };
            if stat == 0 {
                count[n1] = 2;
                count[n2] = 0;
            } else {
                count[n1] = 1;
                count[n2] = 1;
                if stat > 0 {
                    buf.swap(base, base + 1);
                }
            }
            return false;
        }
        let mut n1 = nel / 2;
        let mut n2 = nel - n1;
        let (b1, t1, b2, t2) = (base, temp, base + n1, temp + n1);
        let (mut s1, mut s2, mut ptr, result);
        if self.hanoi(f, buf, b1, n1, t1, count, changed) {
            s2 = if self.hanoi(f, buf, b2, n2, t2, count, changed) {
                t2
            } else {
                b2
            };
            result = false;
            ptr = base;
            s1 = t1;
        } else {
            s2 = if self.hanoi(f, buf, b2, n2, t2, count, changed) {
                t2
            } else {
                b2
            };
            result = true;
            ptr = temp;
            s1 = b1;
        }
        loop {
            let (x, y) = (buf[s1] as usize, buf[s2] as usize);
            let stat = if changed[x] || changed[y] {
                self.compare(f, x, y)
            } else {
                0
            };
            let len1 = count[x] as usize;
            let len2 = count[y] as usize;
            if stat == 0 {
                count[x] = (len1 + len2) as i32;
                count[y] = 0;
                buf.copy_within(s1..s1 + len1, ptr);
                ptr += len1;
                n1 -= len1;
                if n1 == 0 {
                    if ptr != s2 {
                        buf.copy_within(s2..s2 + n2, ptr);
                    }
                    return result;
                }
                s1 += len1;
                buf.copy_within(s2..s2 + len2, ptr);
                ptr += len2;
                n2 -= len2;
                if n2 == 0 {
                    buf.copy_within(s1..s1 + n1, ptr);
                    return result;
                }
                s2 += len2;
            } else if stat < 0 {
                buf.copy_within(s1..s1 + len1, ptr);
                ptr += len1;
                n1 -= len1;
                if n1 == 0 {
                    if ptr != s2 {
                        buf.copy_within(s2..s2 + n2, ptr);
                    }
                    return result;
                }
                s1 += len1;
            } else {
                buf.copy_within(s2..s2 + len2, ptr);
                ptr += len2;
                n2 -= len2;
                if n2 == 0 {
                    buf.copy_within(s1..s1 + n1, ptr);
                    return result;
                }
                s2 += len2;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn refine(
        &mut self,
        f: Functor,
        order: &mut [u32],
        count: &mut [i32],
        activeset: &mut i64,
        next: &mut [i64],
        changed: &mut [bool],
        touched: &mut [bool],
    ) {
        let mut buf: Vec<u32> = std::mem::take(&mut self.scratch);
        while *activeset != -1 {
            let partition = *activeset as usize;
            *activeset = next[partition];
            next[partition] = -2;
            let len = count[partition] as usize;
            let offset = self.atoms[partition].index as usize;
            buf.clear();
            buf.extend_from_slice(&order[offset..offset + len]);
            buf.resize(2 * len, 0);
            let in_temp = self.hanoi(f, &mut buf, 0, len, len, count, changed);
            let src = if in_temp { len } else { 0 };
            order[offset..offset + len].copy_from_slice(&buf[src..src + len]);
            for k in 0..len {
                changed[order[offset + k] as usize] = false;
            }
            let mut index = order[offset] as usize;
            let mut symclass = 0u32;
            let mut i = count[index] as usize;
            while i < len {
                index = order[offset + i] as usize;
                if count[index] != 0 {
                    symclass = (offset + i) as u32;
                }
                self.set_index(index, symclass);
                for &nb in &self.atoms[index].nbr_ids {
                    changed[nb as usize] = true;
                }
                i += 1;
            }
            index = order[offset] as usize;
            let mut i = count[index] as usize;
            while i < len {
                index = order[offset + i] as usize;
                for &nb in &self.atoms[index].nbr_ids {
                    let cls = self.atoms[nb as usize].index as usize;
                    if !touched[cls] {
                        touched[cls] = true;
                        self.touched_list.push(cls);
                    }
                }
                i += 1;
            }
            self.activate_touched(order, count, activeset, next, touched);
        }
        self.scratch = buf;
    }

    /// RDKit walks every class index in ascending order and pushes each
    /// touched class that still holds ties; the touched ones, sorted, are
    /// the same walk.
    fn activate_touched(
        &mut self,
        order: &[u32],
        count: &[i32],
        activeset: &mut i64,
        next: &mut [i64],
        touched: &mut [bool],
    ) {
        self.touched_list.sort_unstable();
        for &ii in &self.touched_list {
            let partition = order[ii] as usize;
            if count[partition] > 1 && next[partition] == -2 {
                next[partition] = *activeset;
                *activeset = partition as i64;
            }
            touched[ii] = false;
        }
        self.touched_list.clear();
    }

    #[allow(clippy::too_many_arguments)]
    fn break_ties(
        &mut self,
        order: &mut [u32],
        count: &mut [i32],
        activeset: &mut i64,
        next: &mut [i64],
        changed: &mut [bool],
        touched: &mut [bool],
    ) {
        let n = self.atoms.len();
        let mut i: usize = 0;
        while i < n {
            let partition = order[i] as usize;
            let old = self.atoms[partition].index;
            while count[partition] > 1 {
                let len = count[partition] as usize;
                let offset = self.atoms[partition].index as usize + len - 1;
                let index = order[offset] as usize;
                self.set_index(index, offset as u32);
                count[partition] = (len - 1) as i32;
                count[index] = 1;
                if self.atoms[index].degree < 1 {
                    continue;
                }
                for k in 0..self.atoms[index].nbr_ids.len() {
                    let nb = self.atoms[index].nbr_ids[k] as usize;
                    let cls = self.atoms[nb].index as usize;
                    if !touched[cls] {
                        touched[cls] = true;
                        self.touched_list.push(cls);
                    }
                    changed[nb] = true;
                }
                self.activate_touched(order, count, activeset, next, touched);
                self.refine(
                    Functor::Atom,
                    order,
                    count,
                    activeset,
                    next,
                    changed,
                    touched,
                );
            }
            if self.atoms[partition].index != old {
                // RDKit decrements its unsigned loop index and lets `++` undo it.
                i = i.wrapping_sub(1);
            }
            i = i.wrapping_add(1);
        }
    }

    fn activate(
        order: &[u32],
        count: &[i32],
        activeset: &mut i64,
        next: &mut [i64],
        changed: &mut [bool],
    ) {
        let n = order.len();
        *activeset = -1;
        next.iter_mut().for_each(|x| *x = -2);
        let mut i = 0usize;
        loop {
            let j = order[i] as usize;
            if count[j] > 1 {
                next[j] = *activeset;
                *activeset = j as i64;
                i += count[j] as usize;
            } else {
                i += 1;
            }
            if i >= n {
                break;
            }
        }
        changed.iter_mut().for_each(|x| *x = true);
    }

    fn rank(mut self) -> Vec<u32> {
        let n = self.atoms.len();
        if n == 0 {
            return Vec::new();
        }
        let mut order: Vec<u32> = (0..n as u32).collect();
        let mut count = vec![0i32; n];
        let mut next = vec![-2i64; n];
        let mut changed = vec![true; n];
        let mut touched = vec![false; n];
        let mut activeset: i64 = -1;
        for a in &mut self.atoms {
            a.index = 0;
        }
        self.epoch += 1;
        count[0] = n as i32;
        self.use_nbrs = true;
        Self::activate(&order, &count, &mut activeset, &mut next, &mut changed);
        self.refine(
            Functor::Atom,
            &mut order,
            &mut count,
            &mut activeset,
            &mut next,
            &mut changed,
            &mut touched,
        );
        if self.inp.special_chirality && count.contains(&0) {
            Self::activate(&order, &count, &mut activeset, &mut next, &mut changed);
            self.refine(
                Functor::SpecialChirality,
                &mut order,
                &mut count,
                &mut activeset,
                &mut next,
                &mut changed,
                &mut touched,
            );
        }
        let nrings = &self.inp.nrings;
        let (mut sym, mut ring) = (0usize, 0usize);
        let mut branching = false;
        let mut ties = false;
        for i in 0..n {
            let o = order[i] as usize;
            if nrings[o] > 0 {
                if count[o] > 2 {
                    sym += count[o] as usize;
                }
                ring += 1;
                if nrings[o] > 1 && count[o] > 1 {
                    branching = true;
                }
            }
            if count[i] == 0 {
                ties = true;
            }
        }
        if ties && ring > 0 && (sym as f32) / (ring as f32) > 0.5 && branching {
            self.ring_atoms_num_neighbors();
            Self::activate(&order, &count, &mut activeset, &mut next, &mut changed);
            self.refine(
                Functor::SpecialSymmetry,
                &mut order,
                &mut count,
                &mut activeset,
                &mut next,
                &mut changed,
                &mut touched,
            );
        }
        if self.inp.break_ties {
            self.break_ties(
                &mut order,
                &mut count,
                &mut activeset,
                &mut next,
                &mut changed,
                &mut touched,
            );
        }
        let mut res = vec![0u32; n];
        for &o in &order {
            res[o as usize] = self.atoms[o as usize].index;
        }
        res
    }
}

