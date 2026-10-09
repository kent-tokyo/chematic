//! RDKit's VF2 substructure enumeration (`Substruct/vf2.hpp`,
//! `vf2_all` without node sorting and without `RDK_VF2_PRUNING`) on the
//! RDKit-model molecule, in RDKit's match order.
//!
//! Graphs are walked as RDKit's `boost::adjacency_list` walks them: a
//! vertex's out-edges in bond-insertion order (`Mol::atom_bonds`), and
//! `boost::edge(u, v)` as the first such bond to `v`.

use super::mol::Mol;

const NULL_NODE: usize = usize::MAX;

/// A VF2 state (`VF2SubState`) matching `g1` (the query) into `g2`.
struct State<'a, VC, EC> {
    g1: &'a Mol,
    g2: &'a Mol,
    vc: &'a VC,
    ec: &'a EC,
    n1: usize,
    n2: usize,
    core_len: usize,
    t1_len: usize,
    t2_len: usize,
    core_1: Vec<usize>,
    core_2: Vec<usize>,
    term_1: Vec<usize>,
    term_2: Vec<usize>,
}

/// `Pair`: the candidate pair iterator of one `MatchAll` frame.
struct Pair {
    n1: usize,
    n2: usize,
    /// The VF2-Plus neighbour iterator: `g2`'s neighbours of the image of
    /// a mapped neighbour of `n1` (a vertex and its position).
    iter: Option<(usize, usize)>,
}

impl<VC, EC> State<'_, VC, EC>
where
    VC: Fn(usize, usize) -> bool,
    EC: Fn(usize, usize) -> bool,
{
    fn next_pair(&self, pair: &mut Pair) -> bool {
        if pair.n1 == NULL_NODE {
            pair.n1 = 0;
        }
        if pair.n2 == NULL_NODE {
            pair.n2 = 0;
        } else {
            pair.n2 += 1;
        }
        if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while pair.n1 < self.n1
                && (self.core_1[pair.n1] != NULL_NODE || self.term_1[pair.n1] == 0)
            {
                pair.n1 += 1;
                pair.n2 = 0;
            }
            if pair.iter.is_none() {
                let anchor = self.g1.atom_bonds[pair.n1]
                    .iter()
                    .map(|&b| self.g1.bonds[b].other(pair.n1))
                    .find(|&nb| self.core_1[nb] != NULL_NODE)
                    .expect("a terminal query atom has a mapped neighbour");
                pair.iter = Some((self.core_1[anchor], 0));
            }
        } else {
            while pair.n1 < self.n1 && self.core_1[pair.n1] != NULL_NODE {
                pair.n1 += 1;
                pair.n2 = 0;
            }
        }
        if let Some((vertex, pos)) = pair.iter.as_mut() {
            let nbrs = &self.g2.atom_bonds[*vertex];
            while *pos < nbrs.len()
                && self.core_2[self.g2.bonds[nbrs[*pos]].other(*vertex)] != NULL_NODE
            {
                *pos += 1;
            }
            if *pos < nbrs.len() {
                pair.n2 = self.g2.bonds[nbrs[*pos]].other(*vertex);
                *pos += 1;
            } else {
                pair.n2 = self.n2;
            }
        } else if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while pair.n2 < self.n2
                && (self.core_2[pair.n2] != NULL_NODE || self.term_2[pair.n2] == 0)
            {
                pair.n2 += 1;
            }
        } else {
            while pair.n2 < self.n2 && self.core_2[pair.n2] != NULL_NODE {
                pair.n2 += 1;
            }
        }
        pair.n1 < self.n1 && pair.n2 < self.n2
    }

    fn is_feasible_pair(&self, node1: usize, node2: usize) -> bool {
        if self.g1.degree(node1) > self.g2.degree(node2) {
            return false;
        }
        if !(self.vc)(node1, node2) {
            return false;
        }
        for &b1 in &self.g1.atom_bonds[node1] {
            let other1 = self.g1.bonds[b1].other(node1);
            if self.core_1[other1] != NULL_NODE {
                let other2 = self.core_1[other1];
                match self.g2.bond_between(node2, other2) {
                    Some(b2) if (self.ec)(b1, b2) => {}
                    _ => return false,
                }
            }
        }
        true
    }

    fn add_pair(&mut self, node1: usize, node2: usize) {
        self.core_len += 1;
        if self.term_1[node1] == 0 {
            self.term_1[node1] = self.core_len;
            self.t1_len += 1;
        }
        if self.term_2[node2] == 0 {
            self.term_2[node2] = self.core_len;
            self.t2_len += 1;
        }
        self.core_1[node1] = node2;
        self.core_2[node2] = node1;
        for &b in &self.g1.atom_bonds[node1] {
            let other = self.g1.bonds[b].other(node1);
            if self.term_1[other] == 0 {
                self.term_1[other] = self.core_len;
                self.t1_len += 1;
            }
        }
        for &b in &self.g2.atom_bonds[node2] {
            let other = self.g2.bonds[b].other(node2);
            if self.term_2[other] == 0 {
                self.term_2[other] = self.core_len;
                self.t2_len += 1;
            }
        }
    }

    fn back_track(&mut self, node1: usize, node2: usize) {
        if self.term_1[node1] == self.core_len {
            self.term_1[node1] = 0;
            self.t1_len -= 1;
        }
        for &b in &self.g1.atom_bonds[node1] {
            let other = self.g1.bonds[b].other(node1);
            if self.term_1[other] == self.core_len {
                self.term_1[other] = 0;
                self.t1_len -= 1;
            }
        }
        if self.term_2[node2] == self.core_len {
            self.term_2[node2] = 0;
            self.t2_len -= 1;
        }
        for &b in &self.g2.atom_bonds[node2] {
            let other = self.g2.bonds[b].other(node2);
            if self.term_2[other] == self.core_len {
                self.term_2[other] = 0;
                self.t2_len -= 1;
            }
        }
        self.core_1[node1] = NULL_NODE;
        self.core_2[node2] = NULL_NODE;
        self.core_len -= 1;
    }

    /// `MatchAll`: true when `lim` matches have been collected.
    fn match_all(&mut self, res: &mut Vec<Vec<(usize, usize)>>, lim: usize) -> bool {
        if self.core_len == self.n1 {
            // GetCoreSet: query atoms in index order (the final check
            // accepts every match without chirality or uniquification).
            res.push(
                (0..self.n1)
                    .filter(|&i| self.core_1[i] != NULL_NODE)
                    .map(|i| (i, self.core_1[i]))
                    .collect(),
            );
            return lim != 0 && res.len() >= lim;
        }
        if self.n1 > self.n2 || self.t1_len > self.t2_len {
            return false;
        }
        let mut pair = Pair {
            n1: NULL_NODE,
            n2: NULL_NODE,
            iter: None,
        };
        while self.next_pair(&mut pair) {
            let (p1, p2) = (pair.n1, pair.n2);
            if self.is_feasible_pair(p1, p2) {
                self.add_pair(p1, p2);
                if self.match_all(res, lim) {
                    return true;
                }
                self.back_track(p1, p2);
            }
        }
        false
    }
}

/// `SubstructMatch(mol, query, params)` with `uniquify=false` and
/// `useChirality=false`: every match of `query` in `mol` as `(query atom,
/// mol atom)` pairs in query-atom order, in RDKit's enumeration order, at
/// most `max_matches` of them. `vc(query atom, mol atom)` and
/// `ec(query bond, mol bond)` are the atom and bond compatibility tests.
pub(crate) fn substruct_matches(
    mol: &Mol,
    query: &Mol,
    vc: &impl Fn(usize, usize) -> bool,
    ec: &impl Fn(usize, usize) -> bool,
    max_matches: usize,
) -> Vec<Vec<(usize, usize)>> {
    let (n1, n2) = (query.atoms.len(), mol.atoms.len());
    let mut res = Vec::new();
    if n2 == 0 || n1 == 0 || n1 > n2 {
        return res;
    }
    let mut state = State {
        g1: query,
        g2: mol,
        vc,
        ec,
        n1,
        n2,
        core_len: 0,
        t1_len: 0,
        t2_len: 0,
        core_1: vec![NULL_NODE; n1],
        core_2: vec![NULL_NODE; n2],
        term_1: vec![0; n1],
        term_2: vec![0; n2],
    };
    state.match_all(&mut res, max_matches);
    res
}
