//! RDKit's `SubstructMatch(mol, query, matches, uniquify=false,
//! recursionPossible=true)` (`maxMatches` 1000): a port of RDKit's VF2
//! (`Code/GraphMol/Substruct/vf2.hpp`, no node sorting, "VF2 Plus"
//! neighbour candidates) so matches come out in RDKit's order, with SMARTS
//! primitives evaluated with RDKit's semantics on an [`RdkitMolView`].

use std::collections::HashMap;

use chematic_core::Element;
use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};
use chematic_smiles::RdkitMolView;

const NULL: usize = usize::MAX;
pub(crate) const MAX_MATCHES: usize = 1000;

/// Per-target-molecule data the predicates need.
pub(crate) struct Target<'a> {
    pub v: &'a RdkitMolView,
    min_ring_size: Vec<usize>,
    num_atom_rings: Vec<usize>,
    ring_bond_count: Vec<usize>,
    num_bond_rings: Vec<usize>,
    /// Recursive-query root sets, keyed by the query's address.
    recursive: std::cell::RefCell<HashMap<usize, Vec<bool>>>,
}

impl<'a> Target<'a> {
    pub(crate) fn new(v: &'a RdkitMolView) -> Self {
        let na = v.num_atoms();
        let nb = v.num_bonds();
        let mut min_ring_size = vec![0; na];
        let mut num_atom_rings = vec![0; na];
        for r in &v.atom_rings {
            for &a in r {
                num_atom_rings[a] += 1;
                if min_ring_size[a] == 0 || r.len() < min_ring_size[a] {
                    min_ring_size[a] = r.len();
                }
            }
        }
        let mut num_bond_rings = vec![0; nb];
        for r in &v.bond_rings {
            for &b in r {
                num_bond_rings[b] += 1;
            }
        }
        let ring_bond_count = (0..na)
            .map(|a| {
                v.atom_bonds[a]
                    .iter()
                    .filter(|&&b| num_bond_rings[b] > 0)
                    .count()
            })
            .collect();
        Target {
            v,
            min_ring_size,
            num_atom_rings,
            ring_bond_count,
            num_bond_rings,
            recursive: Default::default(),
        }
    }

    fn atom_matches(&self, q: &AtomQuery, a: usize) -> bool {
        match q {
            AtomQuery::Primitive(p) => self.primitive(p, a),
            AtomQuery::And(x, y) => self.atom_matches(x, a) && self.atom_matches(y, a),
            AtomQuery::Or(x, y) => self.atom_matches(x, a) || self.atom_matches(y, a),
            AtomQuery::Not(x) => !self.atom_matches(x, a),
        }
    }

    fn primitive(&self, p: &AtomPrimitive, a: usize) -> bool {
        let at = &self.v.atoms[a];
        match p {
            AtomPrimitive::AtomicNum(n) => at.atomic_num == u32::from(*n),
            AtomPrimitive::Symbol(s) => Element::from_symbol(s)
                .is_some_and(|e| u32::from(e.atomic_number()) == at.atomic_num),
            AtomPrimitive::Aromatic(x) => at.is_aromatic == *x,
            AtomPrimitive::Charge(c) => at.formal_charge == i32::from(*c),
            // queryAtomHCount: getTotalNumHs(true)
            AtomPrimitive::HCount(h) => at.total_num_hs_with_neighbors == u32::from(*h),
            // queryAtomImplicitHCount: getTotalNumHs(false)
            AtomPrimitive::ImplicitHCount(h) => at.total_num_hs == u32::from(*h),
            AtomPrimitive::Degree(d) => self.v.degree(a) == usize::from(*d),
            AtomPrimitive::RingMembership(r) => (self.num_atom_rings[a] > 0) == *r,
            AtomPrimitive::RingSize(n) => self.v.is_atom_in_ring_of_size(a, usize::from(*n)),
            AtomPrimitive::MinRingSize(n) => self.min_ring_size[a] == usize::from(*n),
            AtomPrimitive::Wildcard => true,
            AtomPrimitive::Recursive(sub) => self.recursive_root(sub, a),
            AtomPrimitive::RingBondCount(n) => self.ring_bond_count[a] == usize::from(*n),
            // ^n: RDKit's hybridization enum value n + 1 (S = 1, SP = 2, ...)
            AtomPrimitive::Hybridization(h) => u32::from(at.hybridization) == u32::from(*h) + 1,
            // getTotalDegree = degree + getTotalNumHs(false)
            AtomPrimitive::TotalConnectivity(x) => {
                self.v.degree(a) as u32 + at.total_num_hs == u32::from(*x)
            }
            AtomPrimitive::RingCount(n) => self.num_atom_rings[a] == usize::from(*n),
            // not used by the torsion tables; chirality is ignored (useChirality = false)
            AtomPrimitive::Chirality(_) | AtomPrimitive::Isotope(_) => true,
            AtomPrimitive::HeavyDegree(d) => {
                self.v
                    .neighbors(a)
                    .filter(|&n| self.v.atoms[n].atomic_num != 1)
                    .count()
                    == usize::from(*d)
            }
            AtomPrimitive::Valence(_)
            | AtomPrimitive::HeteroNeighborCount(_)
            | AtomPrimitive::AliphaticHeteroNeighborCount(_) => {
                unimplemented!("SMARTS primitive {p:?} is not used by the ETKDG torsion tables")
            }
        }
    }

    fn bond_matches(&self, q: &BondQuery, b: usize) -> bool {
        let bond = &self.v.bonds[b];
        match q {
            // implicit SMARTS bond: single or aromatic
            BondQuery::Any => bond.bond_type == 1 || bond.bond_type == 12,
            BondQuery::And(x, y) => self.bond_matches(x, b) && self.bond_matches(y, b),
            BondQuery::Or(x, y) => self.bond_matches(x, b) || self.bond_matches(y, b),
            BondQuery::Not(x) => !self.bond_matches(x, b),
            BondQuery::Primitive(p) => match p {
                BondPrimitive::Single => bond.bond_type == 1,
                BondPrimitive::Double => bond.bond_type == 2,
                BondPrimitive::Triple => bond.bond_type == 3,
                BondPrimitive::Quadruple => bond.bond_type == 4,
                BondPrimitive::Aromatic => bond.bond_type == 12,
                BondPrimitive::Any => true,
                BondPrimitive::Ring => self.num_bond_rings[b] > 0,
                other => unimplemented!("SMARTS bond primitive {other:?} in a torsion pattern"),
            },
        }
    }

    /// `RecursiveMatcher`: atoms that are the image of query atom 0 in one
    /// of the first `max(maxRecursiveMatches, maxMatches)` matches.
    fn recursive_root(&self, sub: &QueryMolecule, a: usize) -> bool {
        let key = sub as *const QueryMolecule as usize;
        if let Some(set) = self.recursive.borrow().get(&key) {
            return set[a];
        }
        let mut set = vec![false; self.v.num_atoms()];
        for m in substruct_match(self, sub) {
            set[m[0]] = true;
        }
        let r = set[a];
        self.recursive.borrow_mut().insert(key, set);
        r
    }
}

/// Adjacency lists of the query in RDKit's order (bond creation order).
fn query_adj(q: &QueryMolecule) -> Vec<Vec<(usize, usize)>> {
    let mut adj = vec![Vec::new(); q.atoms.len()];
    for (bi, b) in q.bonds.iter().enumerate() {
        adj[b.atom1].push((bi, b.atom2));
        adj[b.atom2].push((bi, b.atom1));
    }
    adj
}

struct State<'q, 't> {
    q: &'q QueryMolecule,
    qadj: Vec<Vec<(usize, usize)>>,
    t: &'t Target<'t>,
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

struct Pair {
    n1: usize,
    n2: usize,
    /// VF2-Plus candidate list (neighbours of a mapped atom) and position.
    nbrs: Option<(Vec<usize>, usize)>,
}

impl State<'_, '_> {
    fn next_pair(&self, p: &mut Pair) -> bool {
        if p.n1 == NULL {
            p.n1 = 0;
        }
        if p.n2 == NULL {
            p.n2 = 0;
        } else {
            p.n2 += 1;
        }
        if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while p.n1 < self.n1 && (self.core_1[p.n1] != NULL || self.term_1[p.n1] == 0) {
                p.n1 += 1;
                p.n2 = 0;
            }
            if p.nbrs.is_none() {
                let mapped = self.qadj[p.n1]
                    .iter()
                    .map(|&(_, o)| o)
                    .find(|&o| self.core_1[o] != NULL)
                    .expect("terminal query atom has a mapped neighbour");
                let list: Vec<usize> = self.t.v.neighbors(self.core_1[mapped]).collect();
                p.nbrs = Some((list, 0));
            }
        } else {
            while p.n1 < self.n1 && self.core_1[p.n1] != NULL {
                p.n1 += 1;
                p.n2 = 0;
            }
        }
        if let Some((list, pos)) = p.nbrs.as_mut() {
            while *pos < list.len() && self.core_2[list[*pos]] != NULL {
                *pos += 1;
            }
            if *pos < list.len() {
                p.n2 = list[*pos];
                *pos += 1;
            } else {
                p.n2 = self.n2;
            }
        } else if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while p.n2 < self.n2 && (self.core_2[p.n2] != NULL || self.term_2[p.n2] == 0) {
                p.n2 += 1;
            }
        } else {
            while p.n2 < self.n2 && self.core_2[p.n2] != NULL {
                p.n2 += 1;
            }
        }
        p.n1 < self.n1 && p.n2 < self.n2
    }

    fn feasible(&self, n1: usize, n2: usize) -> bool {
        if self.qadj[n1].len() > self.t.v.degree(n2) {
            return false;
        }
        if !self.t.atom_matches(&self.q.atoms[n1].query, n2) {
            return false;
        }
        for &(qb, other1) in &self.qadj[n1] {
            if self.core_1[other1] != NULL {
                let other2 = self.core_1[other1];
                match self.t.v.bond_between(n2, other2) {
                    Some(tb) if self.t.bond_matches(&self.q.bonds[qb].query, tb) => {}
                    _ => return false,
                }
            }
        }
        true
    }

    fn add_pair(&mut self, n1: usize, n2: usize) {
        self.core_len += 1;
        if self.term_1[n1] == 0 {
            self.term_1[n1] = self.core_len;
            self.t1_len += 1;
        }
        if self.term_2[n2] == 0 {
            self.term_2[n2] = self.core_len;
            self.t2_len += 1;
        }
        self.core_1[n1] = n2;
        self.core_2[n2] = n1;
        for i in 0..self.qadj[n1].len() {
            let o = self.qadj[n1][i].1;
            if self.term_1[o] == 0 {
                self.term_1[o] = self.core_len;
                self.t1_len += 1;
            }
        }
        for &b in &self.t.v.atom_bonds[n2] {
            let o = self.t.v.other_atom(b, n2);
            if self.term_2[o] == 0 {
                self.term_2[o] = self.core_len;
                self.t2_len += 1;
            }
        }
    }

    fn back_track(&mut self, n1: usize, n2: usize) {
        if self.term_1[n1] == self.core_len {
            self.term_1[n1] = 0;
            self.t1_len -= 1;
        }
        for i in 0..self.qadj[n1].len() {
            let o = self.qadj[n1][i].1;
            if self.term_1[o] == self.core_len {
                self.term_1[o] = 0;
                self.t1_len -= 1;
            }
        }
        if self.term_2[n2] == self.core_len {
            self.term_2[n2] = 0;
            self.t2_len -= 1;
        }
        for &b in &self.t.v.atom_bonds[n2] {
            let o = self.t.v.other_atom(b, n2);
            if self.term_2[o] == self.core_len {
                self.term_2[o] = 0;
                self.t2_len -= 1;
            }
        }
        self.core_1[n1] = NULL;
        self.core_2[n2] = NULL;
        self.core_len -= 1;
    }

    fn match_all(&mut self, res: &mut Vec<Vec<usize>>) -> bool {
        if self.core_len == self.n1 {
            res.push(self.core_1.clone());
            return res.len() >= MAX_MATCHES;
        }
        if self.n1 > self.n2 || self.t1_len > self.t2_len {
            return false;
        }
        let mut p = Pair {
            n1: NULL,
            n2: NULL,
            nbrs: None,
        };
        while self.next_pair(&mut p) {
            if self.feasible(p.n1, p.n2) {
                let (a, b) = (p.n1, p.n2);
                self.add_pair(a, b);
                if self.match_all(res) {
                    return true;
                }
                self.back_track(a, b);
            }
        }
        false
    }
}

/// All matches (query atom index -> target atom index), RDKit's order.
pub(crate) fn substruct_match(t: &Target, q: &QueryMolecule) -> Vec<Vec<usize>> {
    let n1 = q.atoms.len();
    let n2 = t.v.num_atoms();
    if n2 == 0 || n1 == 0 || n1 > n2 {
        return Vec::new();
    }
    let mut s = State {
        q,
        qadj: query_adj(q),
        t,
        n1,
        n2,
        core_len: 0,
        t1_len: 0,
        t2_len: 0,
        core_1: vec![NULL; n1],
        core_2: vec![NULL; n2],
        term_1: vec![0; n1],
        term_2: vec![0; n2],
    };
    let mut res = Vec::new();
    s.match_all(&mut res);
    res
}
