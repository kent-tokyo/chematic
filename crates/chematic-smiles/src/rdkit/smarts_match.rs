//! A SMARTS subset (the patterns RDKit's own code matches internally:
//! MolHash bond flags, the tautomer transforms and scoring terms) and
//! `SubstructMatch(mol, query)` with RDKit's default parameters
//! (`uniquify=true`, `maxMatches=1000`, no chirality) on the RDKit-model
//! molecule, in RDKit's match order.
//!
//! Supported (with RDKit's `smarts.yy`/`smarts.ll` semantics): `*`,
//! organic-subset atoms (`C`, `c`, `N`, ...), bracket atoms with `#n`,
//! element symbols (one or two letters, aromatic `c`, `n`, `se`, ...), `A`,
//! `a`, `Hn`, `Dn`, `Xn`, `vn`, `Rn`, `rn`, `zn`, `Zn`, `xn` (bare forms
//! included) and `{a-b}` ranges on them, charges (`+`, `-`, `+n`, `-0`, ...)
//! and recursive `$(...)`, combined with `!`, implicit `&`, `,` and `;`;
//! bonds `-`, `=`, `#`, `:`, `~`, `@` combined the same way, the default
//! bond (single or aromatic), branches and ring closures.

use super::mol::{Atom, Bond, BondType, Mol};
use super::periodic;
use super::substruct::substruct_matches;

/// An integer atom property RDKit's simple atom queries read.
#[derive(Clone, Copy, Debug)]
enum IntProp {
    /// `queryAtomHCount`: `getTotalNumHs(true)`.
    HCount,
    /// `queryAtomExplicitDegree`.
    Degree,
    /// `queryAtomTotalDegree`.
    TotalDegree,
    /// `queryAtomTotalValence`.
    TotalValence,
    /// `queryIsAtomInNRings`.
    NumRings,
    /// `queryAtomMinRingSize`.
    MinRingSize,
    /// `queryAtomNumHeteroatomNbrs`.
    HeteroNbrs,
    /// `queryAtomNumAliphaticHeteroatomNbrs`.
    AliphaticHeteroNbrs,
    /// `queryAtomRingBondCount`.
    RingBondCount,
    /// `queryAtomFormalCharge`.
    Charge,
    /// `queryAtomIsotope`.
    Isotope,
}

#[derive(Clone, Debug)]
enum AtomQuery {
    Any,
    /// `AtomType`: atomic number and aromaticity.
    Type(u32, bool),
    /// `AtomAtomicNum`.
    Num(u32),
    /// An integer property in an inclusive range.
    Int(IntProp, i32, i32),
    Aliphatic,
    Aromatic,
    /// `queryIsAtomInRing` (`R`, `r`).
    InRing,
    /// `queryAtomHasHeteroatomNbrs` (`z`).
    HasHeteroNbrs,
    /// `queryAtomHasAliphaticHeteroatomNbrs` (`Z`).
    HasAliphaticHeteroNbrs,
    /// `queryAtomHasRingBond` (`x`).
    HasRingBond,
    /// A recursive query (an index into [`Query::recursive`]).
    Recursive(usize),
    And(Vec<AtomQuery>),
    Or(Vec<AtomQuery>),
    Not(Box<AtomQuery>),
}

#[derive(Clone, Debug)]
enum BondQuery {
    Any,
    SingleOrAromatic,
    Order(BondType),
    /// `queryIsBondInRing`.
    InRing,
    And(Vec<BondQuery>),
    Or(Vec<BondQuery>),
    Not(Box<BondQuery>),
}

/// A parsed SMARTS query.
#[derive(Clone, Debug)]
pub(crate) struct Query {
    atoms: Vec<AtomQuery>,
    bonds: Vec<(usize, usize, BondQuery)>,
    recursive: Vec<Query>,
}

impl Query {
    pub(crate) fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// The query atoms bond `b` joins.
    pub(crate) fn bond_atoms(&self, b: usize) -> (usize, usize) {
        (self.bonds[b].0, self.bonds[b].1)
    }
}

/// Parses `smarts` (the subset described in the module documentation).
/// Panics on anything else: the patterns are fixed strings in the port.
pub(crate) fn parse(smarts: &str) -> Query {
    let mut p = Parser {
        s: smarts.as_bytes(),
        pos: 0,
        recursive: Vec::new(),
    };
    let (atoms, bonds) = p.chain();
    assert!(p.pos == p.s.len(), "SMARTS {smarts}: trailing input");
    Query {
        atoms,
        bonds,
        recursive: p.recursive,
    }
}

struct Parser<'a> {
    s: &'a [u8],
    pos: usize,
    recursive: Vec<Query>,
}

type Chain = (Vec<AtomQuery>, Vec<(usize, usize, BondQuery)>);

fn is_bond_char(c: u8) -> bool {
    matches!(c, b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!')
}

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.pos).copied()
    }

    fn peek_at(&self, k: usize) -> Option<u8> {
        self.s.get(self.pos + k).copied()
    }

    fn number(&mut self) -> Option<u32> {
        let start = self.pos;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.pos += 1;
        }
        (self.pos > start).then(|| {
            std::str::from_utf8(&self.s[start..self.pos])
                .expect("ascii")
                .parse()
                .expect("number")
        })
    }

    /// Atoms and bonds up to the end of input or an unmatched `)`. Ring
    /// closure bonds are added after the chain's other bonds, in closing
    /// order (`CloseMolRings`).
    fn chain(&mut self) -> Chain {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut closures = Vec::new();
        let mut open: Vec<(u32, usize, Option<BondQuery>)> = Vec::new();
        let mut prev: Option<usize> = None;
        let mut branch_points: Vec<Option<usize>> = Vec::new();
        let mut pending: Option<BondQuery> = None;
        while let Some(c) = self.peek() {
            match c {
                b'(' => {
                    self.pos += 1;
                    branch_points.push(prev);
                }
                b')' => {
                    if branch_points.is_empty() {
                        break;
                    }
                    self.pos += 1;
                    prev = branch_points.pop().expect("branch");
                }
                c if is_bond_char(c) => pending = Some(self.bond_low_and()),
                c if c.is_ascii_digit() || c == b'%' => {
                    self.pos += 1;
                    let num = if c == b'%' {
                        self.number().expect("ring closure number")
                    } else {
                        u32::from(c - b'0')
                    };
                    let here = prev.expect("ring closure after an atom");
                    let bq = pending.take();
                    if let Some(k) = open.iter().position(|o| o.0 == num) {
                        let (_, other, obq) = open.remove(k);
                        let q = bq.or(obq).unwrap_or(BondQuery::SingleOrAromatic);
                        closures.push((other, here, q));
                    } else {
                        open.push((num, here, bq));
                    }
                }
                _ => {
                    let q = if c == b'[' {
                        self.pos += 1;
                        let q = self.low_and();
                        assert_eq!(self.peek(), Some(b']'), "unterminated bracket atom");
                        self.pos += 1;
                        q
                    } else {
                        self.organic_atom()
                    };
                    let idx = atoms.len();
                    atoms.push(q);
                    if let Some(p) = prev {
                        bonds.push((
                            p,
                            idx,
                            pending.take().unwrap_or(BondQuery::SingleOrAromatic),
                        ));
                    }
                    prev = Some(idx);
                }
            }
        }
        assert!(open.is_empty(), "unclosed ring");
        bonds.extend(closures);
        (atoms, bonds)
    }

    fn bond_low_and(&mut self) -> BondQuery {
        let mut terms = vec![self.bond_or()];
        while self.peek() == Some(b';') {
            self.pos += 1;
            terms.push(self.bond_or());
        }
        combine(terms, BondQuery::And)
    }

    fn bond_or(&mut self) -> BondQuery {
        let mut terms = vec![self.bond_high_and()];
        while self.peek() == Some(b',') {
            self.pos += 1;
            terms.push(self.bond_high_and());
        }
        combine(terms, BondQuery::Or)
    }

    fn bond_high_and(&mut self) -> BondQuery {
        let mut terms = Vec::new();
        while let Some(c) = self.peek() {
            if c == b'&' {
                self.pos += 1;
                continue;
            }
            if !is_bond_char(c) {
                break;
            }
            terms.push(self.bond_unary());
        }
        assert!(!terms.is_empty(), "empty bond expression");
        combine(terms, BondQuery::And)
    }

    fn bond_unary(&mut self) -> BondQuery {
        let c = self.peek().expect("bond");
        self.pos += 1;
        match c {
            b'!' => BondQuery::Not(Box::new(self.bond_unary())),
            b'-' => BondQuery::Order(BondType::Single),
            b'=' => BondQuery::Order(BondType::Double),
            b'#' => BondQuery::Order(BondType::Triple),
            b':' => BondQuery::Order(BondType::Aromatic),
            b'~' => BondQuery::Any,
            b'@' => BondQuery::InRing,
            other => panic!("unsupported SMARTS bond {:?}", other as char),
        }
    }

    fn organic_atom(&mut self) -> AtomQuery {
        let c = self.peek().expect("atom");
        self.pos += 1;
        if c == b'*' {
            return AtomQuery::Any;
        }
        if c == b'a' {
            return AtomQuery::Aromatic;
        }
        if c == b'A' {
            return AtomQuery::Aliphatic;
        }
        // `Cl`, `Br`
        if (c == b'C' && self.peek() == Some(b'l')) || (c == b'B' && self.peek() == Some(b'r')) {
            let sym = [c, self.peek().expect("letter")];
            self.pos += 1;
            let anum = periodic::atomic_number(std::str::from_utf8(&sym).expect("ascii"))
                .expect("element");
            return AtomQuery::Type(anum, false);
        }
        element(c)
    }

    /// `;`-separated terms.
    fn low_and(&mut self) -> AtomQuery {
        let mut terms = vec![self.or()];
        while self.peek() == Some(b';') {
            self.pos += 1;
            terms.push(self.or());
        }
        combine(terms, AtomQuery::And)
    }

    /// `,`-separated terms.
    fn or(&mut self) -> AtomQuery {
        let mut terms = vec![self.high_and()];
        while self.peek() == Some(b',') {
            self.pos += 1;
            terms.push(self.high_and());
        }
        combine(terms, AtomQuery::Or)
    }

    /// Adjacent (or `&`-joined) terms.
    fn high_and(&mut self) -> AtomQuery {
        let mut terms = Vec::new();
        while let Some(c) = self.peek() {
            if matches!(c, b']' | b',' | b';' | b')') {
                break;
            }
            if c == b'&' {
                self.pos += 1;
                continue;
            }
            terms.push(self.unary());
        }
        combine(terms, AtomQuery::And)
    }

    /// An optional `{a-b}`, `{-b}` or `{a-}` range after a primitive.
    fn range(&mut self) -> Option<(i32, i32)> {
        if self.peek() != Some(b'{') {
            return None;
        }
        self.pos += 1;
        let lo = self.number().map_or(i32::MIN, |n| n as i32);
        assert_eq!(self.peek(), Some(b'-'), "range");
        self.pos += 1;
        let hi = self.number().map_or(i32::MAX, |n| n as i32);
        assert_eq!(self.peek(), Some(b'}'), "range");
        self.pos += 1;
        Some((lo, hi))
    }

    /// A counted primitive: `{a-b}`, a number, or the bare form.
    fn counted(&mut self, prop: IntProp, bare: AtomQuery, range_bare: Option<i32>) -> AtomQuery {
        if self.peek() == Some(b'{') && range_bare.is_some() {
            let (lo, hi) = self.range().expect("range");
            return AtomQuery::Int(prop, lo, hi);
        }
        match self.number() {
            Some(n) => AtomQuery::Int(prop, n as i32, n as i32),
            None => bare,
        }
    }

    fn charge(&mut self, sign: i32) -> AtomQuery {
        // `+`, `++`, `+n`
        let sym = if sign > 0 { b'+' } else { b'-' };
        if let Some(n) = self.number() {
            let v = sign * n as i32;
            return AtomQuery::Int(IntProp::Charge, v, v);
        }
        let mut v = sign;
        while self.peek() == Some(sym) {
            self.pos += 1;
            v += sign;
        }
        AtomQuery::Int(IntProp::Charge, v, v)
    }

    fn unary(&mut self) -> AtomQuery {
        if self.peek() == Some(b'!') {
            self.pos += 1;
            return AtomQuery::Not(Box::new(self.unary()));
        }
        let c = self.peek().expect("primitive");
        // Isotope: a leading number.
        if c.is_ascii_digit() {
            let n = self.number().expect("isotope") as i32;
            return AtomQuery::Int(IntProp::Isotope, n, n);
        }
        // Two-letter element symbols (longest match, as the lexer does).
        if c.is_ascii_uppercase()
            && let Some(c2) = self.peek_at(1)
            && c2.is_ascii_lowercase()
        {
            let sym = [c, c2];
            let sym = std::str::from_utf8(&sym).expect("ascii");
            if let Some(anum) = periodic::atomic_number(sym) {
                self.pos += 2;
                return if matches!(sym, "Cl" | "Br") {
                    AtomQuery::Type(anum, false)
                } else {
                    AtomQuery::Num(anum)
                };
            }
        }
        // Aromatic two-letter symbols.
        if let Some(c2) = self.peek_at(1) {
            let anum = match (c, c2) {
                (b's', b'e') => Some(34),
                (b't', b'e') => Some(52),
                (b'a', b's') => Some(33),
                (b's', b'i') => Some(14),
                _ => None,
            };
            if let Some(anum) = anum {
                self.pos += 2;
                return AtomQuery::Type(anum, true);
            }
        }
        self.pos += 1;
        match c {
            b'*' => AtomQuery::Any,
            b'#' => AtomQuery::Num(self.number().expect("atomic number")),
            b'H' => {
                let n = self.number().unwrap_or(1) as i32;
                AtomQuery::Int(IntProp::HCount, n, n)
            }
            b'D' => self.counted(
                IntProp::Degree,
                AtomQuery::Int(IntProp::Degree, 1, 1),
                Some(1),
            ),
            b'X' => self.counted(
                IntProp::TotalDegree,
                AtomQuery::Int(IntProp::TotalDegree, 1, 1),
                Some(1),
            ),
            b'v' => self.counted(
                IntProp::TotalValence,
                AtomQuery::Int(IntProp::TotalValence, 1, 1),
                Some(1),
            ),
            b'R' => self.counted(IntProp::NumRings, AtomQuery::InRing, Some(-1)),
            b'r' => self.counted(IntProp::MinRingSize, AtomQuery::InRing, Some(5)),
            b'z' => self.counted(IntProp::HeteroNbrs, AtomQuery::HasHeteroNbrs, Some(0)),
            b'Z' => self.counted(
                IntProp::AliphaticHeteroNbrs,
                AtomQuery::HasAliphaticHeteroNbrs,
                Some(0),
            ),
            b'x' => self.counted(IntProp::RingBondCount, AtomQuery::HasRingBond, Some(0)),
            b'A' => AtomQuery::Aliphatic,
            b'a' => AtomQuery::Aromatic,
            b'+' => self.charge(1),
            b'-' => self.charge(-1),
            b'$' => {
                assert_eq!(self.peek(), Some(b'('), "recursive SMARTS");
                self.pos += 1;
                let saved = std::mem::take(&mut self.recursive);
                let (atoms, bonds) = self.chain();
                let inner = std::mem::replace(&mut self.recursive, saved);
                assert_eq!(self.peek(), Some(b')'), "unterminated recursive SMARTS");
                self.pos += 1;
                self.recursive.push(Query {
                    atoms,
                    bonds,
                    recursive: inner,
                });
                AtomQuery::Recursive(self.recursive.len() - 1)
            }
            _ => element(c),
        }
    }
}

fn combine<T>(mut terms: Vec<T>, f: fn(Vec<T>) -> T) -> T {
    if terms.len() == 1 {
        terms.pop().expect("one")
    } else {
        f(terms)
    }
}

/// A one-letter element symbol (upper case aliphatic, lower case aromatic).
/// Upper case symbols outside the organic subset match any aromaticity.
fn element(c: u8) -> AtomQuery {
    let sym = (c.to_ascii_uppercase() as char).to_string();
    let anum = periodic::atomic_number(&sym).unwrap_or_else(|| panic!("element {sym}"));
    if c.is_ascii_uppercase() && !matches!(c, b'B' | b'C' | b'N' | b'O' | b'F' | b'P' | b'S' | b'I')
    {
        return AtomQuery::Num(anum);
    }
    AtomQuery::Type(anum, c.is_ascii_lowercase())
}

fn is_hetero(mol: &Mol, a: usize) -> bool {
    let n = mol.atoms[a].anum;
    n != 6 && n != 1
}

fn int_prop(prop: IntProp, mol: &Mol, a: usize) -> i32 {
    let atom = &mol.atoms[a];
    match prop {
        IntProp::HCount => {
            let nbr_hs = mol.nbrs(a).filter(|&x| mol.atoms[x].anum == 1).count() as u32;
            (mol.total_num_hs(a) + nbr_hs) as i32
        }
        IntProp::Degree => mol.degree(a) as i32,
        IntProp::TotalDegree => mol.total_degree(a) as i32,
        IntProp::TotalValence => mol.total_valence(a),
        IntProp::NumRings => mol.num_atom_rings(a) as i32,
        IntProp::MinRingSize => mol.rings.as_ref().map_or(0, |r| {
            r.atom_members(a)
                .iter()
                .map(|&ri| r.atom_rings[ri].len())
                .min()
                .unwrap_or(0) as i32
        }),
        IntProp::HeteroNbrs => mol.nbrs(a).filter(|&x| is_hetero(mol, x)).count() as i32,
        IntProp::AliphaticHeteroNbrs => mol
            .nbrs(a)
            .filter(|&x| is_hetero(mol, x) && !mol.atoms[x].aromatic)
            .count() as i32,
        IntProp::RingBondCount => mol.atom_bonds[a]
            .iter()
            .filter(|&&b| mol.num_bond_rings(b) != 0)
            .count() as i32,
        IntProp::Charge => atom.charge,
        IntProp::Isotope => atom.isotope as i32,
    }
}

fn atom_matches(q: &AtomQuery, mol: &Mol, a: usize, rec: &[Vec<bool>]) -> bool {
    let atom = &mol.atoms[a];
    match q {
        AtomQuery::Any => true,
        AtomQuery::Type(n, arom) => atom.anum == *n && atom.aromatic == *arom,
        AtomQuery::Num(n) => atom.anum == *n,
        AtomQuery::Int(p, lo, hi) => {
            let v = int_prop(*p, mol, a);
            *lo <= v && v <= *hi
        }
        AtomQuery::Aliphatic => !atom.aromatic,
        AtomQuery::Aromatic => atom.aromatic,
        AtomQuery::InRing => mol.num_atom_rings(a) != 0,
        AtomQuery::HasHeteroNbrs => mol.nbrs(a).any(|x| is_hetero(mol, x)),
        AtomQuery::HasAliphaticHeteroNbrs => mol
            .nbrs(a)
            .any(|x| is_hetero(mol, x) && !mol.atoms[x].aromatic),
        AtomQuery::HasRingBond => mol.atom_bonds[a]
            .iter()
            .any(|&b| mol.num_bond_rings(b) != 0),
        AtomQuery::Recursive(i) => rec[*i][a],
        AtomQuery::And(v) => v.iter().all(|x| atom_matches(x, mol, a, rec)),
        AtomQuery::Or(v) => v.iter().any(|x| atom_matches(x, mol, a, rec)),
        AtomQuery::Not(x) => !atom_matches(x, mol, a, rec),
    }
}

fn bond_matches(q: &BondQuery, mol: &Mol, b: usize) -> bool {
    let bt = mol.bonds[b].bt;
    match q {
        BondQuery::Any => true,
        BondQuery::SingleOrAromatic => matches!(bt, BondType::Single | BondType::Aromatic),
        BondQuery::Order(o) => bt == *o,
        BondQuery::InRing => mol.num_bond_rings(b) != 0,
        BondQuery::And(v) => v.iter().all(|x| bond_matches(x, mol, b)),
        BondQuery::Or(v) => v.iter().any(|x| bond_matches(x, mol, b)),
        BondQuery::Not(x) => !bond_matches(x, mol, b),
    }
}

/// Every match (query atom order) in RDKit's enumeration order, without
/// uniquification.
fn all_matches(q: &Query, mol: &Mol, max_matches: usize) -> Vec<Vec<usize>> {
    let rec: Vec<Vec<bool>> = q
        .recursive
        .iter()
        .map(|r| {
            let mut hit = vec![false; mol.atoms.len()];
            for m in all_matches(r, mol, 0) {
                hit[m[0]] = true;
            }
            hit
        })
        .collect();
    let mut topo = Mol::default();
    for _ in &q.atoms {
        topo.add_atom(Atom::new(0));
    }
    for &(i, j, _) in &q.bonds {
        topo.add_bond(Bond::new(i, j, BondType::Single));
    }
    let vc = |qa: usize, ma: usize| atom_matches(&q.atoms[qa], mol, ma, &rec);
    let ec = |qb: usize, mb: usize| bond_matches(&q.bonds[qb].2, mol, mb);
    substruct_matches(mol, &topo, &vc, &ec, max_matches)
        .into_iter()
        .map(|m| m.into_iter().map(|(_, x)| x).collect())
        .collect()
}

/// `SubstructMatch(mol, query)` with the default parameters: matches as
/// mol atom indices in query atom order, the first of each atom set only,
/// at most 1000.
pub(crate) fn substruct_match(q: &Query, mol: &Mol) -> Vec<Vec<usize>> {
    let mut seen = std::collections::HashSet::new();
    let mut res = Vec::new();
    for m in all_matches(q, mol, 0) {
        let mut key = m.clone();
        key.sort_unstable();
        if seen.insert(key) {
            res.push(m);
            if res.len() == 1000 {
                break;
            }
        }
    }
    res
}
