//! A SMARTS subset (the patterns RDKit's own code matches internally) and
//! `SubstructMatch(mol, query)` with RDKit's default parameters
//! (`uniquify=true`, `maxMatches=1000`, no chirality) on the RDKit-model
//! molecule, in RDKit's match order.
//!
//! Supported: `*`, organic-subset atoms (`C`, `c`, `N`, ...), bracket atoms
//! with `#n`, element symbols, `A`, `a`, `Hn`, `Dn`, `+`/`-` charges and
//! recursive `$(...)`, combined with `!`, implicit `&`, `,` and `;`; bonds
//! `-`, `=`, `#`, `:` and `,` lists of them, the default bond (single or
//! aromatic) and branches. Ring closures are not supported.

use super::mol::{Atom, Bond, BondType, Mol};
use super::periodic;
use super::substruct::substruct_matches;

#[derive(Clone, Debug)]
enum AtomQuery {
    Any,
    /// `AtomType`: atomic number and aromaticity.
    Type(u32, bool),
    /// `AtomAtomicNum`.
    Num(u32),
    /// `AtomHCount`: `getTotalNumHs(true)`.
    HCount(u32),
    Charge(i32),
    /// `AtomExplicitDegree`.
    Degree(usize),
    Aliphatic,
    Aromatic,
    /// A recursive query (an index into [`Query::recursive`]).
    Recursive(usize),
    And(Vec<AtomQuery>),
    Or(Vec<AtomQuery>),
    Not(Box<AtomQuery>),
}

#[derive(Clone, Debug)]
enum BondQuery {
    SingleOrAromatic,
    Order(BondType),
    Or(Vec<BondQuery>),
}

/// A parsed SMARTS query.
#[derive(Clone, Debug)]
pub(crate) struct Query {
    atoms: Vec<AtomQuery>,
    bonds: Vec<(usize, usize, BondQuery)>,
    recursive: Vec<Query>,
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

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.pos).copied()
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

    /// Atoms and bonds up to the end of input or an unmatched `)`.
    fn chain(&mut self) -> Chain {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
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
                b'-' | b'=' | b'#' | b':' => pending = Some(self.bond()),
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
        (atoms, bonds)
    }

    fn bond(&mut self) -> BondQuery {
        let mut alts = Vec::new();
        loop {
            let bt = match self.peek() {
                Some(b'-') => BondType::Single,
                Some(b'=') => BondType::Double,
                Some(b'#') => BondType::Triple,
                Some(b':') => BondType::Aromatic,
                other => panic!("unsupported SMARTS bond {other:?}"),
            };
            self.pos += 1;
            alts.push(BondQuery::Order(bt));
            if self.peek() == Some(b',') {
                self.pos += 1;
            } else {
                break;
            }
        }
        if alts.len() == 1 {
            alts.pop().expect("one")
        } else {
            BondQuery::Or(alts)
        }
    }

    fn organic_atom(&mut self) -> AtomQuery {
        let c = self.peek().expect("atom");
        self.pos += 1;
        if c == b'*' {
            return AtomQuery::Any;
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

    fn unary(&mut self) -> AtomQuery {
        if self.peek() == Some(b'!') {
            self.pos += 1;
            return AtomQuery::Not(Box::new(self.unary()));
        }
        let c = self.peek().expect("primitive");
        self.pos += 1;
        match c {
            b'*' => AtomQuery::Any,
            b'#' => AtomQuery::Num(self.number().expect("atomic number")),
            b'H' => AtomQuery::HCount(self.number().unwrap_or(1)),
            b'D' => AtomQuery::Degree(self.number().unwrap_or(1) as usize),
            b'A' => AtomQuery::Aliphatic,
            b'a' => AtomQuery::Aromatic,
            b'+' => AtomQuery::Charge(self.number().map_or(1, |n| n as i32)),
            b'-' => AtomQuery::Charge(-self.number().map_or(1, |n| n as i32)),
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

fn combine(mut terms: Vec<AtomQuery>, f: fn(Vec<AtomQuery>) -> AtomQuery) -> AtomQuery {
    if terms.len() == 1 {
        terms.pop().expect("one")
    } else {
        f(terms)
    }
}

/// A one-letter element symbol (upper case aliphatic, lower case aromatic).
fn element(c: u8) -> AtomQuery {
    let sym = (c.to_ascii_uppercase() as char).to_string();
    let anum = periodic::atomic_number(&sym).unwrap_or_else(|| panic!("element {sym}"));
    AtomQuery::Type(anum, c.is_ascii_lowercase())
}

fn atom_matches(q: &AtomQuery, mol: &Mol, a: usize, rec: &[Vec<bool>]) -> bool {
    let atom = &mol.atoms[a];
    match q {
        AtomQuery::Any => true,
        AtomQuery::Type(n, arom) => atom.anum == *n && atom.aromatic == *arom,
        AtomQuery::Num(n) => atom.anum == *n,
        AtomQuery::HCount(n) => {
            let nbr_hs = mol.nbrs(a).filter(|&x| mol.atoms[x].anum == 1).count() as u32;
            mol.total_num_hs(a) + nbr_hs == *n
        }
        AtomQuery::Charge(c) => atom.charge == *c,
        AtomQuery::Degree(d) => mol.degree(a) == *d,
        AtomQuery::Aliphatic => !atom.aromatic,
        AtomQuery::Aromatic => atom.aromatic,
        AtomQuery::Recursive(i) => rec[*i][a],
        AtomQuery::And(v) => v.iter().all(|x| atom_matches(x, mol, a, rec)),
        AtomQuery::Or(v) => v.iter().any(|x| atom_matches(x, mol, a, rec)),
        AtomQuery::Not(x) => !atom_matches(x, mol, a, rec),
    }
}

fn bond_matches(q: &BondQuery, bt: BondType) -> bool {
    match q {
        BondQuery::SingleOrAromatic => matches!(bt, BondType::Single | BondType::Aromatic),
        BondQuery::Order(o) => bt == *o,
        BondQuery::Or(v) => v.iter().any(|x| bond_matches(x, bt)),
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
    let ec = |qb: usize, mb: usize| bond_matches(&q.bonds[qb].2, mol.bonds[mb].bt);
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
