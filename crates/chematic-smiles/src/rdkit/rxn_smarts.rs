//! RDKit 2026.03.1's SMARTS reader (`smarts.ll`/`smarts.yy`: query trees
//! and the atom properties the grammar sets) and `SmartsWrite` for query
//! molecules, as `Chem.MolToSmarts(Chem.MolFromSmarts(s))` and
//! `rdChemReactions.ReactionToSmarts(ReactionFromSmarts(s))` run them.
//!
//! The query tree of every atom and bond is built as RDKit's grammar builds
//! it (left-associative `;` < `,` < `&`/implicit AND, `expandQuery` with its
//! null-query merging, `!` on point queries), then written by the port of
//! `_recurseGetSmarts`/`_combineChildSmarts`/`getAtomSmartsSimple`. The
//! output order and the written chirality come from the SMILES port's
//! `canonicalizeFragment` on RDKit's molecule of the template (atom indices
//! as ranks, empty ring information, start at the first non-chiral atom).
//! Directional (`/`, `\`) and dative bonds, chirality classes and ranges
//! of `k` are refused.

use super::RdkitSmilesError;
use super::canon::StackElem;
use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Mol};
use super::periodic;

fn unsupported(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(what.into())
}

/// A leaf query (`getDescription()` and value).
#[derive(Clone, Debug, PartialEq)]
enum Leaf {
    /// `AtomType`: atomic number and aromaticity.
    AtomType(u32, bool),
    /// `AtomAtomicNum`.
    AtomicNum(u32),
    /// `AtomNull`.
    Null,
    /// A simple query written as a letter with an optional value (`H`, `D`,
    /// `X`, `v`, `d`, `h`, `x`, `z`, `Z`, `r`, `k`, `R`), with RDKit's
    /// range modifiers.
    Valued(&'static str, Option<i32>, Range),
    /// `AtomIsAromatic` / `AtomIsAliphatic`.
    Aromatic(bool),
    /// `AtomFormalCharge`.
    Charge(i32),
    /// `AtomIsotope`.
    Isotope(i32),
    /// `AtomHybridization`.
    Hybridization(i32),
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Range {
    None,
    /// `{n-}`
    Less(i32),
    /// `{-n}`
    Greater(i32),
    /// `{a-b}`
    Between(i32, i32),
}

#[derive(Clone, Debug)]
enum Node {
    Leaf(Leaf, bool),
    Recursive(Box<QMol>, bool),
    And(Vec<Node>, bool),
    Or(Vec<Node>, bool),
}

impl Node {
    fn is_null(&self) -> bool {
        matches!(self, Node::Leaf(Leaf::Null, _))
    }
    fn negation(&self) -> bool {
        match self {
            Node::Leaf(_, n) | Node::Recursive(_, n) | Node::And(_, n) | Node::Or(_, n) => *n,
        }
    }
    fn set_negation(&mut self, v: bool) {
        match self {
            Node::Leaf(_, n) | Node::Recursive(_, n) | Node::And(_, n) | Node::Or(_, n) => *n = v,
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum How {
    And,
    Or,
}

/// `QueryAtom::expandQuery(what, how, maintainOrder=true)` with
/// `mergeNullQueries`.
fn expand(q: Node, what: Node, how: How) -> Node {
    let (qn, wn) = (q.is_null(), what.is_null());
    if qn && wn {
        let mut q = q;
        match how {
            How::And if !q.negation() && what.negation() => q.set_negation(true),
            How::Or if q.negation() && !what.negation() => q.set_negation(false),
            _ => {}
        }
        return q;
    }
    if qn || wn {
        // `mergeNullQFirst(null, other)`: AND keeps the other query unless
        // the null query is negated; OR keeps the null unless negated.
        let (null, other) = if qn { (q, what) } else { (what, q) };
        let keep_other = match how {
            How::And => !null.negation(),
            How::Or => null.negation(),
        };
        return if keep_other { other } else { null };
    }
    match how {
        How::And => Node::And(vec![q, what], false),
        How::Or => Node::Or(vec![q, what], false),
    }
}

/// `Bond::DATIVEL` / `Bond::DATIVER` (the bond types of `<-` / `->`).
const DATIVEL: i32 = 18;
const DATIVER: i32 = 19;

/// The bond type a bond expression leaves on its `QueryBond`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum PBt {
    /// The bond type and the direction a `/` `\\` token set.
    B(BondType, BondDir),
    DativeR,
    DativeL,
}

type OpenRing = (usize, Option<(Node, PBt)>, usize);

/// `QueryBond::expandQuery(what, how, maintainOrder=true)`: a `BondNull`
/// (`~`) side is merged as `mergeNullQueries` merges atom null queries.
fn bond_expand(q: Node, what: Node, how: How) -> Node {
    let is_null = |n: &Node| matches!(n, Node::Leaf(Leaf::Valued("BondNull", ..), _));
    let (qn, wn) = (is_null(&q), is_null(&what));
    if qn && wn {
        let mut q = q;
        match how {
            How::And if !q.negation() && what.negation() => q.set_negation(true),
            How::Or if q.negation() && !what.negation() => q.set_negation(false),
            _ => {}
        }
        return q;
    }
    if qn || wn {
        let (null, other) = if qn { (q, what) } else { (what, q) };
        let keep_other = match how {
            How::And => !null.negation(),
            How::Or => null.negation(),
        };
        return if keep_other { other } else { null };
    }
    match how {
        How::And => Node::And(vec![q, what], false),
        How::Or => Node::Or(vec![q, what], false),
    }
}

/// `SMARTS_H_MASK` / `SMARTS_CHARGE_MASK`.
const H_MASK: u8 = 1;
const CHARGE_MASK: u8 = 2;

/// A query atom as the grammar leaves it.
#[derive(Clone, Debug)]
struct QAtom {
    query: Node,
    anum: u32,
    aromatic: bool,
    isotope: u32,
    charge: i32,
    num_explicit_hs: u32,
    no_implicit: bool,
    chiral: ChiralTag,
    flags: u8,
    map: Option<u32>,
}

impl QAtom {
    fn new(anum: u32, query: Node) -> QAtom {
        QAtom {
            query,
            anum,
            aromatic: false,
            isotope: 0,
            charge: 0,
            num_explicit_hs: 0,
            no_implicit: false,
            chiral: ChiralTag::Unspecified,
            flags: 0,
            map: None,
        }
    }
    /// `QueryAtom()`: atomic number 0.
    fn with_query(query: Node) -> QAtom {
        QAtom::new(0, query)
    }
    /// `ClearAtomChemicalProps`.
    fn clear_chemical_props(&mut self) {
        self.isotope = 0;
        self.charge = 0;
        self.num_explicit_hs = 0;
    }
}

#[derive(Clone, Debug)]
struct QBond {
    begin: usize,
    end: usize,
    query: Node,
    bt: BondType,
    /// The direction of a `/` `\\` bond (relative to `begin`).
    dir: BondDir,
    /// `_unspecifiedOrder`: an implicit bond.
    unspecified: bool,
}

/// One element of the parsed graph, for the template's SMILES stand-in.
#[derive(Clone, Debug)]
enum Tok {
    Atom(usize),
    Bond(usize),
    Open,
    Close,
    /// Ring-closure digit, its bond, the `/` `\\` `=` written at it.
    Ring(u32, Option<usize>, Option<char>),
    Dot,
}

#[derive(Clone, Debug, Default)]
struct QMol {
    atoms: Vec<QAtom>,
    bonds: Vec<QBond>,
    toks: Vec<Tok>,
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

const TWO_LETTER: &[&str] = &[
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Ar", "Ca", "Sc", "Ti", "Cr", "Mn", "Co", "Fe",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr", "Rb", "Sr", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Fl", "Lv",
];

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }
    fn peek_at(&self, k: usize) -> Option<u8> {
        self.s.get(self.i + k).copied()
    }
    fn err(&self, what: &str) -> RdkitSmilesError {
        unsupported(format!("SMARTS parse error at {}: {what}", self.i))
    }
    fn number(&mut self) -> Option<i32> {
        let start = self.i;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.i += 1;
        }
        if self.i == start {
            return None;
        }
        std::str::from_utf8(&self.s[start..self.i])
            .ok()?
            .parse()
            .ok()
    }

    /// `mol` (one template or a recursive query's molecule).
    fn mol(&mut self, recursion: bool) -> Result<QMol, RdkitSmilesError> {
        let mut m = QMol::default();
        let mut prev: Option<usize> = None;
        let mut branches: Vec<Option<usize>> = Vec::new();
        let mut pending: Option<(Node, PBt)> = None;
        let mut rings: std::collections::HashMap<u32, OpenRing> = std::collections::HashMap::new();
        while let Some(c) = self.peek() {
            match c {
                b')' | b'.' if pending.is_some() => {
                    return Err(self.err("bond without a following atom"));
                }
                b')' if recursion && branches.is_empty() => break,
                b'(' => {
                    self.i += 1;
                    branches.push(prev);
                    m.toks.push(Tok::Open);
                }
                b')' => {
                    self.i += 1;
                    prev = branches
                        .pop()
                        .ok_or_else(|| self.err("extra close parentheses"))?;
                    m.toks.push(Tok::Close);
                }
                b'.' => {
                    self.i += 1;
                    prev = None;
                    m.toks.push(Tok::Dot);
                }
                b'0'..=b'9' | b'%' => {
                    let num = if c == b'%' {
                        self.i += 1;
                        if self.peek() == Some(b'(') {
                            self.i += 1;
                            let n = self.number().ok_or_else(|| self.err("ring number"))?;
                            if self.peek() != Some(b')') {
                                return Err(self.err("ring number"));
                            }
                            self.i += 1;
                            n as u32
                        } else {
                            let a = self
                                .peek()
                                .filter(u8::is_ascii_digit)
                                .ok_or_else(|| self.err("ring"))?;
                            let b = self
                                .peek_at(1)
                                .filter(u8::is_ascii_digit)
                                .ok_or_else(|| self.err("ring"))?;
                            self.i += 2;
                            u32::from(a - b'0') * 10 + u32::from(b - b'0')
                        }
                    } else {
                        self.i += 1;
                        u32::from(c - b'0')
                    };
                    let atom = prev.ok_or_else(|| self.err("ring closure without atom"))?;
                    let bond_here = pending.take();
                    // The direction or double bond written at this digit,
                    // for the stand-in SMILES.
                    let sym = match &bond_here {
                        Some((_, PBt::B(_, BondDir::EndUpRight))) => Some('/'),
                        Some((_, PBt::B(_, BondDir::EndDownRight))) => Some('\\'),
                        Some((_, PBt::B(BondType::Double, _))) => Some('='),
                        _ => None,
                    };
                    if let Some((open_atom, open_bond, tok_pos)) = rings.remove(&num) {
                        let (query, bt, dir, unspecified) = match bond_here.or(open_bond) {
                            Some((_, PBt::DativeR | PBt::DativeL)) => {
                                return Err(unsupported("dative ring-closure bond in SMARTS"));
                            }
                            Some((q, PBt::B(bt, d))) => (q, bt, d, false),
                            None => (
                                Node::Leaf(Leaf::Null, false),
                                BondType::Single,
                                BondDir::None,
                                true,
                            ),
                        };
                        let b = m.bonds.len();
                        m.bonds.push(QBond {
                            begin: open_atom,
                            end: atom,
                            query,
                            bt,
                            dir,
                            unspecified,
                        });
                        if let Tok::Ring(_, slot, _) = &mut m.toks[tok_pos] {
                            *slot = Some(b);
                        }
                        m.toks.push(Tok::Ring(num, Some(b), sym));
                    } else {
                        rings.insert(num, (atom, bond_here, m.toks.len()));
                        m.toks.push(Tok::Ring(num, None, sym));
                    }
                }
                b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!' | b'&' | b',' | b';' | b'$'
                | b'/' | b'\\' | b'<' => {
                    pending = Some(self.bond_expr()?);
                }
                _ => {
                    let atom = self.atom()?;
                    let a = m.atoms.len();
                    m.atoms.push(atom);
                    if let Some(p) = prev {
                        // `setBondType(DATIVE)` replaces the query of a
                        // dative bond expression (`->`: begin at the atom
                        // to the left, `<-`: at the atom to the right).
                        let dative_order = || {
                            Node::Leaf(
                                Leaf::Valued(
                                    "BondOrder",
                                    Some(BondType::Dative as i32),
                                    Range::None,
                                ),
                                false,
                            )
                        };
                        let mut dir = BondDir::None;
                        let (query, bt, unspecified, (begin, end)) = match pending.take() {
                            Some((_, PBt::DativeR)) => {
                                (dative_order(), BondType::Dative, false, (p, a))
                            }
                            Some((_, PBt::DativeL)) => {
                                (dative_order(), BondType::Dative, false, (a, p))
                            }
                            Some((q, PBt::B(bt, d))) => {
                                dir = d;
                                (q, bt, false, (p, a))
                            }
                            None => (
                                Node::Leaf(Leaf::Null, false),
                                BondType::Single,
                                true,
                                (p, a),
                            ),
                        };
                        m.bonds.push(QBond {
                            begin,
                            end,
                            query,
                            bt,
                            dir,
                            unspecified,
                        });
                        m.toks.push(Tok::Bond(m.bonds.len() - 1));
                    } else if pending.is_some() {
                        return Err(self.err("bond without a preceding atom"));
                    }
                    m.toks.push(Tok::Atom(a));
                    prev = Some(a);
                }
            }
        }
        if pending.is_some() {
            return Err(self.err("bond without a following atom"));
        }
        if !rings.is_empty() || !branches.is_empty() {
            return Err(self.err("unclosed ring or branch"));
        }
        // `SetUnspecifiedBondTypes`: implicit bonds between aromatic atoms.
        for b in &mut m.bonds {
            if b.unspecified && m.atoms[b.begin].aromatic && m.atoms[b.end].aromatic {
                b.bt = BondType::Aromatic;
            }
        }
        Ok(m)
    }

    /// `bond_expr`, returning its query and bond type.
    fn bond_expr(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        // ';' < ',' < '&' < implicit AND (bond_query).
        self.bond_semi()
    }
    fn bond_semi(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        let (mut q, bt) = self.bond_or()?;
        while self.peek() == Some(b';') {
            self.i += 1;
            let (r, _) = self.bond_or()?;
            q = bond_expand(q, r, How::And);
        }
        Ok((q, bt))
    }
    fn bond_or(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        let (mut q, bt) = self.bond_and()?;
        while self.peek() == Some(b',') {
            self.i += 1;
            let (r, _) = self.bond_and()?;
            q = bond_expand(q, r, How::Or);
        }
        Ok((q, bt))
    }
    fn bond_and(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        let (mut q, bt) = self.bond_query()?;
        while self.peek() == Some(b'&') {
            self.i += 1;
            let (r, _) = self.bond_query()?;
            q = bond_expand(q, r, How::And);
        }
        Ok((q, bt))
    }
    fn bond_query(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        let (mut q, bt) = self.bondd()?;
        while matches!(
            self.peek(),
            Some(b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'!' | b'$' | b'/' | b'\\')
        ) {
            let (r, _) = self.bondd()?;
            q = bond_expand(q, r, How::And);
        }
        Ok((q, bt))
    }
    fn bondd(&mut self) -> Result<(Node, PBt), RdkitSmilesError> {
        let c = self.peek().ok_or_else(|| self.err("bond"))?;
        self.i += 1;
        let order = |bt: BondType| {
            Node::Leaf(
                Leaf::Valued("BondOrder", Some(bt as i32), Range::None),
                false,
            )
        };
        let dative = |v: i32| Node::Leaf(Leaf::Valued("BondOrder", Some(v), Range::None), false);
        Ok(match c {
            b'-' if self.peek() == Some(b'>') => {
                self.i += 1;
                (dative(DATIVER), PBt::DativeR)
            }
            b'<' if self.peek() == Some(b'-') => {
                self.i += 1;
                (dative(DATIVEL), PBt::DativeL)
            }
            b'-' => (
                order(BondType::Single),
                PBt::B(BondType::Single, BondDir::None),
            ),
            b'=' => (
                order(BondType::Double),
                PBt::B(BondType::Double, BondDir::None),
            ),
            b'#' => (
                order(BondType::Triple),
                PBt::B(BondType::Triple, BondDir::None),
            ),
            b'$' => (
                order(BondType::Quadruple),
                PBt::B(BondType::Quadruple, BondDir::None),
            ),
            b':' => (
                order(BondType::Aromatic),
                PBt::B(BondType::Aromatic, BondDir::None),
            ),
            b'~' => (
                Node::Leaf(Leaf::Valued("BondNull", None, Range::None), false),
                PBt::B(BondType::Single, BondDir::None),
            ),
            b'@' => (
                Node::Leaf(Leaf::Valued("BondInRing", None, Range::None), false),
                PBt::B(BondType::Single, BondDir::None),
            ),
            b'/' => (
                Node::Leaf(Leaf::Null, false),
                PBt::B(BondType::Single, BondDir::EndUpRight),
            ),
            b'\\' => (
                Node::Leaf(Leaf::Null, false),
                PBt::B(BondType::Single, BondDir::EndDownRight),
            ),
            b'!' => {
                let (mut q, bt) = self.bondd()?;
                let n = q.negation();
                q.set_negation(!n);
                (q, bt)
            }
            _ => return Err(self.err("bond")),
        })
    }

    /// `atomd`.
    fn atom(&mut self) -> Result<QAtom, RdkitSmilesError> {
        let c = self.peek().ok_or_else(|| self.err("atom"))?;
        if c != b'[' {
            return self.simple_atom(false).ok_or_else(|| self.err("atom"));
        }
        self.i += 1;
        // `hydrogen_atom`: [H], [H:n], [nH], [nH:n], [H+], [nH+:n] ...
        let save = self.i;
        let iso = self.number();
        if self.peek() == Some(b'H') {
            self.i += 1;
            let charge = self.charge_spec();
            let mut map = None;
            if self.peek() == Some(b':') {
                self.i += 1;
                map = self.number();
            }
            if self.peek() == Some(b']') {
                self.i += 1;
                let mut q = QAtom::new(1, Node::Leaf(Leaf::AtomicNum(1), false));
                if let Some(iso) = iso {
                    q.isotope = iso as u32;
                    q.query = expand(q.query, Node::Leaf(Leaf::Isotope(iso), false), How::And);
                }
                if let Some(ch) = charge {
                    q.charge = ch;
                    q.flags |= CHARGE_MASK;
                    q.query = expand(q.query, Node::Leaf(Leaf::Charge(ch), false), How::And);
                }
                q.map = map.map(|m| m as u32);
                return Ok(q);
            }
        }
        self.i = save;
        let mut q = self.atom_expr_semi()?;
        if self.peek() == Some(b':') {
            self.i += 1;
            q.map = Some(self.number().ok_or_else(|| self.err("atom map"))? as u32);
        }
        if self.peek() != Some(b']') {
            return Err(self.err("expected ']'"));
        }
        self.i += 1;
        Ok(q)
    }

    fn atom_expr_semi(&mut self) -> Result<QAtom, RdkitSmilesError> {
        let mut a = self.atom_expr_or()?;
        while self.peek() == Some(b';') {
            self.i += 1;
            let b = self.atom_expr_or()?;
            a = combine_exprs(a, b, How::And);
        }
        Ok(a)
    }
    fn atom_expr_or(&mut self) -> Result<QAtom, RdkitSmilesError> {
        let mut a = self.atom_expr_and()?;
        while self.peek() == Some(b',') {
            self.i += 1;
            let b = self.atom_expr_and()?;
            a = combine_exprs(a, b, How::Or);
            a.anum = 0;
        }
        Ok(a)
    }
    fn atom_expr_and(&mut self) -> Result<QAtom, RdkitSmilesError> {
        let mut a = self.point_query()?;
        loop {
            match self.peek() {
                Some(b'&') => self.i += 1,
                Some(b';' | b',' | b':' | b']') | None => break,
                _ => {}
            }
            let p = self.point_query()?;
            atom_expr_and_point_query(&mut a, p);
        }
        Ok(a)
    }

    fn point_query(&mut self) -> Result<QAtom, RdkitSmilesError> {
        match self.peek() {
            Some(b'!') => {
                self.i += 1;
                let mut p = self.point_query()?;
                let n = p.query.negation();
                p.query.set_negation(!n);
                p.anum = 0;
                p.clear_chemical_props();
                Ok(p)
            }
            Some(b'$') if self.peek_at(1) == Some(b'(') => {
                self.i += 2;
                let inner = self.mol(true)?;
                if self.peek() != Some(b')') {
                    return Err(self.err("unclosed recursion"));
                }
                self.i += 1;
                Ok(QAtom::with_query(Node::Recursive(Box::new(inner), false)))
            }
            _ => self.atom_query(),
        }
    }

    fn charge_spec(&mut self) -> Option<i32> {
        match self.peek() {
            Some(b'+') => {
                self.i += 1;
                if self.peek() == Some(b'+') {
                    self.i += 1;
                    Some(2)
                } else {
                    Some(self.number().unwrap_or(1))
                }
            }
            Some(b'-') => {
                self.i += 1;
                if self.peek() == Some(b'-') {
                    self.i += 1;
                    Some(-2)
                } else {
                    Some(-self.number().unwrap_or(1))
                }
            }
            _ => None,
        }
    }

    /// `simple_atom` (organic subset, aromatic symbols; `in_atom` adds
    /// `si`, `as`, `se`, `te`) and the `*`, `a`, `A` tokens.
    fn simple_atom(&mut self, in_atom: bool) -> Option<QAtom> {
        let c = self.peek()?;
        let two = self.peek_at(1);
        let organic = |anum: u32, aromatic: bool, len: usize, p: &mut Self| {
            p.i += len;
            let mut q = QAtom::new(anum, Node::Leaf(Leaf::AtomType(anum, aromatic), false));
            q.aromatic = aromatic;
            Some(q)
        };
        if in_atom {
            let pair = [c, two.unwrap_or(0)];
            for (sym, anum) in [(b"si", 14), (b"as", 33), (b"se", 34), (b"te", 52)] {
                if &pair == sym {
                    return organic(anum, true, 2, self);
                }
            }
        }
        match (c, two) {
            (b'C', Some(b'l')) => organic(17, false, 2, self),
            (b'B', Some(b'r')) => organic(35, false, 2, self),
            (b'B', _) => organic(5, false, 1, self),
            (b'C', _) => organic(6, false, 1, self),
            (b'N', _) => organic(7, false, 1, self),
            (b'O', _) => organic(8, false, 1, self),
            (b'F', _) => organic(9, false, 1, self),
            (b'P', _) => organic(15, false, 1, self),
            (b'S', _) => organic(16, false, 1, self),
            (b'I', _) => organic(53, false, 1, self),
            (b'b', _) => organic(5, true, 1, self),
            (b'c', _) => organic(6, true, 1, self),
            (b'n', _) => organic(7, true, 1, self),
            (b'o', _) => organic(8, true, 1, self),
            (b'p', _) => organic(15, true, 1, self),
            (b's', _) => organic(16, true, 1, self),
            (b'*', _) => {
                self.i += 1;
                Some(QAtom::with_query(Node::Leaf(Leaf::Null, false)))
            }
            (b'a', _) => {
                self.i += 1;
                let mut q = QAtom::with_query(Node::Leaf(Leaf::Aromatic(true), false));
                q.aromatic = true;
                Some(q)
            }
            (b'A', _) => {
                self.i += 1;
                Some(QAtom::with_query(Node::Leaf(Leaf::Aromatic(false), false)))
            }
            _ => None,
        }
    }

    /// `atom_query` inside brackets.
    fn atom_query(&mut self) -> Result<QAtom, RdkitSmilesError> {
        let c = self.peek().ok_or_else(|| self.err("atom query"))?;
        // Element tokens: two-letter (and K, V, Y, W, U) before the organic
        // subset, as flex's longest match.
        let iso_start = self.i;
        let iso = self.number();
        if let Some(q) = self.element_token() {
            let mut q = q;
            if let Some(n) = iso {
                q.isotope = n as u32;
                q.query = expand(q.query, Node::Leaf(Leaf::Isotope(n), false), How::And);
            }
            return Ok(q);
        }
        if self.peek() == Some(b'#') {
            self.i += 1;
            let n = self.number().ok_or_else(|| self.err("atomic number"))?;
            let mut q = QAtom::new(n as u32, Node::Leaf(Leaf::AtomicNum(n as u32), false));
            if let Some(iso) = iso {
                q.isotope = iso as u32;
                q.query = expand(q.query, Node::Leaf(Leaf::Isotope(iso), false), How::And);
            }
            return Ok(q);
        }
        if self.peek() == Some(b'H') {
            self.i += 1;
            let n = self.number().unwrap_or(1);
            let h = Node::Leaf(Leaf::Valued("H", Some(n), Range::None), false);
            let mut q = match iso {
                Some(iso) => {
                    let mut q = QAtom::with_query(Node::Leaf(Leaf::Isotope(iso), false));
                    q.isotope = iso as u32;
                    q.query = expand(q.query, h, How::And);
                    q
                }
                None => QAtom::with_query(h),
            };
            q.num_explicit_hs = n as u32;
            q.no_implicit = true;
            q.flags |= H_MASK;
            return Ok(q);
        }
        if let Some(n) = iso {
            // `number` alone: an isotope query.
            let mut q = QAtom::with_query(Node::Leaf(Leaf::Isotope(n), false));
            q.isotope = n as u32;
            return Ok(q);
        }
        self.i = iso_start;
        match c {
            b'+' | b'-' => {
                let ch = self.charge_spec().expect("charge");
                let mut q = QAtom::with_query(Node::Leaf(Leaf::Charge(ch), false));
                q.charge = ch;
                q.flags |= CHARGE_MASK;
                Ok(q)
            }
            b'@' => {
                self.i += 1;
                let chiral = if self.peek() == Some(b'@') {
                    self.i += 1;
                    ChiralTag::Cw
                } else {
                    ChiralTag::Ccw
                };
                let class = [self.peek().unwrap_or(0), self.peek_at(1).unwrap_or(0)];
                if matches!(&class, b"TH" | b"AL" | b"SP" | b"TB" | b"OH") {
                    return Err(unsupported("chirality class in SMARTS"));
                }
                let mut q = QAtom::with_query(Node::Leaf(Leaf::Null, false));
                q.chiral = chiral;
                Ok(q)
            }
            b'^' => {
                self.i += 1;
                let n = self.number().ok_or_else(|| self.err("hybridization"))?;
                Ok(QAtom::with_query(Node::Leaf(Leaf::Hybridization(n), false)))
            }
            b'D' | b'd' | b'X' | b'v' | b'R' | b'x' | b'z' | b'Z' | b'h' | b'r' | b'k' => {
                self.i += 1;
                let sym: &'static str = match c {
                    b'D' => "D",
                    b'd' => "d",
                    b'X' => "X",
                    b'v' => "v",
                    b'R' => "R",
                    b'x' => "x",
                    b'z' => "z",
                    b'Z' => "Z",
                    b'h' => "h",
                    b'r' => "r",
                    _ => "k",
                };
                let n = self.number();
                let mut range = Range::None;
                if self.peek() == Some(b'{') {
                    self.i += 1;
                    if self.peek() == Some(b'-') {
                        self.i += 1;
                        range = Range::Greater(self.number().ok_or_else(|| self.err("range"))?);
                    } else {
                        let lo = self.number().ok_or_else(|| self.err("range"))?;
                        if self.peek() != Some(b'-') {
                            return Err(self.err("range"));
                        }
                        self.i += 1;
                        range = match self.number() {
                            Some(hi) => Range::Between(lo, hi),
                            None => Range::Less(lo),
                        };
                    }
                    if self.peek() != Some(b'}') {
                        return Err(self.err("range"));
                    }
                    self.i += 1;
                }
                // Defaults of the bare tokens.
                let val = match (sym, n) {
                    (_, Some(v)) => Some(v),
                    ("D" | "d" | "X" | "v", None) => Some(1),
                    ("R", None) => Some(-1),
                    _ => None,
                };
                // Bare `r` and `k` are "AtomInRing" (written `R`).
                let sym = match (sym, n, range) {
                    ("r" | "k", None, Range::None) => "InRing",
                    _ => sym,
                };
                Ok(QAtom::with_query(Node::Leaf(
                    Leaf::Valued(sym, val, range),
                    false,
                )))
            }
            _ => self.simple_atom(true).ok_or_else(|| self.err("atom query")),
        }
    }

    /// `ATOM_TOKEN`: element symbols other than the organic subset.
    fn element_token(&mut self) -> Option<QAtom> {
        let c = self.peek()?;
        if let Some(t) = self.peek_at(1) {
            let pair = [c, t];
            let pair = std::str::from_utf8(&pair).ok()?;
            if TWO_LETTER.contains(&pair) {
                // `si`, `as`, `se`, `te` are aromatic tokens (lowercase); the
                // list holds capitalized symbols only.
                self.i += 2;
                let anum = periodic::atomic_number(pair)?;
                return Some(QAtom::new(anum, Node::Leaf(Leaf::AtomicNum(anum), false)));
            }
        }
        let anum = match c {
            b'K' => 19,
            b'V' => 23,
            b'Y' => 39,
            b'W' => 74,
            b'U' => 92,
            _ => return self.simple_atom(true),
        };
        self.i += 1;
        Some(QAtom::new(anum, Node::Leaf(Leaf::AtomicNum(anum), false)))
    }
}

/// `atom_expr AND|OR|SEMI atom_expr`.
fn combine_exprs(mut a: QAtom, b: QAtom, how: How) -> QAtom {
    a.query = expand(a.query, b.query, how);
    if a.chiral == ChiralTag::Unspecified {
        a.chiral = b.chiral;
    }
    a.clear_chemical_props();
    a
}

/// `atom_expr_and_point_query`.
fn atom_expr_and_point_query(e: &mut QAtom, p: QAtom) {
    let q = std::mem::replace(&mut e.query, Node::Leaf(Leaf::Null, false));
    e.query = expand(q, p.query, How::And);
    if e.chiral == ChiralTag::Unspecified {
        e.chiral = p.chiral;
    }
    if p.flags & H_MASK != 0 {
        if e.flags & H_MASK == 0 {
            e.num_explicit_hs = p.num_explicit_hs;
            e.no_implicit = true;
            e.flags |= H_MASK;
        } else if e.num_explicit_hs != p.num_explicit_hs {
            e.num_explicit_hs = 0;
            e.no_implicit = true;
        }
    }
    if p.flags & CHARGE_MASK != 0 {
        if e.flags & CHARGE_MASK == 0 {
            e.charge = p.charge;
            e.flags |= CHARGE_MASK;
        } else if e.charge != p.charge {
            e.charge = 0;
        }
    }
}

// ------------------------------------------------------------------ writer

const HAS_AND: u32 = 1;
const HAS_LOWAND: u32 = 2;
const HAS_OR: u32 = 4;
const HAS_RECURSION: u32 = 8;

/// `_combineChildSmarts`.
fn combine_child_smarts(
    cs1: String,
    f1: u32,
    cs2: String,
    f2: u32,
    how: How,
    features: &mut u32,
) -> Result<String, RdkitSmilesError> {
    let mut res = String::new();
    match how {
        How::Or => {
            if (f1 & HAS_LOWAND != 0 && f1 & HAS_OR != 0)
                || (f2 & HAS_LOWAND != 0 && f2 & HAS_OR != 0)
            {
                return Err(unsupported(
                    "This is a non-smartable query - OR above and below AND in the binary tree",
                ));
            }
            res.push_str(&cs1);
            if !(cs1.is_empty() || cs2.is_empty()) {
                res.push(',');
            }
            res.push_str(&cs2);
            *features |= HAS_OR;
        }
        How::And => {
            let symb = if f1 & HAS_OR != 0 || f2 & HAS_OR != 0 {
                *features |= HAS_LOWAND;
                ';'
            } else {
                *features |= HAS_AND;
                '&'
            };
            res.push_str(&cs1);
            if !(cs1.is_empty() || cs2.is_empty()) {
                res.push(symb);
            }
            res.push_str(&cs2);
        }
    }
    *features |= f1 | f2;
    Ok(res)
}

/// Per-atom state of the writer (`_qatomHasStereoSet`).
struct AtomCtx {
    chiral: ChiralTag,
    stereo_done: bool,
}

/// `getAtomSmartsSimple`.
fn atom_smarts_simple(leaf: &Leaf, ctx: &mut AtomCtx, need_paren: &mut bool) -> String {
    use std::fmt::Write;
    let mut res = String::new();
    let mut val: Option<i32> = None;
    let mut range = Range::None;
    match leaf {
        Leaf::AtomType(anum, aromatic) => {
            let sym = periodic::symbol(*anum);
            if *aromatic {
                res.push(char::from(sym.as_bytes()[0].to_ascii_lowercase()));
                res.push_str(&sym[1..]);
            } else {
                res.push_str(sym);
            }
            if !matches!(anum, 0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53) {
                *need_paren = true;
            }
        }
        Leaf::AtomicNum(n) => {
            let _ = write!(res, "#{n}");
            *need_paren = true;
        }
        Leaf::Null => res.push('*'),
        Leaf::Aromatic(true) => res.push('a'),
        Leaf::Aromatic(false) => res.push('A'),
        Leaf::Charge(v) => {
            res.push(if *v < 0 { '-' } else { '+' });
            if v.abs() != 1 {
                let _ = write!(res, "{}", v.abs());
            }
            *need_paren = true;
        }
        Leaf::Isotope(v) => {
            let _ = write!(res, "{v}*");
            *need_paren = true;
        }
        Leaf::Hybridization(v) => {
            res.push('^');
            if (0..=5).contains(v) {
                let _ = write!(res, "{v}");
            }
            *need_paren = true;
        }
        Leaf::Valued(sym, v, r) => {
            *need_paren = true;
            range = *r;
            match *sym {
                "InRing" => res.push('R'),
                "R" => {
                    res.push('R');
                    if *r == Range::None && v.is_some_and(|x| x >= 0) {
                        val = *v;
                    }
                }
                s => {
                    res.push_str(s);
                    val = *v;
                }
            }
        }
    }
    match range {
        Range::None => {
            if let Some(v) = val {
                let _ = write!(res, "{v}");
            }
        }
        Range::Less(v) => {
            let _ = write!(res, "{{{v}-}}");
        }
        Range::Greater(v) => {
            let _ = write!(res, "{{-{v}}}");
        }
        Range::Between(a, b) => {
            let _ = write!(res, "{{{a}-{b}}}");
        }
    }
    if !ctx.stereo_done && matches!(ctx.chiral, ChiralTag::Cw | ChiralTag::Ccw) {
        ctx.stereo_done = true;
        res.push_str(if ctx.chiral == ChiralTag::Cw {
            "@@"
        } else {
            "@"
        });
        *need_paren = true;
    }
    res
}

/// `_recurseGetSmarts` for a composite node.
fn recurse_atom_smarts(
    node: &Node,
    negate: bool,
    features: &mut u32,
    ctx: &mut AtomCtx,
) -> Result<String, RdkitSmilesError> {
    let (children, mut how) = match node {
        Node::And(c, _) => (c, How::And),
        Node::Or(c, _) => (c, How::Or),
        _ => unreachable!("composite"),
    };
    let child_smarts = |child: &Node,
                        features: &mut u32,
                        child_features: &mut u32,
                        ctx: &mut AtomCtx|
     -> Result<String, RdkitSmilesError> {
        Ok(match child {
            Node::Recursive(m, neg) => {
                *features |= HAS_RECURSION;
                recursive_smarts(m, *neg)?
            }
            Node::Leaf(leaf, neg) => {
                let mut np = false;
                let s = atom_smarts_simple(leaf, ctx, &mut np);
                if negate ^ neg { format!("!{s}") } else { s }
            }
            Node::And(_, neg) | Node::Or(_, neg) => {
                recurse_atom_smarts(child, negate ^ neg, child_features, ctx)?
            }
        })
    };
    let mut f1 = 0;
    let cs1 = child_smarts(&children[0], features, &mut f1, ctx)?;
    if negate {
        how = if how == How::Or { How::And } else { How::Or };
    }
    let mut res = cs1;
    for child in &children[1..] {
        let mut f2 = 0;
        let cs2 = child_smarts(child, features, &mut f2, ctx)?;
        res = combine_child_smarts(res, f1, cs2, f2, how, features)?;
    }
    Ok(res)
}

/// `getRecursiveStructureQuerySmarts`.
fn recursive_smarts(m: &QMol, neg: bool) -> Result<String, RdkitSmilesError> {
    let s = format!("$({})", mol_to_smarts(m)?);
    Ok(if neg { format!("!{s}") } else { s })
}

/// `SmartsWrite::GetAtomSmarts` for a query atom.
fn atom_smarts(atom: &QAtom, chiral: ChiralTag) -> Result<String, RdkitSmilesError> {
    let mut ctx = AtomCtx {
        chiral,
        stereo_done: false,
    };
    let mut need_paren = false;
    let mut res = match &atom.query {
        Node::And(_, neg) | Node::Or(_, neg) => {
            need_paren = true;
            let mut features = 0;
            let r = recurse_atom_smarts(&atom.query, *neg, &mut features, &mut ctx)?;
            if r.len() == 1 {
                need_paren = false;
            }
            r
        }
        Node::Recursive(m, neg) => {
            need_paren = true;
            recursive_smarts(m, *neg)?
        }
        Node::Leaf(leaf, neg) => {
            let mut r = atom_smarts_simple(leaf, &mut ctx, &mut need_paren);
            if *neg {
                r.insert(0, '!');
                need_paren = true;
            }
            r
        }
    };
    if let Some(m) = atom.map {
        need_paren = true;
        res.push_str(&format!(":{m}"));
    }
    if need_paren {
        res = format!("[{res}]");
    }
    Ok(res)
}

/// `getBondSmartsSimple` (no directions); `reverse_dative`: the bond
/// begins at the atom written to its right.
fn bond_smarts_simple(leaf: &Leaf, reverse_dative: bool, dir: BondDir) -> String {
    let dir_char = || match dir {
        BondDir::EndDownRight => Some("\\"),
        BondDir::EndUpRight => Some("/"),
        BondDir::None => None,
    };
    match leaf {
        Leaf::Valued("BondNull", ..) => "~".into(),
        Leaf::Valued("BondInRing", ..) => "@".into(),
        Leaf::Null => dir_char().unwrap_or("").into(), // `SingleOrAromaticBond`
        Leaf::Valued("BondOrder", Some(bt), _) => match *bt {
            x if x == BondType::Single as i32 => dir_char().unwrap_or("-").into(),
            x if x == BondType::Double as i32 => "=".into(),
            x if x == BondType::Triple as i32 => "#".into(),
            x if x == BondType::Quadruple as i32 => "$".into(),
            x if x == BondType::Aromatic as i32 => dir_char().unwrap_or(":").into(),
            x if x == BondType::Dative as i32 => {
                if reverse_dative {
                    "<-".into()
                } else {
                    "->".into()
                }
            }
            _ => String::new(),
        },
        _ => String::new(),
    }
}

/// `_recurseBondSmarts` (including RDKit's handling of a composite second
/// child, which overwrites the first child's text).
fn recurse_bond_smarts(
    node: &Node,
    negate: bool,
    features: &mut u32,
    rev: bool,
    dir: BondDir,
) -> Result<String, RdkitSmilesError> {
    let (children, mut how) = match node {
        Node::And(c, _) => (c, How::And),
        Node::Or(c, _) => (c, How::Or),
        _ => unreachable!("composite"),
    };
    let (mut f1, mut f2) = (0, 0);
    let mut cs1: String;
    let mut cs2 = String::new();
    match &children[0] {
        Node::Leaf(l, neg) => {
            cs1 = bond_smarts_simple(l, rev, dir);
            if negate ^ neg {
                cs1.insert(0, '!');
            }
        }
        c => cs1 = recurse_bond_smarts(c, negate ^ c.negation(), &mut f1, rev, dir)?,
    }
    match &children[1] {
        Node::Leaf(l, neg) => {
            cs2 = bond_smarts_simple(l, rev, dir);
            if negate ^ neg {
                cs2.insert(0, '!');
            }
        }
        c => cs1 = recurse_bond_smarts(c, negate ^ c.negation(), &mut f2, rev, dir)?,
    }
    if negate {
        how = if how == How::Or { How::And } else { How::Or };
    }
    combine_child_smarts(cs1, f1, cs2, f2, how, features)
}

/// `SmartsWrite::GetBondSmarts` for a query bond.
fn bond_smarts(bond: &QBond, left_atom: usize, dir: BondDir) -> Result<String, RdkitSmilesError> {
    let rev = bond.begin != left_atom;
    match &bond.query {
        Node::Leaf(l, neg) => {
            let s = bond_smarts_simple(l, rev, dir);
            Ok(if *neg { format!("!{s}") } else { s })
        }
        n @ (Node::And(_, neg) | Node::Or(_, neg)) => {
            let mut f = 0;
            recurse_bond_smarts(n, *neg, &mut f, rev, dir)
        }
        Node::Recursive(..) => Err(unsupported("recursive bond query")),
    }
}

/// RDKit's molecule of a parsed template: the topology (RDKit's bond
/// order) and chiral tags from the SMILES port's reading of a stand-in
/// SMILES, then the query atoms' own properties. Returns the molecule and,
/// per RDKit bond, the template bond.
fn rdkit_mol_of(m: &QMol) -> Result<(Mol, Vec<usize>), RdkitSmilesError> {
    use std::fmt::Write;
    let directed = m.bonds.iter().any(|b| b.dir != BondDir::None);
    let mut s = String::new();
    for t in &m.toks {
        match *t {
            Tok::Atom(a) => {
                let atom = &m.atoms[a];
                s.push('[');
                if atom.anum == 0 || atom.anum > 118 {
                    s.push('*');
                } else {
                    s.push_str(periodic::symbol(atom.anum));
                }
                match atom.chiral {
                    ChiralTag::Cw => s.push_str("@@"),
                    ChiralTag::Ccw => s.push('@'),
                    _ => {}
                }
                if atom.num_explicit_hs > 0 {
                    let _ = write!(s, "H{}", atom.num_explicit_hs);
                }
                s.push(']');
            }
            // With directions, the stand-in carries them (and the double
            // bonds) so the reader sets each bond's direction as RDKit's
            // grammar does.
            Tok::Bond(b) if directed => s.push(match (m.bonds[b].dir, m.bonds[b].bt) {
                (BondDir::EndUpRight, _) => '/',
                (BondDir::EndDownRight, _) => '\\',
                (_, BondType::Double) => '=',
                _ => '-',
            }),
            Tok::Bond(_) => s.push('-'),
            Tok::Open => s.push('('),
            Tok::Close => s.push(')'),
            Tok::Dot => s.push('.'),
            Tok::Ring(n, _, sym) => {
                if directed && let Some(c) = sym {
                    s.push(c);
                }
                if n < 10 {
                    let _ = write!(s, "{n}");
                } else if n < 100 {
                    let _ = write!(s, "%{n}");
                } else {
                    let _ = write!(s, "%({n})");
                }
            }
        }
    }
    let chem = crate::parse_template(&s)
        .map_err(|e| unsupported(format!("template topology '{s}': {e}")))?;
    if chem.atom_count() != m.atoms.len() || chem.bond_count() != m.bonds.len() {
        return Err(unsupported("template topology mismatch"));
    }
    let mut mol = super::parse::from_chematic(&chem)?;
    for (a, qa) in m.atoms.iter().enumerate() {
        let at = &mut mol.atoms[a];
        at.anum = qa.anum;
        at.aromatic = qa.aromatic;
        at.isotope = qa.isotope;
        at.charge = qa.charge;
        at.num_explicit_hs = qa.num_explicit_hs;
        at.no_implicit = qa.no_implicit;
        at.map = qa.map;
        at.single_h_query = has_single_h_query(&qa.query);
        // `calculateImplicitValence`: 0 next to a bond with a complex bond
        // type query (RDKit computes it; no implicit Hs is the same here).
        if m.bonds
            .iter()
            .any(|b| (b.begin == a || b.end == a) && complex_bond_type_query(&b.query, false))
        {
            at.no_implicit = true;
        }
    }
    let mut which = Vec::with_capacity(mol.bonds.len());
    for b in &mut mol.bonds {
        let k = m
            .bonds
            .iter()
            .position(|q| {
                (q.begin == b.begin && q.end == b.end) || (q.begin == b.end && q.end == b.begin)
            })
            .ok_or_else(|| unsupported("template bond mismatch"))?;
        b.bt = m.bonds[k].bt;
        b.aromatic = b.bt == BondType::Aromatic;
        which.push(k);
    }
    for a in 0..mol.atoms.len() {
        mol.update_atom_property_cache(a, false)?;
    }
    if directed {
        set_bond_stereo_from_directions(&mut mol);
    } else {
        for b in &mut mol.bonds {
            b.dir = BondDir::None;
        }
    }
    Ok((mol, which))
}

/// `Canon::details::hasSingleHQuery`.
fn has_single_h_query(q: &Node) -> bool {
    let Node::And(children, _) = q else {
        return false;
    };
    for c in children {
        match c {
            Node::Leaf(Leaf::Valued("H", v, Range::None), neg) => return !neg && *v == Some(1),
            Node::And(..) if has_single_h_query(c) => return true,
            _ => {}
        }
    }
    false
}

/// `QueryOps::hasComplexBondTypeQuery` (`seen`: a `BondOrder` query was
/// seen higher up or in an earlier sibling).
fn complex_bond_type_query(q: &Node, mut seen: bool) -> bool {
    let (is_order_fn, is_bond_order) = match q {
        Node::Leaf(Leaf::Null, _) => (true, false), // `SingleOrAromaticBond`
        Node::Leaf(Leaf::Valued("BondOrder", ..), _) => (true, true),
        _ => (false, false),
    };
    if is_order_fn && (seen || !is_bond_order || q.negation()) {
        return true;
    }
    if let Node::And(children, _) | Node::Or(children, _) = q {
        for c in children {
            if complex_bond_type_query(c, seen | is_bond_order) {
                return true;
            }
            if matches!(c, Node::Leaf(Leaf::Valued("BondOrder", ..), _)) {
                seen = true;
            }
        }
    }
    false
}

/// `MolOps::setBondStereoFromDirections`: each double bond with a
/// directed (non-double) bond at both ends gets those neighbours as stereo
/// atoms and cis (Z) or trans (E) from the two directions.
fn set_bond_stereo_from_directions(mol: &mut Mol) {
    for d in 0..mol.bonds.len() {
        if mol.bonds[d].bt != BondType::Double || mol.bonds[d].stereo == BondStereo::Any {
            continue;
        }
        let (b0, e0) = (mol.bonds[d].begin, mol.bonds[d].end);
        let directed = |atom: usize| {
            mol.atom_bonds[atom].iter().copied().find(|&nb| {
                mol.bonds[nb].bt != BondType::Double && mol.bonds[nb].dir != BondDir::None
            })
        };
        let (Some(at_begin), Some(at_end)) = (directed(b0), directed(e0)) else {
            continue;
        };
        let begin_atom = mol.bonds[at_begin].other(b0);
        let end_atom = mol.bonds[at_end].other(e0);
        let mut begin_dir = mol.bonds[at_begin].dir;
        if mol.bonds[at_begin].begin == b0 {
            begin_dir = begin_dir.flipped();
        }
        let mut end_dir = mol.bonds[at_end].dir;
        if mol.bonds[at_end].end == e0 {
            end_dir = end_dir.flipped();
        }
        let bond = &mut mol.bonds[d];
        bond.stereo_atoms = vec![begin_atom, end_atom];
        bond.stereo = if begin_dir == end_dir {
            BondStereo::E
        } else {
            BondStereo::Z
        };
    }
}

/// `MolToSmarts(mol)` (isomeric, no root) for a parsed template.
fn mol_to_smarts(m: &QMol) -> Result<String, RdkitSmilesError> {
    if m.atoms.is_empty() {
        return Ok(String::new());
    }
    let (mol, which) = rdkit_mol_of(m)?;
    super::smarts_write::mol_to_smarts_custom(mol, &mut |res, e, sub, atoms, bonds| {
        match e {
            StackElem::Atom(a) => {
                res.push_str(&atom_smarts(&m.atoms[atoms[a]], sub.atoms[a].chiral)?)
            }
            StackElem::Bond(b, left) => res.push_str(&bond_smarts(
                &m.bonds[which[bonds[b]]],
                atoms[left],
                sub.bonds[b].dir,
            )?),
            _ => unreachable!("atoms and bonds only"),
        }
        Ok(())
    })
}

fn parse_template(s: &str) -> Result<QMol, RdkitSmilesError> {
    let mut p = Parser {
        s: s.as_bytes(),
        i: 0,
    };
    let m = p.mol(false)?;
    if p.i != s.len() {
        return Err(p.err("trailing input"));
    }
    Ok(m)
}

/// The number of connected components of a template.
fn num_frags(m: &QMol) -> usize {
    let n = m.atoms.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        p[x] = r;
        r
    }
    for b in &m.bonds {
        let (x, y) = (find(&mut parent, b.begin), find(&mut parent, b.end));
        parent[x] = y;
    }
    (0..n).filter(|&i| find(&mut parent, i) == i).count()
}

/// `Chem.MolToSmarts(Chem.MolFromSmarts(smarts))` (RDKit 2026.03.1).
pub(crate) fn smarts_to_smarts(smarts: &str) -> Result<String, RdkitSmilesError> {
    mol_to_smarts(&parse_template(smarts)?)
}

/// The components of one reaction side: `.`-separated outside brackets,
/// `(A.B)` groups kept together (without the parentheses).
fn side_components(side: &str) -> Result<Vec<&str>, RdkitSmilesError> {
    let b = side.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'(' {
            let (mut depth, mut bracket, mut close) = (0usize, 0usize, None);
            for (k, &c) in b.iter().enumerate().skip(i) {
                match c {
                    b'[' => bracket += 1,
                    b']' => bracket = bracket.saturating_sub(1),
                    b'(' if bracket == 0 => depth += 1,
                    b')' if bracket == 0 => {
                        depth -= 1;
                        if depth == 0 {
                            close = Some(k);
                            break;
                        }
                    }
                    _ => {}
                }
            }
            let close = close.ok_or_else(|| unsupported("unclosed component group"))?;
            out.push(&side[i + 1..close]);
            i = close + 1;
            if i < b.len() && b[i] != b'.' {
                return Err(unsupported("component group must be followed by '.'"));
            }
            i += 1;
            continue;
        }
        let (mut bracket, mut paren) = (0usize, 0usize);
        let mut end = b.len();
        for (k, &c) in b.iter().enumerate().skip(i) {
            match c {
                b'[' => bracket += 1,
                b']' => bracket = bracket.saturating_sub(1),
                b'(' if bracket == 0 => paren += 1,
                b')' if bracket == 0 => paren = paren.saturating_sub(1),
                b'.' if bracket == 0 && paren == 0 => {
                    end = k;
                    break;
                }
                _ => {}
            }
        }
        if end > i {
            out.push(&side[i..end]);
        }
        i = end + 1;
    }
    Ok(out)
}

/// `rdChemReactions.ReactionToSmarts(rdChemReactions.ReactionFromSmarts(s))`
/// (RDKit 2026.03.1): every template written by `MolToSmarts`
/// (non-canonical, isomeric), a template with several components wrapped
/// in parentheses, sides joined by `>`.
pub(crate) fn reaction_to_smarts(rxn: &str) -> Result<String, RdkitSmilesError> {
    // Split as `parseReaction` does: on every '>' not preceded by '-' (a
    // dative bond `->`), brackets or not.
    let b = rxn.as_bytes();
    let mut parts = Vec::new();
    let mut start = 0usize;
    for (k, &c) in b.iter().enumerate() {
        if c == b'>' && (k == 0 || b[k - 1] != b'-') {
            parts.push(&rxn[start..k]);
            start = k + 1;
        }
    }
    parts.push(&rxn[start..]);
    if parts.len() != 3 {
        return Err(unsupported("reaction SMARTS needs two '>'"));
    }
    let mut out = Vec::with_capacity(3);
    for side in parts {
        let mut written = Vec::new();
        for comp in side_components(side)? {
            let m = parse_template(comp)?;
            let s = mol_to_smarts(&m)?;
            written.push(if num_frags(&m) > 1 {
                format!("({s})")
            } else {
                s
            });
        }
        out.push(written.join("."));
    }
    Ok(out.join(">"))
}
