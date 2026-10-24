use super::*;

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

pub(super) fn parse_template(s: &str) -> Result<QMol, RdkitSmilesError> {
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
