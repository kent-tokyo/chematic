//! SMARTS parser: converts a SMARTS string into a `QueryMolecule`.
//!
//! ## Supported grammar (subset)
//!
//! ```text
//! smarts     := chain
//! chain      := atom chain_rest*
//! chain_rest := branch | ring_closure | bond? atom | bond? ring_closure
//! branch     := '(' bond? chain ')'
//! atom       := bracket_atom | organic_atom | '*'
//!
//! bracket_atom := '[' expr ']'
//! expr         := low_and
//! low_and      := or    (';' or)*               // lowest precedence AND
//! or           := high_and (',' high_and)*       // OR
//! high_and     := unary  ('&'? unary)*           // high-prec AND (explicit '&' or juxtaposition)
//! unary        := '!' unary | primitive
//! primitive    := '#' DIGITS                     // atomic number
//!               | 'a'                            // aromatic
//!               | 'A'                            // aliphatic
//!               | '+' DIGITS? | '-' DIGITS?      // charge
//!               | 'H' DIGIT?                     // H count (no digit = H1)
//!               | 'D' DIGIT                      // degree
//!               | 'r' DIGIT                      // ring size (exact in chematic)
//!               | 'k' DIGIT                      // exact ring size (RDKit [kN] compat — alias for [rN])
//!               | 'R'                            // in a ring
//!               | '*'                            // wildcard
//!               | element_symbol                 // 'C', 'Cl', 'n', …
//!
//! bond := '-' | '=' | '#' | ':' | '~' | '@'
//! ```
//!
//! Organic-subset shorthand:
//! - Uppercase letter → `And(Symbol("X"), Aromatic(false))` (aliphatic)
//! - Lowercase letter → `And(Symbol("X"), Aromatic(true))`  (aromatic)
//! - `*` outside brackets → `Wildcard`
//!
//! Implicit bond between two adjacent atoms = `BondQuery::Any`, which matches a
//! single or aromatic bond (Daylight/RDKit semantics; `~` is the explicit any-bond).

mod rings;

use std::collections::HashMap;
use std::fmt;

use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur during SMARTS parsing.
#[derive(Debug, Clone, PartialEq)]
pub enum SmartsError {
    /// Input ended unexpectedly.
    UnexpectedEnd,
    /// An unexpected character was found at the given position.
    UnexpectedChar(char, usize),
    /// A `[` was opened but never closed.
    UnclosedBracket(usize),
    /// A `(` was opened but never closed.
    UnclosedBranch(usize),
    /// A ring closure number was used inconsistently.
    InvalidRingClosure(u8),
    /// Recursive SMARTS `$(…)` nesting exceeds the safety limit.
    RecursionDepthExceeded,
}

impl fmt::Display for SmartsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SmartsError::UnexpectedEnd => write!(f, "unexpected end of SMARTS string"),
            SmartsError::UnexpectedChar(ch, pos) => {
                write!(f, "unexpected character {:?} at position {}", ch, pos)
            }
            SmartsError::UnclosedBracket(pos) => {
                write!(f, "unclosed '[' bracket at position {}", pos)
            }
            SmartsError::UnclosedBranch(pos) => {
                write!(f, "unclosed '(' branch at position {}", pos)
            }
            SmartsError::InvalidRingClosure(num) => {
                write!(f, "invalid ring closure number: {}", num)
            }
            SmartsError::RecursionDepthExceeded => {
                write!(f, "recursive SMARTS nesting depth exceeded safety limit")
            }
        }
    }
}

impl std::error::Error for SmartsError {}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Default maximum nesting depth for recursive SMARTS `$(…)` patterns.
const DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH: usize = 8;
/// Absolute maximum allowed depth (safety limit).
const ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH: usize = 16;

/// Configuration for SMARTS parsing.
#[derive(Debug, Clone)]
pub struct SmartsParserConfig {
    /// Maximum nesting depth for recursive SMARTS `$(…)` patterns.
    /// Default: 8. Max: 16 (safety limit to prevent stack overflow).
    pub max_recursion_depth: usize,
}

impl Default for SmartsParserConfig {
    fn default() -> Self {
        Self {
            max_recursion_depth: DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH,
        }
    }
}

/// Parse a SMARTS pattern string into a `QueryMolecule`.
///
/// Returns `Err(SmartsError)` if the input is not a valid SMARTS pattern.
/// Uses default config (max recursion depth 8).
pub fn parse_smarts(smarts: &str) -> Result<QueryMolecule, SmartsError> {
    parse_smarts_with_config(smarts, &SmartsParserConfig::default())
}

/// Parse a SMARTS pattern string with explicit configuration.
///
/// `config.max_recursion_depth` controls the maximum nesting depth for recursive
/// SMARTS `$(…)` patterns. Clamped to [1, 16].
pub fn parse_smarts_with_config(
    smarts: &str,
    config: &SmartsParserConfig,
) -> Result<QueryMolecule, SmartsError> {
    let max_depth = config
        .max_recursion_depth
        .clamp(1, ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH);
    let mut parser = Parser {
        src: smarts.as_bytes(),
        pos: 0,
        recursion_depth: 0,
        max_recursion_depth: max_depth,
    };
    parser.parse()
}

// ---------------------------------------------------------------------------
// Parser internals
// ---------------------------------------------------------------------------

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    /// Current recursive SMARTS `$(…)` nesting depth.
    recursion_depth: usize,
    /// Maximum allowed recursive SMARTS depth.
    max_recursion_depth: usize,
}

impl<'a> Parser<'a> {
    // -- character helpers ---------------------------------------------------

    #[inline]
    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    #[inline]
    fn advance(&mut self) -> Option<u8> {
        let b = self.src.get(self.pos).copied();
        if b.is_some() {
            self.pos += 1;
        }
        b
    }

    // -- top-level -----------------------------------------------------------

    fn parse(&mut self) -> Result<QueryMolecule, SmartsError> {
        let mut mol = QueryMolecule::new();
        // Preserve ring-opening order for tetrahedral `@`/`@@` parity.
        let mut open_rings: HashMap<u8, (usize, Option<BondQuery>, usize)> = HashMap::new();
        self.parse_chain(&mut mol, None, None, &mut open_rings)?;
        // `.` separates disconnected components (e.g. RDKit QED alert
        // `F.F.F.F`): each component matches independently, atoms distinct.
        while self.peek() == Some(b'.') {
            self.advance(); // consume '.'
            let before = mol.atoms.len();
            self.parse_chain(&mut mol, None, None, &mut open_rings)?;
            if mol.atoms.len() == before {
                return Err(SmartsError::UnexpectedEnd);
            }
        }

        // Any remaining open ring closures are errors.
        if let Some((&num, _)) = open_rings.iter().next() {
            return Err(SmartsError::InvalidRingClosure(num));
        }

        // A direct SMARTS API is deliberately strict: accepting a valid
        // prefix followed by whitespace/junk silently changes a query.  UI
        // callers that want to trim an editor value must do so explicitly.
        if let Some(c) = self.peek() {
            return Err(SmartsError::UnexpectedChar(c as char, self.pos));
        }

        Ok(mol)
    }

    /// Parse one linear chain (atom sequence + branches + ring closures).
    ///
    /// `attach_to`: if `Some(idx)`, connect the first atom of this chain to
    /// atom `idx` via `attach_bond` (or implicit Any if not specified).
    fn parse_chain(
        &mut self,
        mol: &mut QueryMolecule,
        attach_to: Option<usize>,
        attach_bond: Option<BondQuery>,
        open_rings: &mut HashMap<u8, (usize, Option<BondQuery>, usize)>,
    ) -> Result<Option<usize>, SmartsError> {
        // Parse the first atom of this chain.
        let (first_atom, first_map) = match self.try_parse_atom()? {
            Some(a) => a,
            None => return Ok(attach_to),
        };

        let first_idx = mol.add_atom_with_map(first_atom, first_map);

        // Connect to the previous atom if requested.
        if let Some(prev) = attach_to {
            let bond = attach_bond.unwrap_or(BondQuery::Any);
            mol.add_bond(prev, first_idx, bond);
        }

        let mut current = first_idx;

        loop {
            match self.peek() {
                // Branch: `(` opens a new sub-chain branching from `current`.
                Some(b'(') => {
                    let branch_start = self.pos;
                    self.advance(); // consume '('
                    let branch_bond = self.try_parse_bond();
                    self.parse_chain(mol, Some(current), branch_bond, open_rings)?;
                    match self.peek() {
                        Some(b')') => {
                            self.advance(); // consume ')'
                        }
                        _ => return Err(SmartsError::UnclosedBranch(branch_start)),
                    }
                }

                // Ring closure: single digit or `%nn`.
                Some(b'0'..=b'9') | Some(b'%') => {
                    let (ring_num, ring_bond) = self.parse_ring_closure_num(None)?;
                    self.handle_ring_closure(mol, current, ring_num, ring_bond, open_rings)?;
                }

                // End of this chain: `)`, `.`, end-of-input, or other stop character.
                None | Some(b')') | Some(b'.') => break,

                // Explicit bond or next atom.
                _ => {
                    let pending_bond = self.try_parse_bond();

                    match self.peek() {
                        // Ring closure after an explicit bond character.
                        Some(b'0'..=b'9') | Some(b'%') => {
                            let (ring_num, ring_bond) =
                                self.parse_ring_closure_num(pending_bond)?;
                            self.handle_ring_closure(
                                mol, current, ring_num, ring_bond, open_rings,
                            )?;
                        }

                        // Branch after an explicit bond character (e.g. `C=(C)C`).
                        Some(b'(') => {
                            let branch_start = self.pos;
                            self.advance(); // consume '('
                            self.parse_chain(mol, Some(current), pending_bond, open_rings)?;
                            match self.peek() {
                                Some(b')') => {
                                    self.advance();
                                }
                                _ => return Err(SmartsError::UnclosedBranch(branch_start)),
                            }
                        }

                        // End of input — a trailing bond character is an error.
                        None | Some(b')') => {
                            if pending_bond.is_some() {
                                return Err(SmartsError::UnexpectedEnd);
                            }
                            break;
                        }

                        // Next atom in the chain.
                        _ => match self.try_parse_atom()? {
                            Some((next_atom, next_map)) => {
                                let next_idx = mol.add_atom_with_map(next_atom, next_map);
                                let bond = pending_bond.unwrap_or(BondQuery::Any);
                                mol.add_bond(current, next_idx, bond);
                                current = next_idx;
                            }
                            None => {
                                if pending_bond.is_some() {
                                    return Err(SmartsError::UnexpectedEnd);
                                }
                                break;
                            }
                        },
                    }
                }
            }
        }

        Ok(Some(current))
    }

    /// Parse a ring number and its optional, already-consumed bond query.
    fn parse_ring_closure_num(
        &mut self,
        prefix_bond: Option<BondQuery>,
    ) -> Result<(u8, Option<BondQuery>), SmartsError> {
        let ring_num = if self.peek() == Some(b'%') {
            self.advance(); // consume '%'
            let tens = self
                .advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmartsError::UnexpectedEnd)?
                - b'0';
            let units = self
                .advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmartsError::UnexpectedEnd)?
                - b'0';
            tens * 10 + units
        } else {
            // Caller peeked b'0'..=b'9', so advance() is guaranteed to return Some(digit).
            self.advance()
                .expect("ring closure digit guaranteed by caller peek")
                - b'0'
        };
        Ok((ring_num, prefix_bond))
    }

    // -- bond ----------------------------------------------------------------

    /// Try to parse a bond expression and return the corresponding `BondQuery`.
    ///
    /// Handles compound expressions used in PAINS and similar SMARTS catalogs:
    /// - `=,:` → OR(Double, Aromatic)
    /// - `-,:` → OR(Single, Aromatic)
    /// - `=!@` → AND(Double, NOT(Ring))
    /// - `-!@` → AND(Single, NOT(Ring))
    ///
    /// Grammar (simplified):
    /// ```text
    /// bond_expr  := bond_or
    /// bond_or    := bond_and (',' bond_and)*
    /// bond_and   := bond_unary ('&'? bond_unary)*
    /// bond_unary := '!' bond_token | bond_token
    /// bond_token := '-' | '=' | '#' | ':' | '~' | '@'
    /// ```
    fn try_parse_bond(&mut self) -> Option<BondQuery> {
        // Full Daylight precedence: `!` > implicit/`&` AND > `,` OR > `;` AND.
        let first = self.try_parse_bond_or()?;
        Some(self.parse_bond_low_and_tail(first))
    }

    /// `true` when the byte after the current operator starts a bond factor.
    fn next_starts_bond_factor(&self) -> bool {
        self.src
            .get(self.pos + 1)
            .copied()
            .is_some_and(|c| self.is_bond_token_at(self.pos + 1, c) || c == b'!')
    }

    /// `bond_or (';' bond_or)*` — low-precedence AND.
    fn parse_bond_low_and_tail(&mut self, mut left: BondQuery) -> BondQuery {
        while self.peek() == Some(b';') && self.next_starts_bond_factor() {
            self.advance(); // consume ';'
            match self.try_parse_bond_or() {
                Some(right) => left = BondQuery::And(Box::new(left), Box::new(right)),
                None => break,
            }
        }
        left
    }

    /// `bond_and (',' bond_and)*`.
    fn try_parse_bond_or(&mut self) -> Option<BondQuery> {
        let first = self.try_parse_bond_factor()?;
        Some(self.parse_bond_or_tail(first))
    }

    /// Parse one bond factor: `!token` or `token`.
    fn try_parse_bond_factor(&mut self) -> Option<BondQuery> {
        match self.peek()? {
            b'!' => {
                // Only treat `!` as bond negation when a bond factor follows.
                if self.next_starts_bond_factor() {
                    self.advance(); // consume '!'
                    let inner = self.try_parse_bond_factor()?;
                    Some(BondQuery::Not(Box::new(inner)))
                } else {
                    None
                }
            }
            c if self.is_bond_token_at(self.pos, c) => {
                let prim = self.consume_bond_prim()?;
                Some(BondQuery::Primitive(prim))
            }
            _ => None,
        }
    }

    /// Consume a single bond primitive character (caller must verify `peek` is a bond token).
    fn consume_bond_prim(&mut self) -> Option<BondPrimitive> {
        // Dative bonds `->` / `<-` (RDKit extension).
        if self.peek() == Some(b'-') && self.src.get(self.pos + 1) == Some(&b'>') {
            self.pos += 2;
            return Some(BondPrimitive::DativeForward);
        }
        if self.peek() == Some(b'<') {
            if self.src.get(self.pos + 1) == Some(&b'-') {
                self.pos += 2;
                return Some(BondPrimitive::DativeBackward);
            }
            return None;
        }
        let prim = match self.peek()? {
            b'-' => BondPrimitive::Single,
            b'=' => BondPrimitive::Double,
            b'#' => BondPrimitive::Triple,
            b'$' => BondPrimitive::Quadruple,
            b':' => BondPrimitive::Aromatic,
            b'~' => BondPrimitive::Any,
            b'@' => BondPrimitive::Ring,
            b'/' => BondPrimitive::Up,
            b'\\' => BondPrimitive::Down,
            _ => return None,
        };
        self.advance();
        Some(prim)
    }

    #[inline]
    /// `true` when byte `c` at `pos` starts a bond primitive (`<` only as `<-`).
    fn is_bond_token_at(&self, pos: usize, c: u8) -> bool {
        match c {
            b'-' | b'=' | b'#' | b'$' | b':' | b'~' | b'@' | b'/' | b'\\' => true,
            b'<' => self.src.get(pos + 1) == Some(&b'-'),
            _ => false,
        }
    }

    /// Continue parsing bond OR after the first factor.
    fn parse_bond_or_tail(&mut self, left: BondQuery) -> BondQuery {
        // AND binds tighter than ',' (so `=@,-` is `(= & @) , -`).
        let mut acc = self.parse_bond_and_tail(left);
        // ',' → OR (only when followed by a bond token or '!')
        while self.peek() == Some(b',') && self.next_starts_bond_factor() {
            self.advance(); // consume ','
            let Some(right) = self.try_parse_bond_factor() else {
                break;
            };
            let right = self.parse_bond_and_tail(right);
            acc = BondQuery::Or(Box::new(acc), Box::new(right));
        }
        acc
    }

    /// Continue parsing implicit AND after the first factor.
    ///
    /// Recognises `&!<bond>` or a bare `!<bond>` juxtaposed with the previous
    /// factor as a high-precedence AND. A trailing `&` with no operand is
    /// silently consumed (caller treats the result the same as `left`).
    fn parse_bond_and_tail(&mut self, left: BondQuery) -> BondQuery {
        // `&factor` or a directly juxtaposed factor (`=@`, `-!@`, ...) is a
        // high-precedence AND.
        let explicit = self.peek() == Some(b'&') && self.next_starts_bond_factor();
        if explicit {
            self.advance(); // consume '&'
        } else if self.peek() == Some(b'&') {
            // A trailing `&` with no operand is silently consumed.
            self.advance();
            return left;
        }
        let juxtaposed = self.peek().is_some_and(|c| {
            self.is_bond_token_at(self.pos, c) || (c == b'!' && self.next_starts_bond_factor())
        });
        if (explicit || juxtaposed)
            && let Some(right) = self.try_parse_bond_factor()
        {
            let and_expr = BondQuery::And(Box::new(left), Box::new(right));
            return self.parse_bond_and_tail(and_expr);
        }
        left
    }

    // -- atoms ---------------------------------------------------------------

    /// Try to parse one atom (bracket, organic shorthand, or wildcard).
    /// Returns `(AtomQuery, atom_map)` — `atom_map` is `Some(n)` only for
    /// bracket atoms with an explicit `:N` suffix.
    fn try_parse_atom(&mut self) -> Result<Option<(AtomQuery, Option<u16>)>, SmartsError> {
        match self.peek() {
            Some(b'[') => {
                let (q, map) = self.parse_bracket_atom()?;
                Ok(Some((q, map)))
            }
            Some(b'*') => {
                self.advance();
                Ok(Some((AtomQuery::Primitive(AtomPrimitive::Wildcard), None)))
            }
            // `a` — any aromatic atom (standalone shorthand, equivalent to `[a]`).
            Some(b'a') => {
                self.advance();
                Ok(Some((
                    AtomQuery::Primitive(AtomPrimitive::Aromatic(true)),
                    None,
                )))
            }
            // `A` — any aliphatic atom (standalone shorthand, equivalent to `[A]`).
            Some(b'A') => {
                self.advance();
                Ok(Some((
                    AtomQuery::Primitive(AtomPrimitive::Aromatic(false)),
                    None,
                )))
            }
            // Organic-subset atoms (uppercase aliphatic).
            Some(b'B') | Some(b'C') | Some(b'N') | Some(b'O') | Some(b'P') | Some(b'S')
            | Some(b'F') | Some(b'I') => Ok(Some((self.parse_organic_atom(false)?, None))),
            // Aromatic organic-subset atoms (lowercase, includes boron `b`).
            Some(b'b') | Some(b'c') | Some(b'n') | Some(b'o') | Some(b'p') | Some(b's') => {
                Ok(Some((self.parse_organic_atom(true)?, None)))
            }
            _ => Ok(None),
        }
    }

    /// Parse an organic-subset atom shorthand such as `C` or `c`.
    ///
    /// Uppercase → `And(Symbol("X"), Aromatic(false))`.
    /// Lowercase → `And(Symbol("X"), Aromatic(true))`.
    fn parse_organic_atom(&mut self, aromatic: bool) -> Result<AtomQuery, SmartsError> {
        let pos = self.pos;
        let first = self.advance().unwrap() as char;

        // Handle two-character symbols: Cl, Br.
        let symbol = if !aromatic {
            if first == 'C' && self.peek() == Some(b'l') {
                self.advance();
                "Cl".to_string()
            } else if first == 'B' && self.peek() == Some(b'r') {
                self.advance();
                "Br".to_string()
            } else {
                first.to_string()
            }
        } else {
            first.to_ascii_uppercase().to_string()
        };

        // Verify the symbol is a real element.
        if chematic_core::Element::from_symbol(&symbol).is_none() {
            return Err(SmartsError::UnexpectedChar(first, pos));
        }

        let sym_query = AtomQuery::Primitive(AtomPrimitive::Symbol(symbol));
        let arom_query = AtomQuery::Primitive(AtomPrimitive::Aromatic(aromatic));
        Ok(AtomQuery::And(Box::new(sym_query), Box::new(arom_query)))
    }

    /// Parse a bracket atom `[expr]`, including an optional `:N` atom map number.
    /// Returns `(AtomQuery, atom_map)`.
    fn parse_bracket_atom(&mut self) -> Result<(AtomQuery, Option<u16>), SmartsError> {
        let bracket_pos = self.pos;
        self.advance(); // consume '['

        // Daylight/OpenSMARTS exception: a bracket atom whose whole expression
        // is `H` — optionally with an isotope before it and a charge and/or
        // atom map after it (`[H]`, `[2H]`, `[H+]`, `[H:1]`) — is a hydrogen
        // *atom*, not an H-count primitive; `[*H]`, `[H1]`, `[C;H]` keep the
        // count meaning. Issue #734: SMIRKS such as `[c:1][H]>>[c:1]O` rely on
        // this reading.
        let expr = if let Some(expr) = self.try_parse_hydrogen_atom() {
            expr
        } else {
            self.parse_expr()?
        };

        // Optional atom map number: `[O;D1;H0:3]` → atom_map = Some(3).
        // The `:` is metadata only and does not affect matching.
        // A bare `:` with no following digit (e.g. `[C:]`) is an error.
        let atom_map = if self.peek() == Some(b':') {
            let colon_pos = self.pos;
            self.advance(); // consume ':'
            if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                // Parse digits as u16, clamped at u16::MAX via checked arithmetic.
                let mut val: u32 = 0;
                while let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                    self.advance();
                    val = val.saturating_mul(10).saturating_add((d - b'0') as u32);
                }
                Some(val.min(u16::MAX as u32) as u16)
            } else {
                return Err(SmartsError::UnexpectedChar(':', colon_pos));
            }
        } else {
            None
        };

        match self.peek() {
            Some(b']') => {
                self.advance(); // consume ']'
            }
            _ => return Err(SmartsError::UnclosedBracket(bracket_pos)),
        }

        Ok((expr, atom_map))
    }

    /// `[H]`-form hydrogen atom (see `parse_bracket_atom`). Consumes the
    /// isotope, `H` and charge when the form matches, else consumes nothing.
    fn try_parse_hydrogen_atom(&mut self) -> Option<AtomQuery> {
        let mut pos = self.pos;
        let mut isotope: Option<u16> = None;
        let mut mass: u32 = 0;
        while let Some(d) = self.src.get(pos).filter(|b| b.is_ascii_digit()) {
            mass = mass.saturating_mul(10).saturating_add((d - b'0') as u32);
            isotope = Some(mass.min(u16::MAX as u32) as u16);
            pos += 1;
        }
        if self.src.get(pos) != Some(&b'H') {
            return None;
        }
        pos += 1;
        let mut charge: i8 = 0;
        if let Some(&sign @ (b'+' | b'-')) = self.src.get(pos) {
            pos += 1;
            let mut n: i8 = 1;
            if let Some(d) = self.src.get(pos).filter(|b| b.is_ascii_digit()) {
                n = (d - b'0') as i8;
                pos += 1;
            } else {
                while self.src.get(pos) == Some(&sign) && n < i8::MAX {
                    n += 1;
                    pos += 1;
                }
            }
            charge = if sign == b'+' { n } else { -n };
        }
        if !matches!(self.src.get(pos), Some(b']') | Some(b':')) {
            return None;
        }
        self.pos = pos;
        let mut expr = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("H".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
        );
        if let Some(mass) = isotope {
            expr = AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::Isotope(mass))),
                Box::new(expr),
            );
        }
        if charge != 0 {
            expr = AtomQuery::And(
                Box::new(expr),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(charge))),
            );
        }
        Some(expr)
    }

    // -- expression grammar (inside brackets) --------------------------------
    //
    // Precedence (lowest to highest):
    //   low_and  (';')
    //   or       (',')
    //   high_and ('&' or juxtaposition)
    //   unary    ('!')
    //   primitive

    fn parse_expr(&mut self) -> Result<AtomQuery, SmartsError> {
        self.parse_low_and()
    }

    /// Low-precedence AND: `or (';' or)*`
    fn parse_low_and(&mut self) -> Result<AtomQuery, SmartsError> {
        let mut left = self.parse_or()?;
        while self.peek() == Some(b';') {
            self.advance(); // consume ';'
            let right = self.parse_or()?;
            left = AtomQuery::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    /// OR: `high_and (',' high_and)*`
    fn parse_or(&mut self) -> Result<AtomQuery, SmartsError> {
        let mut left = self.parse_high_and()?;
        while self.peek() == Some(b',') {
            self.advance(); // consume ','
            let right = self.parse_high_and()?;
            left = AtomQuery::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    /// High-precedence AND: `unary ('&'? unary)*`
    fn parse_high_and(&mut self) -> Result<AtomQuery, SmartsError> {
        let mut left = self.parse_unary()?;
        loop {
            // Check whether the next character could start another unary expression,
            // either with an explicit '&' or by juxtaposition.
            let explicit_and = self.peek() == Some(b'&');
            if explicit_and {
                self.advance(); // consume '&'
            }

            // Check whether there is a next primitive/unary.
            if self.can_start_primitive() {
                let right = self.parse_unary()?;
                left = AtomQuery::And(Box::new(left), Box::new(right));
            } else if explicit_and {
                // '&' with no following operand is an error.
                return Err(SmartsError::UnexpectedEnd);
            } else {
                break;
            }
        }
        Ok(left)
    }

    /// Returns true if the current position can start a unary/primitive expression.
    fn can_start_primitive(&self) -> bool {
        match self.peek() {
            Some(b'!') => true,
            Some(b'#') => true,
            Some(b'a') | Some(b'A') => true,
            Some(b'+') | Some(b'-') => true,
            Some(b'H') => true,
            Some(b'D') => true,
            Some(b'r') => true,
            Some(b'R') => true,
            // Any-ring size `[kN]` (RDKit PR #9172) — distinct from `[rN]`'s
            // smallest-ring semantics, see the `[rN]`/`[kN]` split below.
            Some(b'k') => true,
            Some(b'*') => true,
            // Recursive SMARTS `$(...)`.
            Some(b'$') => true,
            // Valence `[vN]`, ring-bond count `[xN]`, hybridization `[^N]`.
            Some(b'v') | Some(b'x') | Some(b'^') => true,
            // Isotope mass number `[13C]`, `[2H]`.
            Some(c) if c.is_ascii_digit() => true,
            // Chirality `[@]`, `[@@]`.
            Some(b'@') => true,
            // Uppercase element symbol — check it's not a stop character.
            Some(c) if c.is_ascii_uppercase() => true,
            // Lowercase element symbol (but not 'a' already handled).
            Some(c) if c.is_ascii_lowercase() => true,
            _ => false,
        }
    }

    /// Unary: `'!' unary | primitive`
    fn parse_unary(&mut self) -> Result<AtomQuery, SmartsError> {
        if self.peek() == Some(b'!') {
            self.advance(); // consume '!'
            let inner = self.parse_unary()?;
            return Ok(AtomQuery::Not(Box::new(inner)));
        }
        self.parse_primitive()
    }

    fn parse_recursive_primitive(&mut self, pos: usize) -> Result<AtomQuery, SmartsError> {
        self.advance();
        if self.peek() != Some(b'(') {
            return Err(SmartsError::UnexpectedChar('$', pos));
        }
        self.advance();
        let start = self.pos;
        let mut depth = 1usize;
        let mut end = start;
        while end < self.src.len() {
            match self.src[end] {
                b'(' => depth += 1,
                b')' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                _ => {}
            }
            end += 1;
        }
        if depth != 0 {
            return Err(SmartsError::UnexpectedEnd);
        }
        let inner =
            std::str::from_utf8(&self.src[start..end]).map_err(|_| SmartsError::UnexpectedEnd)?;
        if self.recursion_depth >= self.max_recursion_depth {
            return Err(SmartsError::RecursionDepthExceeded);
        }
        let inner_mol = Parser {
            src: inner.as_bytes(),
            pos: 0,
            recursion_depth: self.recursion_depth + 1,
            max_recursion_depth: self.max_recursion_depth,
        }
        .parse()?;
        self.pos = end + 1;
        Ok(AtomQuery::Primitive(AtomPrimitive::Recursive(Box::new(
            inner_mol,
        ))))
    }

    fn parse_charge_primitive(&mut self, sign: u8) -> Result<AtomQuery, SmartsError> {
        self.advance();
        let direction = if sign == b'+' { 1 } else { -1 };
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(charge_range(lo, hi, direction));
        }
        let magnitude = self.parse_charge_magnitude(sign);
        Ok(AtomQuery::Primitive(AtomPrimitive::Charge(
            direction * magnitude,
        )))
    }

    fn element_follows(&self, first: char) -> bool {
        self.src
            .get(self.pos + 1)
            .copied()
            .filter(u8::is_ascii_lowercase)
            .is_some_and(|second| {
                chematic_core::Element::from_symbol(&format!("{first}{}", second as char)).is_some()
            })
    }

    fn parse_smallest_ring_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        self.advance();
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(range_query(lo, hi, min_ring_size));
        }
        Ok(match self.parse_ring_size_number() {
            Some(n) => AtomQuery::Primitive(min_ring_size(n)),
            None => AtomQuery::Primitive(AtomPrimitive::RingMembership(true)),
        })
    }

    fn parse_any_ring_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        self.advance();
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::RingMembership(true))),
                Box::new(range_query(lo, hi, AtomPrimitive::RingSize)),
            ));
        }
        Ok(match self.parse_count_number() {
            Some(0) => AtomQuery::Primitive(AtomPrimitive::Wildcard),
            Some(n) => AtomQuery::Primitive(AtomPrimitive::RingSize(n)),
            None => AtomQuery::Primitive(AtomPrimitive::RingMembership(true)),
        })
    }

    fn parse_ring_membership_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        if self.element_follows('R') {
            return self.parse_element_primitive();
        }
        self.advance();
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(range_query(lo, hi, AtomPrimitive::RingCount));
        }
        Ok(match self.parse_count_number() {
            Some(n) => AtomQuery::Primitive(AtomPrimitive::RingCount(n)),
            None => AtomQuery::Primitive(AtomPrimitive::RingMembership(true)),
        })
    }

    fn parse_isotope_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        let mut mass: u16 = 0;
        while let Some(d) = self.peek().filter(u8::is_ascii_digit) {
            self.advance();
            mass = mass
                .checked_mul(10)
                .and_then(|value| value.checked_add((d - b'0') as u16))
                .ok_or_else(|| SmartsError::UnexpectedChar(d as char, self.pos - 1))?;
        }
        Ok(AtomQuery::Primitive(AtomPrimitive::Isotope(mass)))
    }

    fn parse_chirality_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        self.advance();
        let kind = if self.peek() == Some(b'@') {
            self.advance();
            2u8
        } else {
            1u8
        };
        if kind == 1 {
            let rest = &self.src[self.pos..];
            if let Some(class) = [b"TH", b"AL", b"SP", b"TB", b"OH"]
                .into_iter()
                .find(|class| rest.starts_with(*class))
            {
                self.pos += 2;
                let n = self.parse_digits_u8().ok_or(SmartsError::UnexpectedEnd)?;
                let max = match class {
                    b"TH" | b"AL" => 2,
                    b"SP" => 3,
                    b"TB" => 20,
                    _ => 30,
                };
                if n == 0 || n > max {
                    return Err(SmartsError::UnexpectedChar(
                        (b'0' + n.min(9)) as char,
                        self.pos - 1,
                    ));
                }
                return Ok(AtomQuery::Primitive(AtomPrimitive::Chirality(
                    if class == b"TH" { n } else { 3 },
                )));
            }
        }
        Ok(AtomQuery::Primitive(AtomPrimitive::Chirality(kind)))
    }

    fn parse_hybridization_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        self.advance();
        let n = self
            .parse_single_digit()
            .ok_or(SmartsError::UnexpectedEnd)?;
        if n > 5 {
            return Err(SmartsError::UnexpectedChar(
                (b'0' + n) as char,
                self.pos - 1,
            ));
        }
        Ok(AtomQuery::Primitive(AtomPrimitive::Hybridization(n)))
    }

    /// Parse one primitive expression.
    fn parse_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        let pos = self.pos;
        match self.peek() {
            Some(b'*') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Wildcard))
            }
            Some(b'$') => self.parse_recursive_primitive(pos),
            // Keep aromatic/aliphatic primitives before element parsing.
            Some(b'a') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(true)))
            }
            Some(b'A') => {
                if self.element_follows('A') {
                    return self.parse_element_primitive();
                }
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(false)))
            }
            Some(b'#') => {
                self.advance();
                let n = self.parse_digits_u8().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)))
            }
            Some(b'+') => self.parse_charge_primitive(b'+'),
            Some(b'-') => self.parse_charge_primitive(b'-'),
            Some(b'H') => {
                if self.element_follows('H') {
                    return self.parse_element_primitive();
                }
                self.advance();
                self.count_primitive(Some(1), AtomPrimitive::HCount)
            }
            Some(b'h') => {
                self.advance();
                self.neighbor_count_primitive(AtomPrimitive::ImplicitHCount)
            }
            Some(b'D') => {
                if self.element_follows('D') {
                    return self.parse_element_primitive();
                }
                self.advance();
                self.count_primitive(Some(1), AtomPrimitive::Degree)
            }
            Some(b'r') => self.parse_smallest_ring_primitive(),
            Some(b'k') => self.parse_any_ring_primitive(),
            Some(b'R') => self.parse_ring_membership_primitive(),
            Some(b'v') => {
                self.advance();
                self.count_primitive(Some(1), AtomPrimitive::Valence)
            }
            Some(b'x') => {
                self.advance();
                self.neighbor_count_primitive(AtomPrimitive::RingBondCount)
            }
            Some(b'^') => self.parse_hybridization_primitive(),
            // `X` is connectivity unless it starts xenon.
            Some(b'X') if self.src.get(self.pos + 1) != Some(&b'e') => {
                self.advance();
                self.count_primitive(Some(1), AtomPrimitive::TotalConnectivity)
            }
            Some(b'z') => {
                self.advance();
                self.neighbor_count_primitive(AtomPrimitive::HeteroNeighborCount)
            }
            Some(b'd') => {
                self.advance();
                self.count_primitive(Some(1), AtomPrimitive::HeavyDegree)
            }
            Some(b'Z')
                if !self
                    .src
                    .get(self.pos + 1)
                    .is_some_and(|c| matches!(c, b'n' | b'r')) =>
            {
                self.advance();
                self.neighbor_count_primitive(AtomPrimitive::AliphaticHeteroNeighborCount)
            }
            Some(c) if c.is_ascii_digit() => self.parse_isotope_primitive(),
            Some(b'@') => self.parse_chirality_primitive(),
            Some(c) if c.is_ascii_alphabetic() => self.parse_element_primitive(),
            Some(c) => Err(SmartsError::UnexpectedChar(c as char, pos)),
            None => Err(SmartsError::UnexpectedEnd),
        }
    }

    /// Parse an element symbol as a primitive inside a bracket atom.
    ///
    /// Uppercase (e.g. `[O]`) → aliphatic atom of that element.
    /// Lowercase (e.g. `[o]`) → aromatic atom of that element (adds `Aromatic(true)` constraint).
    fn parse_element_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        let pos = self.pos;
        let first = self.advance().unwrap() as char;
        let aromatic = first.is_ascii_lowercase();
        let upper_first = first.to_ascii_uppercase();

        // Two-character symbols. An uppercase first letter takes the next
        // lowercase letter when they spell an element (`Cl`, `Na`). A
        // lowercase first letter spells only the aromatic `se`, `as`, `te`
        // and `si`, as in RDKit; otherwise it is the one-letter aromatic atom
        // and the next letter is the next primitive: `[ca]` is aromatic
        // carbon and aromatic (not calcium), `[nr6]` aromatic nitrogen in a
        // six-membered ring, `[cl]` an error.
        if let Some(second) = self.peek()
            && second.is_ascii_lowercase()
        {
            let candidate = format!("{upper_first}{}", second as char);
            if !aromatic && chematic_core::Element::from_symbol(&candidate).is_some() {
                self.advance();
                return Ok(AtomQuery::Primitive(AtomPrimitive::Symbol(candidate)));
            }
            if aromatic && matches!(candidate.as_str(), "Se" | "As" | "Te" | "Si") {
                self.advance();
                return Ok(AtomQuery::And(
                    Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol(candidate))),
                    Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
                ));
            }
        }

        // Single-character symbol. Lowercase (aromatic) is only `b c n o p s`,
        // as in Daylight and RDKit: `[i]`, `[f]` are errors, not atoms that
        // can never match.
        let sym = upper_first.to_string();
        if aromatic && !matches!(first, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
            return Err(SmartsError::UnexpectedChar(first, pos));
        }
        if chematic_core::Element::from_symbol(&sym).is_some() {
            let sym_q = AtomQuery::Primitive(AtomPrimitive::Symbol(sym));
            // Lowercase (e.g. `[o]`, `[n]`) → aromatic; uppercase (e.g. `[O]`, `[N]`) → aliphatic.
            Ok(AtomQuery::And(
                Box::new(sym_q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(aromatic))),
            ))
        } else {
            Err(SmartsError::UnexpectedChar(first, pos))
        }
    }

    // -- digit helpers -------------------------------------------------------

    /// Parse a sequence of ASCII digits and return as `u8` (or `None` if no digits).
    fn parse_digits_u8(&mut self) -> Option<u8> {
        if !self.peek().map(|c| c.is_ascii_digit()).unwrap_or(false) {
            return None;
        }
        let mut val: u16 = 0;
        while let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
            self.advance();
            val = val.saturating_mul(10).saturating_add((d - b'0') as u16);
        }
        if val > 255 { None } else { Some(val as u8) }
    }

    /// Parse a single ASCII digit, if present. Returns the digit value 0–9.
    /// Parse a 1-3 digit ring size after `r`, or `None` if no digit follows.
    fn parse_ring_size_number(&mut self) -> Option<u8> {
        let mut value: u32 = 0;
        let mut digits = 0;
        while digits < 3 {
            match self.peek() {
                Some(c @ b'0'..=b'9') => {
                    value = value * 10 + u32::from(c - b'0');
                    self.advance();
                    digits += 1;
                }
                _ => break,
            }
        }
        (digits > 0).then(|| value.min(u32::from(u8::MAX)) as u8)
    }

    /// A multi-digit count (RDKit reads `[D12]`, `[H12]`, `[+10]`), saturating
    /// at 255: no atom reaches that count, so larger values match nothing.
    fn parse_count_number(&mut self) -> Option<u8> {
        let mut value: Option<u16> = None;
        while let Some(d) = self.peek().filter(u8::is_ascii_digit) {
            self.advance();
            let v = value.unwrap_or(0);
            value = Some(
                v.saturating_mul(10)
                    .saturating_add(u16::from(d - b'0'))
                    .min(255),
            );
        }
        value.map(|v| v as u8)
    }

    fn parse_single_digit(&mut self) -> Option<u8> {
        match self.peek() {
            Some(d) if d.is_ascii_digit() => {
                self.advance();
                Some(d - b'0')
            }
            _ => None,
        }
    }

    /// An RDKit range `{a-b}`, `{a-}` or `{-b}` at the current position, if
    /// any: `(a, b)` with `b = None` for an open upper bound.
    fn parse_range(&mut self) -> Result<Option<(u8, Option<u8>)>, SmartsError> {
        if self.peek() != Some(b'{') {
            return Ok(None);
        }
        let start = self.pos;
        self.advance(); // consume '{'
        let lo = self.parse_count_number();
        if self.peek() != Some(b'-') {
            return Err(SmartsError::UnexpectedChar(
                self.peek().map_or('}', char::from),
                self.pos,
            ));
        }
        self.advance(); // consume '-'
        let hi = self.parse_count_number();
        if self.peek() != Some(b'}') || (lo.is_none() && hi.is_none()) {
            return Err(SmartsError::UnexpectedChar('{', start));
        }
        self.advance(); // consume '}'
        Ok(Some((lo.unwrap_or(0), hi)))
    }

    /// A count primitive: a range, a number (multi-digit), or the bare form
    /// (`default`; an error when there is none).
    fn count_primitive(
        &mut self,
        default: Option<u8>,
        make: fn(u8) -> AtomPrimitive,
    ) -> Result<AtomQuery, SmartsError> {
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(range_query(lo, hi, make));
        }
        let n = match self.parse_count_number() {
            Some(n) => n,
            None => default.ok_or(SmartsError::UnexpectedEnd)?,
        };
        Ok(AtomQuery::Primitive(make(n)))
    }

    /// `z`/`Z`: a range, a number, or bare (at least one neighbour).
    fn neighbor_count_primitive(
        &mut self,
        make: fn(u8) -> AtomPrimitive,
    ) -> Result<AtomQuery, SmartsError> {
        if let Some((lo, hi)) = self.parse_range()? {
            return Ok(range_query(lo, hi, make));
        }
        Ok(match self.parse_count_number() {
            Some(n) => AtomQuery::Primitive(make(n)),
            None => AtomQuery::Not(Box::new(AtomQuery::Primitive(make(0)))),
        })
    }

    /// Magnitude of a bracket charge whose first `sign` has just been consumed:
    /// a following digit gives the magnitude (`+2`, `+0`), otherwise one plus
    /// the number of immediately repeated signs (`+` = 1, `++` = 2, `+++` = 3).
    fn parse_charge_magnitude(&mut self, sign: u8) -> i8 {
        if let Some(n) = self.parse_count_number() {
            return n.min(i8::MAX as u8) as i8;
        }
        let mut n: i8 = 1;
        while self.peek() == Some(sign) && n < i8::MAX {
            self.advance();
            n += 1;
        }
        n
    }
}

/// Largest charge magnitude an open charge range (`+{1-}`) enumerates.
const MAX_RANGE_CHARGE: u8 = 15;

/// A charge range `+{lo-hi}` (`sign` 1) or `-{lo-hi}` (`sign` -1); bounds
/// beyond the `i8` charge range cannot match any atom.
fn charge_range(lo: u8, hi: Option<u8>, sign: i8) -> AtomQuery {
    let max = i8::MAX as u8;
    let hi = hi.unwrap_or(MAX_RANGE_CHARGE.max(lo)).min(max);
    if lo > max {
        return AtomQuery::Not(Box::new(AtomQuery::Primitive(AtomPrimitive::Wildcard)));
    }
    range_query(lo, Some(hi), |c| AtomPrimitive::Charge(sign * c as i8))
}

/// `rN`: smallest ring of size N; `r0` is an acyclic atom (RDKit).
fn min_ring_size(n: u8) -> AtomPrimitive {
    if n == 0 {
        AtomPrimitive::RingMembership(false)
    } else {
        AtomPrimitive::MinRingSize(n)
    }
}

/// Desugar an RDKit range `{lo-hi}` over count primitive `make` into existing
/// primitives: `make(lo) OR … OR make(hi)`; an open upper bound is
/// `NOT (make(0) OR … OR make(lo - 1))`; an empty range matches nothing.
fn range_query(lo: u8, hi: Option<u8>, make: impl Fn(u8) -> AtomPrimitive) -> AtomQuery {
    let any_of = |from: u8, to: u8| -> AtomQuery {
        let mut q = AtomQuery::Primitive(make(from));
        for v in from.saturating_add(1)..=to {
            q = AtomQuery::Or(Box::new(q), Box::new(AtomQuery::Primitive(make(v))));
        }
        q
    };
    let nothing = || AtomQuery::Not(Box::new(AtomQuery::Primitive(AtomPrimitive::Wildcard)));
    match hi {
        Some(hi) if lo > hi => nothing(),
        Some(hi) => any_of(lo, hi),
        None if lo == 0 => AtomQuery::Primitive(AtomPrimitive::Wildcard),
        None => AtomQuery::Not(Box::new(any_of(0, lo - 1))),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;
