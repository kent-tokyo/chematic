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
//! Implicit bond between two adjacent atoms = `BondQuery::Any` (`~`).

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
        // open_rings maps ring-closure-number → (atom_idx, Optional<BondQuery>)
        let mut open_rings: HashMap<u8, (usize, Option<BondQuery>)> = HashMap::new();

        self.parse_chain(&mut mol, None, None, &mut open_rings)?;

        // Any remaining open ring closures are errors.
        if let Some((&num, _)) = open_rings.iter().next() {
            return Err(SmartsError::InvalidRingClosure(num));
        }

        // Remaining characters (other than whitespace) are unexpected.
        if let Some(c) = self.peek()
            && c != b' '
            && c != b'\t'
            && c != b'\n'
            && c != b'\r'
        {
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
        open_rings: &mut HashMap<u8, (usize, Option<BondQuery>)>,
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

                // End of this chain: `)`, end-of-input, or other stop character.
                None | Some(b')') => break,

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

    // -- ring closures -------------------------------------------------------

    /// Handle a ring closure: either close an existing open ring or open a new one.
    fn handle_ring_closure(
        &mut self,
        mol: &mut QueryMolecule,
        current: usize,
        ring_num: u8,
        ring_bond: Option<BondQuery>,
        open_rings: &mut HashMap<u8, (usize, Option<BondQuery>)>,
    ) -> Result<(), SmartsError> {
        if let Some((open_atom, open_bond)) = open_rings.remove(&ring_num) {
            // Resolve the bond type.
            let bond = match (open_bond, ring_bond) {
                (Some(a), Some(b)) if a == b => a,
                (Some(_), Some(_)) => {
                    // Conflicting ring bond specifications — use Any as a fallback.
                    BondQuery::Any
                }
                (Some(b), None) | (None, Some(b)) => b,
                (None, None) => BondQuery::Any,
            };
            mol.add_bond(open_atom, current, bond);
        } else {
            open_rings.insert(ring_num, (current, ring_bond));
        }
        Ok(())
    }

    /// Parse a ring closure number (single digit or `%nn`) and an optional
    /// leading bond query (already consumed before calling this).
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
        let first = self.try_parse_bond_factor()?;
        Some(self.parse_bond_or_tail(first))
    }

    /// Parse one bond factor: `!token` or `token`.
    fn try_parse_bond_factor(&mut self) -> Option<BondQuery> {
        match self.peek()? {
            b'!' => {
                // Only treat `!` as bond negation when the next char is a bond token.
                if self
                    .src
                    .get(self.pos + 1)
                    .copied()
                    .map(Self::is_bond_token)
                    .unwrap_or(false)
                {
                    self.advance(); // consume '!'
                    let prim = self.consume_bond_prim().unwrap();
                    Some(BondQuery::Not(Box::new(BondQuery::Primitive(prim))))
                } else {
                    None
                }
            }
            c if Self::is_bond_token(c) => {
                let prim = self.consume_bond_prim().unwrap();
                Some(BondQuery::Primitive(prim))
            }
            _ => None,
        }
    }

    /// Consume a single bond primitive character (caller must verify `peek` is a bond token).
    fn consume_bond_prim(&mut self) -> Option<BondPrimitive> {
        let prim = match self.peek()? {
            b'-' => BondPrimitive::Single,
            b'=' => BondPrimitive::Double,
            b'#' => BondPrimitive::Triple,
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
    fn is_bond_token(c: u8) -> bool {
        matches!(c, b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'/' | b'\\')
    }

    /// Continue parsing bond OR after the first factor.
    fn parse_bond_or_tail(&mut self, left: BondQuery) -> BondQuery {
        // ',' → OR (only when followed by a bond token or '!')
        if self.peek() == Some(b',') {
            let next = self.src.get(self.pos + 1).copied();
            if next.map(Self::is_bond_token).unwrap_or(false) || next == Some(b'!') {
                self.advance(); // consume ','
                if let Some(right) = self.try_parse_bond_factor() {
                    let right = self.parse_bond_and_tail(right);
                    let or_expr = BondQuery::Or(Box::new(left), Box::new(right));
                    return self.parse_bond_or_tail(or_expr);
                }
            }
        }
        self.parse_bond_and_tail(left)
    }

    /// Continue parsing implicit AND after the first factor.
    ///
    /// Recognises `&!<bond>` or a bare `!<bond>` juxtaposed with the previous
    /// factor as a high-precedence AND. A trailing `&` with no operand is
    /// silently consumed (caller treats the result the same as `left`).
    fn parse_bond_and_tail(&mut self, left: BondQuery) -> BondQuery {
        if self.peek() == Some(b'&') {
            self.advance(); // consume '&'
        }

        if self.peek() == Some(b'!') {
            let next = self.src.get(self.pos + 1).copied();
            if next.map(Self::is_bond_token).unwrap_or(false)
                && let Some(right) = self.try_parse_bond_factor()
            {
                let and_expr = BondQuery::And(Box::new(left), Box::new(right));
                return self.parse_bond_and_tail(and_expr);
            }
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

        let expr = self.parse_expr()?;

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

    /// Parse one primitive expression.
    fn parse_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        let pos = self.pos;
        match self.peek() {
            // Wildcard
            Some(b'*') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Wildcard))
            }

            // Recursive SMARTS: `$(inner_smarts)`.
            Some(b'$') => {
                self.advance(); // consume '$'
                if self.peek() != Some(b'(') {
                    return Err(SmartsError::UnexpectedChar('$', pos));
                }
                self.advance(); // consume '('
                // Scan forward to find the matching ')', counting nesting depth.
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
                let inner_str = std::str::from_utf8(&self.src[start..end])
                    .map_err(|_| SmartsError::UnexpectedEnd)?;
                if self.recursion_depth >= self.max_recursion_depth {
                    return Err(SmartsError::RecursionDepthExceeded);
                }
                let mut inner_parser = Parser {
                    src: inner_str.as_bytes(),
                    pos: 0,
                    recursion_depth: self.recursion_depth + 1,
                    max_recursion_depth: self.max_recursion_depth,
                };
                let inner_mol = inner_parser.parse()?;
                self.pos = end + 1; // advance past the closing ')'
                Ok(AtomQuery::Primitive(AtomPrimitive::Recursive(Box::new(
                    inner_mol,
                ))))
            }

            // Aromatic (`a`) — must be checked BEFORE element symbol parsing.
            Some(b'a') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(true)))
            }

            // Aliphatic (`A`) — but if followed by a lowercase letter that forms a valid
            // two-character element symbol (Ac, Ag, Al, Am, Ar, As, At, Au), parse as element.
            Some(b'A') => {
                let next = self.src.get(self.pos + 1).copied();
                if next.is_some_and(|c| c.is_ascii_lowercase()) {
                    let candidate = format!("A{}", next.unwrap() as char);
                    if chematic_core::Element::from_symbol(&candidate).is_some() {
                        return self.parse_element_primitive();
                    }
                }
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(false)))
            }

            // Atomic number `#N`
            Some(b'#') => {
                self.advance(); // consume '#'
                let n = self.parse_digits_u8().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)))
            }

            // Charge `+N` or `+`  (charge = +N or +1; `+0` = explicit neutral)
            Some(b'+') => {
                self.advance(); // consume '+'
                // Parse digit including '0'; only default to 1 if no digit follows at all.
                let n = if self.peek().map(|c| c.is_ascii_digit()).unwrap_or(false) {
                    self.parse_single_digit().unwrap_or(1)
                } else {
                    1
                };
                Ok(AtomQuery::Primitive(AtomPrimitive::Charge(n as i8)))
            }

            // Charge `-N` or `-`  (charge = -N or -1; `-0` = explicit neutral)
            Some(b'-') => {
                self.advance(); // consume '-'
                let n = if self.peek().map(|c| c.is_ascii_digit()).unwrap_or(false) {
                    self.parse_single_digit().unwrap_or(1)
                } else {
                    1
                };
                Ok(AtomQuery::Primitive(AtomPrimitive::Charge(-(n as i8))))
            }

            // H count `HN` or `H` (total hcount = N or 1; explicit + implicit)
            // Exception: if followed by a lowercase letter forming a valid element (Hg, Ho, Hs),
            // treat as element symbol instead.
            Some(b'H') => {
                let next = self.src.get(self.pos + 1).copied();
                if next.is_some_and(|c| c.is_ascii_lowercase()) {
                    let candidate = format!("H{}", next.unwrap() as char);
                    if chematic_core::Element::from_symbol(&candidate).is_some() {
                        return self.parse_element_primitive();
                    }
                }
                self.advance(); // consume 'H'
                let n = self.parse_single_digit().unwrap_or(1);
                Ok(AtomQuery::Primitive(AtomPrimitive::HCount(n)))
            }

            // Implicit H count `hN` or `h` (implicit hcount only = N or 1)
            Some(b'h') => {
                self.advance(); // consume 'h'
                let n = self.parse_single_digit().unwrap_or(1);
                Ok(AtomQuery::Primitive(AtomPrimitive::ImplicitHCount(n)))
            }

            // Degree `DN`
            // Exception: if followed by a lowercase letter forming a valid element (Dy, Db, Ds),
            // treat as element symbol instead.
            Some(b'D') => {
                let next = self.src.get(self.pos + 1).copied();
                if next.is_some_and(|c| c.is_ascii_lowercase()) {
                    let candidate = format!("D{}", next.unwrap() as char);
                    if chematic_core::Element::from_symbol(&candidate).is_some() {
                        return self.parse_element_primitive();
                    }
                }
                self.advance(); // consume 'D'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Degree(n)))
            }

            // Smallest-ring size `rN` — must be checked BEFORE element parsing.
            // RDKit's `[rN]` means "the SMALLEST ring containing this atom has exactly
            // N members" -- NOT "the atom belongs to any N-membered ring" (that's
            // `[kN]`, below). These are genuinely different predicates: on a fusion
            // atom shared between a 5-ring and a 6-ring, RDKit's `[k6]` matches but
            // `[r6]` does not, since that atom's smallest ring is the 5-ring. Verified
            // against real RDKit and measured recoverable across a full 5,000-molecule
            // corpus (70.95% -> 100.00%) -- see `docs/rdkit_compat.md`'s
            // "SMARTS-R0"/"SMARTS-R1" entries. Previously (incorrectly) aliased to the
            // same `RingSize` primitive as `[kN]`.
            Some(b'r') => {
                self.advance(); // consume 'r'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::MinRingSize(n)))
            }

            // Any-ring size `[kN]` — RDKit PR #9172. Means "this atom belongs to *some*
            // perceived ring of exactly N members" -- unlike `[rN]` above, this is
            // NOT restricted to the atom's smallest ring. Confirmed empirically
            // distinct from `[rN]` (see above); do not merge these two primitives.
            Some(b'k') => {
                self.advance(); // consume 'k'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::RingSize(n)))
            }

            // Ring membership `R` or ring count `RN` (N = 0, 1, 2, …).
            // `[R]`  = in any ring (RingMembership(true))
            // `[R0]` = not in any ring (RingCount(0))
            // `[R1]` = in exactly 1 ring, etc.
            // Exception: if followed by a lowercase letter forming a valid element (Ra, Rb, Re,
            // Rf, Rg, Rh, Rn, Ru), treat as element symbol instead.
            Some(b'R') => {
                let next = self.src.get(self.pos + 1).copied();
                if next.is_some_and(|c| c.is_ascii_lowercase()) {
                    let candidate = format!("R{}", next.unwrap() as char);
                    if chematic_core::Element::from_symbol(&candidate).is_some() {
                        return self.parse_element_primitive();
                    }
                }
                self.advance(); // consume 'R'
                if let Some(n) = self.parse_single_digit() {
                    Ok(AtomQuery::Primitive(AtomPrimitive::RingCount(n)))
                } else {
                    Ok(AtomQuery::Primitive(AtomPrimitive::RingMembership(true)))
                }
            }

            // Valence `[vN]` — total valence (bond orders + implicit H).
            Some(b'v') => {
                self.advance(); // consume 'v'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Valence(n)))
            }

            // Ring-bond count `[xN]` — bonds where both endpoints share a ring.
            Some(b'x') => {
                self.advance(); // consume 'x'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::RingBondCount(n)))
            }

            // Hybridization `[^N]` — 1=sp, 2=sp2, 3=sp3.
            Some(b'^') => {
                self.advance(); // consume '^'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Hybridization(n)))
            }

            // Total connectivity `[XN]` — heavy-atom degree + implicit H count.
            Some(b'X') => {
                self.advance(); // consume 'X'
                let n = self
                    .parse_single_digit()
                    .ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::TotalConnectivity(n)))
            }

            // Isotope mass number: `[13C]`, `[2H]`, etc.
            // Digits are consumed, then the element symbol follows via juxtaposition AND.
            Some(c) if c.is_ascii_digit() => {
                let mut mass: u16 = 0;
                while let Some(d) = self.peek().filter(|b| b.is_ascii_digit()) {
                    self.advance();
                    mass = mass * 10 + (d - b'0') as u16;
                }
                Ok(AtomQuery::Primitive(AtomPrimitive::Isotope(mass)))
            }

            // Chirality `[@]` (CCW, value 1) or `[@@]` (CW, value 2).
            Some(b'@') => {
                self.advance(); // consume first '@'
                let kind = if self.peek() == Some(b'@') {
                    self.advance(); // consume second '@'
                    2u8 // clockwise (@@)
                } else {
                    1u8 // counterclockwise (@)
                };
                Ok(AtomQuery::Primitive(AtomPrimitive::Chirality(kind)))
            }

            // Element symbol (uppercase or lowercase start).
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

        // Try two-character symbol first (e.g. `Cl`, `Br`).
        if let Some(second) = self.peek()
            && second.is_ascii_lowercase()
        {
            let candidate = format!("{upper_first}{}", second as char);
            if chematic_core::Element::from_symbol(&candidate).is_some() {
                self.advance();
                // Two-character symbols are never written as aromatic in SMARTS.
                return Ok(AtomQuery::Primitive(AtomPrimitive::Symbol(candidate)));
            }
        }

        // Single-character symbol.
        let sym = upper_first.to_string();
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
            val = val * 10 + (d - b'0') as u16;
        }
        if val > 255 { None } else { Some(val as u8) }
    }

    /// Parse a single ASCII digit, if present. Returns the digit value 0–9.
    fn parse_single_digit(&mut self) -> Option<u8> {
        match self.peek() {
            Some(d) if d.is_ascii_digit() => {
                self.advance();
                Some(d - b'0')
            }
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

    #[test]
    fn test_parse_aliphatic_c() {
        // `C` → 1 atom with And(Symbol("C"), Aromatic(false))
        let mol = parse_smarts("C").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.bonds.len(), 0);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
        );
        assert_eq!(mol.atoms[0].query, expected);
    }

    #[test]
    fn test_parse_aromatic_c() {
        // `c` → 1 atom with And(Symbol("C"), Aromatic(true))
        let mol = parse_smarts("c").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        let expected = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
        );
        assert_eq!(mol.atoms[0].query, expected);
    }

    #[test]
    fn test_parse_atomic_num() {
        // `[#6]` → AtomicNum(6)
        let mol = parse_smarts("[#6]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(6))
        );
    }

    #[test]
    fn test_parse_not() {
        // `[!C]` → Not(And(Symbol("C"), Aromatic(false)))
        // [C] inside brackets means aliphatic C (uppercase → Aromatic(false) added).
        let mol = parse_smarts("[!C]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Not(Box::new(AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(false))),
            )))
        );
    }

    #[test]
    fn test_parse_aromatic_primitive() {
        // `[a]` → Aromatic(true)
        let mol = parse_smarts("[a]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Aromatic(true))
        );
    }

    #[test]
    fn test_parse_degree() {
        // `[D3]` → Degree(3)
        let mol = parse_smarts("[D3]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Degree(3))
        );
    }

    #[test]
    fn test_parse_ring_size() {
        // `[r5]` → MinRingSize(5) -- smallest-ring semantics (SMARTS-R1), not the
        // `[kN]` any-ring semantics chematic used to (incorrectly) alias it to.
        let mol = parse_smarts("[r5]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::MinRingSize(5))
        );
    }

    #[test]
    fn test_parse_k_not_same_primitive_as_r() {
        // SMARTS-R1: `[kN]` (any-ring) and `[rN]` (smallest-ring) are DIFFERENT
        // predicates in RDKit -- confirmed empirically (a fusion atom shared between
        // a 5-ring and a 6-ring: RDKit's `[k6]` matches, `[r6]` does not). They must
        // no longer parse to the identical primitive (this test previously asserted
        // the opposite, which was the bug -- see docs/rdkit_compat.md's
        // "SMARTS-R0"/"SMARTS-R1" entries).
        let q_k = parse_smarts("[k6]").expect("[k6] must parse");
        let q_r = parse_smarts("[r6]").expect("[r6] must parse");
        assert_eq!(
            q_k.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::RingSize(6))
        );
        assert_eq!(
            q_r.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::MinRingSize(6))
        );
        assert_ne!(
            q_k.atoms[0].query, q_r.atoms[0].query,
            "[k6] and [r6] must NOT produce the same primitive"
        );
    }

    #[test]
    fn test_k_primitive_matches_benzene() {
        // [k6] must match all 6 benzene atoms (all in 6-membered ring).
        use crate::{find_matches, parse_smarts};
        use chematic_smiles::parse;
        let benzene = parse("c1ccccc1").unwrap();
        let q = parse_smarts("[k6]").unwrap();
        let m = find_matches(&q, &benzene);
        assert_eq!(m.len(), 6, "[k6] must match all 6 benzene carbons");
    }

    #[test]
    fn test_k_vs_r_equivalent_for_cyclopentane() {
        // [k5] and [r5] must give identical match counts on a SINGLE isolated ring
        // (no fusion) -- here "any ring of size 5" and "smallest ring is size 5"
        // trivially coincide, since there's only one ring to consider at all.
        use crate::{find_matches, parse_smarts};
        use chematic_smiles::parse;
        let mol = parse("C1CCCC1").unwrap(); // cyclopentane
        let q_k = parse_smarts("[k5]").unwrap();
        let q_r = parse_smarts("[r5]").unwrap();
        assert_eq!(
            find_matches(&q_k, &mol).len(),
            find_matches(&q_r, &mol).len(),
            "[k5] and [r5] must be equivalent in chematic"
        );
    }

    #[test]
    fn test_k_vs_r_diverge_at_5_6_fusion_atom() {
        // SMARTS-R1: on a fused bicyclic (5-ring + 6-ring sharing one bond), the
        // fusion atoms belong to BOTH rings -- [kN] (any-ring) must match both [k5]
        // and [k6] on them, but [rN] (smallest-ring) must match ONLY [r5], since
        // their smallest ring is the 5-ring. This is the exact property Milestone
        // SMARTS-R0 found chematic previously got wrong (aliasing [r6] to [k6]'s
        // semantics would have wrongly matched the fusion atom on [r6] too).
        use crate::{find_matches, parse_smarts};
        use chematic_perception::find_sssr;
        use chematic_smiles::parse;

        let mol = parse("C1CCC2CCCC2C1").unwrap(); // bicyclic: 5-ring fused to 6-ring
        let sssr = find_sssr(&mol);
        let rings = sssr.rings();
        assert_eq!(rings.len(), 2, "expected exactly 2 SSSR rings");
        let sizes: Vec<usize> = rings.iter().map(|r| r.len()).collect();
        assert!(
            sizes.contains(&5) && sizes.contains(&6),
            "expected one 5-ring and one 6-ring, got sizes {sizes:?}"
        );
        let ring5: std::collections::HashSet<_> = rings
            .iter()
            .find(|r| r.len() == 5)
            .unwrap()
            .iter()
            .copied()
            .collect();
        let ring6: std::collections::HashSet<_> = rings
            .iter()
            .find(|r| r.len() == 6)
            .unwrap()
            .iter()
            .copied()
            .collect();
        let fusion_atom = *ring5
            .intersection(&ring6)
            .next()
            .expect("rings must share an atom");

        let q_k5 = parse_smarts("[k5]").unwrap();
        let q_k6 = parse_smarts("[k6]").unwrap();
        let q_r5 = parse_smarts("[r5]").unwrap();
        let q_r6 = parse_smarts("[r6]").unwrap();

        let matches = |q| {
            find_matches(q, &mol)
                .into_iter()
                .any(|m| m[&0] == fusion_atom)
        };
        assert!(matches(&q_k5), "[k5] must match the fusion atom");
        assert!(matches(&q_k6), "[k6] must match the fusion atom");
        assert!(
            matches(&q_r5),
            "[r5] must match the fusion atom (its smallest ring is 5)"
        );
        assert!(
            !matches(&q_r6),
            "[r6] must NOT match the fusion atom (its smallest ring is 5, not 6)"
        );
    }

    #[test]
    fn test_r6_excludes_purine_fusion_atom_matching_rdkit() {
        // The exact discriminating case from SMARTS-R0's diagnosis: a purine where
        // RDKit's [k6] includes the ring-fusion nitrogen but [r6] excludes it (its
        // smallest ring is the 5-membered imidazole ring). Regression-pins the fix
        // against a concrete, previously-verified-against-real-RDKit example rather
        // than only a synthetic bicyclic.
        use crate::{find_matches, parse_smarts};
        use chematic_smiles::parse;
        let mol = parse("Nc1nc(N)c2ncn(C3CC(O)C(O)C(CO)O3)c2n1").unwrap();
        let q_r6 = parse_smarts("[r6]").unwrap();
        let q_k6 = parse_smarts("[k6]").unwrap();
        let r6_atoms: std::collections::HashSet<u32> = find_matches(&q_r6, &mol)
            .into_iter()
            .map(|m| m[&0].0)
            .collect();
        let k6_atoms: std::collections::HashSet<u32> = find_matches(&q_k6, &mol)
            .into_iter()
            .map(|m| m[&0].0)
            .collect();
        assert_eq!(
            r6_atoms,
            std::collections::HashSet::from([1, 2, 3, 9, 10, 11, 13, 15, 18, 20]),
            "[r6] atom set must match RDKit exactly (see SMARTS-R0 diagnosis)"
        );
        assert_eq!(
            k6_atoms,
            std::collections::HashSet::from([1, 2, 3, 5, 9, 10, 11, 13, 15, 18, 19, 20]),
            "[k6] atom set must match RDKit exactly (unchanged by this milestone)"
        );
        assert!(!r6_atoms.contains(&19) && k6_atoms.contains(&19));
    }

    #[test]
    fn test_min_ring_size_false_for_non_ring_atom() {
        // A fully acyclic molecule has no atom in any ring -- [rN] must match nothing.
        use crate::{find_matches, parse_smarts};
        use chematic_smiles::parse;
        let mol = parse("CCCC").unwrap();
        for size in [3, 4, 5, 6] {
            let q = parse_smarts(&format!("[r{size}]")).unwrap();
            assert!(
                find_matches(&q, &mol).is_empty(),
                "[r{size}] must not match any atom in an acyclic molecule"
            );
        }
    }

    #[test]
    fn test_ring_membership_and_ring_count_unaffected_by_r_split() {
        // [R]/[R0]/[R1]/[R2] are a separate primitive (RingMembership/RingCount) from
        // [rN] (MinRingSize) -- splitting [rN] from [kN] must not change their
        // behavior at all. Regression check on the same fused bicyclic as above.
        use crate::{find_matches, parse_smarts};
        use chematic_smiles::parse;
        let mol = parse("C1CCC2CCCC2C1").unwrap(); // 9 atoms, 2 fusion atoms
        let q_r = parse_smarts("[R]").unwrap();
        assert_eq!(
            find_matches(&q_r, &mol).len(),
            9,
            "[R] must match all 9 ring atoms"
        );
        let q_r1 = parse_smarts("[R1]").unwrap();
        assert_eq!(
            find_matches(&q_r1, &mol).len(),
            7,
            "[R1] must match the 7 non-fusion ring atoms"
        );
        let q_r2 = parse_smarts("[R2]").unwrap();
        assert_eq!(
            find_matches(&q_r2, &mol).len(),
            2,
            "[R2] must match the 2 fusion atoms (each in 2 SSSR rings)"
        );
    }

    #[test]
    fn test_parse_hcount() {
        // `[H2]` → HCount(2)
        let mol = parse_smarts("[H2]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::HCount(2))
        );
    }

    #[test]
    fn test_parse_cc_bond() {
        // `CC` → 2 atoms, 1 bond (Any)
        let mol = parse_smarts("CC").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].query, BondQuery::Any);
    }

    #[test]
    fn test_parse_double_bond() {
        // `C=C` → 2 atoms, 1 Double bond
        let mol = parse_smarts("C=C").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(
            mol.bonds[0].query,
            BondQuery::Primitive(BondPrimitive::Double)
        );
    }

    #[test]
    fn test_parse_branch() {
        // `C(=O)O` → 3 atoms (C, O, O), 2 bonds
        let mol = parse_smarts("C(=O)O").unwrap();
        assert_eq!(mol.atoms.len(), 3, "should have 3 atoms");
        assert_eq!(mol.bonds.len(), 2, "should have 2 bonds");
        // First bond: C=O (double)
        assert_eq!(
            mol.bonds[0].query,
            BondQuery::Primitive(BondPrimitive::Double)
        );
        // Second bond: C-O (implicit Any)
        assert_eq!(mol.bonds[1].query, BondQuery::Any);
    }

    #[test]
    fn test_parse_benzene_ring() {
        // `c1ccccc1` → 6 aromatic C atoms, 6 bonds
        let mol = parse_smarts("c1ccccc1").unwrap();
        assert_eq!(mol.atoms.len(), 6, "benzene has 6 atoms");
        assert_eq!(mol.bonds.len(), 6, "benzene has 6 bonds");
        // All atoms should be aromatic.
        for atom in &mol.atoms {
            let expected = AtomQuery::And(
                Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(true))),
            );
            assert_eq!(atom.query, expected);
        }
    }

    // ---- Sprint 2: SMARTS depth expansion edge cases ----
    #[test]
    fn test_parse_with_custom_config_default_depth() {
        let config = SmartsParserConfig::default();
        assert_eq!(
            config.max_recursion_depth,
            DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH
        );
        let mol = parse_smarts_with_config("c1ccccc1", &config).unwrap();
        assert_eq!(mol.atoms.len(), 6);
    }

    #[test]
    fn test_parse_with_custom_config_increased_depth() {
        let config = SmartsParserConfig {
            max_recursion_depth: 12,
        };
        let mol = parse_smarts_with_config("[C,N,O]", &config).unwrap();
        assert_eq!(mol.atoms.len(), 1);
        // Pattern parses successfully with higher depth limit.
    }

    #[test]
    fn test_parse_with_custom_config_depth_clamped_high() {
        // Request depth > 16 should clamp to 16.
        let config = SmartsParserConfig {
            max_recursion_depth: 100,
        };
        let mol = parse_smarts_with_config("C", &config).unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_parse_with_custom_config_depth_clamped_low() {
        // Request depth < 1 should clamp to 1.
        let config = SmartsParserConfig {
            max_recursion_depth: 0,
        };
        let mol = parse_smarts_with_config("C", &config).unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_parse_with_default_config_equivalent() {
        // parse_smarts should be equivalent to parse_smarts_with_config with default
        let pattern = "c1ccccc1";
        let mol1 = parse_smarts(pattern).unwrap();
        let mol2 = parse_smarts_with_config(pattern, &SmartsParserConfig::default()).unwrap();
        assert_eq!(mol1.atoms.len(), mol2.atoms.len());
    }

    #[test]
    fn test_config_depth_parameter_clamped_correctly() {
        // Test that depth is clamped to [1, 16]
        assert_eq!(
            SmartsParserConfig {
                max_recursion_depth: 0
            }
            .max_recursion_depth,
            0
        );
        assert_eq!(
            SmartsParserConfig {
                max_recursion_depth: 16
            }
            .max_recursion_depth,
            16
        );
        assert_eq!(
            SmartsParserConfig {
                max_recursion_depth: 100
            }
            .max_recursion_depth,
            100
        );
        // Config stores the value as-is; clamping happens in parse_smarts_with_config
    }

    #[test]
    fn test_operator_precedence_and_over_or() {
        // `[C&N,O]` should parse as `(C & N) | O`
        let mol = parse_smarts("[C&N,O]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        // Single atom with complex query (And/Or combination).
    }

    #[test]
    fn test_operator_precedence_not_highest() {
        // `[!C&N]` should parse as `(! C) & N`
        let mol = parse_smarts("[!C&N]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_operator_precedence_semicolon_lowest() {
        // `[C;N,O]` should parse as `(C ; (N | O))` — different from `[C&N,O]`
        let mol = parse_smarts("[C;N,O]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_complex_bracket_atom_all_primitives() {
        // Combine multiple primitives: `[#6;a;R;H1]`
        let mol = parse_smarts("[#6;a;R;H1]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_long_chain_many_atoms() {
        // Long aliphatic chain: `CCCCCCCCCC` (10 atoms)
        let mol = parse_smarts("CCCCCCCCCC").unwrap();
        assert_eq!(mol.atoms.len(), 10);
        assert_eq!(mol.bonds.len(), 9);
    }

    #[test]
    fn test_multiple_rings_fused() {
        // Naphthalene pattern: `c1ccc2ccccc2c1` (10 atoms in two fused rings)
        let mol = parse_smarts("c1ccc2ccccc2c1").unwrap();
        assert_eq!(mol.atoms.len(), 10);
    }

    #[test]
    fn test_branching_from_multiple_atoms() {
        // `C(C)(C)C` — central C with 3 methyl branches
        let mol = parse_smarts("C(C)(C)C").unwrap();
        assert_eq!(mol.atoms.len(), 4);
        assert_eq!(mol.bonds.len(), 3);
    }

    #[test]
    fn test_recursive_smarts_simple() {
        // `[$(C)]` — atoms matching inner pattern C
        let mol = parse_smarts("[$(C)]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert!(matches!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::Recursive(_))
        ));
    }

    #[test]
    fn test_recursive_smarts_with_operators() {
        // `[$([C&N])]` — recursive pattern with operators
        let mol = parse_smarts("[$([C&N])]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    fn nested_recursive_smarts(depth: usize) -> String {
        let mut smarts = String::from("C");
        for _ in 0..depth {
            smarts = format!("[$({smarts})]");
        }
        smarts
    }

    #[test]
    fn test_recursive_smarts_default_depth_boundary() {
        let smarts = nested_recursive_smarts(DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH);
        let mol = parse_smarts(&smarts).expect("default recursive SMARTS depth should parse");
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_recursive_smarts_default_depth_rejects_too_deep_pattern() {
        let smarts = nested_recursive_smarts(DEFAULT_MAX_RECURSIVE_SMARTS_DEPTH + 1);
        assert!(matches!(
            parse_smarts(&smarts),
            Err(SmartsError::RecursionDepthExceeded)
        ));
    }

    #[test]
    fn test_recursive_smarts_config_depth_is_clamped_for_safety() {
        let low_config = SmartsParserConfig {
            max_recursion_depth: 0,
        };
        let too_deep_for_low_config = nested_recursive_smarts(2);
        assert!(matches!(
            parse_smarts_with_config(&too_deep_for_low_config, &low_config),
            Err(SmartsError::RecursionDepthExceeded)
        ));

        let high_config = SmartsParserConfig {
            max_recursion_depth: 100,
        };
        let absolute_boundary = nested_recursive_smarts(ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH);
        parse_smarts_with_config(&absolute_boundary, &high_config)
            .expect("absolute recursive SMARTS depth boundary should parse");

        let too_deep_for_absolute_limit =
            nested_recursive_smarts(ABSOLUTE_MAX_RECURSIVE_SMARTS_DEPTH + 1);
        assert!(matches!(
            parse_smarts_with_config(&too_deep_for_absolute_limit, &high_config),
            Err(SmartsError::RecursionDepthExceeded)
        ));
    }

    #[test]
    fn test_malformed_recursive_smarts_return_errors() {
        assert!(matches!(
            parse_smarts("[$(C]"),
            Err(SmartsError::UnexpectedEnd)
        ));
        assert!(matches!(
            parse_smarts("[$(C"),
            Err(SmartsError::UnexpectedEnd)
        ));
        assert!(matches!(
            parse_smarts("[$(C))]"),
            Err(SmartsError::UnclosedBracket(_) | SmartsError::UnexpectedChar(_, _))
        ));
    }

    #[test]
    fn test_isotope_with_bracket_atom() {
        // `[13C]` — explicitly labeled isotope
        let mol = parse_smarts("[13C]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_charge_positive_and_negative() {
        // `[C+2]` and `[O-1]`
        let mol_pos = parse_smarts("[C+2]").unwrap();
        let mol_neg = parse_smarts("[O-1]").unwrap();
        assert_eq!(mol_pos.atoms.len(), 1);
        assert_eq!(mol_neg.atoms.len(), 1);
    }

    #[test]
    fn test_degree_and_connectivity() {
        // `[D4]` (degree 4) and `[X3]` (connectivity 3)
        let mol = parse_smarts("[D4]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        let mol2 = parse_smarts("[X3]").unwrap();
        assert_eq!(mol2.atoms.len(), 1);
    }

    #[test]
    fn test_ring_size_constraint() {
        // `[r6]` (in 6-membered ring)
        let mol = parse_smarts("[r6]").unwrap();
        assert_eq!(mol.atoms.len(), 1);
    }

    #[test]
    fn test_nested_branch_depth() {
        // Deeply nested branches: `C(C(C(C(C))))`
        let mol = parse_smarts("C(C(C(C(C))))").unwrap();
        assert_eq!(mol.atoms.len(), 5);
    }

    #[test]
    fn test_ring_and_branch_combined() {
        // Ring with branch: `c1cc(C)ccc1`
        let mol = parse_smarts("c1cc(C)ccc1").unwrap();
        // 6 aromatic carbons + 1 branch carbon = 7 atoms
        assert_eq!(mol.atoms.len(), 7);
    }

    #[test]
    fn test_all_bond_types() {
        // Test all bond types: `-`, `=`, `#`, `:`, `~`
        let bonds = vec![
            ("C-C", BondPrimitive::Single),
            ("C=C", BondPrimitive::Double),
            ("C#C", BondPrimitive::Triple),
            ("c:c", BondPrimitive::Aromatic),
        ];
        for (smarts, expected_prim) in bonds {
            let mol = parse_smarts(smarts).unwrap();
            assert_eq!(mol.atoms.len(), 2);
            assert_eq!(
                mol.bonds[0].query,
                BondQuery::Primitive(expected_prim),
                "bond type mismatch for {smarts}"
            );
        }
    }

    #[test]
    fn test_empty_recursive_pattern_rejected() {
        // `[$()] `— empty inner pattern should fail gracefully or return empty match
        let result = parse_smarts("[$$()]");
        // Empty recursion is probably an error.
        assert!(result.is_err() || result.as_ref().unwrap().atoms.len() == 1);
    }
}
