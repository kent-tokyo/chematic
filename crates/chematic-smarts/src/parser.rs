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
//!               | 'r' DIGIT                      // ring size
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

use crate::query::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule,
};

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
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Parse a SMARTS pattern string into a `QueryMolecule`.
///
/// Returns `Err(SmartsError)` if the input is not a valid SMARTS pattern.
pub fn parse_smarts(smarts: &str) -> Result<QueryMolecule, SmartsError> {
    let mut parser = Parser::new(smarts.as_bytes());
    parser.parse()
}

// ---------------------------------------------------------------------------
// Parser internals
// ---------------------------------------------------------------------------

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self { src, pos: 0 }
    }

    // -- character helpers ---------------------------------------------------

    #[inline]
    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    #[inline]
    #[allow(dead_code)]
    fn peek_at(&self, offset: usize) -> Option<u8> {
        self.src.get(self.pos + offset).copied()
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
        if let Some(c) = self.peek() {
            if c != b' ' && c != b'\t' && c != b'\n' && c != b'\r' {
                return Err(SmartsError::UnexpectedChar(c as char, self.pos));
            }
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
        let first_atom = match self.try_parse_atom()? {
            Some(a) => a,
            None => return Ok(attach_to),
        };

        let first_idx = mol.add_atom(first_atom);

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
                            self.handle_ring_closure(mol, current, ring_num, ring_bond, open_rings)?;
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
                        _ => {
                            match self.try_parse_atom()? {
                                Some(next_atom) => {
                                    let next_idx = mol.add_atom(next_atom);
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
                            }
                        }
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
            self.advance().unwrap() - b'0'
        };
        Ok((ring_num, prefix_bond))
    }

    // -- bond ----------------------------------------------------------------

    /// Try to parse one bond character and return the corresponding `BondQuery`.
    /// Returns `None` if the next character is not a bond character.
    fn try_parse_bond(&mut self) -> Option<BondQuery> {
        let prim = match self.peek()? {
            b'-' => BondPrimitive::Single,
            b'=' => BondPrimitive::Double,
            b'#' => {
                // '#' could be part of a bracket `[#6]` — but at the chain level
                // (outside brackets) '#' is unambiguously a bond character.
                BondPrimitive::Triple
            }
            b':' => BondPrimitive::Aromatic,
            b'~' => BondPrimitive::Any,
            b'@' => BondPrimitive::Ring,
            _ => return None,
        };
        self.advance();
        Some(BondQuery::Primitive(prim))
    }

    // -- atoms ---------------------------------------------------------------

    /// Try to parse one atom (bracket, organic shorthand, or wildcard).
    fn try_parse_atom(&mut self) -> Result<Option<AtomQuery>, SmartsError> {
        match self.peek() {
            Some(b'[') => Ok(Some(self.parse_bracket_atom()?)),
            Some(b'*') => {
                self.advance();
                Ok(Some(AtomQuery::Primitive(AtomPrimitive::Wildcard)))
            }
            // Organic-subset atoms (uppercase aliphatic).
            Some(b'B') | Some(b'C') | Some(b'N') | Some(b'O')
            | Some(b'P') | Some(b'S') | Some(b'F') | Some(b'I') => {
                Ok(Some(self.parse_organic_atom(false)?))
            }
            // Aromatic organic-subset atoms (lowercase).
            Some(b'c') | Some(b'n') | Some(b'o') | Some(b'p') | Some(b's') => {
                Ok(Some(self.parse_organic_atom(true)?))
            }
            // 'b' could be boron lowercase — handle similarly.
            Some(b'b') => Ok(Some(self.parse_organic_atom(true)?)),
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

    /// Parse a bracket atom `[expr]`.
    fn parse_bracket_atom(&mut self) -> Result<AtomQuery, SmartsError> {
        let bracket_pos = self.pos;
        self.advance(); // consume '['

        let expr = self.parse_expr()?;

        match self.peek() {
            Some(b']') => {
                self.advance(); // consume ']'
            }
            _ => return Err(SmartsError::UnclosedBracket(bracket_pos)),
        }

        Ok(expr)
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
            Some(b'*') => true,
            // Recursive SMARTS `$(...)`.
            Some(b'$') => true,
            // Valence `[vN]`, ring-bond count `[xN]`, hybridization `[^N]`.
            Some(b'v') | Some(b'x') | Some(b'^') => true,
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
                let inner_mol = parse_smarts(inner_str)?;
                self.pos = end + 1; // advance past the closing ')'
                Ok(AtomQuery::Primitive(AtomPrimitive::Recursive(Box::new(inner_mol))))
            }

            // Aromatic (`a`) — must be checked BEFORE element symbol parsing.
            Some(b'a') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(true)))
            }

            // Aliphatic (`A`) — checked BEFORE element parsing (there is no element 'A').
            Some(b'A') => {
                self.advance();
                Ok(AtomQuery::Primitive(AtomPrimitive::Aromatic(false)))
            }

            // Atomic number `#N`
            Some(b'#') => {
                self.advance(); // consume '#'
                let n = self
                    .parse_digits_u8()
                    .ok_or(SmartsError::UnexpectedEnd)?;
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

            // H count `HN` or `H` (hcount = N or 1)
            Some(b'H') => {
                self.advance(); // consume 'H'
                let n = self.parse_single_digit().unwrap_or(1);
                Ok(AtomQuery::Primitive(AtomPrimitive::HCount(n)))
            }

            // Degree `DN`
            Some(b'D') => {
                self.advance(); // consume 'D'
                let n = self.parse_single_digit().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Degree(n)))
            }

            // Ring size `rN` — must be checked BEFORE element parsing.
            Some(b'r') => {
                self.advance(); // consume 'r'
                let n = self.parse_single_digit().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::RingSize(n)))
            }

            // Ring membership `R`
            Some(b'R') => {
                self.advance(); // consume 'R'
                Ok(AtomQuery::Primitive(AtomPrimitive::RingMembership(true)))
            }

            // Valence `[vN]` — total valence (bond orders + implicit H).
            Some(b'v') => {
                self.advance(); // consume 'v'
                let n = self.parse_single_digit().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Valence(n)))
            }

            // Ring-bond count `[xN]` — bonds where both endpoints share a ring.
            Some(b'x') => {
                self.advance(); // consume 'x'
                let n = self.parse_single_digit().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::RingBondCount(n)))
            }

            // Hybridization `[^N]` — 1=sp, 2=sp2, 3=sp3.
            Some(b'^') => {
                self.advance(); // consume '^'
                let n = self.parse_single_digit().ok_or(SmartsError::UnexpectedEnd)?;
                Ok(AtomQuery::Primitive(AtomPrimitive::Hybridization(n)))
            }

            // Element symbol (uppercase or lowercase start).
            Some(c) if c.is_ascii_alphabetic() => {
                self.parse_element_primitive()
            }

            Some(c) => Err(SmartsError::UnexpectedChar(c as char, pos)),
            None => Err(SmartsError::UnexpectedEnd),
        }
    }

    /// Parse an element symbol as a primitive inside a bracket atom.
    ///
    /// Accepts both uppercase (aliphatic) and lowercase (aromatic) symbols.
    fn parse_element_primitive(&mut self) -> Result<AtomQuery, SmartsError> {
        let pos = self.pos;
        let first = self.advance().unwrap() as char;
        let _aromatic = first.is_ascii_lowercase();
        let upper_first = first.to_ascii_uppercase();

        // Try two-character symbol first (e.g. `Cl`, `Br`).
        if let Some(second) = self.peek() {
            if second.is_ascii_lowercase() {
                let candidate = format!("{upper_first}{}", second as char);
                if chematic_core::Element::from_symbol(&candidate).is_some() {
                    self.advance();
                    return Ok(AtomQuery::Primitive(AtomPrimitive::Symbol(candidate)));
                }
            }
        }

        // Single-character symbol.
        let sym = upper_first.to_string();
        if chematic_core::Element::from_symbol(&sym).is_some() {
            Ok(AtomQuery::Primitive(AtomPrimitive::Symbol(sym)))
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
        Some(val as u8)
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
        // `[!C]` → Not(Symbol("C"))
        let mol = parse_smarts("[!C]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Not(Box::new(AtomQuery::Primitive(AtomPrimitive::Symbol("C".to_string()))))
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
        // `[r5]` → RingSize(5)
        let mol = parse_smarts("[r5]").unwrap();
        assert_eq!(
            mol.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::RingSize(5))
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
}
