//! OpenSMILES recursive-descent parser.
//!
//! Grammar (simplified):
//! ```text
//! smiles     := chain ('.' chain)*
//! chain      := atom chain_rest*
//! chain_rest := branch | ring_closure | bond? atom | bond? ring_closure
//! branch     := '(' bond? chain ')'
//! ```
//!
//! Reference: OpenSMILES specification <http://opensmiles.org/opensmiles-spec.html>

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder};

use crate::error::SmilesError;

pub use chematic_core::Molecule;

/// Parse an OpenSMILES string into a [`Molecule`].
pub fn parse(input: &str) -> Result<Molecule, SmilesError> {
    if input.trim().is_empty() {
        return Err(SmilesError::EmptyInput);
    }
    let bytes = input.as_bytes();
    let mut p = Parser::new(bytes);
    p.parse_smiles()
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self { src, pos: 0 }
    }

    #[inline]
    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    #[inline]
    fn peek_at(&self, offset: usize) -> Option<u8> {
        self.src.get(self.pos + offset).copied()
    }

    #[inline]
    fn advance(&mut self) -> Option<u8> {
        let b = self.src.get(self.pos).copied();
        if b.is_some() { self.pos += 1; }
        b
    }

    fn parse_smiles(&mut self) -> Result<Molecule, SmilesError> {
        let mut mol = MoleculeBuilder::new();
        let mut open_rings: HashMap<u8, (AtomIdx, Option<BondOrder>)> = HashMap::new();

        // Parse the first fragment
        self.parse_chain(&mut mol, None, None, &mut open_rings)?;

        // Parse additional disconnected fragments separated by '.'
        while self.peek() == Some(b'.') {
            self.advance(); // consume '.'
            self.parse_chain(&mut mol, None, None, &mut open_rings)?;
        }

        // Trailing unmatched ring closures are errors
        if let Some((&num, _)) = open_rings.iter().next() {
            return Err(SmilesError::UnmatchedRingClosure { ring_num: num, pos: self.pos });
        }

        Ok(mol.build())
    }

    // Parse one atom-chain, attaching it to `attach_to` via `attach_bond`.
    // Returns the last atom index in this chain (or `attach_to` if nothing parsed).
    fn parse_chain(
        &mut self,
        mol: &mut MoleculeBuilder,
        attach_to: Option<AtomIdx>,
        attach_bond: Option<BondOrder>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<BondOrder>)>,
    ) -> Result<Option<AtomIdx>, SmilesError> {
        // Parse the first atom of this chain
        let first_atom = match self.try_parse_atom()? {
            Some(a) => a,
            None => return Ok(attach_to),
        };

        let first_idx = mol.add_atom(first_atom);

        // Connect to the preceding atom if requested
        if let Some(prev) = attach_to {
            let bond = attach_bond.unwrap_or_else(|| implicit_bond(mol, prev, first_idx));
            mol.add_bond(prev, first_idx, bond)
                .map_err(|_| SmilesError::InvalidBracketAtom {
                    detail: "duplicate bond".to_string(),
                    pos: self.pos,
                })?;
        }

        let mut current = first_idx;

        // Process the rest of the chain
        loop {
            match self.peek() {
                Some(b'(') => self.parse_branch(mol, current, None, open_rings)?,

                // Ring closure (digit or %nn) — no preceding bond char.
                Some(b'0'..=b'9') | Some(b'%') => {
                    let (ring_num, ring_bond) = self.parse_ring_num(None)?;
                    self.close_or_open_ring(mol, current, ring_num, ring_bond, open_rings)?;
                }

                // End of this chain: ')' closes a branch, '.' starts new fragment, or EOF.
                None | Some(b')') | Some(b'.') => break,

                // Explicit bond or next atom.
                _ => {
                    let pending_bond = self.try_parse_bond();
                    match self.peek() {
                        Some(b'0'..=b'9') | Some(b'%') => {
                            let (ring_num, ring_bond) = self.parse_ring_num(pending_bond)?;
                            self.close_or_open_ring(mol, current, ring_num, ring_bond, open_rings)?;
                        }
                        // Branch after explicit bond (unusual but valid: e.g. C=(C)C).
                        Some(b'(') => {
                            self.parse_branch(mol, current, pending_bond, open_rings)?;
                        }
                        // Disconnected or end — explicit bond with nothing after is an error.
                        None | Some(b')') | Some(b'.') => {
                            if pending_bond.is_some() {
                                return Err(SmilesError::UnexpectedEnd { pos: self.pos });
                            }
                            break;
                        }
                        _ => match self.try_parse_atom()? {
                            Some(next_atom) => {
                                let next_idx = mol.add_atom(next_atom);
                                let bond = pending_bond
                                    .unwrap_or_else(|| implicit_bond(mol, current, next_idx));
                                mol.add_bond(current, next_idx, bond).map_err(|_| {
                                    SmilesError::InvalidBracketAtom {
                                        detail: "duplicate bond".to_string(),
                                        pos: self.pos,
                                    }
                                })?;
                                current = next_idx;
                            }
                            None => {
                                if pending_bond.is_some() {
                                    return Err(SmilesError::UnexpectedEnd { pos: self.pos });
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

    /// Parse a `(...)` branch: consume `(`, parse the inner chain, then `)`.
    fn parse_branch(
        &mut self,
        mol: &mut MoleculeBuilder,
        attach_to: AtomIdx,
        explicit_bond: Option<BondOrder>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<BondOrder>)>,
    ) -> Result<(), SmilesError> {
        self.advance(); // consume '('
        let bond = explicit_bond.or_else(|| self.try_parse_bond());
        self.parse_chain(mol, Some(attach_to), bond, open_rings)?;
        if self.peek() != Some(b')') {
            return Err(SmilesError::MismatchedParentheses { pos: self.pos });
        }
        self.advance(); // consume ')'
        Ok(())
    }

    /// Handle ring closure: close an existing open ring, or register a new one.
    fn close_or_open_ring(
        &self,
        mol: &mut MoleculeBuilder,
        current: AtomIdx,
        ring_num: u8,
        ring_bond: Option<BondOrder>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<BondOrder>)>,
    ) -> Result<(), SmilesError> {
        if let Some((open_atom, open_bond)) = open_rings.remove(&ring_num) {
            // Resolve the bond type (both ends may specify one; they must agree)
            let bond = match (open_bond, ring_bond) {
                (Some(a), Some(b)) if a == b => a,
                (Some(_), Some(_)) => {
                    return Err(SmilesError::ConflictingRingBond { ring_num, pos: self.pos });
                }
                (Some(b), None) | (None, Some(b)) => b,
                (None, None) => implicit_bond(mol, open_atom, current),
            };
            mol.add_bond(open_atom, current, bond)
                .map_err(|_| SmilesError::InvalidBracketAtom {
                    detail: format!("duplicate ring bond {ring_num}"),
                    pos: self.pos,
                })?;
        } else {
            open_rings.insert(ring_num, (current, ring_bond));
        }
        Ok(())
    }

    /// Parse a ring closure number (single digit or `%nn`), together with an
    /// optional `prefix_bond` that was already consumed before this call.
    fn parse_ring_num(
        &mut self,
        prefix_bond: Option<BondOrder>,
    ) -> Result<(u8, Option<BondOrder>), SmilesError> {
        let ring_num = if self.peek() == Some(b'%') {
            self.advance(); // consume '%'
            let tens = self.advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmilesError::UnexpectedEnd { pos: self.pos })?
                - b'0';
            let units = self.advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmilesError::UnexpectedEnd { pos: self.pos })?
                - b'0';
            tens * 10 + units
        } else {
            // Caller peeked b'0'..=b'9', so advance() is guaranteed to return Some(digit).
            self.advance().expect("ring closure digit guaranteed by caller peek") - b'0'
        };
        Ok((ring_num, prefix_bond))
    }

    /// Consume a bond character and return the corresponding order, or `None`.
    fn try_parse_bond(&mut self) -> Option<BondOrder> {
        let order = match self.peek()? {
            b'-' => BondOrder::Single,
            b'=' => BondOrder::Double,
            b'#' => BondOrder::Triple,
            b'$' => BondOrder::Quadruple,
            b':' => BondOrder::Aromatic,
            b'/' => BondOrder::Up,
            b'\\' => BondOrder::Down,
            _ => return None,
        };
        self.advance();
        Some(order)
    }

    fn try_parse_atom(&mut self) -> Result<Option<Atom>, SmilesError> {
        match self.peek() {
            Some(b'[') => Ok(Some(self.parse_bracket_atom()?)),
            Some(b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I') => {
                Ok(Some(self.parse_organic_atom()?))
            }
            Some(b'b' | b'c' | b'n' | b'o' | b'p' | b's') => {
                Ok(Some(self.parse_aromatic_organic()?))
            }
            _ => Ok(None),
        }
    }

    fn parse_organic_atom(&mut self) -> Result<Atom, SmilesError> {
        let pos = self.pos;
        let first = self.advance().unwrap() as char;

        let symbol = if first == 'C' && self.peek() == Some(b'l') {
            self.advance();
            "Cl".to_string()
        } else if first == 'B' && self.peek() == Some(b'r') {
            self.advance();
            "Br".to_string()
        } else {
            first.to_string()
        };

        let element = Element::from_symbol(&symbol)
            .ok_or_else(|| SmilesError::UnknownElement { symbol: symbol.clone(), pos })?;

        let chirality = self.parse_chirality();
        let mut atom = Atom::organic(element);
        atom.chirality = chirality;
        Ok(atom)
    }

    fn parse_aromatic_organic(&mut self) -> Result<Atom, SmilesError> {
        let pos = self.pos;
        let first = self.advance().unwrap() as char;

        // Handle `se` and `as` written without brackets (rare but valid per spec)
        let (symbol, _multi) = if first == 's' && self.peek() == Some(b'e') {
            self.advance();
            ("Se".to_string(), true)
        } else if first == 'a' && self.peek() == Some(b's') {
            self.advance();
            ("As".to_string(), true)
        } else {
            (first.to_ascii_uppercase().to_string(), false)
        };

        let element = Element::from_symbol(&symbol)
            .ok_or_else(|| SmilesError::UnknownElement { symbol: symbol.clone(), pos })?;

        let chirality = self.parse_chirality();
        let mut atom = Atom::aromatic(element);
        atom.chirality = chirality;
        Ok(atom)
    }

    fn parse_bracket_atom(&mut self) -> Result<Atom, SmilesError> {
        let start_pos = self.pos;
        self.advance(); // consume '['

        // Optional isotope
        let isotope = self.parse_leading_digits_u16();

        // Element symbol (required)
        let (symbol, aromatic) = self.parse_bracket_symbol()
            .ok_or_else(|| SmilesError::InvalidBracketAtom {
                detail: "missing element symbol".to_string(),
                pos: self.pos,
            })?;

        // Handle wildcard [*] — return immediately with a dedicated wildcard atom.
        if symbol == "*" {
            let chirality = self.parse_chirality();
            let _hcount   = self.parse_hcount();
            let _charge   = self.parse_charge();
            if self.peek() == Some(b':') {
                self.advance();
                let _ = self.parse_leading_digits_u16(); // skip atom map
            }
            if self.peek() != Some(b']') {
                return Err(SmilesError::InvalidBracketAtom {
                    detail: "missing ']'".to_string(),
                    pos: self.pos,
                });
            }
            self.advance();
            let mut wc = Atom::wildcard();
            wc.chirality = chirality;
            return Ok(wc);
        }

        let element = Element::from_symbol(&symbol)
            .ok_or_else(|| SmilesError::UnknownElement { symbol: symbol.clone(), pos: start_pos })?;

        let chirality = self.parse_chirality();
        let hcount    = self.parse_hcount();
        let charge    = self.parse_charge();

        let atom_map = if self.peek() == Some(b':') {
            self.advance();
            self.parse_leading_digits_u16()
        } else {
            None
        };

        if self.peek() != Some(b']') {
            return Err(SmilesError::InvalidBracketAtom {
                detail: "missing ']'".to_string(),
                pos: self.pos,
            });
        }
        self.advance(); // consume ']'

        let mut atom = Atom::bracket(element, isotope, chirality, hcount, charge, atom_map);
        atom.aromatic = aromatic;
        Ok(atom)
    }

    /// Parse element symbol inside `[...]`. Returns (canonical_symbol, is_aromatic).
    fn parse_bracket_symbol(&mut self) -> Option<(String, bool)> {
        let first = self.peek()?;

        if first == b'*' {
            self.advance();
            return Some(("*".to_string(), false));
        }

        let aromatic = first.is_ascii_lowercase();
        let upper_first = first.to_ascii_uppercase() as char;

        // Try two-character symbol first
        if let Some(second) = self.peek_at(1) {
            if second.is_ascii_lowercase() {
                let candidate = format!("{upper_first}{}", second as char);
                if Element::from_symbol(&candidate).is_some() {
                    self.advance();
                    self.advance();
                    return Some((candidate, aromatic));
                }
            }
        }

        // Single-character symbol
        let sym = upper_first.to_string();
        if Element::from_symbol(&sym).is_some() {
            self.advance();
            Some((sym, aromatic))
        } else {
            None
        }
    }

    fn parse_chirality(&mut self) -> Chirality {
        if self.peek() == Some(b'@') {
            self.advance();
            if self.peek() == Some(b'@') {
                self.advance();
                Chirality::Clockwise
            } else {
                Chirality::CounterClockwise
            }
        } else {
            Chirality::None
        }
    }

    fn parse_hcount(&mut self) -> u8 {
        if self.peek() == Some(b'H') {
            self.advance();
            match self.peek().filter(|c| c.is_ascii_digit()) {
                Some(d) => { self.advance(); d - b'0' }
                None => 1,
            }
        } else {
            0
        }
    }

    fn parse_charge(&mut self) -> i8 {
        match self.peek() {
            Some(b'+') => {
                self.advance();
                if self.peek() == Some(b'+') {
                    self.advance();
                    return 2;
                }
                if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                    self.advance();
                    return (d - b'0') as i8;
                }
                1
            }
            Some(b'-') => {
                self.advance();
                if self.peek() == Some(b'-') {
                    self.advance();
                    return -2;
                }
                if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
                    self.advance();
                    return -((d - b'0') as i8);
                }
                -1
            }
            _ => 0,
        }
    }

    fn parse_leading_digits_u16(&mut self) -> Option<u16> {
        if !self.peek().is_some_and(|c| c.is_ascii_digit()) {
            return None;
        }
        let mut val: u16 = 0;
        while let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
            self.advance();
            val = val * 10 + (d - b'0') as u16;
        }
        Some(val)
    }
}

/// Determine the implicit bond between two adjacent atoms already in the builder.
///
/// Rule: if both atoms are aromatic → Aromatic bond; otherwise → Single bond.
fn implicit_bond(mol: &MoleculeBuilder, a: AtomIdx, b: AtomIdx) -> BondOrder {
    if mol.atom_at(a).aromatic && mol.atom_at(b).aromatic {
        BondOrder::Aromatic
    } else {
        BondOrder::Single
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::AtomIdx;

    #[test]
    fn test_parse_methane() {
        let mol = parse("C").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn test_parse_ethane() {
        let mol = parse("CC").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_parse_propane() {
        let mol = parse("CCC").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn test_parse_isobutane() {
        // CC(C)C — branched structure
        let mol = parse("CC(C)C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
    }

    #[test]
    fn test_parse_double_bond() {
        let mol = parse("C=C").unwrap();
        assert_eq!(mol.bond_count(), 1);
        let (_, bond) = mol.bonds().next().unwrap();
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn test_parse_triple_bond() {
        let mol = parse("C#N").unwrap();
        let (_, bond) = mol.bonds().next().unwrap();
        assert_eq!(bond.order, BondOrder::Triple);
    }

    #[test]
    fn test_parse_benzene_kekulized() {
        let mol = parse("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn test_parse_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        for (_, atom) in mol.atoms() {
            assert!(atom.aromatic);
        }
        for (_, bond) in mol.bonds() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_parse_pyridine() {
        let mol = parse("c1ccncc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn test_parse_naphthalene() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        assert_eq!(mol.atom_count(), 10);
        assert_eq!(mol.bond_count(), 11);
    }

    #[test]
    fn test_parse_bracket_water() {
        let mol = parse("[OH2]").unwrap();
        let atom = mol.atom(AtomIdx(0));
        assert_eq!(atom.element, Element::O);
        assert_eq!(atom.hydrogen_count, Some(2));
    }

    #[test]
    fn test_parse_ammonium() {
        let mol = parse("[NH4+]").unwrap();
        let atom = mol.atom(AtomIdx(0));
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.hydrogen_count, Some(4));
    }

    #[test]
    fn test_parse_13c() {
        let mol = parse("[13C]").unwrap();
        let atom = mol.atom(AtomIdx(0));
        assert_eq!(atom.isotope, Some(13));
    }

    #[test]
    fn test_parse_ethanol() {
        let mol = parse("CCO").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn test_parse_disconnected() {
        // Salt: sodium chloride
        let mol = parse("[Na+].[Cl-]").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn test_parse_acetic_acid() {
        // CC(=O)O
        let mol = parse("CC(=O)O").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
    }

    #[test]
    fn test_empty_smiles_error() {
        assert!(matches!(parse(""), Err(SmilesError::EmptyInput)));
    }

    #[test]
    fn test_parse_cyclohexane() {
        let mol = parse("C1CCCCC1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn test_parse_percent_ring() {
        // %10 ring closure
        let mol = parse("C%10CCCCC%10").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn test_parse_chlorobenzene() {
        let mol = parse("c1ccccc1Cl").unwrap();
        assert_eq!(mol.atom_count(), 7);
        assert_eq!(mol.bond_count(), 7);
    }
}
