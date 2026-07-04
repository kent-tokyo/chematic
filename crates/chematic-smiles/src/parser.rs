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

use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder, STEREO_H_SENTINEL,
};

use crate::error::SmilesError;

pub use chematic_core::Molecule;

/// Parse an OpenSMILES string into a [`Molecule`].
pub fn parse(input: &str) -> Result<Molecule, SmilesError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(SmilesError::EmptyInput);
    }
    let bytes = input.as_bytes();
    let mut p = Parser::new(bytes);
    p.parse_smiles()
}

const MAX_BRANCH_DEPTH: usize = 500;

/// Maximum number of atoms allowed in a SMILES molecule (prevents memory exhaustion).
const MAX_ATOMS: usize = 100_000;

/// An entry in the stereo neighbor sequence accumulated during parsing.
#[derive(Clone, Copy)]
enum StereoEntry {
    Atom(AtomIdx),
    ImplicitH,
    /// Ring opened at the chiral atom; will be resolved to `Atom` when the ring closes.
    PendingRing(u8),
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    depth: usize,
    /// Completed stereo records (may still contain `PendingRing` until `resolve_rings`).
    stereo_records: Vec<(AtomIdx, Vec<StereoEntry>)>,
    /// Stereo record being built for the current chiral atom.
    current_stereo: Option<(AtomIdx, Vec<StereoEntry>)>,
    /// ring_num → index in `stereo_records` for records with an unresolved `PendingRing(n)`.
    pending_ring_stereo: HashMap<u8, usize>,
    /// ring_num → close_atom, populated when rings close, for final resolution.
    ring_close_partners: HashMap<u8, AtomIdx>,
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            pos: 0,
            depth: 0,
            stereo_records: Vec::new(),
            current_stereo: None,
            pending_ring_stereo: HashMap::new(),
            ring_close_partners: HashMap::new(),
        }
    }

    /// Push an entry to the active stereo record (if one is open for the given atom).
    fn stereo_push(&mut self, current: AtomIdx, entry: StereoEntry) {
        if let Some((tracked, ref mut entries)) = self.current_stereo
            && tracked == current
        {
            entries.push(entry);
        }
    }

    /// Finalise the current stereo record: move it to `stereo_records` and register
    /// any `PendingRing` entries so they can be resolved when the ring closes.
    fn finalize_current_stereo(&mut self) {
        if let Some((atom_idx, entries)) = self.current_stereo.take() {
            let record_idx = self.stereo_records.len();
            for e in &entries {
                if let StereoEntry::PendingRing(rn) = e {
                    // Last-writer-wins for reused ring numbers (safe: ring must be closed
                    // before the same number can be reused).
                    self.pending_ring_stereo.insert(*rn, record_idx);
                }
            }
            self.stereo_records.push((atom_idx, entries));
        }
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
        if b.is_some() {
            self.pos += 1;
        }
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
            return Err(SmilesError::UnmatchedRingClosure {
                ring_num: num,
                pos: self.pos,
            });
        }

        // Anything left unconsumed is not valid SMILES (e.g. an unrecognised
        // atom symbol) — `parse_chain` treats an unparseable token as "end of
        // chain" rather than an error, so the check has to happen here.
        if self.pos < self.src.len() {
            return Err(SmilesError::UnexpectedCharacter { pos: self.pos });
        }

        // Finalise any active stereo record (shouldn't normally be needed).
        self.finalize_current_stereo();

        // Resolve any remaining PendingRing entries using recorded close partners.
        for (_, entries) in &mut self.stereo_records {
            for entry in entries.iter_mut() {
                if let StereoEntry::PendingRing(rn) = entry
                    && let Some(&partner) = self.ring_close_partners.get(rn)
                {
                    *entry = StereoEntry::Atom(partner);
                }
            }
        }

        // Store stereo neighbor orders in the molecule builder.
        let records = std::mem::take(&mut self.stereo_records);
        for (atom_idx, entries) in records {
            let order: Vec<u32> = entries
                .iter()
                .map(|e| match e {
                    StereoEntry::Atom(a) => a.0,
                    StereoEntry::ImplicitH => STEREO_H_SENTINEL,
                    StereoEntry::PendingRing(_) => STEREO_H_SENTINEL, // unresolved → skip
                })
                .collect();
            mol.set_stereo_neighbor_order(atom_idx, order);
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

        if mol.atom_count() >= MAX_ATOMS {
            return Err(SmilesError::InvalidBracketAtom {
                detail: format!("molecule exceeds maximum atom count {}", MAX_ATOMS),
                pos: self.pos,
            });
        }

        let first_idx = mol.add_atom(first_atom.clone());

        // Connect to the preceding atom if requested
        if let Some(prev) = attach_to {
            let bond = attach_bond.unwrap_or_else(|| implicit_bond(mol, prev, first_idx));
            mol.add_bond(prev, first_idx, bond)
                .map_err(|_| SmilesError::InvalidBracketAtom {
                    detail: "duplicate bond".to_string(),
                    pos: self.pos,
                })?;
        }

        // Save any parent stereo record (branches interrupt the parent's tracking).
        let saved_stereo = self.current_stereo.take();

        // Begin stereo tracking if the first atom is chiral.
        Self::begin_stereo_if_chiral(&first_atom, first_idx, attach_to, &mut self.current_stereo);

        // `current` is the atom we're currently processing.
        // `current_from` is the from-atom for `current` (needed to start stereo for next atoms).
        let mut current = first_idx;
        let mut current_from: Option<AtomIdx> = attach_to;

        // Process the rest of the chain
        loop {
            match self.peek() {
                Some(b'(') => {
                    let first_branch_atom = self.parse_branch(mol, current, None, open_rings)?;
                    if let Some(fa) = first_branch_atom {
                        self.stereo_push(current, StereoEntry::Atom(fa));
                    }
                }

                // Ring closure (digit or %nn) — no preceding bond char.
                Some(b'0'..=b'9') | Some(b'%') => {
                    let (ring_num, ring_bond) = self.parse_ring_num(None)?;
                    let ring_entry =
                        self.close_or_open_ring(mol, current, ring_num, ring_bond, open_rings)?;
                    self.stereo_push(current, ring_entry);
                }

                // End of this chain: ')' closes a branch, '.' starts new fragment, or EOF.
                None | Some(b')') | Some(b'.') => break,

                // Explicit bond or next atom.
                _ => {
                    let pending_bond = self.try_parse_bond();
                    match self.peek() {
                        Some(b'0'..=b'9') | Some(b'%') => {
                            let (ring_num, ring_bond) = self.parse_ring_num(pending_bond)?;
                            let ring_entry = self.close_or_open_ring(
                                mol, current, ring_num, ring_bond, open_rings,
                            )?;
                            self.stereo_push(current, ring_entry);
                        }
                        // Branch after explicit bond (unusual but valid: e.g. C=(C)C).
                        Some(b'(') => {
                            let first_branch_atom =
                                self.parse_branch(mol, current, pending_bond, open_rings)?;
                            if let Some(fa) = first_branch_atom {
                                self.stereo_push(current, StereoEntry::Atom(fa));
                            }
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
                                if mol.atom_count() >= MAX_ATOMS {
                                    return Err(SmilesError::InvalidBracketAtom {
                                        detail: format!(
                                            "molecule exceeds maximum atom count {}",
                                            MAX_ATOMS
                                        ),
                                        pos: self.pos,
                                    });
                                }
                                let next_idx = mol.add_atom(next_atom.clone());
                                // `/` and `\` between two aromatic atoms specify geometry of an
                                // adjacent double bond, not a stereo single bond. Aromatic atoms
                                // must remain connected by Aromatic bonds so SMARTS `:a` queries
                                // match correctly (e.g. Crippen `[c](:a)(:a)=[C,N,O]`). The
                                // original direction is stashed on the side (`bond_directions`)
                                // so an exocyclic E/Z double bond anchored on this ring bond
                                // survives into the canonical writer instead of being lost.
                                let mut stashed_direction = None;
                                let bond = match pending_bond {
                                    Some(dir @ (BondOrder::Up | BondOrder::Down))
                                        if mol.atom_at(current).aromatic
                                            && mol.atom_at(next_idx).aromatic =>
                                    {
                                        stashed_direction = Some(dir);
                                        BondOrder::Aromatic
                                    }
                                    Some(bo) => bo,
                                    None => implicit_bond(mol, current, next_idx),
                                };
                                let new_bond_idx =
                                    mol.add_bond(current, next_idx, bond).map_err(|_| {
                                        SmilesError::InvalidBracketAtom {
                                            detail: "duplicate bond".to_string(),
                                            pos: self.pos,
                                        }
                                    })?;
                                if let Some(dir) = stashed_direction {
                                    mol.set_bond_direction(new_bond_idx, dir);
                                }
                                // next_idx is the last stereo entry for `current`.
                                self.stereo_push(current, StereoEntry::Atom(next_idx));
                                // Finalise stereo for `current` before advancing.
                                self.finalize_current_stereo();
                                // Advance: next_idx is the new current, old current is its from-atom.
                                let old_current = current;
                                current = next_idx;
                                current_from = Some(old_current);
                                // Begin stereo for next_idx if it's chiral.
                                Self::begin_stereo_if_chiral(
                                    &next_atom,
                                    current,
                                    current_from,
                                    &mut self.current_stereo,
                                );
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

        // Finalise stereo for the last atom in this chain segment.
        self.finalize_current_stereo();
        // Restore the parent's stereo record (the branch is done).
        if saved_stereo.is_some() {
            self.current_stereo = saved_stereo;
        }

        let _ = current_from; // suppress unused warning
        Ok(Some(current))
    }

    /// Start a stereo record for `atom_idx` if `atom` has chirality.
    fn begin_stereo_if_chiral(
        atom: &Atom,
        atom_idx: AtomIdx,
        from_atom: Option<AtomIdx>,
        current_stereo: &mut Option<(AtomIdx, Vec<StereoEntry>)>,
    ) {
        if atom.chirality == Chirality::None {
            return;
        }
        let mut entries: Vec<StereoEntry> = Vec::new();
        if let Some(prev) = from_atom {
            entries.push(StereoEntry::Atom(prev));
        }
        if atom.hydrogen_count.is_some_and(|h| h > 0) {
            entries.push(StereoEntry::ImplicitH);
        }
        *current_stereo = Some((atom_idx, entries));
    }

    /// Parse a `(...)` branch: consume `(`, parse the inner chain, then `)`.
    /// Returns the first atom added inside the branch, or `None` if empty.
    fn parse_branch(
        &mut self,
        mol: &mut MoleculeBuilder,
        attach_to: AtomIdx,
        explicit_bond: Option<BondOrder>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<BondOrder>)>,
    ) -> Result<Option<AtomIdx>, SmilesError> {
        if self.depth >= MAX_BRANCH_DEPTH {
            return Err(SmilesError::NestingTooDeep { pos: self.pos });
        }
        self.depth += 1;
        self.advance(); // consume '('
        let bond = explicit_bond.or_else(|| self.try_parse_bond());
        let count_before = mol.atom_count();
        self.parse_chain(mol, Some(attach_to), bond, open_rings)?;
        self.depth -= 1;
        if self.peek() != Some(b')') {
            return Err(SmilesError::MismatchedParentheses { pos: self.pos });
        }
        self.advance(); // consume ')'
        let first_branch_atom = if mol.atom_count() > count_before {
            Some(AtomIdx(count_before as u32))
        } else {
            None
        };
        Ok(first_branch_atom)
    }

    /// Handle ring closure: close an existing open ring, or register a new one.
    /// Returns the stereo entry to push for the current atom's stereo sequence:
    /// - `Atom(open_atom)` when the ring is closed (partner = open_atom)
    /// - `PendingRing(ring_num)` when the ring is opened (partner unknown until close)
    fn close_or_open_ring(
        &mut self,
        mol: &mut MoleculeBuilder,
        current: AtomIdx,
        ring_num: u8,
        ring_bond: Option<BondOrder>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<BondOrder>)>,
    ) -> Result<StereoEntry, SmilesError> {
        if let Some((open_atom, open_bond)) = open_rings.remove(&ring_num) {
            // Resolve the bond type (both ends may specify one; they must agree)
            let bond = match (open_bond, ring_bond) {
                (Some(a), Some(b)) if a == b => a,
                (Some(_), Some(_)) => {
                    return Err(SmilesError::ConflictingRingBond {
                        ring_num,
                        pos: self.pos,
                    });
                }
                (Some(b), None) | (None, Some(b)) => b,
                (None, None) => implicit_bond(mol, open_atom, current),
            };
            mol.add_bond(open_atom, current, bond).map_err(|_| {
                SmilesError::InvalidBracketAtom {
                    detail: format!("duplicate ring bond {ring_num}"),
                    pos: self.pos,
                }
            })?;
            // Record the close partner for final PendingRing resolution.
            self.ring_close_partners.insert(ring_num, current);
            // Also resolve any PendingRing(ring_num) in already-finalized stereo records.
            if let Some(rec_idx) = self.pending_ring_stereo.remove(&ring_num)
                && let Some((_, entries)) = self.stereo_records.get_mut(rec_idx)
            {
                for entry in entries.iter_mut() {
                    if matches!(entry, StereoEntry::PendingRing(n) if *n == ring_num) {
                        *entry = StereoEntry::Atom(current);
                        break;
                    }
                }
            }
            // Return: the stereo entry for `current` is the open atom.
            Ok(StereoEntry::Atom(open_atom))
        } else {
            open_rings.insert(ring_num, (current, ring_bond));
            // Return: we opened a ring; partner not yet known.
            Ok(StereoEntry::PendingRing(ring_num))
        }
    }

    /// Parse a ring closure number (single digit or `%nn`), together with an
    /// optional `prefix_bond` that was already consumed before this call.
    fn parse_ring_num(
        &mut self,
        prefix_bond: Option<BondOrder>,
    ) -> Result<(u8, Option<BondOrder>), SmilesError> {
        let ring_num = if self.peek() == Some(b'%') {
            self.advance(); // consume '%'
            let tens = self
                .advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmilesError::UnexpectedEnd { pos: self.pos })?
                - b'0';
            let units = self
                .advance()
                .filter(|c| c.is_ascii_digit())
                .ok_or(SmilesError::UnexpectedEnd { pos: self.pos })?
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

    /// Consume a bond character and return the corresponding order, or `None`.
    fn try_parse_bond(&mut self) -> Option<BondOrder> {
        if self.peek() == Some(b'-') && self.peek_at(1) == Some(b'>') {
            self.advance();
            self.advance();
            return Some(BondOrder::Dative);
        }
        if self.peek() == Some(b'<') && self.peek_at(1) == Some(b'-') {
            self.advance();
            self.advance();
            return Some(BondOrder::Dative);
        }
        let order = match self.peek()? {
            b'-' => BondOrder::Single,
            b'=' => BondOrder::Double,
            b'#' => BondOrder::Triple,
            b'$' => BondOrder::Quadruple,
            b':' => BondOrder::Aromatic,
            b'/' => BondOrder::Up,
            b'\\' => BondOrder::Down,
            b'~' => BondOrder::QueryAny,
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

        let element = Element::from_symbol(&symbol).ok_or_else(|| SmilesError::UnknownElement {
            symbol: symbol.clone(),
            pos,
        })?;

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

        let element = Element::from_symbol(&symbol).ok_or_else(|| SmilesError::UnknownElement {
            symbol: symbol.clone(),
            pos,
        })?;

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
        let (symbol, aromatic) =
            self.parse_bracket_symbol()
                .ok_or_else(|| SmilesError::InvalidBracketAtom {
                    detail: "missing element symbol".to_string(),
                    pos: self.pos,
                })?;

        // Handle wildcard [*] — return immediately with a dedicated wildcard atom.
        if symbol == "*" {
            let chirality = self.parse_chirality();
            let _hcount = self.parse_hcount();
            let _charge = self.parse_charge();
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

        let element = Element::from_symbol(&symbol).ok_or_else(|| SmilesError::UnknownElement {
            symbol: symbol.clone(),
            pos: start_pos,
        })?;

        let chirality = self.parse_chirality();
        let hcount = self.parse_hcount();
        let charge = self.parse_charge();

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
        if let Some(second) = self.peek_at(1)
            && second.is_ascii_lowercase()
        {
            let candidate = format!("{upper_first}{}", second as char);
            if Element::from_symbol(&candidate).is_some() {
                self.advance();
                self.advance();
                return Some((candidate, aromatic));
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
                Some(d) => {
                    self.advance();
                    d - b'0'
                }
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
            val = val.saturating_mul(10).saturating_add((d - b'0') as u16);
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

    fn nested_branch_smiles(depth: usize) -> String {
        let mut smiles = String::from("C");
        for _ in 0..depth {
            smiles.push_str("(C");
        }
        for _ in 0..depth {
            smiles.push(')');
        }
        smiles
    }

    #[test]
    fn test_branch_depth_limit_accepts_boundary() {
        let smiles = nested_branch_smiles(MAX_BRANCH_DEPTH);
        let mol = parse(&smiles).expect("maximum configured branch depth should parse");
        assert_eq!(mol.atom_count(), MAX_BRANCH_DEPTH + 1);
    }

    #[test]
    fn test_branch_depth_limit_rejects_too_deep_input() {
        let smiles = nested_branch_smiles(MAX_BRANCH_DEPTH + 1);
        assert!(matches!(
            parse(&smiles),
            Err(SmilesError::NestingTooDeep { .. })
        ));
    }

    #[test]
    fn test_malformed_ring_closures_return_errors() {
        assert!(matches!(
            parse("C1"),
            Err(SmilesError::UnmatchedRingClosure { ring_num: 1, .. })
        ));
        assert!(matches!(
            parse("C%"),
            Err(SmilesError::UnexpectedEnd { .. })
        ));
        assert!(matches!(
            parse("C%1"),
            Err(SmilesError::UnexpectedEnd { .. })
        ));
        assert!(matches!(
            parse("C%AA"),
            Err(SmilesError::UnexpectedEnd { .. })
        ));
        assert!(matches!(
            parse("C=1CC-1"),
            Err(SmilesError::ConflictingRingBond { ring_num: 1, .. })
        ));
    }

    #[test]
    fn test_stereo_marker_between_aromatic_atoms_stays_aromatic() {
        // SMILES like `c1\c(O)c(O)\c1=N` use `/`/`\` to specify geometry of an
        // exocyclic double bond; the ring bonds themselves must remain Aromatic.
        // Regression: previously these were stored as Down/Up, breaking SMARTS (:a).
        use chematic_core::BondOrder;
        let mol = parse(r"N=c1\c(O)c(O)\c1=N").unwrap();
        let aromatic_ring_bonds: Vec<_> = mol
            .bonds()
            .filter(|(_, bond)| {
                mol.atom(bond.atom1).aromatic
                    && mol.atom(bond.atom2).aromatic
                    && bond.order == BondOrder::Aromatic
            })
            .collect();
        assert_eq!(
            aromatic_ring_bonds.len(),
            4,
            "all 4 ring bonds should be Aromatic"
        );
    }

    #[test]
    fn test_malformed_bracket_atoms_return_errors() {
        assert!(matches!(
            parse("[C"),
            Err(SmilesError::InvalidBracketAtom { .. })
        ));
        assert!(matches!(
            parse("[CH4"),
            Err(SmilesError::InvalidBracketAtom { .. })
        ));
        assert!(matches!(
            parse("[Q]"),
            Err(SmilesError::InvalidBracketAtom { .. })
        ));
        assert!(matches!(
            parse("[C@"),
            Err(SmilesError::InvalidBracketAtom { .. })
        ));
    }
}
