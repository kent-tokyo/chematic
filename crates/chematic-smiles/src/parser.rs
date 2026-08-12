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
    SquarePlanarPermutation,
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
    /// Ring opened at the chiral atom; will be resolved to `Atom` when the ring
    /// closes. Keyed by a unique per-occurrence slot id (see `next_ring_slot`),
    /// NOT the raw ring digit -- ring digits are reused within a single SMILES
    /// (closed, then reused for an unrelated ring later), and resolving by
    /// digit let a later, unrelated reuse of the same digit silently
    /// overwrite/steal an earlier occurrence's still-pending resolution.
    PendingRing(u32),
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    depth: usize,
    /// Completed stereo records (may still contain `PendingRing` until `resolve_rings`).
    stereo_records: Vec<(AtomIdx, Vec<StereoEntry>)>,
    /// Stereo record being built for the current chiral atom.
    current_stereo: Option<(AtomIdx, Vec<StereoEntry>)>,
    /// Next fresh id to hand out for a ring-opening occurrence (see
    /// `PendingRing`). Monotonically increasing, never reused, so it uniquely
    /// identifies one specific open/close pair regardless of ring-digit reuse.
    next_ring_slot: u32,
    /// slot id → index in `stereo_records` for records with an unresolved `PendingRing(slot)`.
    pending_ring_stereo: HashMap<u32, usize>,
    /// slot id → close_atom, populated when rings close, for final resolution.
    ring_close_partners: HashMap<u32, AtomIdx>,
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            pos: 0,
            depth: 0,
            stereo_records: Vec::new(),
            current_stereo: None,
            next_ring_slot: 0,
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
                if let StereoEntry::PendingRing(slot) = e {
                    // Safe now: each slot id is unique to one open/close pair
                    // (see `next_ring_slot`), so there is no reuse to race.
                    self.pending_ring_stereo.insert(*slot, record_idx);
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
        let mut open_rings: HashMap<u8, (AtomIdx, Option<ParsedBond>, u32)> = HashMap::new();

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
        attach_bond: Option<ParsedBond>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<ParsedBond>, u32)>,
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
            let bond = attach_bond
                .unwrap_or_else(|| ParsedBond::plain(implicit_bond(mol, prev, first_idx)));
            // `<-` stores the bond donor→acceptor, i.e. with the endpoints
            // swapped relative to reading order.
            let (a1, a2) = bond.endpoints(prev, first_idx);
            let (bond, stash) = resolve_aromatic_direction_stash(mol, a1, a2, bond.order);
            let new_bond_idx =
                mol.add_bond(a1, a2, bond)
                    .map_err(|_| SmilesError::InvalidBracketAtom {
                        detail: "duplicate bond".to_string(),
                        pos: self.pos,
                    })?;
            if let Some(dir) = stash {
                mol.set_bond_direction(new_bond_idx, dir);
            }
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
                                let bond = match pending_bond {
                                    Some(bo) => bo,
                                    None => {
                                        ParsedBond::plain(implicit_bond(mol, current, next_idx))
                                    }
                                };
                                // `<-` stores the bond donor→acceptor, i.e.
                                // with the endpoints swapped relative to
                                // reading order.
                                let (a1, a2) = bond.endpoints(current, next_idx);
                                // See `resolve_aromatic_direction_stash`: a `/`/`\` between two
                                // aromatic atoms specifies geometry of an adjacent double bond,
                                // not a stereo single bond on this edge -- the original direction
                                // is stashed on the side (`bond_directions`) so it survives into
                                // the canonical writer instead of being lost.
                                let (bond, stashed_direction) =
                                    resolve_aromatic_direction_stash(mol, a1, a2, bond.order);
                                let new_bond_idx = mol.add_bond(a1, a2, bond).map_err(|_| {
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
        explicit_bond: Option<ParsedBond>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<ParsedBond>, u32)>,
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
    /// - `PendingRing(slot)` when the ring is opened (partner unknown until close)
    fn close_or_open_ring(
        &mut self,
        mol: &mut MoleculeBuilder,
        current: AtomIdx,
        ring_num: u8,
        ring_bond: Option<ParsedBond>,
        open_rings: &mut HashMap<u8, (AtomIdx, Option<ParsedBond>, u32)>,
    ) -> Result<StereoEntry, SmilesError> {
        if let Some((open_atom, open_bond, slot)) = open_rings.remove(&ring_num) {
            // A directional marker (`/`, `\`) is read "toward" the ring digit
            // from wherever it's written. At the OPENING occurrence (e.g.
            // "C/1..."), that's already the open->close direction, matching
            // the `mol.add_bond(open_atom, current, ..)` convention below. At
            // the CLOSING occurrence (e.g. "...C/1"), the marker instead reads
            // close->open (from the current atom back to its ring partner) --
            // the opposite traversal direction over the same bond -- so it
            // must be flipped before it can be compared against `open_bond` or
            // stored as the open->close direction.
            let ring_bond = flip_direction(ring_bond);
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
                (None, None) => ParsedBond::plain(implicit_bond(mol, open_atom, current)),
            };
            // `<-` stores the bond donor→acceptor, i.e. with the endpoints
            // swapped relative to the open→close reading order above.
            let (a1, a2) = bond.endpoints(open_atom, current);
            // See `resolve_aromatic_direction_stash`: without this guard, a
            // `/`/`\` ring-closure marker between two aromatic atoms becomes
            // a literal Up/Down bond between them -- inconsistent (aromatic
            // atoms must stay connected by Aromatic bonds), and downstream
            // E/Z perception reads the stray marker as a genuine stereo
            // descriptor that was never in the input, since the canonical
            // writer can route a stashed direction through a ring-closure
            // digit when that bond is chosen as the canonical back-edge
            // (see `dfs_mark`'s own stashed-direction handling above).
            let (bond, stash) = resolve_aromatic_direction_stash(mol, a1, a2, bond.order);
            let new_bond_idx =
                mol.add_bond(a1, a2, bond)
                    .map_err(|_| SmilesError::InvalidBracketAtom {
                        detail: format!("duplicate ring bond {ring_num}"),
                        pos: self.pos,
                    })?;
            if let Some(dir) = stash {
                mol.set_bond_direction(new_bond_idx, dir);
            }
            // Record the close partner for final PendingRing resolution, keyed
            // by this occurrence's unique slot -- NOT the ring digit, which
            // may be reused by an unrelated ring later in the same SMILES.
            self.ring_close_partners.insert(slot, current);
            // Also resolve PendingRing(slot) if the opener's stereo record was
            // already finalized (e.g. the closer is NOT nested inside the
            // opener's own branch). If the opener's record is still open (the
            // closer sits inside the opener's branch subtree, as in
            // `[C@]1(...[closes 1 here]...)`), this is a no-op here and the
            // entry is instead patched by the final resolution pass in
            // `parse_smiles`, which is safe now that `slot` can never collide
            // with a later, unrelated reuse of the same ring digit.
            if let Some(rec_idx) = self.pending_ring_stereo.remove(&slot)
                && let Some((_, entries)) = self.stereo_records.get_mut(rec_idx)
            {
                for entry in entries.iter_mut() {
                    if matches!(entry, StereoEntry::PendingRing(s) if *s == slot) {
                        *entry = StereoEntry::Atom(current);
                        break;
                    }
                }
            }
            // Return: the stereo entry for `current` is the open atom.
            Ok(StereoEntry::Atom(open_atom))
        } else {
            let slot = self.next_ring_slot;
            self.next_ring_slot += 1;
            open_rings.insert(ring_num, (current, ring_bond, slot));
            // Return: we opened a ring; partner not yet known.
            Ok(StereoEntry::PendingRing(slot))
        }
    }

    /// Parse a ring closure number (single digit or `%nn`), together with an
    /// optional `prefix_bond` that was already consumed before this call.
    fn parse_ring_num(
        &mut self,
        prefix_bond: Option<ParsedBond>,
    ) -> Result<(u8, Option<ParsedBond>), SmilesError> {
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

    /// Consume a bond token and return it, or `None`.
    fn try_parse_bond(&mut self) -> Option<ParsedBond> {
        if self.peek() == Some(b'-') && self.peek_at(1) == Some(b'>') {
            self.advance();
            self.advance();
            return Some(ParsedBond::plain(BondOrder::Dative));
        }
        if self.peek() == Some(b'<') && self.peek_at(1) == Some(b'-') {
            self.advance();
            self.advance();
            // `A<-B`: B is the donor, so the bond must be stored B→A.
            return Some(ParsedBond {
                order: BondOrder::Dative,
                reversed: true,
            });
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
        Some(ParsedBond::plain(order))
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

        let chirality = self.parse_chirality(false)?;
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

        let chirality = self.parse_chirality(false)?;
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
            let chirality = self.parse_chirality(true)?;
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

        let chirality = self.parse_chirality(true)?;
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

    /// `extended`: whether `@`-followed-by-a-class-tag (`@SP1`, `@TB4`, ...) is
    /// recognized at all. `true` only from bracket-atom contexts (`parse_bracket_atom`) --
    /// organic-subset atoms (`C@SP1...`) keep today's behavior byte-for-byte: `@`
    /// there is always plain tetrahedral, and a following `SP1` is left unconsumed
    /// for whatever error the caller's own subsequent parsing produces. Without this
    /// gating, `C@SP1(F)(Cl)Br` would silently become a "square-planar carbon" --
    /// meaningless outside a bracket atom, and a real behavior regression.
    fn parse_chirality(&mut self, extended: bool) -> Result<Chirality, SmilesError> {
        if self.peek() != Some(b'@') {
            return Ok(Chirality::None);
        }
        let at_pos = self.pos;
        self.advance();
        if self.peek() == Some(b'@') {
            self.advance();
            return Ok(Chirality::Clockwise);
        }
        if extended && let Some(class) = self.peek_chirality_class() {
            // Consume the whole token (class + digits) even when rejecting, so the
            // error points at this token instead of cascading into an unrelated
            // "missing ']'" a few characters later.
            let digits = self.parse_leading_digits_u16().unwrap_or(0);
            return match (class, digits) {
                ("SP", 1) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP1)),
                ("SP", 2) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP2)),
                ("SP", 3) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP3)),
                _ => Err(SmilesError::UnsupportedChiralityClass {
                    class: format!("{class}{digits}"),
                    pos: at_pos,
                }),
            };
        }
        Ok(Chirality::CounterClockwise)
    }

    /// If the next 2 bytes are exactly one of the OpenSMILES extended chirality
    /// class tags (`TH`/`AL`/`SP`/`TB`/`OH`), consumes them and returns the tag;
    /// otherwise leaves the position untouched (so a bare `@H` hydrogen-count
    /// marker, `@]`, `@+`, etc. are unaffected).
    fn peek_chirality_class(&mut self) -> Option<&'static str> {
        const CLASSES: [&str; 5] = ["TH", "AL", "SP", "TB", "OH"];
        let a = self.peek()?;
        let b = self.peek_at(1)?;
        let candidate = [a, b];
        for class in CLASSES {
            if class.as_bytes() == candidate {
                self.advance();
                self.advance();
                return Some(class);
            }
        }
        None
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

/// Flip a directional bond marker (`/` <-> `\`) to reverse its reading
/// direction; other bond orders pass through unchanged. Used to normalize a
/// ring-closure bond symbol captured at the closing occurrence (read
/// close->open) into the open->close sense used for storage and for
/// agreement-checking against the opening occurrence's symbol.
fn flip_direction(bond: Option<ParsedBond>) -> Option<ParsedBond> {
    bond.map(|b| match b.order {
        BondOrder::Up => ParsedBond::plain(BondOrder::Down),
        BondOrder::Down => ParsedBond::plain(BondOrder::Up),
        // A dative arrow is directional in exactly the same way, but its
        // direction lives in `reversed` rather than in the order itself.
        BondOrder::Dative => ParsedBond {
            order: BondOrder::Dative,
            reversed: !b.reversed,
        },
        _ => b,
    })
}

/// A bond token as written in the input, before it is attached to atoms.
///
/// `reversed` is `true` only for a `<-` dative arrow. `BondOrder` has a
/// single `Dative` variant whose stored `atom1 → atom2` order *is* the
/// donor → acceptor direction, so the arrow that points backwards relative
/// to the reading order cannot be encoded in the order itself — it has to
/// swap the two atoms when the bond is added. Without this, `A<-B` and
/// `A->B` produce byte-identical molecules and every `<-` in the input is
/// silently read as its own mirror image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ParsedBond {
    order: BondOrder,
    reversed: bool,
}

impl ParsedBond {
    fn plain(order: BondOrder) -> Self {
        Self {
            order,
            reversed: false,
        }
    }

    /// `(from, to)` in reading order, returned as the `(atom1, atom2)` pair
    /// the bond should actually be stored with.
    fn endpoints(self, from: AtomIdx, to: AtomIdx) -> (AtomIdx, AtomIdx) {
        if self.reversed {
            (to, from)
        } else {
            (from, to)
        }
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

/// A `/`/`\` bond order resolved between two atoms that are BOTH aromatic
/// describes the geometry of an adjacent (exocyclic) double bond, not a
/// stereo single bond on this edge itself. Aromatic atoms must stay
/// connected via `Aromatic` bonds (so SMARTS `:a` queries and re-perception
/// both work); the true direction is stashed on the side (`bond_directions`)
/// instead, so a stereo double bond anchored on this bond survives into the
/// canonical writer rather than being lost -- or, if this guard is skipped,
/// misapplied to this bond's own order on a later re-parse (which is
/// syntactically indistinguishable from a genuine directional single bond,
/// and can fabricate a stereo descriptor that was never in the input).
///
/// This is the single source of truth for that rule; every call site that
/// resolves a bond order from a possibly-directional pending symbol
/// (a plain chain bond, a branch-attachment bond, or a ring-closure bond)
/// must route through this so the guard can't drift out of sync between
/// them -- see Canonical-Stereo-D0.
fn resolve_aromatic_direction_stash(
    mol: &MoleculeBuilder,
    a: AtomIdx,
    b: AtomIdx,
    order: BondOrder,
) -> (BondOrder, Option<BondOrder>) {
    match order {
        dir @ (BondOrder::Up | BondOrder::Down)
            if mol.atom_at(a).aromatic && mol.atom_at(b).aromatic =>
        {
            (BondOrder::Aromatic, Some(dir))
        }
        other => (other, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::canonical_smiles;
    use chematic_core::{AtomIdx, BondIdx};

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

    // -----------------------------------------------------------------
    // Canonical-Stereo-D0: `resolve_aromatic_direction_stash` must be
    // applied consistently at all three sites that can create a bond from a
    // pending `/`/`\` marker -- the normal chain-edge path (already had the
    // guard), the ring-closure path, and the branch-attachment path (the
    // first bond of a `(...)` group) -- both of which silently stored the
    // marker as the bond's own literal `Up`/`Down` order instead of
    // stashing it. That inconsistency meant the SAME physical bond could be
    // correctly stashed on one parse (when it happened to fall on a chain
    // edge) and incorrectly stored literally on a later parse of chematic's
    // own canonical output (when the SAME bond became a ring-closure or
    // branch-attachment edge instead) -- silently changing which CIP
    // elements were visible to `assign_ez` between the two parses, purely
    // as an artifact of which of the three paths happened to run.
    //
    // These tests check the structural representation (bond order + stash)
    // is correct and stable across a round trip for each of the three
    // paths. They deliberately do NOT assert `assign_cip`'s E/Z output --
    // `assign_ez`/`substituent_is_up` do not read `bond_direction` at all
    // (a separate, pre-existing, real gap: this direction genuinely encodes
    // legitimate exocyclic-imine E/Z per an independent RDKit check, not
    // fabricated stereo -- see EZ-A0/EZ-S1). Pinning `[]` here as the
    // expected CIP output would calcify that gap as intended behavior and
    // get in the way of the follow-up that fixes it.
    // -----------------------------------------------------------------

    /// Assert `bidx` is `Aromatic` with a stashed direction, on a molecule
    /// where both endpoints are aromatic -- the correct representation
    /// regardless of which parser path created the bond.
    #[track_caller]
    fn assert_aromatic_with_stash(mol: &Molecule, bidx: BondIdx, path_label: &str) {
        let bond = mol.bond(bidx);
        assert!(
            mol.atom(bond.atom1).aromatic && mol.atom(bond.atom2).aromatic,
            "{path_label}: test setup sanity -- both endpoints must be aromatic"
        );
        assert_eq!(
            bond.order,
            BondOrder::Aromatic,
            "{path_label}: bond order must stay Aromatic, not the literal Up/Down marker"
        );
        assert!(
            mol.bond_direction(bidx).is_some(),
            "{path_label}: the direction must be stashed on the side channel"
        );
    }

    /// Re-parse `canonical_smiles(mol)` and assert the same structural
    /// invariant holds again: some aromatic-aromatic bond is `Aromatic`
    /// order with a stashed direction (not necessarily the same `BondIdx`,
    /// since canonicalization renumbers atoms/bonds).
    #[track_caller]
    fn assert_round_trip_preserves_stash_representation(mol: &Molecule, path_label: &str) {
        let c1 = canonical_smiles(mol);
        let mol2 = parse(&c1).unwrap_or_else(|e| panic!("{path_label}: re-parse '{c1}': {e}"));
        let has_stashed_aromatic_bond = (0..mol2.bond_count()).any(|i| {
            let bidx = BondIdx(i as u32);
            let bond = mol2.bond(bidx);
            mol2.atom(bond.atom1).aromatic
                && mol2.atom(bond.atom2).aromatic
                && bond.order == BondOrder::Aromatic
                && mol2.bond_direction(bidx).is_some()
        });
        assert!(
            has_stashed_aromatic_bond,
            "{path_label}: round-tripped molecule '{c1}' must still represent the direction \
             as a stash on an Aromatic bond, not a literal Up/Down order"
        );
        let c2 = canonical_smiles(&mol2);
        assert_eq!(
            c1, c2,
            "{path_label}: canonical_smiles must be stable across a round trip"
        );
    }

    #[test]
    fn direction_stash_normal_chain_edge() {
        // The `/`/`\` sits between the ring-opening atom and the very next
        // chain atom -- the plain tree-edge path (already had the guard
        // before this fix; kept as a path-1 regression pin, not a new fix).
        let mol = parse(r"N=c1\c(O)c(O)\c1=N").unwrap();
        assert_aromatic_with_stash(&mol, BondIdx(1), "path1(chain-edge)");
        assert_round_trip_preserves_stash_representation(&mol, "path1(chain-edge)");
    }

    #[test]
    fn direction_stash_ring_closure_edge() {
        // The `/` sits directly on the ring-closing bond itself (`...c/1`),
        // routed through `close_or_open_ring` -- previously stored as a
        // literal Up/Down order between two aromatic atoms.
        let mol = parse(r"C/N=c1ccccc/1").unwrap();
        assert_aromatic_with_stash(&mol, BondIdx(7), "path2(ring-closure)");
        assert_round_trip_preserves_stash_representation(&mol, "path2(ring-closure)");
    }

    #[test]
    fn direction_stash_branch_attachment_edge() {
        // The `/` is the first bond inside a `(...)` branch, attaching the
        // branch's first atom to its parent -- routed through
        // `parse_chain`'s `attach_to` handling, not the plain chain-edge
        // loop. Previously stored as a literal Up/Down order between two
        // aromatic atoms (this is the real-world corpus mechanism: a
        // canonical DFS can route what was a chain edge on the original
        // parse through a branch attachment on the next).
        let mol = parse(r"Cc1ccc(/c1)N").unwrap();
        assert_aromatic_with_stash(&mol, BondIdx(4), "path3(branch-attachment)");
        assert_round_trip_preserves_stash_representation(&mol, "path3(branch-attachment)");
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
