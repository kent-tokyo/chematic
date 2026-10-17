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
    STEREO_VACANT_SENTINEL, SquarePlanarPermutation,
};

use crate::error::SmilesError;

mod stereo;

pub use chematic_core::Molecule;

/// Parse an OpenSMILES string into a [`Molecule`].
pub fn parse(input: &str) -> Result<Molecule, SmilesError> {
    parse_with_limits(input, &SmilesParseLimits::default())
}

/// Resource limits for parsing untrusted SMILES.
#[derive(Debug, Clone, Copy)]
pub struct SmilesParseLimits {
    /// Maximum UTF-8 input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum atoms in the parsed molecule.
    pub max_atoms: usize,
    /// Maximum bonds in the parsed molecule.
    pub max_bonds: usize,
}

impl Default for SmilesParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_atoms: MAX_ATOMS,
            max_bonds: 200_000,
        }
    }
}

/// Parse an OpenSMILES string while enforcing input and graph-size limits.
pub fn parse_with_limits(input: &str, limits: &SmilesParseLimits) -> Result<Molecule, SmilesError> {
    parse_checked(input, limits, true)
}

/// Parse one side of a reaction template (SMIRKS or reaction SMARTS read as
/// SMILES) into its template graph.
///
/// The same grammar and limits as [`parse`], without the check that rejects
/// a neutral oxygen or fluorine with too many bonds
/// ([`SmilesError::InvalidValence`], #769): a template atom is a pattern, not
/// an atom of a molecule (RDKit does not sanitize templates either), and a
/// product it would write with an impossible valence is refused when the
/// template is applied.
pub fn parse_template(input: &str) -> Result<Molecule, SmilesError> {
    parse_checked(input, &SmilesParseLimits::default(), false)
}

fn parse_checked(
    input: &str,
    limits: &SmilesParseLimits,
    check_valence: bool,
) -> Result<Molecule, SmilesError> {
    crate::rdkit::register_rdkit_model_hook();
    if input.len() > limits.max_input_bytes {
        return Err(SmilesError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let mol = parse_unbounded(input)?;
    if mol.atom_count() > limits.max_atoms {
        return Err(SmilesError::ResourceLimit {
            resource: "atoms",
            actual: mol.atom_count(),
            limit: limits.max_atoms,
        });
    }
    if mol.bond_count() > limits.max_bonds {
        return Err(SmilesError::ResourceLimit {
            resource: "bonds",
            actual: mol.bond_count(),
            limit: limits.max_bonds,
        });
    }
    if check_valence {
        check_neutral_valence(&mol)?;
        // Every RDKit-compatible descriptor/fingerprint asks the same
        // question before deciding whether it must reconstruct RDKit's
        // sanitized molecule.  A parsed SMILES is immutable on return and
        // mutations clear derived data, so answer and memoize that shared
        // gate once while the input is still hot in cache.  This is not
        // eager aromaticity or descriptor perception: it is the same linear
        // safety predicate those calls would otherwise repeat on first use.
        let _ = chematic_perception::rdkit_model_may_disagree(&mol);
        // Likewise, many hot compatibility APIs need to know whether an
        // acyclic parsed graph already is RDKit's aromatic representation.
        // The parser has already proved and seeded the no-ring case, so this
        // remains a cheap verdict.  Cyclic inputs deliberately stay lazy:
        // deciding them may perform ring perception and does not belong in
        // parsing merely to improve a downstream benchmark.
        if mol.smiles_ring_closure_bond_count() == 0 {
            let _ = chematic_perception::rdkit_parity_view_is_identity(&mol);
        }
    }
    Ok(mol)
}

fn parse_unbounded(input: &str) -> Result<Molecule, SmilesError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(SmilesError::EmptyInput);
    }
    let bytes = input.as_bytes();
    let mut p = Parser::new(bytes);
    p.parse_smiles()
}

/// Reject a neutral oxygen with explicit valence above 2 or a neutral
/// fluorine above 1 (#769). Neither element has a hypervalent neutral state,
/// so such a graph (`O=O1C=CC=C1`, `CO(C)C`, `F(C)C`) cannot be a molecule;
/// accepting it let aromaticity perception aromatize the ring through the
/// invalid oxygen. The valence is the sum of bond orders plus explicit
/// hydrogens. An aromatic bond counts 1 (an aromatic oxygen or furan-type
/// atom has two sigma bonds), a dative bond counts only at its acceptor, and
/// zero-order and query bonds count nothing. Charged atoms and the
/// hypervalent states of S, P, N and the heavier halogens are left to the
/// existing contracts.
fn check_neutral_valence(mol: &Molecule) -> Result<(), SmilesError> {
    for (idx, atom) in mol.atoms() {
        if atom.charge != 0 || atom.wildcard {
            continue;
        }
        let max = match atom.element {
            Element::O => 2,
            Element::F => 1,
            _ => continue,
        };
        let mut valence = u32::from(atom.hydrogen_count.unwrap_or(0));
        for (nbr, bond) in mol.neighbors(idx) {
            let entry = mol.bond(bond);
            // A single bond to a metal is a coordination bond (RDKit's
            // `cleanUpOrganometallics` makes it dative from the ligand), as
            // in aqua and alkoxide complexes: it adds nothing here.
            let other = mol.atom(nbr);
            if entry.order == BondOrder::Single
                && !other.wildcard
                && is_metal(other.element.atomic_number())
            {
                continue;
            }
            valence += match entry.order {
                BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Aromatic => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
                BondOrder::Dative => u32::from(entry.atom2 == idx),
                _ => 0,
            };
        }
        if valence > max {
            return Err(SmilesError::InvalidValence {
                element: atom.element.symbol(),
                atom: idx.0 as usize,
                valence,
                max,
            });
        }
    }
    Ok(())
}

/// RDKit's `QueryOps::isMetal`: everything but H, the noble gases and the
/// non-metals B, C, N, O, F, Si, P, S, Cl, As, Se, Br, Te, I, At.
fn is_metal(anum: u8) -> bool {
    !matches!(
        anum,
        0 | 1 | 2 | 5..=10 | 14..=18 | 33..=36 | 52..=54 | 85 | 86
    )
}

const MAX_BRANCH_DEPTH: usize = 500;

/// Maximum number of atoms allowed in a SMILES molecule (prevents memory exhaustion).
const MAX_ATOMS: usize = 100_000;

// OpenSMILES ring tokens are a digit or exactly two digits after `%`. SMILES+
// additionally uses `%(n)` for labels >= 100. Keeping only currently-open
// labels makes the extended form bounded by the parsed graph rather than by a
// label value supplied by untrusted input.
type OpenRings = HashMap<u32, (AtomIdx, Option<ParsedBond>, u32)>;
const MAX_OPEN_RING_LABELS: usize = MAX_ATOMS;

/// An entry in the stereo neighbor sequence accumulated during parsing.
#[derive(Clone, Copy)]
enum StereoEntry {
    Atom(AtomIdx),
    ImplicitH,
    Vacant,
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
    ring_close_partners: smallvec::SmallVec<[Option<AtomIdx>; 8]>,
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
            ring_close_partners: smallvec::SmallVec::new(),
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
        let mut open_rings = OpenRings::new();

        // Parse the first fragment
        self.parse_chain(&mut mol, None, None, &mut open_rings)?;

        // Parse additional disconnected fragments separated by '.'
        while self.peek() == Some(b'.') {
            self.advance(); // consume '.'
            self.parse_chain(&mut mol, None, None, &mut open_rings)?;
        }

        // Trailing unmatched ring closures are errors
        // Keep the legacy diagnostic deterministic: the old fixed-size table
        // reported the lowest unmatched label. HashMap iteration is random.
        if let Some(&num) = open_rings.keys().min() {
            return Err(ring_closure_unmatched_error(num, self.pos));
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
                    && let Some(&Some(partner)) = self.ring_close_partners.get(*rn as usize)
                {
                    *entry = StereoEntry::Atom(partner);
                }
            }
        }
        // Store stereo neighbor orders in the molecule builder.
        let records = std::mem::take(&mut self.stereo_records);
        let tagged_atoms: Vec<AtomIdx> = records.iter().map(|(atom, _)| *atom).collect();
        for (atom_idx, entries) in records {
            let order: Vec<u32> = entries
                .iter()
                .map(|e| match e {
                    StereoEntry::Atom(a) => a.0,
                    StereoEntry::ImplicitH => STEREO_H_SENTINEL,
                    StereoEntry::Vacant => STEREO_VACANT_SENTINEL,
                    StereoEntry::PendingRing(_) => STEREO_H_SENTINEL, // unresolved → skip
                })
                .collect();
            mol.set_stereo_neighbor_order(atom_idx, order);
        }

        let mut molecule = mol.build();
        normalize_ring_wildcards(&mut molecule);
        stereo::retain_vacant_sites(&mut molecule, &tagged_atoms);
        crate::allene::normalize_allene_stereo(&mut molecule, &tagged_atoms);
        // A SMILES graph with no ring-closure token is a forest by
        // construction. Publish that already-known result so the first
        // fingerprint or descriptor call does not need a full bridge DFS to
        // rediscover it. This is parser provenance, not eager ring
        // perception; edits invalidate the derived cache as usual.
        if molecule.smiles_ring_closure_bond_count() == 0 {
            molecule.seed_derived(
                chematic_core::DerivedSlot::RingBondFlags,
                std::sync::Arc::new(vec![false; molecule.bond_count()]),
            );
        }
        Ok(molecule)
    }

    // Parse one atom-chain, attaching it to `attach_to` via `attach_bond`.
    // Returns the last atom index in this chain (or `attach_to` if nothing parsed).
    fn parse_chain(
        &mut self,
        mol: &mut MoleculeBuilder,
        attach_to: Option<AtomIdx>,
        attach_bond: Option<ParsedBond>,
        open_rings: &mut OpenRings,
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
                mol.set_bond_direction_anchor(new_bond_idx, a1);
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
                                    mol.set_bond_direction_anchor(new_bond_idx, a1);
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

    /// Parse a `(...)` branch: consume `(`, parse the inner chain, then `)`.
    /// Returns the first atom added inside the branch, or `None` if empty.
    fn parse_branch(
        &mut self,
        mol: &mut MoleculeBuilder,
        attach_to: AtomIdx,
        explicit_bond: Option<ParsedBond>,
        open_rings: &mut OpenRings,
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
        ring_num: u32,
        ring_bond: Option<ParsedBond>,
        open_rings: &mut OpenRings,
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
                    return Err(ring_closure_conflict_error(ring_num, self.pos));
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
                mol.set_bond_direction_anchor(new_bond_idx, a1);
            }
            mol.set_smiles_ring_closure(new_bond_idx, ring_num, slot);
            // RDKit's `CloseMolRings` keeps the opening partial bond (begin =
            // opening atom) only when it carried an explicit, non-directional
            // bond symbol; dative arrows set their own orientation.
            if bond != BondOrder::Dative {
                let at_open = matches!(open_bond, Some(b) if !matches!(b.order, BondOrder::Up | BondOrder::Down));
                mol.set_smiles_ring_closure_begins_at_open(new_bond_idx, at_open);
            }
            // Record the close partner for final PendingRing resolution, keyed
            // by this occurrence's unique slot -- NOT the ring digit, which
            // may be reused by an unrelated ring later in the same SMILES.
            self.ring_close_partners[slot as usize] = Some(current);
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
            self.ring_close_partners.push(None);
            if open_rings.len() >= MAX_OPEN_RING_LABELS {
                return Err(SmilesError::ResourceLimit {
                    resource: "open ring labels",
                    actual: open_rings.len() + 1,
                    limit: MAX_OPEN_RING_LABELS,
                });
            }
            open_rings.insert(ring_num, (current, ring_bond, slot));
            // Return: we opened a ring; partner not yet known.
            Ok(StereoEntry::PendingRing(slot))
        }
    }

    /// Parse a ring closure number (single digit, `%nn`, or SMILES+ `%(n)`), together with an
    /// optional `prefix_bond` that was already consumed before this call.
    fn parse_ring_num(
        &mut self,
        prefix_bond: Option<ParsedBond>,
    ) -> Result<(u32, Option<ParsedBond>), SmilesError> {
        let ring_num = if self.peek() == Some(b'%') {
            self.advance(); // consume '%'
            if self.peek() == Some(b'(') {
                self.advance(); // consume '('
                let start = self.pos;
                let mut value = 0_u32;
                let mut digits = 0_usize;
                while let Some(byte) = self.peek() {
                    if byte == b')' {
                        break;
                    }
                    if !byte.is_ascii_digit() {
                        return Err(SmilesError::UnexpectedCharacter { pos: self.pos });
                    }
                    digits += 1;
                    value = value
                        .checked_mul(10)
                        .and_then(|prior| prior.checked_add((byte - b'0') as u32))
                        .ok_or(SmilesError::ResourceLimit {
                            resource: "ring label value",
                            actual: digits,
                            limit: u32::MAX.ilog10() as usize + 1,
                        })?;
                    self.advance();
                }
                if digits == 0 || self.advance() != Some(b')') {
                    return Err(SmilesError::UnexpectedEnd { pos: self.pos });
                }
                if value < 100 {
                    return Err(SmilesError::UnexpectedCharacter { pos: start });
                }
                value
            } else {
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
                u32::from(tens) * 10 + u32::from(units)
            }
        } else {
            // Caller peeked b'0'..=b'9', so advance() is guaranteed to return Some(digit).
            u32::from(
                self.advance()
                    .expect("ring closure digit guaranteed by caller peek")
                    - b'0',
            )
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
            // OpenSMILES: an unbracketed `*` is a wildcard atom with no
            // hydrogens, the same atom as `[*]` (RDKit writes it bare).
            Some(b'*') => {
                self.advance();
                let mut atom = Atom::wildcard();
                atom.chirality = self.parse_chirality(false)?;
                Ok(Some(atom))
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
            self.advance();
            let mut wc = Atom::wildcard();
            wc.chirality = chirality;
            // Isotope, H count and charge are kept as RDKit keeps them on a
            // dummy atom (`[*-]`, `[13*]`); they used to be dropped.
            wc.isotope = isotope;
            wc.charge = charge;
            if hcount > 0 {
                wc.hydrogen_count = Some(hcount);
            }
            // Wildcard atoms participate in atom-mapped reactions and in
            // RDKit CXSMILES attachment labels (`[*:7] |$_AP1;$|`). Dropping
            // the map here made the label's atom identity untraceable after a
            // CheMatic round trip, unlike every non-wildcard bracket atom.
            wc.atom_map = atom_map;
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
            // RDKit's grammar: a class without a number is permutation 0
            // (kept as "unknown"), an explicit 0 is an error; `@TH1`/`@TH`
            // are `@`, `@TH2` is `@@`. Allene classes are normalized after parsing.
            let number = self.parse_leading_digits_u16();
            let digits = number.unwrap_or(0);
            return match (class, number) {
                ("TH", None | Some(1)) => Ok(Chirality::CounterClockwise),
                ("TH", Some(2)) => Ok(Chirality::Clockwise),
                ("AL", None | Some(1)) => Ok(Chirality::Allene1),
                ("AL", Some(2)) => Ok(Chirality::Allene2),
                ("SP", None) => Ok(Chirality::SquarePlanarUnnumbered),
                ("SP", Some(1)) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP1)),
                ("SP", Some(2)) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP2)),
                ("SP", Some(3)) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP3)),
                ("TB", None) => Ok(Chirality::TrigonalBipyramidal(0)),
                ("TB", Some(n @ 1..=20)) => Ok(Chirality::TrigonalBipyramidal(n as u8)),
                ("OH", None) => Ok(Chirality::Octahedral(0)),
                ("OH", Some(n @ 1..=30)) => Ok(Chirality::Octahedral(n as u8)),
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

    /// Up to two digits, as RDKit reads `[CH10]` and `[C+10]` (OpenSMILES
    /// allows two-digit charges; a one-digit reading rejected the rest of
    /// the bracket atom).
    fn parse_two_digits(&mut self) -> Option<u8> {
        let first = self.peek().filter(|c| c.is_ascii_digit())?;
        self.advance();
        let mut value = first - b'0';
        if let Some(d) = self.peek().filter(|c| c.is_ascii_digit()) {
            self.advance();
            value = value * 10 + (d - b'0');
        }
        Some(value)
    }

    fn parse_hcount(&mut self) -> u8 {
        if self.peek() == Some(b'H') {
            self.advance();
            self.parse_two_digits().unwrap_or(1)
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
                if let Some(n) = self.parse_two_digits() {
                    return n as i8;
                }
                1
            }
            Some(b'-') => {
                self.advance();
                if self.peek() == Some(b'-') {
                    self.advance();
                    return -2;
                }
                if let Some(n) = self.parse_two_digits() {
                    return -(n as i8);
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
fn ring_closure_unmatched_error(ring_num: u32, pos: usize) -> SmilesError {
    match u8::try_from(ring_num) {
        Ok(ring_num) => SmilesError::UnmatchedRingClosure { ring_num, pos },
        Err(_) => SmilesError::UnmatchedExtendedRingClosure { ring_num, pos },
    }
}

fn ring_closure_conflict_error(ring_num: u32, pos: usize) -> SmilesError {
    match u8::try_from(ring_num) {
        Ok(ring_num) => SmilesError::ConflictingRingBond { ring_num, pos },
        Err(_) => SmilesError::ConflictingExtendedRingBond { ring_num, pos },
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

/// A wildcard between aromatic atoms on a 4n+2 ring is aromatic in RDKit's
/// SMILES model, including when its two written bonds were explicit single
/// bonds. The ring-size gate avoids promoting seven- and eight-membered rings
/// that RDKit de-aromatizes. A branch or inter-ring bridge stays aliphatic.
fn normalize_ring_wildcards(mol: &mut Molecule) {
    if mol.smiles_ring_closure_bond_count() == 0 || !mol.atoms().any(|(_, atom)| atom.wildcard) {
        return;
    }
    // Keep perception's type-erased caches off the returned parser graph.
    // Downstream callers may prepare that graph with a different perception
    // build (for example a crate's test dependency graph).
    let ring_view = mol.clone();
    let ring_bonds = chematic_perception::ring_bond_flags_shared(&ring_view);
    let rings = chematic_perception::find_sssr_shared(&ring_view);
    let mut changes = Vec::new();
    for (idx, atom) in mol.atoms() {
        if !atom.wildcard
            || atom.aromatic
            || !rings
                .rings()
                .iter()
                .any(|ring| ring.len() % 4 == 2 && ring.contains(&idx))
        {
            continue;
        }
        let aromatic_ring_bonds: Vec<_> = mol
            .neighbors(idx)
            .filter(|&(neighbor, bond)| ring_bonds[bond.0 as usize] && mol.atom(neighbor).aromatic)
            .map(|(_, bond)| bond)
            .collect();
        if aromatic_ring_bonds.len() >= 2 {
            changes.push((idx, aromatic_ring_bonds));
        }
    }
    for (idx, bonds) in changes {
        mol.set_atom_aromatic(idx, true);
        for bond in bonds {
            mol.set_bond_order(bond, BondOrder::Aromatic);
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
#[path = "parser/tests.rs"]
mod tests;

#[cfg(test)]
mod parse_error_boundary_contract_tests {
    use super::*;

    #[test]
    fn actual_parse_failures_keep_positions_and_relevant_labels() {
        for (source, context) in [
            ("C%", "unexpected end"),
            ("[Xx]", "missing element symbol"),
            ("C%(256)", "extended ring closure 256"),
            (
                "C/%(256)CC/%(256)",
                "conflicting bond types for extended ring closure 256",
            ),
            ("C=1CC-1", "conflicting bond types for ring closure 1"),
            ("[Pt@SP4](F)(Cl)(Br)I", "unsupported chirality class '@SP4'"),
        ] {
            let error = parse(source).err().unwrap();
            let message = error.to_string();
            assert!(message.contains(context), "{source}: {message}");
            assert!(message.contains("position"), "{source}: {message}");
        }
        let empty = parse("").err().unwrap();
        assert_eq!(empty, SmilesError::EmptyInput);
        assert!(empty.to_string().contains("input is empty"));
        let depth = MAX_BRANCH_DEPTH + 1;
        let source = format!("C{}{}", "(C".repeat(depth), ")".repeat(depth));
        let error = parse(&source).err().unwrap();
        assert!(matches!(error, SmilesError::NestingTooDeep { .. }));
        assert!(
            error
                .to_string()
                .contains("branch nesting too deep at position")
        );
    }
}
