//! SMILES writer: serialize a Molecule to an OpenSMILES string via DFS traversal.
//!
//! Produces valid (non-canonical) SMILES. Canonical SMILES (Morgan-rank ordering)
//! is a planned future milestone.

use std::collections::{HashMap, HashSet};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL, STEREO_VACANT_SENTINEL,
    remap_square_planar_tag, remap_tetrahedral_parity,
};

/// Write the `H`/`Hn` token for a bracket atom's hydrogen count.
///
/// Uses `implicit_hcount` rather than `atom.hydrogen_count` directly so that
/// atoms forced into bracket notation by isotope/charge/atom-map (but with no
/// explicit H count recorded) still get their inferred hydrogens written,
/// e.g. `[NH4+]` rather than `[N+]`.
pub(crate) fn emit_bracket_hydrogens(out: &mut String, mol: &Molecule, idx: AtomIdx) {
    let atom = mol.atom(idx);
    let h = if atom.wildcard {
        atom.hydrogen_count.unwrap_or(0)
    } else {
        chematic_core::implicit_hcount(mol, idx)
    };
    if h > 0 {
        out.push('H');
        if h > 1 {
            out.push_str(&h.to_string());
        }
    }
}

/// True when bond `bidx` flanks a double bond whose *other* end also has a
/// directional substituent bond, so a `/`/`\\` token on `bidx` takes part in
/// an E/Z specification. A lone marker -- e.g. a MOL wedge drawn from a
/// stereocentre to an alkene carbon -- specifies nothing, and RDKit drops it
/// on parse.
pub(crate) fn flanks_marked_double_bond(mol: &Molecule, bidx: BondIdx) -> bool {
    let bond = mol.bond(bidx);
    [bond.atom1, bond.atom2].into_iter().any(|end| {
        mol.neighbors(end).any(|(far, db)| {
            db != bidx
                && mol.bond(db).order == BondOrder::Double
                && mol
                    .neighbors(far)
                    .any(|(_, sb)| sb != db && sb != bidx && raw_bond_direction(mol, sb).is_some())
                && !in_ring_smaller_than_eight(mol, db)
        })
    })
}

/// Whether `bond` lies on a ring of fewer than eight atoms, where a double
/// bond has no E/Z (RDKit drops such markers).
pub(crate) fn in_ring_smaller_than_eight(mol: &Molecule, bond: BondIdx) -> bool {
    let (start, goal) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    let mut seen = vec![false; mol.atom_count()];
    seen[start.0 as usize] = true;
    let mut frontier = vec![start];
    for _ in 0..6 {
        let mut next = Vec::new();
        for &a in &frontier {
            for (nb, bi) in mol.neighbors(a) {
                if bi == bond {
                    continue;
                }
                if nb == goal {
                    return true;
                }
                if !seen[nb.0 as usize] {
                    seen[nb.0 as usize] = true;
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    false
}

/// Demote a `BondOrder::Up`/`Down` bond to `Single` for writer-emission
/// purposes when it does not flank a double bond.
///
/// `Up`/`Down` is overloaded in `chematic_core` for two unrelated concepts:
/// 2D wedge/hash depiction (tetrahedral stereocenter drawing, set by
/// MDL/CDXML/MRV/KET-style readers) and the OpenSMILES `/`/`\` E/Z
/// directional-bond marker (only meaningful adjacent to a double bond). A
/// wedge/hash bond with no adjacent double bond has nothing to mark; writing
/// `/`/`\` for it is a meaningless token that any spec-compliant SMILES
/// parser (e.g. RDKit) silently drops on re-parse. This only changes what
/// the writer decides to print — the bond's real `order` field on the
/// `Molecule` is never touched, so wedge/hash depiction data survives.
///
/// Gates on the bond's *stored* order (`mol.bond(bidx).order`), not the
/// `order` parameter, so a SMILES-parser E/Z direction stashed on an
/// `Aromatic`-order bond (`Molecule::bond_direction`, a separate mechanism
/// entirely unrelated to 2D wedge readers) is never affected by this
/// suppression — its stored order is `Aromatic`, never a literal
/// `Up`/`Down`, so it can never match the guard below.
pub(crate) fn suppress_standalone_wedge(
    mol: &Molecule,
    bidx: BondIdx,
    order: BondOrder,
) -> BondOrder {
    if matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
        && !flanks_marked_double_bond(mol, bidx)
    {
        BondOrder::Single
    } else {
        order
    }
}

/// The "raw" (atom1→atom2-relative, unreoriented) directional marker for
/// bond `bidx`, if any: a literal `Up`/`Down` bond order, or -- when the
/// bond's real order was overwritten to something else (e.g. `Aromatic`, for
/// a ring bond stashing an adjacent exocyclic double bond's direction, or a
/// reader-perceived 2D E/Z direction stashed on a plain `Single` bond) --
/// whatever [`chematic_core::Molecule::bond_direction`] holds for it.
///
/// Mirrors `crate::canonical::CanonicalWriter::raw_input_direction` exactly
/// -- both the plain and canonical writers must read the same effective
/// direction (docs/rfcs/stereo2d_reader_integration_rfc.md), so this crate keeps
/// exactly one copy of the rule rather than two that could silently drift.
pub(crate) fn raw_bond_direction(mol: &Molecule, bidx: BondIdx) -> Option<BondOrder> {
    // A stored direction wins over a literal `Up`/`Down` order: a MOL
    // reader keeps a wedge's tetrahedral meaning in `order` and may put the
    // E/Z direction of an adjacent double bond on the same bond.
    if let Some(direction) = mol.bond_direction(bidx) {
        return Some(direction);
    }
    let order = mol.bond(bidx).order;
    matches!(order, BondOrder::Up | BondOrder::Down).then_some(order)
}

/// Re-orient a raw (atom1→atom2-relative) directional marker for reading
/// "from `from_atom` toward the bond's other endpoint".
///
/// OpenSMILES `/`/`\` is relative to *written text order*, not a bond's
/// internal atom1/atom2 storage order: a bond visited in the direction
/// opposite to how its stored marker was oriented must flip, or the printed
/// character encodes the wrong geometry. A DFS write can visit either
/// endpoint first regardless of which one happens to be `atom1`, so this
/// must be applied at every emission site, not just assumed consistent.
/// Mirrors the per-site match blocks in `crate::canonical`'s
/// `dfs_mark`/`write_chain`.
pub(crate) fn direction_from(
    mol: &Molecule,
    bidx: BondIdx,
    dir: BondOrder,
    from_atom: AtomIdx,
) -> BondOrder {
    let bond = mol.bond(bidx);
    // A parser-side aromatic direction stash is relative to the endpoint
    // where the `/` or `\\` token was read, not necessarily to the bond's
    // internal `atom1`.  Use the preserved anchor whenever one exists so a
    // canonical DFS cannot reinterpret the same carrier merely because the
    // molecule was rebuilt with the bond endpoints reversed.
    let direction_atom = mol.bond_direction_anchor(bidx).unwrap_or(bond.atom1);
    match dir {
        BondOrder::Up => {
            if direction_atom == from_atom {
                BondOrder::Up
            } else {
                BondOrder::Down
            }
        }
        BondOrder::Down => {
            if direction_atom == from_atom {
                BondOrder::Down
            } else {
                BondOrder::Up
            }
        }
        other => other,
    }
}

/// The SMILES token to print for bond `bidx` when the DFS is writing it
/// *starting from* `from_atom`.
///
/// `Dative` is the only bond order whose token is longer than one character
/// (`"->"`), and the only one whose *text* encodes an asymmetric fact about
/// its two endpoints: `BondOrder::Dative` stores donor→acceptor as
/// `atom1`→`atom2`, while `->`/`<-` are read left-to-right in written order.
/// A DFS can reach either endpoint first, so writing `"->"` unconditionally
/// would assert "the atom I just wrote is the donor" — false, and actively
/// misleading, whenever the traversal arrived at `atom2` first. Every other
/// order's token is direction-free, so this is exactly
/// [`BondOrder::smiles_token`] for them.
///
/// This is the dative analogue of [`direction_from`], which does the same
/// job for `/`/`\`. It has to live here rather than in `direction_from`
/// because there is only one `Dative` variant — the flip has no `BondOrder`
/// to be expressed as, only a different token string.
pub(crate) fn bond_token_from(
    mol: &Molecule,
    bidx: BondIdx,
    order: BondOrder,
    from_atom: AtomIdx,
) -> &'static str {
    if order == BondOrder::Dative && mol.bond(bidx).atom1 != from_atom {
        "<-"
    } else {
        order.smiles_token()
    }
}

/// SMILES text token for a square-planar permutation tag (`@SP1`/`@SP2`/`@SP3`,
/// including the leading `@`).
pub(crate) fn square_planar_token(p: chematic_core::SquarePlanarPermutation) -> &'static str {
    use chematic_core::SquarePlanarPermutation::*;
    match p {
        SP1 => "@SP1",
        SP2 => "@SP2",
        SP3 => "@SP3",
    }
}

/// Re-express a square-planar tag with one vacant ligand position. RDKit's
/// three-neighbor SMILES convention puts the vacancy in its own slot; applying
/// the four-neighbor table to only the three explicit neighbors changes the
/// configuration. This table is the 18-case RDKit 2025.09.3 permutation oracle
/// for three distinct ligands, indexed by the written order relative to the
/// stored order. Refuse missing or duplicate ligand ids.
pub(crate) fn remap_three_ligand_square_planar(
    tag: chematic_core::SquarePlanarPermutation,
    original: [u32; 3],
    written: [u32; 3],
) -> Option<chematic_core::SquarePlanarPermutation> {
    use chematic_core::SquarePlanarPermutation::{SP1, SP2, SP3};
    if original[0] == original[1] || original[0] == original[2] || original[1] == original[2] {
        return None;
    }
    let mut permutation = [0usize; 3];
    for (slot, ligand) in written.into_iter().enumerate() {
        permutation[slot] = original.iter().position(|&id| id == ligand)?;
    }
    if permutation[0] == permutation[1]
        || permutation[0] == permutation[2]
        || permutation[1] == permutation[2]
    {
        return None;
    }
    let mapped = match permutation {
        [0, 1, 2] => [SP1, SP2, SP3],
        [0, 2, 1] => [SP3, SP2, SP1],
        [1, 0, 2] => [SP1, SP3, SP2],
        [1, 2, 0] => [SP3, SP1, SP2],
        [2, 0, 1] => [SP2, SP3, SP1],
        [2, 1, 0] => [SP2, SP1, SP3],
        _ => unreachable!("validated three-element permutation"),
    };
    Some(
        mapped[match tag {
            SP1 => 0,
            SP2 => 1,
            SP3 => 2,
        }],
    )
}

/// Whether the recorded coordination has exactly two unoccupied positions.
pub(crate) fn has_two_vacant_sites(order: &[u32]) -> bool {
    order
        .iter()
        .filter(|&&id| id == STEREO_VACANT_SENTINEL)
        .count()
        == 2
}

/// Carry a non-tetrahedral tag from the stored ligand order to the written
/// order with RDKit's permutation tables. Refuse mismatched ligand sets.
pub(crate) fn remap_nontetrahedral(
    stored: chematic_core::Chirality,
    original: &[u32],
    written: &[u32],
) -> chematic_core::Chirality {
    let Some((class, perm)) = stored.nontetrahedral() else {
        return chematic_core::Chirality::None;
    };
    let two_vacancies = has_two_vacant_sites(original) && has_two_vacant_sites(written);
    let root_changes = (original.first() == Some(&STEREO_VACANT_SENTINEL))
        != (written.first() == Some(&STEREO_VACANT_SENTINEL));
    // The two empty positions are indistinguishable in the molecular graph,
    // but RDKit's SMILES permutation convention exchanges their slots when
    // a ligand moves between the incoming and first-child positions.
    let order = |ids: &[u32], reverse_vacancies: bool| {
        let mut vacancy = 0;
        ids.iter()
            .map(|&id| {
                if two_vacancies && id == STEREO_VACANT_SENTINEL {
                    vacancy += 1;
                    if (vacancy == 1) != reverse_vacancies {
                        None
                    } else {
                        Some(id as usize)
                    }
                } else {
                    Some(id as usize)
                }
            })
            .collect::<Vec<Option<usize>>>()
    };
    let o = order(original, false);
    let w = order(written, two_vacancies && root_changes);
    chematic_core::nontetrahedral::permute(class, perm, &o, &w)
        .and_then(|p| chematic_core::Chirality::from_nontetrahedral(class, p))
        .unwrap_or(chematic_core::Chirality::None)
}

/// Write a [`Molecule`] to a SMILES string.
///
/// Disconnected fragments are joined with `.`.
/// Aromatic atoms are written in lowercase.
pub fn write(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }
    SmilesWriter::new(mol).write_all().0
}

/// Write a [`Molecule`] to SMILES and report the order in which its atoms
/// appear in the string.
///
/// The string is exactly [`write`]'s. `order[k]` is the index in `mol` of the
/// `k`-th atom written, which is also the index that atom receives when the
/// string is parsed back with [`crate::parse`]: for every `k`,
/// `parse(&smiles)?.atom(AtomIdx(k))` corresponds to `mol.atom(order[k])`.
/// Every atom of `mol` is written exactly once, so `order` is a permutation
/// of `0..mol.atom_count()`. This is the DFS visit order (RDKit's
/// `_smilesAtomOutputOrder`), not a canonical rank.
pub fn write_with_atom_order(mol: &Molecule) -> (String, Vec<AtomIdx>) {
    if mol.atom_count() == 0 {
        return (String::new(), Vec::new());
    }
    SmilesWriter::new(mol).write_all()
}

struct SmilesWriter<'a> {
    mol: &'a Molecule,
    /// Bonds that are back-edges in the DFS tree (ring closures).
    ring_bonds: HashSet<BondIdx>,
    /// ring number(s) each atom must write when serialized.
    /// Both the "open" ancestor and "close" descendant of a ring store the same number.
    /// `BondIdx` is kept alongside the order so the emission site can ask
    /// [`bond_token_from`] for a direction-correct token (dative arrows).
    atom_ring_nums: HashMap<AtomIdx, Vec<(u16, BondOrder, BondIdx)>>,
    /// Whether each atom has been serialized in phase 2.
    written: Vec<bool>,
    next_ring: u16,
    out: String,
    /// Atoms in the order they were emitted.
    order: Vec<AtomIdx>,
}

impl<'a> SmilesWriter<'a> {
    fn new(mol: &'a Molecule) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            ring_bonds: HashSet::new(),
            atom_ring_nums: HashMap::new(),
            written: vec![false; n],
            next_ring: 1,
            out: String::new(),
            order: Vec::with_capacity(n),
        }
    }

    fn write_all(mut self) -> (String, Vec<AtomIdx>) {
        // Phase 1: find all back-edges and assign ring-closure numbers.
        self.find_ring_closures();

        // Phase 2: DFS serialization, one fragment at a time.
        let mut first = true;
        for i in 0..self.mol.atom_count() {
            if !self.written[i] {
                if !first {
                    self.out.push('.');
                }
                first = false;
                self.write_chain(AtomIdx(i as u32), None, None);
            }
        }

        (self.out, self.order)
    }

    fn find_ring_closures(&mut self) {
        let n = self.mol.atom_count();
        let mut visited = vec![false; n];
        let mut in_stack = vec![false; n];

        for start in 0..n {
            if !visited[start] {
                self.dfs_mark(AtomIdx(start as u32), None, &mut visited, &mut in_stack);
            }
        }
    }

    /// DFS that marks back-edges and assigns ring-closure numbers.
    ///
    /// `from_bond`: the bond used to arrive at `atom` (skip it to avoid re-visiting the parent).
    fn dfs_mark(
        &mut self,
        atom: AtomIdx,
        from_bond: Option<BondIdx>,
        visited: &mut Vec<bool>,
        in_stack: &mut Vec<bool>,
    ) {
        visited[atom.0 as usize] = true;
        in_stack[atom.0 as usize] = true;

        for (neighbor, bidx) in self.mol.neighbors(atom) {
            // Skip the edge we came from (undirected graph: would look like a back-edge otherwise).
            if Some(bidx) == from_bond {
                continue;
            }
            // Skip bonds already classified.
            if self.ring_bonds.contains(&bidx) {
                continue;
            }

            if !visited[neighbor.0 as usize] {
                // Tree edge: recurse.
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                // Back-edge: `atom` (descendant) closes a ring back to `neighbor` (ancestor).
                self.ring_bonds.insert(bidx);
                let rn = self.next_ring;
                self.next_ring += 1;
                // Resolve the effective (literal-or-stashed) direction once,
                // then re-orient it separately for each endpoint -- the two
                // sides of a ring closure are visited from opposite ends, so
                // a single un-reoriented value would print the same
                // character on both, backwards on whichever side isn't
                // `bond.atom1`.
                let raw = raw_bond_direction(self.mol, bidx).unwrap_or(self.mol.bond(bidx).order);
                let order_at_open = suppress_standalone_wedge(
                    self.mol,
                    bidx,
                    direction_from(self.mol, bidx, raw, neighbor),
                );
                let order_at_close = suppress_standalone_wedge(
                    self.mol,
                    bidx,
                    direction_from(self.mol, bidx, raw, atom),
                );
                // Both endpoints need to emit this ring number when serialized.
                self.atom_ring_nums
                    .entry(neighbor)
                    .or_default()
                    .push((rn, order_at_open, bidx)); // open
                self.atom_ring_nums
                    .entry(atom)
                    .or_default()
                    .push((rn, order_at_close, bidx)); // close
            }
        }

        in_stack[atom.0 as usize] = false;
    }

    /// Write `atom` and then recurse into its unvisited tree-edge neighbors.
    /// `incoming_bond`: the already-oriented token for the edge leading to
    /// this atom (None for the root, or when the bond is implicit). It is a
    /// token rather than a `BondOrder` because a dative arrow's direction
    /// depends on which endpoint is written first — see [`bond_token_from`].
    fn write_chain(
        &mut self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        incoming_bond: Option<&'static str>,
    ) {
        self.written[atom.0 as usize] = true;

        // Write the incoming bond (if explicit / non-default).
        if let Some(token) = incoming_bond {
            self.out.push_str(token);
        }

        // Write the atom symbol.
        self.order.push(atom);
        self.emit_atom(atom, from_atom);

        // Write ring-closure digits for this atom (both open and close digits).
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, bond_order, bidx) in rings {
                // Write bond type unless the parser would infer it from the
                // two endpoints: a bare digit reads back as aromatic between
                // two aromatic atoms and as single otherwise (the same rule
                // as tree edges, and as the canonical writer since #395).
                let bond = self.mol.bond(bidx);
                let partner = if bond.atom1 == atom {
                    bond.atom2
                } else {
                    bond.atom1
                };
                let both_aromatic = self.mol.atom(atom).aromatic && self.mol.atom(partner).aromatic;
                let implicit = match bond_order {
                    BondOrder::Single => !both_aromatic,
                    BondOrder::Aromatic => both_aromatic,
                    _ => false,
                };
                if !implicit {
                    // Oriented from the atom being written right now: the two
                    // ends of a dative ring closure print opposite arrows
                    // (`->` at the donor, `<-` at the acceptor), which is the
                    // same bond read from opposite directions -- exactly how
                    // `/`/`\` ring-closure markers already behave here.
                    self.out
                        .push_str(bond_token_from(self.mol, bidx, bond_order, atom));
                }
                // Ring number: single digit for 1-9, `%NN` form for 10-99,
                // and SMILES+ `%(n)` for 100+.
                if rn >= 100 {
                    self.out.push_str("%(");
                    self.out.push_str(&rn.to_string());
                    self.out.push(')');
                } else if rn >= 10 {
                    self.out.push('%');
                    self.out
                        .push(char::from_digit((rn / 10) as u32, 10).unwrap());
                    self.out
                        .push(char::from_digit((rn % 10) as u32, 10).unwrap());
                } else {
                    self.out.push(char::from_digit(rn as u32, 10).unwrap());
                }
            }
        }

        // Collect tree-edge children (unvisited, non-ring-closure bonds).
        let children: Vec<(AtomIdx, BondIdx, BondOrder)> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, bidx)| {
                // Direction seen from `atom` going toward `nb` -- re-oriented
                // from the bond's raw (atom1→atom2-relative) marker so a
                // literal/stashed direction prints correctly regardless of
                // which endpoint the DFS happens to write first.
                let raw = raw_bond_direction(self.mol, bidx).unwrap_or(self.mol.bond(bidx).order);
                let oriented = direction_from(self.mol, bidx, raw, atom);
                (
                    nb,
                    bidx,
                    suppress_standalone_wedge(self.mol, bidx, oriented),
                )
            })
            .collect();

        // Write children: all but the last one are branches (wrapped in parentheses).
        let n = children.len();
        for (i, (child, bidx, bond_order)) in children.into_iter().enumerate() {
            let is_last = i == n - 1;

            // Determine whether the bond should be written explicitly.
            let parent_arom = self.mol.atom(atom).aromatic;
            let child_arom = self.mol.atom(child).aromatic;
            let implicit = match bond_order {
                BondOrder::Single => !(parent_arom && child_arom), // single is implicit
                BondOrder::Aromatic => parent_arom && child_arom,  // aromatic is implicit
                _ => false,
            };
            let written_bond = if implicit {
                None
            } else {
                Some(bond_token_from(self.mol, bidx, bond_order, atom))
            };

            if !is_last {
                self.out.push('(');
                self.write_chain(child, Some(atom), written_bond);
                self.out.push(')');
            } else {
                self.write_chain(child, Some(atom), written_bond);
            }
        }
    }

    /// The atom's stereo tag re-expressed for the neighbour order this
    /// writer emits: the incoming atom, an implicit/bracket H, ring-closure
    /// partners in digit order, then the tree children. A stored `@`/`@@` is
    /// relative to [`Molecule::stereo_neighbor_order`]; when the atom has
    /// none, the tag is written as stored.
    fn written_chirality(&self, idx: AtomIdx, from_atom: Option<AtomIdx>) -> Chirality {
        let stored = self.mol.atom(idx).chirality;
        if stored == Chirality::None {
            return stored;
        }
        let Some(original) = self.mol.stereo_neighbor_order(idx) else {
            return stored;
        };
        // The non-canonical traversal preserves the parse-order atom layout.
        // Canonical serialization performs the full virtual-terminal parity
        // remap; retaining the tag here is safer than dropping axial stereo.
        if self.mol.is_allene_stereo(idx) {
            return stored;
        }
        let mut written: Vec<u32> = Vec::with_capacity(original.len());
        if let Some(prev) = from_atom {
            written.push(prev.0);
        }
        written.extend(
            original
                .iter()
                .filter(|&&id| id == STEREO_VACANT_SENTINEL)
                .copied(),
        );
        // Three-coordinate sulfoxides use the H sentinel for the stereogenic
        // lone pair, even though the atom has no implicit hydrogen.
        if original.contains(&STEREO_H_SENTINEL)
            || chematic_core::implicit_hcount(self.mol, idx) > 0
        {
            written.push(STEREO_H_SENTINEL);
        }
        if let Some(rings) = self.atom_ring_nums.get(&idx) {
            for &(_, _, bidx) in rings {
                let bond = self.mol.bond(bidx);
                let partner = if bond.atom1 == idx {
                    bond.atom2
                } else {
                    bond.atom1
                };
                written.push(partner.0);
            }
        }
        for (nb, bidx) in self.mol.neighbors(idx) {
            if Some(nb) != from_atom
                && !self.written[nb.0 as usize]
                && !self.ring_bonds.contains(&bidx)
            {
                written.push(nb.0);
            }
        }
        if written.len() != original.len() {
            return if stored.is_tetrahedral() {
                stored
            } else {
                Chirality::None
            };
        }
        let original_arr = <[u32; 4]>::try_from(original).ok();
        let written_arr = <[u32; 4]>::try_from(written.as_slice()).ok();
        match stored {
            Chirality::CounterClockwise | Chirality::Clockwise => {
                let odd = match (original_arr, written_arr) {
                    (Some(o), Some(w)) => remap_tetrahedral_parity(o, w).unwrap_or(false),
                    _ => permutation_is_odd(original, &written),
                };
                match (odd, stored) {
                    (false, tag) => tag,
                    (true, Chirality::CounterClockwise) => Chirality::Clockwise,
                    (true, _) => Chirality::CounterClockwise,
                }
            }
            Chirality::SquarePlanar(_) if has_two_vacant_sites(original) => {
                remap_nontetrahedral(stored, original, &written)
            }
            Chirality::SquarePlanar(tag) => match (original_arr, written_arr) {
                (Some(o), Some(w)) => remap_square_planar_tag(tag, o, w)
                    .map(Chirality::SquarePlanar)
                    .unwrap_or(Chirality::None),
                _ => match (
                    <[u32; 3]>::try_from(original).ok(),
                    <[u32; 3]>::try_from(written.as_slice()).ok(),
                ) {
                    (Some(o), Some(w)) => remap_three_ligand_square_planar(tag, o, w)
                        .map(Chirality::SquarePlanar)
                        .unwrap_or(Chirality::None),
                    _ => Chirality::None,
                },
            },
            Chirality::None => Chirality::None,
            other => remap_nontetrahedral(other, original, &written),
        }
    }

    fn emit_atom(&mut self, idx: AtomIdx, from_atom: Option<AtomIdx>) {
        let atom = self.mol.atom(idx);
        let chirality = self.written_chirality(idx, from_atom);

        // An atom needs bracket notation when:
        //  - it has an isotope, charge, explicit H count, atom map, or
        //    a stereo tag, or
        //  - it is not in the organic subset (cannot rely on implicit-H rules).
        let needs_bracket = atom.wildcard
            || atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogen_count.is_some()
            || chirality != Chirality::None
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some();

        if needs_bracket {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            let sym = if atom.wildcard {
                "*".to_string()
            } else if atom.aromatic {
                atom.element.symbol().to_lowercase()
            } else {
                atom.element.symbol().to_string()
            };
            self.out.push_str(&sym);

            match chirality {
                chematic_core::Chirality::CounterClockwise if self.mol.is_allene_stereo(idx) => {
                    self.out.push_str("@AL1")
                }
                chematic_core::Chirality::Clockwise if self.mol.is_allene_stereo(idx) => {
                    self.out.push_str("@AL2")
                }
                chematic_core::Chirality::CounterClockwise => self.out.push('@'),
                chematic_core::Chirality::Clockwise => self.out.push_str("@@"),
                chematic_core::Chirality::None => {}
                chematic_core::Chirality::SquarePlanar(p) => {
                    self.out.push_str(square_planar_token(p))
                }
                other => self.out.push_str(&other.smiles_token()),
            }

            emit_bracket_hydrogens(&mut self.out, self.mol, idx);

            match atom.charge {
                0 => {}
                1 => self.out.push('+'),
                -1 => self.out.push('-'),
                c if c > 0 => self.out.push_str(&format!("+{c}")),
                c => self.out.push_str(&c.to_string()),
            }

            if let Some(m) = atom.atom_map {
                self.out.push(':');
                self.out.push_str(&m.to_string());
            }

            self.out.push(']');
        } else if atom.aromatic {
            self.out.push_str(&atom.element.symbol().to_lowercase());
        } else {
            self.out.push_str(atom.element.symbol());
        }
    }
}

/// Whether reordering `original` into `written` is an odd permutation (both
/// name the same distinct ids).
/// Not a permutation (a different id set or a repeated id): `false`, the
/// pass-through used for every unverifiable order.
fn permutation_is_odd(original: &[u32], written: &[u32]) -> bool {
    let Some(mut perm) = written
        .iter()
        .map(|w| original.iter().position(|o| o == w))
        .collect::<Option<Vec<usize>>>()
    else {
        return false;
    };
    let mut seen = vec![false; perm.len()];
    for &p in &perm {
        if p >= seen.len() || std::mem::replace(&mut seen[p], true) {
            return false;
        }
    }
    let mut swaps = 0;
    for i in 0..perm.len() {
        while perm[i] != i {
            let j = perm[i];
            perm.swap(i, j);
            swaps += 1;
        }
    }
    swaps % 2 == 1
}

#[cfg(test)]
mod two_vacant_tests;

#[cfg(test)]
#[path = "writer/tests.rs"]
mod tests;
