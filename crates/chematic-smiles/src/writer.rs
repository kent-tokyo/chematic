//! SMILES writer: serialize a Molecule to an OpenSMILES string via DFS traversal.
//!
//! Produces valid (non-canonical) SMILES. Canonical SMILES (Morgan-rank ordering)
//! is a planned future milestone.

use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// Write the `H`/`Hn` token for a bracket atom's hydrogen count.
///
/// Uses `implicit_hcount` rather than `atom.hydrogen_count` directly so that
/// atoms forced into bracket notation by isotope/charge/atom-map (but with no
/// explicit H count recorded) still get their inferred hydrogens written,
/// e.g. `[NH4+]` rather than `[N+]`.
pub(crate) fn emit_bracket_hydrogens(out: &mut String, mol: &Molecule, idx: AtomIdx) {
    let h = chematic_core::implicit_hcount(mol, idx);
    if h > 0 {
        out.push('H');
        if h > 1 {
            out.push_str(&h.to_string());
        }
    }
}

/// True when bond `bidx` shares an atom with some *other* `BondOrder::Double`
/// bond — i.e. whether a `/`/`\` token on this bond could ever carry genuine
/// OpenSMILES E/Z meaning. Mirrors the double-bond-adjacency test
/// `CanonicalWriter::build_ez_groups` uses to collect real E/Z side bonds,
/// just queried per-bond instead of built as a whole-molecule group.
pub(crate) fn flanks_double_bond(mol: &Molecule, bidx: BondIdx) -> bool {
    let bond = mol.bond(bidx);
    [bond.atom1, bond.atom2].into_iter().any(|endpoint| {
        mol.neighbors(endpoint)
            .any(|(_, nb)| nb != bidx && mol.bond(nb).order == BondOrder::Double)
    })
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
        && !flanks_double_bond(mol, bidx)
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
    let order = mol.bond(bidx).order;
    if matches!(order, BondOrder::Up | BondOrder::Down) {
        return Some(order);
    }
    mol.bond_direction(bidx)
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
    match dir {
        BondOrder::Up => {
            if bond.atom1 == from_atom {
                BondOrder::Up
            } else {
                BondOrder::Down
            }
        }
        BondOrder::Down => {
            if bond.atom1 == from_atom {
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

/// Write a [`Molecule`] to a SMILES string.
///
/// Disconnected fragments are joined with `.`.
/// Aromatic atoms are written in lowercase.
pub fn write(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
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
        }
    }

    fn write_all(mut self) -> String {
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

        self.out
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
        self.emit_atom(atom);

        // Write ring-closure digits for this atom (both open and close digits).
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, bond_order, bidx) in rings {
                // Write bond type unless it is implicit.
                let atom_aromatic = self.mol.atom(atom).aromatic;
                // For ring closures we can't know the other atom's aromaticity here,
                // so we emit the bond type unless it is a plain aromatic ring bond.
                if !(bond_order == BondOrder::Aromatic && atom_aromatic)
                    && bond_order != BondOrder::Single
                {
                    // Oriented from the atom being written right now: the two
                    // ends of a dative ring closure print opposite arrows
                    // (`->` at the donor, `<-` at the acceptor), which is the
                    // same bond read from opposite directions -- exactly how
                    // `/`/`\` ring-closure markers already behave here.
                    self.out
                        .push_str(bond_token_from(self.mol, bidx, bond_order, atom));
                }
                // Ring number: single digit for 1-9, `%NN` form for 10-99, `%NNN` for 100+.
                if rn >= 100 {
                    self.out.push('%');
                    for ch in rn.to_string().chars() {
                        self.out.push(ch);
                    }
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

    fn emit_atom(&mut self, idx: AtomIdx) {
        let atom = self.mol.atom(idx);

        // An atom needs bracket notation when:
        //  - it has an isotope, charge, explicit H count, atom map, or
        //  - it is not in the organic subset (cannot rely on implicit-H rules).
        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogen_count.is_some()
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some();

        if needs_bracket {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            let sym = if atom.aromatic {
                atom.element.symbol().to_lowercase()
            } else {
                atom.element.symbol().to_string()
            };
            self.out.push_str(&sym);

            match atom.chirality {
                chematic_core::Chirality::CounterClockwise => self.out.push('@'),
                chematic_core::Chirality::Clockwise => self.out.push_str("@@"),
                chematic_core::Chirality::None => {}
                chematic_core::Chirality::SquarePlanar(p) => {
                    self.out.push_str(square_planar_token(p))
                }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Parse → write → re-parse → verify atom/bond counts are preserved.
    fn roundtrip(smiles: &str) {
        let mol1 = parse(smiles).expect(smiles);
        let out = write(&mol1);
        let mol2 = parse(&out).unwrap_or_else(|e| {
            panic!(
                "roundtrip failed for '{}': wrote '{}', error: {e}",
                smiles, out
            )
        });
        assert_eq!(
            mol1.atom_count(),
            mol2.atom_count(),
            "atom count mismatch: input='{}' output='{}'",
            smiles,
            out
        );
        assert_eq!(
            mol1.bond_count(),
            mol2.bond_count(),
            "bond count mismatch: input='{}' output='{}'",
            smiles,
            out
        );
    }

    #[test]
    fn test_write_methane() {
        assert_eq!(write(&parse("C").unwrap()), "C");
    }
    #[test]
    fn test_write_ethane() {
        assert_eq!(write(&parse("CC").unwrap()), "CC");
    }

    #[test]
    fn test_roundtrip_propane() {
        roundtrip("CCC");
    }
    #[test]
    fn test_roundtrip_isobutane() {
        roundtrip("CC(C)C");
    }
    #[test]
    fn test_roundtrip_ethanol() {
        roundtrip("CCO");
    }
    #[test]
    fn test_roundtrip_acetic_acid() {
        roundtrip("CC(=O)O");
    }
    #[test]
    fn test_roundtrip_cyclohexane() {
        roundtrip("C1CCCCC1");
    }
    #[test]
    fn test_roundtrip_benzene_kekule() {
        roundtrip("C1=CC=CC=C1");
    }
    #[test]
    fn test_roundtrip_benzene_arom() {
        roundtrip("c1ccccc1");
    }
    #[test]
    fn test_roundtrip_pyridine() {
        roundtrip("c1ccncc1");
    }
    #[test]
    fn test_roundtrip_naphthalene() {
        roundtrip("c1ccc2ccccc2c1");
    }
    #[test]
    fn test_roundtrip_chlorobenzene() {
        roundtrip("c1ccccc1Cl");
    }
    #[test]
    fn test_roundtrip_13c() {
        roundtrip("[13C]");
    }
    #[test]
    fn test_roundtrip_ammonium() {
        roundtrip("[NH4+]");
    }
    #[test]
    fn test_roundtrip_disconnected() {
        roundtrip("[Na+].[Cl-]");
    }
    #[test]
    fn test_roundtrip_aspirin() {
        roundtrip("CC(=O)Oc1ccccc1C(=O)O");
    }
    #[test]
    fn test_roundtrip_caffeine() {
        roundtrip("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
    }

    // Bracket atoms with `hydrogen_count: None` (implicit H left uninferred by the
    // caller, e.g. after a programmatic charge/isotope/atom-map edit) must still get
    // their implicit hydrogens written — regression tests for the bracket-H bug found
    // via MRV oracle validation (isotope_0/2/3, charge_0/3, atom_map_0/1/2, disconnected_3
    // in validation/mrv_io_parity_summary.json).
    use chematic_core::{Atom, Element, MoleculeBuilder};

    #[test]
    fn test_bracket_implicit_h_ammonium_charge_only() {
        let mut b = MoleculeBuilder::new();
        let mut n = Atom::new(Element::N);
        n.charge = 1;
        b.add_atom(n);
        assert_eq!(write(&b.build()), "[NH4+]");
    }

    #[test]
    fn test_bracket_implicit_h_isotope_only() {
        let mut b = MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.isotope = Some(13);
        b.add_atom(c);
        assert_eq!(write(&b.build()), "[13CH4]");
    }

    #[test]
    fn test_bracket_implicit_h_atom_map_only() {
        let mut b = MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(write(&b.build()), "[CH3:7]C");
    }

    #[test]
    fn test_bracket_implicit_h_isotope_and_atom_map() {
        let mut b = MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.isotope = Some(13);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(write(&b.build()), "[13CH3:7]C");
    }

    #[test]
    fn test_bracket_implicit_h_hydroxide_charge_only() {
        let mut b = MoleculeBuilder::new();
        let mut o = Atom::new(Element::O);
        o.charge = -1;
        b.add_atom(o);
        assert_eq!(write(&b.build()), "[OH-]");
    }

    // ── Standalone wedge/hash bond must not emit a meaningless SMILES
    // directional token (docs/rfcs/stereo2d_reader_integration_rfc.md §3/§7,
    // docs/rfcs/stereo2d_local_parity_calibration.md "Scope note") ────────────────
    //
    // `BondOrder::Up`/`Down` is overloaded for two unrelated concepts: 2D
    // wedge/hash depiction (set by MDL/CDXML/MRV/KET-style readers on a
    // stereocenter's substituent bond) and the OpenSMILES `/`/`\` E/Z
    // directional-bond marker (only meaningful adjacent to a double bond).
    // A wedge bond with no adjacent double bond has nothing to mark; writing
    // `/`/`\` for it is self-consistent only within chematic's own
    // round-trip and is silently dropped by any spec-compliant SMILES
    // parser (e.g. RDKit) on re-parse.

    /// Minimal repro from the RFC's own end-to-end trace (§3): a MOL V2000
    /// wedge bond (`BondOrder::Up`) on a tetrahedral center's substituent,
    /// no double bond anywhere in the molecule -- the RFC's own trace shows
    /// naive `write()` producing `"C(F)(Cl)/Br"`, a token any RDKit re-parse
    /// silently discards.
    #[test]
    fn test_standalone_solid_wedge_not_written_as_slash() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone solid wedge (no adjacent double bond) must not be written \
             as a directional token: got '{out}'"
        );

        // The fix only changes what gets printed -- the bond's real order in
        // memory must be untouched (this is 2D depiction metadata, not
        // something the writer should mutate).
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);

        // Round-trip: connectivity survives re-parsing the token-free output.
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// Same shape as above but with a hash wedge (`BondOrder::Down`).
    #[test]
    fn test_standalone_hash_wedge_not_written_as_backslash() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Down).unwrap();
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone hash wedge (no adjacent double bond) must not be written \
             as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Down);

        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// A genuinely 4-heavy-atom stereocenter (no implicit/explicit H at all)
    /// with a standalone wedge -- PR #130's `tetrahedral_4heavy_no_h`
    /// fixture shape, added there specifically so it isn't conflated with
    /// the 3-heavy+H case above.
    #[test]
    fn test_standalone_wedge_four_heavy_neighbors_no_h() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "4-heavy-neighbor standalone wedge must not be written as a \
             directional token: got '{out}'"
        );
    }

    /// A wedge bond that lands on a ring-closure edge (not a plain tree
    /// edge) must be suppressed identically -- the fix touches both
    /// emission sites in `write_chain`/`dfs_mark`.
    #[test]
    fn test_standalone_wedge_on_ring_closure_bond() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        // Whichever of these two becomes the DFS back-edge, it carries a
        // wedge with no adjacent double bond anywhere in this molecule.
        b.add_bond(c1, c3, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a ring-closure bond must not be written as \
             a directional token: got '{out}'"
        );
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 3);
    }

    /// The fix must not disturb `Atom.chirality`/`Molecule::stereo_neighbor_order`
    /// (the P1-S1a-core local-parity metadata) while suppressing a wedge
    /// bond's spurious directional token on the same stereocenter -- these
    /// are two unrelated mechanisms (bond-token emission vs `emit_atom`'s
    /// chirality symbol) and must keep working independently.
    #[test]
    fn test_standalone_wedge_does_not_disturb_stereocenter_chirality() {
        let mut center = Atom::new(Element::C);
        center.chirality = chematic_core::Chirality::Clockwise;
        // Force bracket notation so the chirality symbol is actually
        // printed (`needs_bracket` doesn't key off chirality alone).
        center.hydrogen_count = Some(0);
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(center);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let iodine = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, iodine, BondOrder::Up).unwrap();
        b.set_stereo_neighbor_order(c, vec![f.0, cl.0, br.0, iodine.0]);
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a stereocenter's substituent must still be \
             suppressed: got '{out}'"
        );
        assert!(
            out.contains('@'),
            "the stereocenter's own chirality symbol must still be printed \
             alongside the (now-suppressed) wedge bond: got '{out}'"
        );

        // The fix must not have touched the stereocenter metadata itself --
        // only what gets printed for the wedge bond's token.
        assert_eq!(mol.atom(c).chirality, chematic_core::Chirality::Clockwise);
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
        assert_eq!(
            mol.stereo_neighbor_order(c),
            Some([f.0, cl.0, br.0, iodine.0].as_slice())
        );
    }

    /// Round-trip check with charge and isotope present, per the required
    /// test list: connectivity/charge/isotope must survive re-parsing the
    /// (now token-free) output of a molecule carrying a standalone wedge.
    #[test]
    fn test_standalone_wedge_roundtrip_preserves_charge_and_isotope() {
        let mut b = MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.isotope = Some(13);
        let c = b.add_atom(c);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let mut o = Atom::new(Element::O);
        o.charge = -1;
        let o = b.add_atom(o);
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, o, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = write(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge must not be written as a directional token: got '{out}'"
        );

        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
        assert!(
            (0..mol2.atom_count()).any(|i| mol2.atom(AtomIdx(i as u32)).isotope == Some(13)),
            "isotope must survive the round trip: '{out}'"
        );
        assert!(
            (0..mol2.atom_count()).any(|i| mol2.atom(AtomIdx(i as u32)).charge == -1),
            "charge must survive the round trip: '{out}'"
        );
        // The original wedge bond's order in memory is still untouched.
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
    }

    /// Regression guard: a genuine E/Z directional marker (real double
    /// bond, real geometry) must still be emitted -- this fix must not
    /// suppress legitimate cases.
    #[test]
    fn test_genuine_ez_directional_bond_still_written() {
        let mol = parse("F/C=C/F").unwrap();
        let out = write(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine double-bond-adjacent directional marker must survive \
             plain write(): got '{out}'"
        );
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        let has_directional = (0..mol2.bond_count()).any(|i| {
            matches!(
                mol2.bond(BondIdx(i as u32)).order,
                BondOrder::Up | BondOrder::Down
            )
        });
        assert!(
            has_directional,
            "re-parsed molecule lost its E/Z directional bonds: '{out}'"
        );
    }

    // ── plain write() must read the bond_direction side channel too
    // (previously canonical.rs-only; see docs/rfcs/stereo2d_reader_integration_rfc.md
    // "both the non-canonical and canonical writer read the same effective
    // direction") ─────────────────────────────────────────────────────────

    /// F-C=C-F with the direction stashed via `Molecule::bond_direction`
    /// instead of a literal `BondOrder::Up`/`Down` -- exactly the shape a
    /// 2D-coordinate-derived E/Z direction (Track B) or the pre-existing
    /// aromatic-bond stash produces. Before this fix, plain `write()` never
    /// consulted the stash at all, so it silently emitted no `/`/`\` here.
    #[test]
    fn test_plain_write_reads_bond_direction_stash() {
        let mut b = MoleculeBuilder::new();
        let f1 = b.add_atom(Atom::new(Element::F));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let f2 = b.add_atom(Atom::new(Element::F));
        let b1 = b.add_bond(f1, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Double).unwrap();
        let b2 = b.add_bond(c2, f2, BondOrder::Single).unwrap();
        let mut mol = b.build();
        // Stash directions instead of literal Up/Down -- mirrors how the new
        // 2D E/Z perception stage writes (never mutating `bond.order`).
        mol.set_bond_direction(b1, BondOrder::Up);
        mol.set_bond_direction(b2, BondOrder::Up);

        let out = write(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "plain write() must surface a stashed bond_direction: got '{out}'"
        );
        // Round-trip through chematic's own parser must preserve real E/Z.
        let mol2 = parse(&out).unwrap();
        let has_directional = (0..mol2.bond_count()).any(|i| {
            matches!(
                mol2.bond(BondIdx(i as u32)).order,
                BondOrder::Up | BondOrder::Down
            )
        });
        assert!(has_directional, "stashed E/Z lost on round-trip: '{out}'");
    }

    /// A ring-closure bond is emitted from BOTH endpoints (the "open" atom
    /// and, later, the "close" atom) but is one physical bond with one
    /// `atom1`/`atom2` pair -- each side must independently re-derive its own
    /// correct character rather than sharing one (necessarily backwards on
    /// one side) value. This directly asserts the shared helper both
    /// `dfs_mark` emission sites depend on: a bond declared `atom1 = c2`
    /// must read as `Up` from `c2`'s side and `Down` from `c1`'s side.
    #[test]
    fn test_ring_closure_directional_bond_reoriented_per_endpoint() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let bidx = b.add_bond(c2, c1, BondOrder::Up).unwrap(); // atom1=c2, atom2=c1
        let mol = b.build();

        // Seen from c2 (== atom1): Up stays Up.
        assert_eq!(direction_from(&mol, bidx, BondOrder::Up, c2), BondOrder::Up);
        // Seen from c1 (!= atom1): Up flips to Down.
        assert_eq!(
            direction_from(&mol, bidx, BondOrder::Up, c1),
            BondOrder::Down
        );
    }

    /// A wedge whose direction is ALREADY suppressed (no adjacent double
    /// bond) must stay suppressed even after being resolved through
    /// `raw_bond_direction`/`direction_from` -- the reorientation step must
    /// not accidentally resurrect a token `suppress_standalone_wedge` would
    /// otherwise drop.
    #[test]
    fn test_standalone_wedge_still_suppressed_after_reorientation_helpers() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Up).unwrap();
        let mol = b.build();
        let out = write(&mol);
        assert!(!out.contains('/') && !out.contains('\\'), "got '{out}'");
    }
}
