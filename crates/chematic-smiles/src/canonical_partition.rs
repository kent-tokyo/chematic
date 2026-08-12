//! Exact colored-graph coloring + partition refinement for the
//! automorphism-orbit-aware canonical search (`canonical_search.rs`).
//!
//! This module is a *parallel* structure to `canonical.rs`'s existing
//! `morgan_ranks`/`refine_ranks` (the hash-based Morgan refinement whose
//! public meaning this PR does not change -- see `docs/
//! canonical_automorphism_pruning.md`). It is consulted only to decide
//! which individualize-refine branches are *provably* redundant; the rank
//! vectors actually handed to `CanonicalWriter` always come from the
//! original `individualize` + `refine_ranks` pipeline, unperturbed.
//!
//! Two things distinguish this from `refine_ranks`:
//! 1. The vertex/edge coloring audits every attribute `CanonicalWriter`
//!    (and the plain, non-canonical writer) actually distinguish in their
//!    output -- not just the coarser set `initial_invariant` uses for
//!    Morgan-rank search-order heuristics. See the design doc for the full
//!    writer-visible-attribute audit.
//! 2. Refinement here NEVER merges two cells on hash collision alone: cell
//!    membership is decided by full signature *equality* (`Eq`/`Ord`
//!    comparison of the whole signature struct). A `u64`/hash may still be
//!    used internally by a `HashMap`/sort for fast lookup, but Rust's
//!    `HashMap`/`sort`/`dedup` already resolve collisions via real `Eq`/`Ord`
//!    comparison -- this is standard hash-table behavior, not the forbidden
//!    "hash collision as proof of equivalence" pattern.
//!
//!    This is true of `exact_refine`'s *own* iterations. Its *starting
//!    point* (`initial_partition`'s `ranks` component) is a different story:
//!    `ranks` comes from `crate::canonical`'s pre-existing, unchanged
//!    `individualize`/`refine_ranks`, and the latter's `normalize_ranks`
//!    step does group by raw FNV-1a hash-value equality. See
//!    `canonical_search::exact_orbit_representatives`'s doc comment for the
//!    full account of what that means for this module's callers (in short:
//!    a hypothetical hash collision there is a pre-existing, crate-wide,
//!    practically-unreachable risk already relied on by
//!    `equivalent_atom_classes`/`are_atoms_equivalent`, not one this PR
//!    introduces -- but this PR does change its potential consequence from
//!    "redundant exploration" to "a silently skipped branch").

use std::collections::HashSet;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, SquarePlanarPermutation, implicit_hcount,
};

pub(crate) type CellId = u32;

/// Every molecule attribute the SMILES writer (`crate::canonical::
/// CanonicalWriter` and the plain `crate::writer`) distinguishes in its
/// output, for one atom. Two atoms with different `VertexColor` values can
/// NEVER be automorphic for canonicalization purposes.
///
/// `stereo_unique`, when `Some(atom_index)`, forces this atom's color to be
/// globally unique (bakes its own index in) -- the deliberate conservative
/// simplification this PR uses for every atom whose stereo meaning
/// (tetrahedral parity, or E/Z direction) is not cheaply provable to be
/// preserved by a candidate mapping. See the design doc, "judgment call:
/// stereo-bearing atoms are never merged". A false negative here (failing to
/// prune two genuinely-automorphic stereo atoms) only costs performance; it
/// can never cause a false merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct VertexColor {
    wildcard: bool,
    atomic_number: u8,
    isotope: u32,
    charge: i8,
    aromatic: bool,
    h_state: u16,
    atom_map: u32,
    chirality: u8,
    stereo_unique: Option<u32>,
}

/// Every bond attribute the writer distinguishes, queried from one
/// endpoint's perspective (`from`) so directional meaning (dative
/// donor/acceptor) is preserved correctly by an automorphism check that
/// always compares `edge_color(u, ..)` against `edge_color(phi(u), ..)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct EdgeColor {
    order_class: u8,
    /// `true` only for a `Dative` bond whose donor (`bond.atom1`, per
    /// `chematic_core::BondOrder::Dative`'s own doc) is the queried-from
    /// atom. Always `false` for non-dative bonds.
    from_is_donor: bool,
}

fn order_class(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Quadruple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Up => 5,
        BondOrder::Down => 6,
        BondOrder::Zero => 7,
        BondOrder::Dative => 8,
        BondOrder::QueryAny => 9,
        BondOrder::QuerySingleOrDouble => 10,
        BondOrder::QuerySingleOrAromatic => 11,
        BondOrder::QueryDoubleOrAromatic => 12,
    }
}

/// `true` when `bidx` carries any real-or-potential stereo direction
/// information: a literal `Up`/`Down` bond order, or a stashed direction
/// (`Molecule::bond_direction`, used e.g. for a ring bond next to an
/// exocyclic stereo double bond). Mirrors the raw-direction check
/// `crate::canonical::CanonicalWriter::raw_input_direction`/`crate::writer::
/// raw_bond_direction` use, minus the `resolve_ez_markers` carrier-choice
/// overlay (which depends on fully-resolved ranks, not available yet at
/// vertex-color-construction time -- irrelevant here since this function
/// only needs to decide "is this atom stereo-sensitive at all", not which
/// specific bond ends up carrying the marker).
fn bond_has_direction_info(mol: &Molecule, bidx: BondIdx) -> bool {
    matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
        || mol.bond_direction(bidx).is_some()
}

/// Every atom whose stereo meaning is not handled by this PR's exact
/// automorphism machinery, so it (and, critically, its direct neighbors)
/// must never be merged with any other atom -- see `VertexColor::
/// stereo_unique` and the design doc's "judgment call" section.
///
/// **Why 1-hop, not just the stereo atom itself**: a tetrahedral center's
/// chirality tag (`@`/`@@`) is defined relative to the *order* of its direct
/// neighbors (`Molecule::stereo_neighbor_order`); transposing two of those
/// neighbors with each other inverts the encoded configuration even though
/// the stereocenter atom itself never moves. Pinning the stereocenter alone
/// (mapping it only to itself) does NOT stop an automorphism from
/// transposing two of *its own neighbors* with each other -- caught by this
/// PR's own idempotence regression test on a stereocenter inside a ring
/// (`ring_digit_reuse_inside_stereocenter_branch_minimal`): a candidate
/// automorphism swapping the stereocenter's two ring-neighbor CH2 atoms
/// passed every per-atom color/edge check yet silently inverted the
/// encoded chirality. Pinning every direct neighbor too makes any accepted
/// automorphism fix the stereocenter's whole neighbor list pointwise, which
/// is sufficient (chirality parity depends only on that direct neighbor
/// list, never on anything further away). The same 1-hop rule is applied
/// uniformly to E/Z-direction-bearing bonds for the analogous reason (a
/// stereogenic alkene end's *un-marked* substituent, implied via
/// substituent-count rather than an explicit bond marker, must also stay
/// pinned) rather than deriving a separate, narrower rule for each case.
///
/// False negative (conservatively pinning an atom that, with more careful
/// analysis, could safely have been pruned) only costs performance. See
/// section 8/19: this is the sanctioned trade.
fn stereo_sensitive_atoms(mol: &Molecule) -> HashSet<AtomIdx> {
    let mut base: HashSet<AtomIdx> = HashSet::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        if mol.atom(idx).chirality != Chirality::None {
            base.insert(idx);
        }
    }
    for bidx in 0..mol.bond_count() {
        let bidx = BondIdx(bidx as u32);
        if bond_has_direction_info(mol, bidx) {
            let bond = mol.bond(bidx);
            base.insert(bond.atom1);
            base.insert(bond.atom2);
        }
    }

    let mut expanded = base.clone();
    for &idx in &base {
        for (nb, _) in mol.neighbors(idx) {
            expanded.insert(nb);
        }
    }
    expanded
}

/// `canonical_fidelity` gates every device in this function that exists
/// only to make canonical-SMILES search pruning safe, never to describe
/// real molecular symmetry -- see `CanonicalColoredGraph::new_topological`'s
/// doc comment for the full reasoning. Three such devices, not one:
/// - the raw `chirality` tag + `stereo_unique` (stereo, as before),
/// - `h_state` reading the *raw bracket-spelled* `hydrogen_count` rather
///   than the *effective* count (`implicit_hcount`). The raw form exists so
///   this module's coloring matches literally everything `CanonicalWriter`
///   can distinguish in its output (this module's own top-of-file doc);
///   `crate::canonical::initial_invariant` already made the opposite call
///   for the exact same reason issue #205 required it to: bracket
///   *spelling* alone (`[Cl]` vs organic-subset `Cl`, both H=0) must never
///   change a semantic answer. Using the raw form here for topological
///   classes would reintroduce exactly that bug one layer up -- two
///   spellings of one molecule getting different equivalence classes; and
/// - `atom_map`. Reaction atom-map numbers (`[CH3:1]`) are bookkeeping a
///   caller attached to individual atoms, not a molecular structural
///   property -- two atoms that differ *only* by map number are still
///   really the same topological class (`crate::canonical::
///   initial_invariant` already agrees: it never reads `atom.atom_map`
///   either). `CanonicalWriter` must still preserve map numbers exactly
///   when writing, so the canonicalization-fidelity coloring keeps folding
///   it in.
fn vertex_color(
    mol: &Molecule,
    idx: AtomIdx,
    force_unique: bool,
    canonical_fidelity: bool,
) -> VertexColor {
    let atom = mol.atom(idx);
    let stereo_unique = if force_unique { Some(idx.0) } else { None };
    let h_state = if canonical_fidelity {
        atom.hydrogen_count.map(|h| h as u16 + 1).unwrap_or(0)
    } else {
        implicit_hcount(mol, idx) as u16
    };
    let atom_map = if canonical_fidelity {
        atom.atom_map.map(|m| m as u32 + 1).unwrap_or(0)
    } else {
        0
    };
    VertexColor {
        wildcard: atom.wildcard,
        atomic_number: if atom.wildcard {
            0
        } else {
            atom.element.atomic_number()
        },
        isotope: atom.isotope.map(|i| i as u32 + 1).unwrap_or(0),
        charge: atom.charge,
        aromatic: atom.aromatic,
        h_state,
        atom_map,
        chirality: if canonical_fidelity {
            match atom.chirality {
                Chirality::None => 0,
                Chirality::CounterClockwise => 1,
                Chirality::Clockwise => 2,
                Chirality::SquarePlanar(SquarePlanarPermutation::SP1) => 3,
                Chirality::SquarePlanar(SquarePlanarPermutation::SP2) => 4,
                Chirality::SquarePlanar(SquarePlanarPermutation::SP3) => 5,
            }
        } else {
            0
        },
        stereo_unique,
    }
}

/// The molecular graph plus its (precomputed, static) writer-visible
/// coloring, shared read-only across one `canonical_smiles` call's whole
/// search.
pub(crate) struct CanonicalColoredGraph<'a> {
    mol: &'a Molecule,
    vcolor: Vec<VertexColor>,
    /// See `new`/`new_topological`. Also gates `edge_color`: an `Up`/`Down`
    /// bond order is itself an E/Z *direction* marker on what is
    /// chemically a single bond (`bond_has_direction_info`'s doc comment),
    /// so it is exactly as canonicalization-only as the vertex-side stereo
    /// devices `vertex_color` documents -- omitted here for the same
    /// reason. Caught empirically, not by inspection alone: an earlier
    /// version of this PR left `edge_color` reading raw `Up`/`Down`
    /// unconditionally, and `F/C=C\F` (cis-1,2-difluoroethene, a real
    /// mirror-symmetric molecule -- swap the two `=CF` ends) came back as
    /// 4 distinct singleton classes instead of the correct 2 merged pairs.
    canonical_fidelity: bool,
}

impl<'a> CanonicalColoredGraph<'a> {
    /// Canonical-SMILES-generation coloring: stereo-bearing atoms (and their
    /// direct neighbors) are pinned globally unique via `stereo_unique`, and
    /// `h_state` reads the raw bracket-spelled H count -- see
    /// `stereo_sensitive_atoms`'s and `vertex_color`'s doc comments. This is
    /// the *only* coloring production `canonical_search.rs` ever uses.
    pub(crate) fn new(mol: &'a Molecule) -> Self {
        Self::build(mol, true)
    }

    /// Coloring for [`crate::topological_equivalence_classes`]: identical to
    /// `new` except every canonicalization-only device `vertex_color` (and
    /// this struct's own `edge_color`) documents is turned off (see
    /// `topological_equivalence_classes`'s doc comment for the full
    /// reasoning -- the short version: `new`'s stereo and raw-H-spelling
    /// handling exist to make canonical-SMILES tie-breaking
    /// conservative/safe, a different goal from reporting real topological
    /// symmetry; applying any of them here would wrongly split atoms that a
    /// real symmetry query must merge -- a meso compound's two
    /// mirror-equivalent stereocenters via stereo, `[Cl]` vs `Cl` via raw
    /// H-spelling, or a cis-alkene's two mirror-equivalent ends via
    /// `Up`/`Down` bond-direction markers).
    pub(crate) fn new_topological(mol: &'a Molecule) -> Self {
        Self::build(mol, false)
    }

    fn build(mol: &'a Molecule, canonical_fidelity: bool) -> Self {
        let sensitive = if canonical_fidelity {
            stereo_sensitive_atoms(mol)
        } else {
            HashSet::new()
        };
        let vcolor = (0..mol.atom_count())
            .map(|i| {
                let idx = AtomIdx(i as u32);
                vertex_color(mol, idx, sensitive.contains(&idx), canonical_fidelity)
            })
            .collect();
        Self {
            mol,
            vcolor,
            canonical_fidelity,
        }
    }

    pub(crate) fn mol(&self) -> &'a Molecule {
        self.mol
    }

    pub(crate) fn n(&self) -> usize {
        self.vcolor.len()
    }

    pub(crate) fn vertex_color(&self, a: AtomIdx) -> VertexColor {
        self.vcolor[a.0 as usize]
    }

    pub(crate) fn edge_color(&self, from: AtomIdx, bidx: BondIdx) -> EdgeColor {
        let bond = self.mol.bond(bidx);
        let from_is_donor = bond.order == BondOrder::Dative && bond.atom1 == from;
        // `Up`/`Down` are E/Z direction markers on an otherwise-single
        // bond; in topological mode they collapse to plain `Single`, same
        // as every other stereo device this struct excludes there.
        let order =
            if !self.canonical_fidelity && matches!(bond.order, BondOrder::Up | BondOrder::Down) {
                BondOrder::Single
            } else {
                bond.order
            };
        EdgeColor {
            order_class: order_class(order),
            from_is_donor,
        }
    }

    pub(crate) fn neighbors(&self, a: AtomIdx) -> impl Iterator<Item = (AtomIdx, BondIdx)> + '_ {
        self.mol.neighbors(a)
    }
}

/// A partition of atom indices into cells, `cell_of[i]` = the cell of atom
/// `i`. Cell ids are gap-free ordinals; a lower cell id does not carry any
/// meaning beyond "distinct from every other cell id" for this module's own
/// purposes (unlike `crate::canonical`'s `ranks`, which additionally encodes
/// a *search-order* preference -- this partition is never itself handed to
/// `CanonicalWriter`).
#[derive(Debug, Clone)]
pub(crate) struct Partition {
    pub(crate) cell_of: Vec<CellId>,
}

impl Partition {
    /// Group atom indices by cell, cell ids ascending; atoms within a cell
    /// in ascending atom-index order. Test-only utility (production code
    /// only ever needs `cell_of[atom_index]` lookups, see
    /// `canonical_search.rs`).
    #[cfg(test)]
    pub(crate) fn cell_members(&self) -> Vec<Vec<AtomIdx>> {
        let max_cell = self
            .cell_of
            .iter()
            .copied()
            .max()
            .map(|m| m as usize + 1)
            .unwrap_or(0);
        let mut cells: Vec<Vec<AtomIdx>> = vec![Vec::new(); max_cell];
        for (i, &c) in self.cell_of.iter().enumerate() {
            cells[c as usize].push(AtomIdx(i as u32));
        }
        cells
    }
}

fn assign_cell_ids<K: Ord + Clone>(keys: &[K]) -> Vec<CellId> {
    let mut distinct: Vec<K> = keys.to_vec();
    distinct.sort();
    distinct.dedup();
    keys.iter()
        .map(|k| {
            distinct
                .binary_search(k)
                .expect("key present in distinct set") as CellId
        })
        .collect()
}

fn count_distinct_cells(cell_of: &[CellId]) -> usize {
    let mut seen: Vec<CellId> = cell_of.to_vec();
    seen.sort_unstable();
    seen.dedup();
    seen.len()
}

/// Build the initial (pre-refinement) partition for the current search node:
/// combine the current Morgan-derived `ranks` (already gap-free ordinals,
/// encoding every individualization done so far) with the full exact
/// `VertexColor` (which additionally distinguishes atom-map/stereo/isotope/
/// charge differences the hash-based `initial_invariant` does not carry) so
/// a subsequent exact refinement pass can never wrongly merge two atoms
/// `ranks` alone would conflate.
pub(crate) fn initial_partition(graph: &CanonicalColoredGraph, ranks: &[u64]) -> Partition {
    let keys: Vec<(u64, VertexColor)> = (0..ranks.len())
        .map(|i| (ranks[i], graph.vertex_color(AtomIdx(i as u32))))
        .collect();
    Partition {
        cell_of: assign_cell_ids(&keys),
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone)]
struct RefineSignature {
    own: CellId,
    neighbors: Vec<(u8, bool, CellId)>,
}

/// Refine `partition` to a fixpoint using an exact structural signature
/// (own cell + sorted multiset of (edge color, neighbor cell)), comparing
/// full signatures via real equality -- never a hash-collision shortcut.
/// Mirrors `crate::canonical::refine_ranks`'s neighbor-aggregation shape.
pub(crate) fn exact_refine(graph: &CanonicalColoredGraph, mut partition: Partition) -> Partition {
    let n = partition.cell_of.len();
    let max_iter = n + 2;
    for _ in 0..max_iter {
        let old_cells = count_distinct_cells(&partition.cell_of);

        let sigs: Vec<RefineSignature> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut neighbors: Vec<(u8, bool, CellId)> = graph
                    .neighbors(idx)
                    .map(|(nb, bidx)| {
                        let ec = graph.edge_color(idx, bidx);
                        (
                            ec.order_class,
                            ec.from_is_donor,
                            partition.cell_of[nb.0 as usize],
                        )
                    })
                    .collect();
                neighbors.sort_unstable();
                RefineSignature {
                    own: partition.cell_of[i],
                    neighbors,
                }
            })
            .collect();

        let new_cell_of = assign_cell_ids(&sigs);
        let new_cells = count_distinct_cells(&new_cell_of);
        partition.cell_of = new_cell_of;

        if new_cells <= old_cells {
            break;
        }
    }
    partition
}

/// Real (hash-collision-free) topological equivalence classes for every
/// atom of `mol`: `classes[i] == classes[j]` iff atoms `i` and `j` are
/// indistinguishable by exact structural color refinement -- element,
/// isotope, charge, aromaticity, effective H-count, and bond
/// order/connectivity, propagated to a fixpoint by `exact_refine`.
/// Deliberately excluded: stereo and reaction atom-map numbers -- see
/// "Canonicalization-only devices" below. Class ids are gap-free ordinals
/// with no meaning beyond "same id = same class" (mirrors `Partition`'s own
/// `cell_of`, which this function returns almost verbatim).
///
/// Conceptually the same query as RDKit's `Chem.CanonicalRankAtoms(mol,
/// breakTies=False)` restricted to constitutional (non-stereo) symmetry;
/// this doc comment states chematic's own exclusion set explicitly rather
/// than asserting parity with any particular RDKit parameter defaults,
/// which this PR did not verify.
///
/// `crate::canonical::equivalent_atom_classes` already exposes a similarly-
/// shaped result (also chirality-blind -- `initial_invariant` never reads
/// `atom.chirality`); the two agree on orbit *structure* on every case this
/// PR's test suite checked, including the meso-compound case. The
/// difference this function is for: `equivalent_atom_classes` is built on
/// `morgan_ranks`/`refine_ranks`, whose own `normalize_ranks` step groups by
/// raw FNV-1a hash-value equality (see this module's top-of-file doc
/// comment) -- a real, if practically negligible, hash-collision risk. This
/// function is built on `exact_refine`, which never merges two cells on
/// anything but full signature `Eq` -- the "real (hash-collision-free)"
/// guarantee above.
///
/// # Known limitation: this is a sound over-approximation, not exact orbits
///
/// "Same class" does **not** strictly prove "provably interchangeable by a
/// real automorphism" -- only the converse is guaranteed (different class
/// always means provably NOT interchangeable). Like any Morgan/1-WL-style
/// color refinement, `exact_refine` never merges two atoms that differ in
/// any local, propagatable attribute, but on a pathological input (a
/// WL-indistinguishable disjoint union of differently-sized regular
/// subgraphs -- see
/// `canonical_automorphism::tests::triangle_and_square_same_wl_cell_different_orbits`
/// for a constructed witness) it can under-split: report two atoms in one
/// class that a full automorphism-group search would separate. No real
/// molecule in this crate's test corpus has ever exhibited this; it is
/// documented, not fixed, exactly as it already is for
/// `morgan_ranks`/`equivalent_atom_classes` (same underlying algorithm
/// shape, same limitation, longstanding and never a reported problem).
///
/// # Canonicalization-only devices are deliberately excluded, not an oversight
///
/// This crate's other automorphism-adjacent machinery
/// (`CanonicalColoredGraph::new`, used by canonical-SMILES generation) folds
/// in devices that exist purely to make canonical-SMILES search pruning
/// *safe*, not to describe real molecular symmetry -- see `vertex_color`'s
/// doc comment for the full list. Four matter most:
/// - **Stereo (atom parity)**: `new`'s coloring reads an atom's raw `@`/`@@` parity tag
///   and pins every stereo-bearing atom (and its direct neighbors) globally
///   unique via `stereo_unique`, because verifying that a candidate mapping
///   *actually* preserves tetrahedral/E-Z parity (which depends on neighbor
///   *order*, not just neighbor identity) is out of scope for that module --
///   see `stereo_sensitive_atoms`'s doc comment. A false negative there
///   (failing to merge two atoms a full CIP-aware analysis would prove
///   equivalent) only costs canonical-SMILES search performance, never
///   correctness, so it is the right conservative default *for that
///   caller*. It is the wrong default here: this function's whole purpose
///   is to expose *real* symmetry (e.g. for meso-compound detection, a
///   stereocenter mapping onto another under the molecule's own internal
///   mirror symmetry). A meso pair's raw `@`/`@@` tags legitimately
///   *differ* -- that textual difference is exactly what the mirror
///   relationship looks like once written down as SMILES
///   neighbor-order-relative parity -- so folding either stereo device in
///   here would put a meso pair in two different classes, silently
///   downgrading this from a "topological symmetry" answer to a "canonical
///   tie-break" answer and defeating the caller's whole reason for asking.
/// - **H-count spelling**: `new`'s coloring reads the *raw bracket-spelled*
///   hydrogen count, so `[Cl]` and organic-subset `Cl` (same molecule, same
///   effective H=0, different spelling) get different colors there. That
///   is intentional for that caller (writer-output-exact auditing) but
///   would be exactly the issue #205 bug one layer up if reused here: two
///   spellings of the same molecule getting different equivalence classes.
/// - **Atom-map numbers**: reaction bookkeeping (`[CH3:1]`) a caller
///   attached to a specific atom, not a molecular structural property. Two
///   atoms differing *only* by map number are still the same real
///   topological class; `crate::canonical::initial_invariant` already
///   agrees (it never reads `atom.atom_map` either). Verified empirically,
///   not just asserted: on `[CH3:1]c1ccc(cc1)[CH3:2]` the two chemically
///   identical methyls get different classes under `new`'s
///   canonicalization-fidelity coloring (map numbers must round-trip
///   through the writer) but the same class here (see
///   `atom_map_number_alone_does_not_split_the_class` below).
/// - **Stereo (bond direction)**: `new`'s `edge_color` reads an `Up`/`Down`
///   bond order literally -- the E/Z direction marker (`/`/`\`) on what is
///   chemically a single bond -- as a distinct edge color from plain
///   `Single`. Same category of device as atom parity above, same fix:
///   `CanonicalColoredGraph::new_topological` collapses `Up`/`Down` to
///   `Single` there. Verified empirically: `F/C=C\F` (cis-1,2-
///   difluoroethene, real mirror symmetry across the two `=CF` ends) came
///   back as 4 singleton classes before this collapse was added, 2 merged
///   pairs after (see `ez_bond_direction_marker_alone_does_not_split_the_class`
///   below). Note `edge_color` only ever reads `bond.order` -- the
///   separately-stashed `Molecule::bond_direction` (`bond_has_direction_info`'s
///   doc comment: used for a ring bond next to an exocyclic stereo double
///   bond) is not part of `EdgeColor` in either mode, so there is nothing
///   to exclude there today; if a future change ever folds it in, it needs
///   the same `canonical_fidelity` gate.
///
/// `CanonicalColoredGraph::new_topological` therefore omits all four,
/// using the *effective* H count (`implicit_hcount`, matching
/// `crate::canonical::initial_invariant`'s post-#205 choice), no atom-map
/// signal, and no stereo signal (atom parity or bond direction) at all --
/// matching `morgan_ranks`'s own long-standing `initial_invariant`, which
/// has never read `atom.chirality` or `atom.atom_map`, and whose own
/// `bond_order_value` collapses `Up`/`Down` to the same value as `Single`.
pub fn topological_equivalence_classes(mol: &Molecule) -> Vec<usize> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }
    let graph = CanonicalColoredGraph::new_topological(mol);
    // No pre-individualization: an all-equal starting `ranks` collapses
    // `initial_partition` to grouping by raw vertex color alone, the
    // correct from-scratch starting point for a whole-molecule query (as
    // opposed to `canonical_search.rs`'s use, which seeds from the current
    // search node's already-partially-individualized `ranks`). This also
    // sidesteps `crate::canonical::normalize_ranks`'s raw-FNV-1a-hash-value
    // grouping entirely (see this module's own top-of-file doc comment) --
    // `exact_refine`'s own iterations never rely on anything but real
    // `Eq`/`Ord` signature comparison.
    let ranks = vec![0u64; n];
    let partition = exact_refine(&graph, initial_partition(&graph, &ranks));
    partition.cell_of.into_iter().map(|c| c as usize).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn benzene_all_one_cell_before_refine_after_refine() {
        let mol = parse("c1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let ranks = crate::canonical::morgan_ranks(&mol);
        let p = exact_refine(&graph, initial_partition(&graph, &ranks));
        let cells = p.cell_members();
        let nonempty: Vec<_> = cells.into_iter().filter(|c| !c.is_empty()).collect();
        assert_eq!(nonempty.len(), 1, "plain benzene: all 6 atoms in one cell");
        assert_eq!(nonempty[0].len(), 6);
    }

    #[test]
    fn isotope_difference_splits_cell() {
        // One ring carbon isotope-labeled -- must land in its own cell even
        // though morgan_ranks (initial_invariant has no isotope bit at the
        // Morgan-search-heuristic layer... actually it does include iso, but
        // this exercises the exact-refine path independently of that).
        let mol = parse("[13cH]1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let ranks = crate::canonical::morgan_ranks(&mol);
        let p = exact_refine(&graph, initial_partition(&graph, &ranks));
        let cells = p.cell_members();
        let iso_cell = cells.iter().find(|c| c.contains(&AtomIdx(0))).unwrap();
        assert_eq!(
            iso_cell.len(),
            1,
            "isotope-labeled atom must be alone in its cell"
        );
    }

    // --- topological_equivalence_classes ------------------------------

    #[test]
    fn bracket_h_spelling_does_not_split_the_class() {
        // [cH]1ccccc1 is still just benzene: atom 0's H is spelled out
        // explicitly in the bracket, atoms 1..5 get theirs implicitly, but
        // the *effective* H count (1 each) is identical. issue #205's own
        // bug shape one layer up: raw hydrogen_count bracket-spelling must
        // never leak into a semantic answer.
        let mol = parse("[cH]1ccccc1").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 6);
        assert!(
            classes.iter().all(|&c| c == classes[0]),
            "bracket-spelled H must not split an otherwise-symmetric ring: {classes:?}"
        );
    }

    #[test]
    fn wildcard_never_collides_with_a_zero_h_real_atom() {
        // [*]C(C)(C)C: the wildcard's *effective* H is 0 (implicit_hcount
        // returns 0 for wildcards, per chematic-core), matching the central
        // quaternary carbon's own effective H (0, bonded to 4 heavy atoms).
        // Dropping the canonical branch's "+1, None->0" offset in the
        // topological branch (this function reads plain `implicit_hcount`)
        // must not let the two collide on `h_state` alone -- `wildcard`
        // and `atomic_number` (forced to 0 only for wildcards, and no real
        // element has atomic_number 0) still separate them.
        let mol = parse("[*]C(C)(C)C").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_ne!(
            classes[0], classes[1],
            "wildcard must never share a class with a real 0-H atom: {classes:?}"
        );
    }

    #[test]
    fn benzene_all_six_carbons_one_class() {
        let mol = parse("c1ccccc1").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 6);
        assert!(
            classes.iter().all(|&c| c == classes[0]),
            "all 6 benzene carbons must share one class, got {classes:?}"
        );
    }

    #[test]
    fn toluene_five_classes_ortho_and_meta_pairs_merge() {
        // Cc1ccccc1: 0=methyl, 1=ipso (nbrs 0,2,6), 2&6=ortho, 3&5=meta,
        // 4=para (confirmed via atom.neighbors(), not assumed).
        let mol = parse("Cc1ccccc1").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 7);

        // Ortho pair and meta pair each merge; everything else is a
        // singleton class, distinct from every other singleton.
        assert_eq!(classes[2], classes[6], "ortho carbons must be equivalent");
        assert_eq!(classes[3], classes[5], "meta carbons must be equivalent");
        assert_ne!(classes[0], classes[1], "methyl != ipso");
        assert_ne!(classes[0], classes[4], "methyl != para");
        assert_ne!(classes[1], classes[4], "ipso != para");
        assert_ne!(classes[0], classes[2], "methyl != ortho");
        assert_ne!(classes[1], classes[2], "ipso != ortho");
        assert_ne!(classes[4], classes[2], "para != ortho");
        assert_ne!(classes[0], classes[3], "methyl != meta");
        assert_ne!(classes[1], classes[3], "ipso != meta");
        assert_ne!(classes[4], classes[3], "para != meta");
        assert_ne!(classes[2], classes[3], "ortho != meta");

        let mut distinct: Vec<usize> = classes.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(
            distinct.len(),
            5,
            "toluene: methyl, ipso, ortho-pair, meta-pair, para = 5 distinct classes, got {classes:?}"
        );
    }

    #[test]
    fn asymmetric_chain_every_atom_its_own_class() {
        // F-C-C-C-C-Cl: no two atoms are topologically interchangeable
        // (the two ends differ, so no rank-refinement plateau survives).
        let mol = parse("FCCCCCl").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 6);
        let mut distinct = classes.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(
            distinct.len(),
            6,
            "no-symmetry chain: every atom must be its own class, got {classes:?}"
        );
    }

    #[test]
    fn meso_tartaric_acid_stereocenters_share_a_class_despite_opposite_tags() {
        // OC(=O)[C@H](O)[C@@H](O)C(=O)O -- meso-tartaric acid. Atom map
        // (confirmed via atom.neighbors()/chirality, not assumed):
        //   0=O(H)  1=C(=O)  2=O(=)      3=C@H (stereocenter)  4=O(H)
        //   5=C@@H (stereocenter)  6=O(H)  7=C(=O)  8=O(=)  9=O(H)
        // The molecule has an internal mirror relating the two halves:
        // {0,1,2,3,4} <-> {9,7,8,5,6}. Atoms 3 and 5 are the two
        // stereocenters; their raw parity tags are opposite
        // (CounterClockwise vs Clockwise) -- that is the meso relationship,
        // not a sign of inequivalence -- so with stereo correctly excluded
        // they must land in the SAME topological class.
        let mol = parse("OC(=O)[C@H](O)[C@@H](O)C(=O)O").unwrap();
        assert_eq!(mol.atom(AtomIdx(3)).chirality, Chirality::CounterClockwise);
        assert_eq!(mol.atom(AtomIdx(5)).chirality, Chirality::Clockwise);

        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 10);
        assert_eq!(
            classes[3], classes[5],
            "meso stereocenters must share a class despite opposite @ / @@ tags: {classes:?}"
        );
        // The rest of the mirror mapping too, for good measure.
        assert_eq!(classes[0], classes[9], "the two carboxyl OH oxygens");
        assert_eq!(classes[1], classes[7], "the two carboxyl carbons");
        assert_eq!(classes[2], classes[8], "the two carbonyl (=O) oxygens");
        assert_eq!(classes[4], classes[6], "the two stereocenter hydroxyls");
    }

    #[test]
    fn ez_bond_direction_marker_alone_does_not_split_the_class() {
        // cis-1,2-difluoroethene: F/C=C\F. atoms: 0=F, 1=C, 2=C, 3=F. Real
        // mirror symmetry swaps the two =CF ends (0<->3, 1<->2), even
        // though one C-F bond is marked Up and the other Down -- an E/Z
        // direction marker on what is chemically a single bond, exactly as
        // canonicalization-only as tetrahedral parity. Caught empirically:
        // an earlier version of this function read raw Up/Down through
        // `edge_color` unconditionally and returned 4 singleton classes
        // here instead of 2 merged pairs.
        let mol = parse(r"F/C=C\F").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 4);
        assert_eq!(
            classes[0], classes[3],
            "the two mirror-equivalent F: {classes:?}"
        );
        assert_eq!(
            classes[1], classes[2],
            "the two mirror-equivalent C: {classes:?}"
        );
    }

    #[test]
    fn atom_map_number_alone_does_not_split_the_class() {
        // p-xylene, the two methyls given different reaction atom-map
        // numbers (:1 / :2). Map numbers are caller bookkeeping, not
        // molecular structure -- the two methyls (and, by the same mirror
        // symmetry, the two ortho/meta ring-carbon pairs either side) must
        // still land in one class each, unlike `CanonicalColoredGraph::new`
        // (canonical-SMILES's own coloring), which must keep map numbers
        // distinguishable so the writer round-trips them exactly.
        let mol = parse("[CH3:1]c1ccc(cc1)[CH3:2]").unwrap();
        let classes = topological_equivalence_classes(&mol);
        assert_eq!(classes.len(), 8);
        assert_eq!(
            classes[0], classes[7],
            "differently-mapped-but-otherwise-identical methyls must share a class: {classes:?}"
        );
    }
}
