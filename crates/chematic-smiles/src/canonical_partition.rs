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

use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, Molecule};

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

fn vertex_color(mol: &Molecule, idx: AtomIdx, force_unique: bool) -> VertexColor {
    let atom = mol.atom(idx);
    let stereo_unique = if force_unique { Some(idx.0) } else { None };
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
        h_state: atom.hydrogen_count.map(|h| h as u16 + 1).unwrap_or(0),
        atom_map: atom.atom_map.map(|m| m as u32 + 1).unwrap_or(0),
        chirality: match atom.chirality {
            Chirality::None => 0,
            Chirality::CounterClockwise => 1,
            Chirality::Clockwise => 2,
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
}

impl<'a> CanonicalColoredGraph<'a> {
    pub(crate) fn new(mol: &'a Molecule) -> Self {
        let sensitive = stereo_sensitive_atoms(mol);
        let vcolor = (0..mol.atom_count())
            .map(|i| {
                let idx = AtomIdx(i as u32);
                vertex_color(mol, idx, sensitive.contains(&idx))
            })
            .collect();
        Self { mol, vcolor }
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
        EdgeColor {
            order_class: order_class(bond.order),
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
}
