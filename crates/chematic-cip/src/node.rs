//! Digraph node types.

use chematic_core::AtomIdx;

use crate::edge::EdgeId;
use crate::rational::AtomicNumberKey;

/// Index of a [`CipNode`] in a [`crate::digraph::CipDigraph`]'s arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NodeId(pub u32);

/// What a digraph node represents.
///
/// CIP's hierarchical digraph is built from a molecule's real atoms, plus explicit
/// *duplicate* nodes standing in for the extra bonds of a multiple bond or the closure
/// of a ring. Representing these as distinct, provenance-carrying nodes (rather than
/// pooling atom keys into a shell multiset, as the existing fast/approximate engine in
/// `chematic-chem` does) is the entire point of this crate -- see
/// `docs/rfcs/cip_accurate_rfc.md`. A duplicate node's provenance (which real atom it stands
/// in for, and via which bond) is exactly the information a shell-multiset comparison
/// discards, and discarding it is what made a locally-correct per-bond-type duplication
/// rule unsafe to add to that engine (proven: a reverted triple-bond fix went net
/// negative there).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipNodeKind {
    /// A real atom in the molecule.
    Atom { atom_idx: AtomIdx },
    /// A phantom duplicate contributed by a multiple bond. `source_atom` is the atom
    /// whose substituent list this duplicate appears in; `duplicated_atom` is the real
    /// atom it stands in for; `bond_order` is the bond's multiplicity
    /// (`BondOrder::order_int()`: 2 for double, 3 for triple). A bond of multiplicity
    /// `k` contributes `k - 1` duplicate leaves on each of its two atoms.
    MultipleBondDuplicate {
        source_atom: AtomIdx,
        duplicated_atom: AtomIdx,
        bond_order: u8,
    },
    /// A phantom duplicate marking a ring closure: `source_atom` is the atom being
    /// expanded, `closure_atom` is an ancestor of `source_atom` on the current
    /// root-to-node path that `source_atom` is bonded back to. Emitted instead of
    /// re-expanding `closure_atom` as a real subtree, which is what guarantees every
    /// root-to-leaf path is finite (see [`crate::digraph`] module docs).
    RingDuplicate {
        source_atom: AtomIdx,
        closure_atom: AtomIdx,
    },
    /// The implicit hydrogen of a chiral bracket atom (`[C@H]`), when present.
    ImplicitHydrogen,
}

/// A single node in a [`crate::digraph::CipDigraph`].
#[derive(Debug, Clone)]
pub struct CipNode {
    pub id: NodeId,
    pub kind: CipNodeKind,
    pub parent: Option<NodeId>,
    pub depth: u32,
    pub incoming_edge: Option<EdgeId>,
    /// The atomic number this node compares by, computed once at construction time.
    /// Always `Integral` for `Atom`/`RingDuplicate`/`ImplicitHydrogen` nodes -- MANCUDE
    /// treatment only ever applies to a `MultipleBondDuplicate` whose *owner*
    /// (`source_atom`) sits in a resonance component (see `crate::mancude`'s module
    /// docs for why the owner, not `duplicated_atom`, is the correct key). **Not yet
    /// read by the comparator** (`crate::compare`) -- Milestone 3B-1a computes and
    /// tests this field in isolation; wiring it into ranking is Milestone 3B-1b.
    pub atomic_number: AtomicNumberKey,
}
