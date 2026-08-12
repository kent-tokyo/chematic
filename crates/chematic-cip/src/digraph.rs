//! The hierarchical digraph itself.
//!
//! A [`CipDigraph`] is rooted at a candidate stereocenter and lazily expands into its
//! substituent branches. Two structural rules do the real work, and are what this
//! crate exists to get right (see `docs/rfcs/cip_accurate_rfc.md`):
//!
//! - **Multiple-bond duplication**: a bond of order `k` between atoms `a` and `b`
//!   contributes `k - 1` [`CipNodeKind::MultipleBondDuplicate`] leaves to *each* of
//!   `a`'s and `b`'s own child lists (a double bond duplicates its partner into both
//!   ends; a triple bond duplicates it twice into both ends). Deliberately implemented
//!   as two explicit halves in [`CipDigraph::expand_children`] -- an "arrival side" (a
//!   node whose incoming edge was a multiple bond gets `k - 1` duplicates of its own
//!   parent) and a "departure side" (iterating a node's neighbors, a fresh neighbor
//!   found via a multiple bond contributes `k - 1` duplicates too) -- because forgetting
//!   either half is exactly the bug `d0e726b` fixed in the older, approximate engine in
//!   `chematic-chem`: a real atom "sees" its multiply-bonded partner twice, from *both*
//!   atoms' own local perspective, regardless of which direction a tree traversal
//!   happened to cross the bond.
//! - **Ring-closure termination**: when expanding atom `a`'s neighbors, a neighbor `b`
//!   that is already an *ancestor of `a` on the current root-to-`a` path* becomes a
//!   [`CipNodeKind::RingDuplicate`] leaf instead of a real, re-expanded subtree. This
//!   is the single property that guarantees every root-to-leaf path is finite: an atom
//!   can appear as a real node at most once per path (it would have been caught as an
//!   ancestor on any later occurrence), so path length is bounded by
//!   [`Molecule::atom_count`]. [`crate::budget::CipBudget`] is a *separate* backstop
//!   against total node-count blow-up on richly fused/bridged/cage systems (many
//!   individually-finite branches can still add up to a lot of nodes) -- it is not
//!   what makes expansion terminate.
//!
//! Both rules mirror the *already-correct* mechanism in
//! `crates/chematic-chem/src/cip.rs`'s `cip_branch_spheres` (its ancestor-tracking
//! `visited: HashSet<AtomIdx>`, cloned and extended per branch), just expressed as
//! explicit, provenance-carrying nodes instead of atom keys pooled into a shell
//! multiset. That pooling is exactly what this crate is designed not to repeat.

use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};

use crate::CipError;
use crate::budget::CipBudget;
use crate::edge::{CipEdge, EdgeId};
use crate::mancude::MancudeContext;
use crate::node::{CipNode, CipNodeKind, NodeId};
use crate::rational::AtomicNumberKey;

/// A lazily-expanding hierarchical digraph rooted at one atom (typically a candidate
/// stereocenter). See the module docs for the two structural rules that make this a
/// real replacement for shell-multiset pooling.
pub struct CipDigraph<'m> {
    mol: &'m Molecule,
    /// `None` for every existing call site (`Self::new`) -- Milestone 3B-1a's MANCUDE
    /// fractional-atomic-number treatment is only active via [`Self::new_with_mancude`],
    /// a separate, not-yet-wired-into-`assign_cip_accurate_experimental` entry point (see
    /// `crate::mancude`'s module docs). Attaching a `MancudeContext` computed for a
    /// *different* molecule than `mol` would silently misattribute fractional values --
    /// callers are responsible for computing it from the exact same (Kekulé-form) `mol`
    /// passed in here.
    mancude: Option<&'m MancudeContext>,
    /// `None` for every existing call site (`Self::new`/`Self::new_with_mancude`) --
    /// [`Self::new_with_artificial_ancestor`] is a separate, Milestone 4B-1 entry point
    /// that pre-seeds the root's own ancestor set with one extra atom, so a neighbor
    /// equal to it terminates as a [`CipNodeKind::RingDuplicate`] leaf immediately
    /// instead of expanding as a real subtree. See that constructor's doc comment for
    /// why: it is what makes an *embedded* stereocenter's true auxiliary descriptor
    /// (as seen from one specific branch of an outer, still-unresolved stereocenter)
    /// computable at all -- rooting fresh at that embedded atom with no artificial
    /// ancestor reproduces its ordinary *global/molecular* digraph instead, which is a
    /// different (and, for a Rule-4-tied atom, differently ambiguous) computation.
    artificial_ancestor: Option<AtomIdx>,
    nodes: Vec<CipNode>,
    edges: Vec<CipEdge>,
    /// Parallel to `nodes`: `None` until a node's children have been computed once.
    children_cache: Vec<Option<Vec<NodeId>>>,
    root: NodeId,
    budget: CipBudget,
    node_count: usize,
    expansion_count: usize,
    /// Count of constructed `MultipleBondDuplicate` nodes whose `atomic_number` is
    /// `AtomicNumberKey::Rational` -- i.e. the MANCUDE fractional-atomic-number path was
    /// actually exercised, not just wired in. A pure "path reached" counter: it says
    /// nothing about whether the fraction ever *decided* a comparison (see
    /// `CompareContext::fractional_decisions` for that). Exists because Milestone 3B's
    /// fractional machinery was kept-but-measured-inert on the available corpus (see
    /// `docs/rfcs/cip_accurate_rfc.md`'s Milestone 3B closeout entry) -- "correct and never
    /// fired" needs a test that can tell the difference from "correct and never even
    /// reached."
    fractional_nodes_emitted: u64,
}

impl<'m> CipDigraph<'m> {
    /// Create a digraph rooted at `root_atom`. Materializes only the root node --
    /// children are computed on first access via [`Self::expand_children`] or the
    /// [`crate::DigraphExpander`] trait.
    ///
    /// Every `MultipleBondDuplicate` node's `atomic_number` is a plain integer (the
    /// represented atom's real atomic number) -- exactly today's existing behavior. Use
    /// [`Self::new_with_mancude`] to additionally apply MANCUDE fractional atomic numbers.
    pub fn new(mol: &'m Molecule, root_atom: AtomIdx, budget: CipBudget) -> Result<Self, CipError> {
        Self::new_impl(mol, root_atom, budget, None, None)
    }

    /// Stable identity for the attached [`MancudeContext`], if any -- the raw pointer
    /// address of the reference this digraph was constructed with. Used only as a
    /// same-resolution memoization cache-key discriminator (issue #107's
    /// `compare_ligands`/`rank_children` pairwise cache in `crate::compare`); never
    /// compared across different digraphs or persisted beyond one resolution's
    /// lifetime, so pointer identity (not content equality) is exactly the right
    /// notion -- two digraphs attaching *equal-content* but distinct `MancudeContext`
    /// values must never be treated as interchangeable by a cache, and pointer
    /// identity guarantees that.
    pub(crate) fn mancude_identity(&self) -> Option<usize> {
        self.mancude.map(|m| m as *const MancudeContext as usize)
    }

    /// This digraph's own construction-time budget -- part of the same cache-key
    /// discriminator as [`Self::mancude_identity`] (see that method's doc).
    pub(crate) fn budget(&self) -> CipBudget {
        self.budget
    }

    /// Like [`Self::new`], but `artificial_ancestor` is treated as already visited from
    /// the start -- a neighbor of any node in this digraph equal to `artificial_ancestor`
    /// becomes a [`CipNodeKind::RingDuplicate`] leaf immediately, instead of a real,
    /// re-expanded subtree (exactly the ordinary ring-closure rule in the module docs,
    /// just pre-seeded rather than discovered by walking up the tree).
    ///
    /// **Milestone 4B-1.** Motivation: when an outer stereocenter's two branches tie
    /// under Rules 1a/2 (Rule 4 territory), resolving the tie needs each branch's
    /// *auxiliary* descriptor for whatever embedded stereocenter it reaches -- not that
    /// atom's ordinary *molecular* descriptor. Rooting a fresh, independent digraph at
    /// the embedded atom (`Self::new`, no artificial ancestor) computes the *molecular*
    /// descriptor: all of that atom's real neighbors compete as ordinary children,
    /// including the one leading back toward the outer stereocenter -- which, for a
    /// constitutionally symmetric branch pair, is typically itself ambiguous (that
    /// neighbor's own subtree mirrors the embedded atom's other ring-continuation
    /// child). Pre-seeding that one specific neighbor as an artificial ancestor instead
    /// terminates it as a childless duplicate leaf immediately, matching what CIP's
    /// auxiliary-descriptor procedure actually asks for: this atom's local priority
    /// order *as seen from one specific incoming direction*, not its independent
    /// global ranking. `root_atom` is the embedded stereocenter; `artificial_ancestor`
    /// is the real neighbor the branch was reached through (i.e. the atom immediately
    /// preceding `root_atom` on the outer stereocenter's own root-to-`root_atom` path).
    pub fn new_with_artificial_ancestor(
        mol: &'m Molecule,
        root_atom: AtomIdx,
        artificial_ancestor: AtomIdx,
        budget: CipBudget,
    ) -> Result<Self, CipError> {
        Self::new_impl(mol, root_atom, budget, None, Some(artificial_ancestor))
    }

    /// Like [`Self::new`], but `MultipleBondDuplicate` nodes whose owner atom
    /// (`source_atom`) has a fractional atomic number in `mancude` get
    /// `AtomicNumberKey::Rational` instead of the plain integer of `duplicated_atom` --
    /// see `crate::mancude`'s module docs for why the owner, not the represented atom, is
    /// the correct key. `mol` **must** be the same (Kekulé-form) molecule `mancude` was
    /// computed from -- `MancudeContext` indexes by `AtomIdx` into that specific molecule.
    pub fn new_with_mancude(
        mol: &'m Molecule,
        root_atom: AtomIdx,
        budget: CipBudget,
        mancude: &'m MancudeContext,
    ) -> Result<Self, CipError> {
        Self::new_impl(mol, root_atom, budget, Some(mancude), None)
    }

    fn new_impl(
        mol: &'m Molecule,
        root_atom: AtomIdx,
        budget: CipBudget,
        mancude: Option<&'m MancudeContext>,
        artificial_ancestor: Option<AtomIdx>,
    ) -> Result<Self, CipError> {
        let mut g = Self {
            mol,
            mancude,
            artificial_ancestor,
            nodes: Vec::new(),
            edges: Vec::new(),
            children_cache: Vec::new(),
            root: NodeId(0),
            budget,
            node_count: 0,
            expansion_count: 0,
            fractional_nodes_emitted: 0,
        };
        // Root has no incoming edge; the bond order passed here is unused (no edge is
        // created when `parent` is `None`).
        let root = g.push_node(
            CipNodeKind::Atom {
                atom_idx: root_atom,
            },
            None,
            0,
            BondOrder::Single,
        )?;
        g.root = root;
        Ok(g)
    }

    pub fn root(&self) -> NodeId {
        self.root
    }

    pub fn node(&self, id: NodeId) -> &CipNode {
        &self.nodes[id.0 as usize]
    }

    pub fn nodes(&self) -> &[CipNode] {
        &self.nodes
    }

    pub fn edges(&self) -> &[CipEdge] {
        &self.edges
    }

    pub fn molecule(&self) -> &'m Molecule {
        self.mol
    }

    /// How many constructed `MultipleBondDuplicate` nodes carry a real MANCUDE fraction
    /// (`AtomicNumberKey::Rational`), i.e. exercised the path Milestone 3B-1a/1b added --
    /// see the field doc for what this does and does not prove.
    pub fn fractional_nodes_emitted(&self) -> u64 {
        self.fractional_nodes_emitted
    }

    /// Compute (or return the cached) children of `node`. Duplicate and implicit-H
    /// nodes are always leaves: calling this on one returns `Ok(&[])` without touching
    /// the budget, since no expansion work happens.
    pub fn expand_children(&mut self, node_id: NodeId) -> Result<Vec<NodeId>, CipError> {
        if let Some(cached) = &self.children_cache[node_id.0 as usize] {
            return Ok(cached.clone());
        }

        let node = &self.nodes[node_id.0 as usize];
        let atom_idx = match node.kind {
            CipNodeKind::Atom { atom_idx } => atom_idx,
            // Duplicate/implicit-H nodes are leaves by design (see module docs) --
            // they represent a terminated branch, not a real atom to expand further.
            _ => {
                self.children_cache[node_id.0 as usize] = Some(Vec::new());
                return Ok(Vec::new());
            }
        };
        let depth = node.depth;
        let child_depth = depth + 1;

        self.expansion_count += 1;
        if self.expansion_count > self.budget.max_expansions {
            return Err(CipError::BudgetExceeded {
                reason: format!(
                    "expansion count exceeded max_expansions={}",
                    self.budget.max_expansions
                ),
            });
        }

        let parent_atom = node
            .parent
            .and_then(|p| match self.nodes[p.0 as usize].kind {
                CipNodeKind::Atom { atom_idx } => Some(atom_idx),
                _ => None,
            });
        let incoming_bond_order = node
            .incoming_edge
            .map(|e| self.edges[e.0 as usize].bond_order);
        let ancestors = self.ancestor_atoms(node_id);

        let mut children = Vec::new();

        // Arrival-side duplication: a multiple bond duplicates its partner into BOTH
        // atoms' own substituent lists, not just the one doing the iterating below
        // ("departure side"). This node's own list needs k-1 duplicates of its parent
        // if the edge that reached it was itself a multiple bond -- symmetric with the
        // departure-side duplication a few lines down. Getting this wrong (doing only
        // one side) is exactly the bug `d0e726b` fixed in the older, approximate engine
        // (crates/chematic-chem/src/cip.rs): a real atom "sees" its multiply-bonded
        // partner twice, from *both* atoms' own local perspective, regardless of which
        // direction a tree traversal happened to cross the bond.
        if let (Some(parent), Some(bond_order)) = (parent_atom, incoming_bond_order) {
            let multiplicity = bond_order.order_int();
            for _ in 0..multiplicity.saturating_sub(1) {
                let dup = self.push_node(
                    CipNodeKind::MultipleBondDuplicate {
                        source_atom: atom_idx,
                        duplicated_atom: parent,
                        bond_order: multiplicity,
                    },
                    Some(node_id),
                    child_depth,
                    bond_order,
                )?;
                children.push(dup);
            }
        }

        for (nb, bond_idx) in self.mol.neighbors(atom_idx) {
            if Some(nb) == parent_atom {
                // The bond we arrived through: not re-descended, not a ring closure --
                // just the single tree edge that already exists as `node`'s own
                // `incoming_edge`.
                continue;
            }
            let bond_order = self.mol.bond(bond_idx).order;
            let multiplicity = bond_order.order_int();

            // Departure-side duplication (the other half of the symmetric rule above):
            // k-1 extra leaves in *this* node's own list for each fresh neighbor found
            // via a multiple bond, regardless of whether `nb` turns out to be a fresh
            // atom or a ring closure (both can co-occur, e.g. a
            // ring-closing double bond in a bridged bicyclic alkene).
            for _ in 0..multiplicity.saturating_sub(1) {
                let dup = self.push_node(
                    CipNodeKind::MultipleBondDuplicate {
                        source_atom: atom_idx,
                        duplicated_atom: nb,
                        bond_order: multiplicity,
                    },
                    Some(node_id),
                    child_depth,
                    bond_order,
                )?;
                children.push(dup);
            }

            if ancestors.contains(&nb) {
                let dup = self.push_node(
                    CipNodeKind::RingDuplicate {
                        source_atom: atom_idx,
                        closure_atom: nb,
                    },
                    Some(node_id),
                    child_depth,
                    bond_order,
                )?;
                children.push(dup);
            } else {
                let real = self.push_node(
                    CipNodeKind::Atom { atom_idx: nb },
                    Some(node_id),
                    child_depth,
                    bond_order,
                )?;
                children.push(real);
            }
        }

        let h_count = implicit_hcount(self.mol, atom_idx);
        for _ in 0..h_count {
            let h = self.push_node(
                CipNodeKind::ImplicitHydrogen,
                Some(node_id),
                child_depth,
                BondOrder::Single,
            )?;
            children.push(h);
        }

        self.children_cache[node_id.0 as usize] = Some(children.clone());
        Ok(children)
    }

    /// The set of atoms on the root-to-`node` path, including `node`'s own atom (if
    /// it's an `Atom` node), plus [`Self::artificial_ancestor`] if one was seeded at
    /// construction (see [`Self::new_with_artificial_ancestor`]) -- unioned in
    /// unconditionally since it must terminate a matching neighbor at *any* depth in
    /// this digraph, not just at the root. Duplicate/implicit-H nodes never appear as
    /// ancestors -- they're leaves and can't have descendants to be an ancestor of.
    fn ancestor_atoms(&self, node_id: NodeId) -> HashSet<AtomIdx> {
        let mut set = HashSet::new();
        if let Some(seed) = self.artificial_ancestor {
            set.insert(seed);
        }
        let mut current = Some(node_id);
        while let Some(id) = current {
            let node = &self.nodes[id.0 as usize];
            if let CipNodeKind::Atom { atom_idx } = node.kind {
                set.insert(atom_idx);
            }
            current = node.parent;
        }
        set
    }

    fn push_node(
        &mut self,
        kind: CipNodeKind,
        parent: Option<NodeId>,
        depth: u32,
        edge_bond_order: BondOrder,
    ) -> Result<NodeId, CipError> {
        let prospective_count = self.node_count + 1;
        if prospective_count > self.budget.max_nodes {
            return Err(CipError::BudgetExceeded {
                reason: format!("node count exceeded max_nodes={}", self.budget.max_nodes),
            });
        }
        if depth as usize > self.budget.max_depth {
            return Err(CipError::BudgetExceeded {
                reason: format!("depth exceeded max_depth={}", self.budget.max_depth),
            });
        }
        self.node_count = prospective_count;

        let node_id = NodeId(self.nodes.len() as u32);
        let incoming_edge = parent.map(|p| {
            let edge_id = EdgeId(self.edges.len() as u32);
            self.edges.push(CipEdge {
                id: edge_id,
                parent: p,
                child: node_id,
                bond_order: edge_bond_order,
            });
            edge_id
        });
        let atomic_number = self.atomic_number_for(kind);
        if matches!(atomic_number, AtomicNumberKey::Rational(_)) {
            self.fractional_nodes_emitted += 1;
        }
        self.nodes.push(CipNode {
            id: node_id,
            kind,
            parent,
            depth,
            incoming_edge,
            atomic_number,
        });
        self.children_cache.push(None);
        Ok(node_id)
    }

    /// The `AtomicNumberKey` a freshly-constructed node of `kind` should carry. Real
    /// `Atom`, `RingDuplicate`, and `ImplicitHydrogen` nodes are always `Integral` --
    /// MANCUDE never touches them (see `crate::mancude`'s module docs: applying a
    /// fractional value to a real atom's own identity, rather than a duplicate that
    /// represents one of its resonance-averaged partners, would be incoherent). A
    /// `MultipleBondDuplicate` is `Rational` when its *owner* (`source_atom`) has a
    /// fractional atomic number in `self.mancude`, else it falls back to today's existing
    /// plain-integer behavior (the represented atom's real atomic number).
    fn atomic_number_for(&self, kind: CipNodeKind) -> AtomicNumberKey {
        match kind {
            CipNodeKind::Atom { atom_idx } => {
                AtomicNumberKey::Integral(self.mol.atom(atom_idx).element.atomic_number())
            }
            CipNodeKind::MultipleBondDuplicate {
                source_atom,
                duplicated_atom,
                ..
            } => self
                .mancude
                .and_then(|ctx| ctx.fractional_atomic_number(source_atom))
                .map(AtomicNumberKey::Rational)
                .unwrap_or_else(|| {
                    AtomicNumberKey::Integral(
                        self.mol.atom(duplicated_atom).element.atomic_number(),
                    )
                }),
            CipNodeKind::RingDuplicate { closure_atom, .. } => {
                AtomicNumberKey::Integral(self.mol.atom(closure_atom).element.atomic_number())
            }
            CipNodeKind::ImplicitHydrogen => AtomicNumberKey::Integral(1),
        }
    }

    /// Recursively expand every reachable node from `node` (typically the root),
    /// forcing full materialization. Milestone 1 has no ranking that would otherwise
    /// drive expansion, so tests that need to exercise the whole digraph (determinism,
    /// the residual-corpus smoke test) call this explicitly.
    pub fn expand_all(&mut self, node: NodeId) -> Result<(), CipError> {
        let children = self.expand_children(node)?;
        for child in children {
            self.expand_all(child)?;
        }
        Ok(())
    }

    /// A structural signature for the subtree rooted at `node`: an
    /// Aho-Hopcroft-Ullman-style recursive hash of `(node kind + relevant atom info,
    /// sorted child signatures)`. Two subtrees have equal signatures iff they have the
    /// same shape up to atom identity -- invariant to child traversal order (sorted
    /// before hashing) and, transitively, to the `AtomIdx`/`NodeId` numbering of either
    /// input. Used by the atom-renumbering and SMILES-representation invariance tests;
    /// forces full expansion of `node`'s subtree as a side effect.
    pub fn branch_signature(&mut self, node: NodeId) -> Result<u64, CipError> {
        let children = self.expand_children(node)?;
        let mut child_sigs = Vec::with_capacity(children.len());
        for child in children {
            child_sigs.push(self.branch_signature(child)?);
        }
        child_sigs.sort_unstable();

        let kind = self.nodes[node.0 as usize].kind;
        let mut hasher = DefaultHasher::new();
        match kind {
            CipNodeKind::Atom { atom_idx } => {
                0u8.hash(&mut hasher);
                let atom = self.mol.atom(atom_idx);
                atom.element.atomic_number().hash(&mut hasher);
                atom.isotope.hash(&mut hasher);
            }
            CipNodeKind::MultipleBondDuplicate {
                duplicated_atom,
                bond_order,
                ..
            } => {
                1u8.hash(&mut hasher);
                let atom = self.mol.atom(duplicated_atom);
                atom.element.atomic_number().hash(&mut hasher);
                atom.isotope.hash(&mut hasher);
                bond_order.hash(&mut hasher);
            }
            CipNodeKind::RingDuplicate { closure_atom, .. } => {
                2u8.hash(&mut hasher);
                let atom = self.mol.atom(closure_atom);
                atom.element.atomic_number().hash(&mut hasher);
                atom.isotope.hash(&mut hasher);
            }
            CipNodeKind::ImplicitHydrogen => {
                3u8.hash(&mut hasher);
            }
        }
        child_sigs.hash(&mut hasher);
        Ok(hasher.finish())
    }
}

/// Interface for lazily fetching a digraph node's children.
///
// ponytail: single implementor (`CipDigraph`) and no consumer yet -- this trait exists
// only so Milestone 2's ranking comparator can depend on an interface instead of a
// concrete type, per docs/rfcs/cip_accurate_rfc.md. If that never materializes, fold this
// back into an inherent method.
pub trait DigraphExpander {
    fn children(&mut self, node: NodeId) -> Result<Vec<NodeId>, CipError>;
}

impl<'m> DigraphExpander for CipDigraph<'m> {
    fn children(&mut self, node: NodeId) -> Result<Vec<NodeId>, CipError> {
        self.expand_children(node)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mancude::prepare_kekule_form;
    use crate::rational::RationalAtomicNumber;

    /// The concrete-value assertion that guards the owner-vs-represented-atom design
    /// decision (see `crate::mancude`'s module docs): quinoline's ring carbon directly
    /// bonded to N must produce a `MultipleBondDuplicate` node valued `Rational(13/2)` --
    /// the *owner*'s (this carbon's own) fractional atomic number -- never
    /// `Integral(6)`/`Integral(7)` (the represented atom's real, per-Kekulé-form value)
    /// or `Rational(19/3)`/`Rational(20/3)` (the oracle's global-Kekulé-enumeration mean,
    /// a different formula entirely -- see the divergence table in `crate::mancude`'s docs).
    #[test]
    fn new_with_mancude_quinoline_n_adjacent_carbon_duplicate_is_owner_fraction() {
        let mol = chematic_smiles::parse("n1ccc2ccccc2c1").unwrap();
        let (kekule_mol, mancude) = prepare_kekule_form(&mol).unwrap();

        // atom1 is N-adjacent (per quinoline's atom order, matching mancude.rs's own
        // divergence-table test). Root the digraph there so its own duplicate (from the
        // ring double bond touching it, in whichever direction apply_kekule resolved it)
        // is a *root child*, not requiring deeper expansion to reach.
        let atom1 = AtomIdx(1);
        let budget = CipBudget::default_budget();
        let mut graph = CipDigraph::new_with_mancude(&kekule_mol, atom1, budget, &mancude).unwrap();
        let root = graph.root();
        let children = graph.expand_children(root).unwrap();

        let duplicates: Vec<_> = children
            .iter()
            .filter_map(|&id| {
                let node = graph.node(id);
                matches!(node.kind, CipNodeKind::MultipleBondDuplicate { .. }).then_some(node)
            })
            .collect();
        assert_eq!(
            duplicates.len(),
            1,
            "atom1 has exactly one double bond in any Kekulé form, so exactly 1 duplicate"
        );
        let expected = RationalAtomicNumber::mean(&[6, 7]); // 13/2
        assert_eq!(
            duplicates[0].atomic_number,
            AtomicNumberKey::Rational(expected),
            "must be the OWNER's fraction (13/2), not the represented atom's real value \
             or the oracle's global-Kekule-enumeration mean"
        );
    }

    /// Without a `MancudeContext` (today's existing `CipDigraph::new` path, still what
    /// `assign_cip_accurate_experimental` calls), the same duplicate stays a plain
    /// integer -- confirms `new_with_mancude` is strictly additive, not a silent change
    /// to `new`'s existing behavior.
    #[test]
    fn new_without_mancude_keeps_plain_integer_duplicates() {
        let mol = chematic_smiles::parse("n1ccc2ccccc2c1").unwrap();
        let (kekule_mol, _mancude) = prepare_kekule_form(&mol).unwrap();
        let atom1 = AtomIdx(1);
        let budget = CipBudget::default_budget();
        let mut graph = CipDigraph::new(&kekule_mol, atom1, budget).unwrap();
        let root = graph.root();
        let children = graph.expand_children(root).unwrap();
        let duplicate = children
            .iter()
            .map(|&id| graph.node(id))
            .find(|node| matches!(node.kind, CipNodeKind::MultipleBondDuplicate { .. }))
            .unwrap();
        assert!(
            matches!(duplicate.atomic_number, AtomicNumberKey::Integral(_)),
            "no MancudeContext attached -- must stay the existing plain-integer behavior"
        );
    }

    /// Firing test for Milestone 3B's kept-but-measured-inert fractional path (see
    /// `docs/rfcs/cip_accurate_rfc.md`'s Milestone 3B closeout entry): "correct and never
    /// fired" needs a test that can tell the difference from "correct and never even
    /// reached." Uses a curated `aromatic_mancude`-bucket corpus molecule (atom 13 of
    /// `validation/cip_label_corpus.jsonl`'s first `aromatic_mancude` row) whose
    /// stereocenter is directly bonded to a MANCUDE-typed phenol ring, run through the
    /// same live path (`prepare_kekule_form` + `new_with_mancude` + `rank_children`) as
    /// `assign_cip_accurate_experimental`.
    #[test]
    fn live_path_fires_fractional_nodes_and_comparisons_on_curated_corpus_molecule() {
        use crate::compare::{CompareContext, rank_children};

        let mol =
            chematic_smiles::parse("C=CCCC[C@H](c1ccc(O)cc1)[C@@](C)(CC)c1ccc(O)cc1").unwrap();
        let (kekule_mol, mancude) = prepare_kekule_form(&mol).unwrap();
        let atom13 = AtomIdx(13);
        let budget = CipBudget::default_budget();
        let mut graph =
            CipDigraph::new_with_mancude(&kekule_mol, atom13, budget, &mancude).unwrap();
        let root = graph.root();
        let root_children = graph.expand_children(root).unwrap();

        let mut ctx = CompareContext::new();
        rank_children(&mut graph, &root_children, &mut ctx).unwrap();

        assert!(
            graph.fractional_nodes_emitted() > 0,
            "expected at least one MultipleBondDuplicate node with a MANCUDE fraction on \
             this molecule's phenol ring"
        );
        assert!(
            ctx.fractional_comparisons > 0,
            "expected at least one Rule-1a/2 comparison to reach a fractional key"
        );
        // fractional_decisions is NOT asserted > 0 here: this molecule's phenol rings are
        // pure-carbon-neighborhood positions (Rational(6/1), denominator 1), so the
        // fraction never actually decides anything -- 0 is the correct, expected value,
        // consistent with Milestone 3B-1b's own attribution finding. See
        // `CompareContext::fractional_decisions`'s doc comment.
    }

    /// [`CipDigraph::new_with_artificial_ancestor`]'s whole point: the seeded atom must
    /// terminate as a childless [`CipNodeKind::RingDuplicate`] leaf the moment it's
    /// reached as a neighbor, not expand as a real subtree -- checked concretely on one
    /// of Milestone 4B-0's own quinic-acid residual rows (atom4 is a ring neighbor of
    /// atom3; rooting at atom4 with atom3 as the artificial ancestor must produce a
    /// `RingDuplicate{closure_atom: atom3}` child, never an `Atom{atom_idx: atom3}` one).
    #[test]
    fn new_with_artificial_ancestor_terminates_seeded_neighbor_immediately() {
        let mol =
            chematic_smiles::parse("O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1")
                .unwrap();
        let root_atom = AtomIdx(4);
        let artificial_ancestor = AtomIdx(3);
        let budget = CipBudget::default_budget();
        let mut graph =
            CipDigraph::new_with_artificial_ancestor(&mol, root_atom, artificial_ancestor, budget)
                .unwrap();
        let root = graph.root();
        let children = graph.expand_children(root).unwrap();

        let seeded_child = children
            .iter()
            .find(|&&id| match graph.node(id).kind {
                CipNodeKind::Atom { atom_idx } => atom_idx == artificial_ancestor,
                CipNodeKind::RingDuplicate { closure_atom, .. } => {
                    closure_atom == artificial_ancestor
                }
                _ => false,
            })
            .copied()
            .expect("artificial_ancestor atom must appear among root's children");

        assert!(
            matches!(
                graph.node(seeded_child).kind,
                CipNodeKind::RingDuplicate { .. }
            ),
            "artificial_ancestor must terminate as a RingDuplicate leaf immediately, not \
             expand as a real Atom subtree"
        );
        assert_eq!(
            graph.expand_children(seeded_child).unwrap(),
            Vec::new(),
            "a RingDuplicate node is always childless"
        );
    }

    /// `new`'s existing behavior (no artificial ancestor) is unaffected: the same root
    /// atom's real neighbor set is unchanged when `new_with_artificial_ancestor` is not
    /// used -- confirms the new field is strictly additive, mirroring the
    /// `new_without_mancude_keeps_plain_integer_duplicates` precedent above.
    #[test]
    fn new_without_artificial_ancestor_keeps_ordinary_neighbor_expansion() {
        let mol =
            chematic_smiles::parse("O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1")
                .unwrap();
        let root_atom = AtomIdx(4);
        let budget = CipBudget::default_budget();
        let mut graph = CipDigraph::new(&mol, root_atom, budget).unwrap();
        let root = graph.root();
        let children = graph.expand_children(root).unwrap();

        let atom3_child = children
            .iter()
            .find(|&&id| matches!(graph.node(id).kind, CipNodeKind::Atom { atom_idx } if atom_idx == AtomIdx(3)))
            .copied();
        assert!(
            atom3_child.is_some(),
            "without an artificial ancestor, atom3 must appear as an ordinary real Atom \
             child, exactly today's existing behavior"
        );
    }
}
