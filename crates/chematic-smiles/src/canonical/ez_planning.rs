impl<'a> CanonicalWriter<'a> {
    pub(crate) fn new(mol: &'a Molecule, ranks: &'a [u64]) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            ranks,
            written: vec![false; n],
            ring_bonds: vec![false; mol.bond_count()],
            atom_ring_nums: vec![Vec::new(); n],
            ring_marker_on_close: HashSet::new(),
            small_ring_double: std::cell::RefCell::new(HashMap::new()),
            next_ring: 1,
            out: String::with_capacity(n.saturating_mul(4) + mol.bond_count().saturating_mul(2)),
            order: Vec::with_capacity(n),
            ez_group: HashMap::new(),
            ez_flip: HashMap::new(),
            forced_ez_flip: None,
            ez_marker: HashMap::new(),
            ez_geometry_facts: Vec::new(),
            #[cfg(test)]
            ez_shared_bond_abstains: Vec::new(),
            #[cfg(test)]
            traversal_bond_preference: None,
            direct_ez_orders: None,
            canonical_tree_parent: if cfg!(test)
                || mol.atoms().any(|(idx, _)| mol.is_allene_stereo(idx))
            {
                vec![None; n]
            } else {
                Vec::new()
            },
        }
    }

    /// The raw parse-time direction of `bidx`, if any: a literal `Up`/`Down`
    /// bond order, or (when both endpoints are aromatic) the stashed
    /// direction of an adjacent exocyclic double bond. This is the
    /// *unresolved* reading — [`Self::effective_order`] is what emission
    /// code should use instead, since it additionally applies
    /// `resolve_ez_markers`'s carrier choice on top of this.
    ///
    /// Delegates to `crate::writer::raw_bond_direction` -- the exact same
    /// rule the plain (non-canonical) writer uses, kept as one shared
    /// implementation on purpose so the two writers can never silently
    /// disagree on what a bond's stored/stashed direction means.
    fn raw_input_direction(&self, bidx: BondIdx) -> Option<BondOrder> {
        crate::writer::raw_bond_direction(self.mol, bidx)
    }

    /// Strip directionality from a bond order, leaving its "plain" chemical
    /// order untouched: `Up`/`Down` becomes `Single` (the literal-marker
    /// case), anything else (notably `Aromatic`, the stash case) is
    /// returned as-is. Used to demote a substituent bond that must no
    /// longer carry the `/`/`\` marker in the output.
    fn plain_order(order: BondOrder) -> BondOrder {
        if matches!(order, BondOrder::Up | BondOrder::Down) {
            BondOrder::Single
        } else {
            order
        }
    }

    /// The effective direction to use when emitting `bidx`: `resolve_ez_markers`'s
    /// resolved choice if this bond participates in a tri-/tetra-substituted
    /// stereo alkene end, otherwise the unresolved raw parse-time direction
    /// (literal order or aromatic-bond-direction stash), otherwise the
    /// bond's own real chemical order. Every write-time site that needs "the
    /// order to show for this bond" (ring-closure emission and tree-edge
    /// emission alike) must go through this, not read `bond_direction`/
    /// `order` directly, or the two sites can disagree on which bond a
    /// moved E/Z marker landed on.
    /// Drop a `/`/`\\` token that takes part in no E/Z specification: a
    /// standalone wedge, or a marker flanking only double bonds whose other
    /// end will carry no marker in this output (RDKit drops such tokens, so
    /// writing one would make the canonical form unstable).
    fn suppress_lone_marker(&self, bidx: BondIdx, order: BondOrder) -> BondOrder {
        let order = suppress_standalone_wedge(self.mol, bidx, order);
        if !matches!(order, BondOrder::Up | BondOrder::Down) {
            return order;
        }
        let bond = self.mol.bond(bidx);
        let small_ring_double = |db: BondIdx| {
            let mut cache = self.small_ring_double.borrow_mut();
            *cache
                .entry(db)
                .or_insert_with(|| crate::writer::in_ring_smaller_than_eight(self.mol, db))
        };
        let specifies = [bond.atom1, bond.atom2].into_iter().any(|end| {
            self.mol.neighbors(end).any(|(far, db)| {
                db != bidx
                    && self.mol.bond(db).order == BondOrder::Double
                    && self.mol.neighbors(far).any(|(_, sb)| {
                        sb != db
                            && sb != bidx
                            && matches!(self.effective_order(sb), BondOrder::Up | BondOrder::Down)
                    })
                    && !small_ring_double(db)
            })
        });
        if specifies {
            order
        } else {
            Self::plain_order(self.mol.bond(bidx).order)
        }
    }

    fn effective_order(&self, bidx: BondIdx) -> BondOrder {
        if let Some(orders) = &self.direct_ez_orders
            && let Some(&order) = orders.get(&bidx)
        {
            return order;
        }
        if let Some(&resolved) = self.ez_marker.get(&bidx) {
            return resolved;
        }
        self.raw_input_direction(bidx)
            .unwrap_or(self.mol.bond(bidx).order)
    }

    /// Return `true` if `dir` encodes "up" (`/`, read atom1→atom2) as seen
    /// from `alkene_end`'s side of `bond` — the same atom1/atom2-relative
    /// convention `chematic_chem::cip::substituent_is_up` uses to read a
    /// marker back out, so a value written here round-trips identically.
    fn direction_is_up(dir: BondOrder, bond_atom1: AtomIdx, alkene_end: AtomIdx) -> bool {
        match dir {
            BondOrder::Up => bond_atom1 == alkene_end,
            BondOrder::Down => bond_atom1 != alkene_end,
            _ => false, // never called with a non-directional `dir`
        }
    }

    /// The inverse of `direction_is_up`: which `Up`/`Down` order to store on
    /// `bond` so that, read from `alkene_end`'s side, it encodes `want_up`.
    fn direction_for_up(bond_atom1: AtomIdx, alkene_end: AtomIdx, want_up: bool) -> BondOrder {
        if (bond_atom1 == alkene_end) == want_up {
            BondOrder::Up
        } else {
            BondOrder::Down
        }
    }

    /// Return the endpoint from which a raw direction was originally read.
    /// Literal `Up`/`Down` orders are stored relative to the bond's atom1;
    /// aromatic direction stashes retain the parser-side source endpoint
    /// separately. Carrier selection must use the latter or two equivalent
    /// spellings can encode the same geometry with opposite internal bond
    /// orientation.
    fn raw_direction_anchor(&self, bidx: BondIdx) -> AtomIdx {
        self.mol
            .bond_direction_anchor(bidx)
            .unwrap_or(self.mol.bond(bidx).atom1)
    }

    /// Maximum coupled-component size [`Self::resolve_component_jointly`]
    /// will attempt to jointly resolve via bounded enumeration
    /// (`2^size` candidate assignments). Measured, not guessed, and fully
    /// reproducible from a corpus committed to this repo (no external
    /// download required):
    ///
    /// ```text
    /// cargo run -p chematic-smiles --release --example ez_shared_carrier_component_audit -- \
    ///     scripts/fixtures/descriptor_census_corpus.smi
    /// ```
    ///
    /// Expected output (re-run to confirm, not assumed): 18/18 pinned
    /// fixtures each produce exactly one size-2 component; the 5,000-line
    /// committed corpus produces 31 coupling components total, still every
    /// one of them size exactly 2, 0 cycles. (An independent run against
    /// the larger, non-committed ChEMBL corpus used for this crate's own
    /// full-corpus residual measurements — see
    /// `docs/canonical-ez-plan.md` — found the same: every
    /// component observed has been size exactly 2.) Every node has at most
    /// 2 candidate substituent bonds, so every component is a simple path
    /// or cycle, never a general graph. This cap gives an 8x margin over
    /// the measured maximum. Worst case is bounded at `2^16` = 65,536
    /// candidate-assignment evaluations, each O(size) with small constants
    /// (no independent wall-clock measurement is claimed here — every
    /// corpus scanned so far never exceeds 4 evaluations, i.e. component
    /// size 2, so the `2^16` bound is a safety margin, not an observed
    /// cost). A component that exceeds the cap abstains as a whole --
    /// never partially applied.
    const MAX_JOINT_COMPONENT_SIZE: usize = 16;

    /// Resolve, for every tri-/tetra-substituted stereo alkene end, which ONE
    /// of its ≥2 non-double-bond substituent bonds should carry the `/`/`\`
    /// marker in the output — deterministically, from the molecule's own
    /// (already fully individualized) canonical ranks, rather than just
    /// inheriting whichever bond the original parse happened to mark.
    ///
    /// **The bug this fixes** (canonical-residual RFC, Root cause 1): SMILES
    /// only requires ONE of an alkene carbon's ≥2 substituent bonds to carry
    /// an explicit marker (the geometry of the other is implied — at a
    /// trigonal alkene carbon with two substituents, they are always on
    /// opposite sides). *Which* substituent gets marked is a free choice at
    /// write time, but the writer previously just checked "does this
    /// specific bond already carry a direction", inherited straight from
    /// parse time — so two RDKit-valid respellings of the identical molecule
    /// that mark different substituents produced two different (but
    /// chemically identical) canonical outputs. Resolving the carrier here,
    /// once, from topology + rank alone, makes the choice depend only on the
    /// molecule itself, not on which substituent the input happened to mark.
    ///
    /// Two independently stereogenic double bonds can additionally share one
    /// physical candidate substituent bond between them (issue #149 —
    /// typically a ring-closure bond connecting two adjacent stereocenters
    /// each bearing their own exocyclic stereo double bond). Resolving one
    /// end's carrier in isolation from such a coupled partner risks moving
    /// or demoting a mark the *other* system's own resolution simultaneously
    /// relies on for the same physical bond. [`Self::coupling_components`]
    /// partitions every ambiguous end into maximal connected clusters joined
    /// by exactly this kind of shared-bond edge; [`Self::resolve_component_
    /// jointly`] solves each cluster as one unit (a lone, uncoupled end is
    /// simply a cluster of size 1, resolved by the same code path at
    /// negligible extra cost — see that method's own doc comment).
    ///
    /// Runs before `build_ez_groups`/ring discovery/DFS write, so it depends
    /// only on molecule topology, never on write order.
    fn resolve_ez_markers(&mut self) {
        if self.ranks.is_empty() {
            return;
        }
        self.ez_geometry_facts = self.extract_ez_geometry_facts();
        let stereo_alkene_ends = Self::compute_stereo_alkene_ends(self.mol);
        for component in Self::coupling_components(self.mol, &stereo_alkene_ends) {
            self.resolve_component_jointly(&component);
        }
    }

    fn end_has_substituent(mol: &Molecule, end: AtomIdx) -> bool {
        mol.neighbors(end)
            .any(|(_, b)| mol.bond(b).order != BondOrder::Double)
    }

    #[inline]
    fn substituent_count(mol: &Molecule, end: AtomIdx) -> usize {
        mol.neighbors(end)
            .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
            .count()
    }

    /// Non-double-bond neighbors of `end` — an alkene carbon's up-to-two
    /// sigma substituents. Filtered by bond order, not by comparing against
    /// a specific double-bond `BondIdx`, so an allene/cumulene terminus
    /// correctly excludes *every* double bond at `end`, matching
    /// `chematic_chem::cip::assign_ez`'s own substituent-collection
    /// convention.
    fn substituents(mol: &Molecule, end: AtomIdx) -> Vec<(AtomIdx, BondIdx)> {
        mol.neighbors(end)
            .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
            .collect()
    }

    /// Whether the double bond between `a` and `b` is endocyclic in a ring
    /// smaller than 8 atoms — i.e. some SSSR ring contains both atoms and
    /// has fewer than 8 members. Such a bond's real-world geometry is fixed
    /// by the ring itself, not a free stereochemical choice (issue #149's
    /// ring-constrained residual, `docs/rfcs/ez_ring_constrained_residual_
    /// audit.md`): RDKit independently agrees the bond is not a real
    /// stereocenter (`STEREONONE`/absent from `FindPotentialStereo`) with
    /// 0/1,387 row-level disagreements on a 5,000-molecule corpus at exactly
    /// this threshold. Deliberately checks ring membership of the BOND (both
    /// endpoints in the same ring), not merely of the atom — a ring carbon
    /// with a genuinely free *exocyclic* double bond (its own C=X pointing
    /// outward) must not be excluded, and the audit found the naive
    /// atom-membership version over-excludes by ~9 percentage points of a
    /// real corpus for exactly this reason.
    fn double_bond_endocyclic_in_small_ring(
        rings: &[Vec<AtomIdx>],
        a: AtomIdx,
        b: AtomIdx,
    ) -> bool {
        rings
            .iter()
            .any(|r| r.len() < 8 && r.contains(&a) && r.contains(&b))
    }

    /// Every atom that is a stereogenic double bond's terminus with exactly
    /// 2 candidate substituents — the ambiguity precondition a marker-carrier
    /// choice exists for at all (0 or 1 substituent has nothing to choose
    /// between; a real double bond's end never has more than 2). A double
    /// bond only counts when BOTH ends have at least one substituent
    /// (matching `chematic_chem::cip::assign_ez`'s own
    /// `subs_a1.is_empty() || subs_a2.is_empty()` guard — e.g. a
    /// ketone/aldehyde `C=O` never does, the O side has none), so the carbon
    /// side of such a bond is never treated as ambiguous purely because it
    /// happens to have 2 ring-bond substituents of its own.
    ///
    /// Also excludes both ends of a double bond that is itself endocyclic in
    /// a ring smaller than 8 atoms (see
    /// [`Self::double_bond_endocyclic_in_small_ring`]) — such a bond has no
    /// real stereochemical freedom, so treating it as ambiguous only
    /// destabilizes marker-carrier selection for a genuinely stereogenic
    /// bond it happens to share a candidate bond with (issue #149).
    ///
    /// This runs on every `canonical_smiles()` call via `resolve_ez_markers`,
    /// so `chematic_perception::find_sssr` (Horton: O(V·E) candidates + GF(2)
    /// elimination) is only ever computed when there is at least one
    /// double-bond candidate that could need the ring check — most
    /// molecules (no double bonds, or only non-candidate ones like a
    /// ketone's `C=O`) skip it entirely.
    fn compute_stereo_alkene_ends(mol: &Molecule) -> HashSet<AtomIdx> {
        let candidates: Vec<_> = mol
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .filter(|(_, bond)| {
                let count1 = Self::substituent_count(mol, bond.atom1);
                let count2 = Self::substituent_count(mol, bond.atom2);
                count1 > 0 && count2 > 0 && (count1 == 2 || count2 == 2)
            })
            .collect();
        if candidates.is_empty() {
            return HashSet::new();
        }

        // The small-ring exclusion below can only apply when the candidate
        // double bond is itself a cycle edge.  Most stereo double bonds are
        // exocyclic or acyclic; avoid running the substantially more
        // expensive whole-molecule SSSR algorithm for those.  Removing an
        // undirected edge and checking endpoint reachability is an exact
        // cycle-membership test.
        let needs_ring_check = candidates.iter().any(|&(candidate_idx, bond)| {
            let mut seen = vec![false; mol.atom_count()];
            let mut stack = vec![bond.atom1];
            seen[bond.atom1.0 as usize] = true;
            while let Some(atom) = stack.pop() {
                for (neighbor, edge) in mol.neighbors(atom) {
                    if edge == candidate_idx || seen[neighbor.0 as usize] {
                        continue;
                    }
                    if neighbor == bond.atom2 {
                        return true;
                    }
                    seen[neighbor.0 as usize] = true;
                    stack.push(neighbor);
                }
            }
            false
        });
        let ring_set = needs_ring_check.then(|| chematic_perception::find_sssr(mol));
        let rings = ring_set.as_ref().map(|set| set.rings());
        // A double bond with two equivalent substituents at one end (=CF2,
        // =C(C#N)2, =CD2) is not stereogenic; its ends must not take part in
        // carrier resolution, or a marker merely adjacent to it would steer
        // the choice for a real stereo bond next to it.
        // The equivalence classes are computed only when two substituents
        // could be equivalent at all (same element, degree and H count).
        let classes: std::cell::OnceCell<std::sync::Arc<Vec<usize>>> = std::cell::OnceCell::new();
        // Without any `/`/`\\` marker there is no carrier to resolve.
        let has_markers = mol
            .bonds()
            .any(|(b, _)| crate::writer::raw_bond_direction(mol, b).is_some());
        let symmetric_end = |end: AtomIdx, other: AtomIdx| {
            if !has_markers {
                return false;
            }
            let mut subs = mol
                .neighbors(end)
                .filter(|&(nb, _)| nb != other)
                .map(|(nb, _)| nb);
            let (Some(a), Some(b), None) = (subs.next(), subs.next(), subs.next()) else {
                return false;
            };
            let (x, y) = (mol.atom(a), mol.atom(b));
            if x.element != y.element
                || x.isotope != y.isotope
                || x.charge != y.charge
                || mol.degree(a) != mol.degree(b)
                || chematic_core::implicit_hcount(mol, a) != chematic_core::implicit_hcount(mol, b)
            {
                return false;
            }
            let classes = classes.get_or_init(|| {
                mol.derived(chematic_core::DerivedSlot::TopologicalClasses, || {
                    crate::canonical_partition::topological_equivalence_classes(mol)
                })
            });
            classes[a.0 as usize] == classes[b.0 as usize]
        };
        let mut ends = HashSet::new();
        for (_, bond) in candidates {
            if rings.is_some_and(|rings| {
                Self::double_bond_endocyclic_in_small_ring(rings, bond.atom1, bond.atom2)
            }) {
                continue;
            }
            if symmetric_end(bond.atom1, bond.atom2) || symmetric_end(bond.atom2, bond.atom1) {
                continue;
            }
            for end in [bond.atom1, bond.atom2] {
                if Self::substituent_count(mol, end) == 2 {
                    ends.insert(end);
                }
            }
        }
        ends
    }

    /// Partition every stereo-alkene-end atom into maximal connected
    /// components of the shared-candidate-bond coupling graph: nodes are the
    /// ends themselves, edges connect two ends directly bonded via one of
    /// their own (≤2) candidate substituent bonds — the only way one
    /// physical bond can be a marker-carrier candidate for two different
    /// double bonds' own resolution at once (a candidate bond is by
    /// definition incident to its owning end, so it can only ever be
    /// "shared" with the atom on its *other* side). Since every node has ≤2
    /// candidate bonds, every component is a simple path or cycle, never a
    /// general graph (measured directly — see [`Self::MAX_JOINT_COMPONENT_
    /// SIZE`]'s doc comment). A molecule with no coupling at all yields one
    /// singleton component per ambiguous end (the overwhelming common
    /// case). Component *membership* is a pure function of molecule
    /// topology; the order components (and atoms within them) are returned
    /// in has no effect on the result — [`Self::resolve_component_jointly`]
    /// re-sorts by canonical rank internally before making any decision.
    ///
    /// Traversal order (which node BFS starts from, and the order its
    /// neighbors are pushed) is made deterministic by sorting on raw
    /// `AtomIdx` — not because that order is canonical (it isn't), but
    /// because `ends` is a `HashSet` and iterating it directly would leak
    /// process-random hash-seed order into `component`'s element order.
    /// [`Self::resolve_component_jointly`]'s own rank-based sort washes
    /// that out for every rank-*distinct* member regardless, but ties in
    /// `self.ranks` are only possible for genuinely automorphic ends, and
    /// leaving *that* residual order-dependence to an unseeded HashSet
    /// would be process-random input to an otherwise-deterministic
    /// decision — so it is removed here at the source instead.
    fn coupling_components(mol: &Molecule, ends: &HashSet<AtomIdx>) -> Vec<Vec<AtomIdx>> {
        let mut starts: Vec<AtomIdx> = ends.iter().copied().collect();
        starts.sort_by_key(|a| a.0);
        let mut visited: HashSet<AtomIdx> = HashSet::new();
        let mut components = Vec::new();
        for start in starts {
            if visited.contains(&start) {
                continue;
            }
            let mut queue = vec![start];
            let mut component = Vec::new();
            visited.insert(start);
            while let Some(cur) = queue.pop() {
                component.push(cur);
                let mut nbs: Vec<AtomIdx> = Self::substituents(mol, cur)
                    .into_iter()
                    .filter(|(nb, _)| ends.contains(nb))
                    .map(|(nb, _)| nb)
                    .collect();
                nbs.sort_by_key(|a| a.0);
                for nb in nbs {
                    if visited.insert(nb) {
                        queue.push(nb);
                    }
                }
            }
            components.push(component);
        }
        components
    }

    /// Whether a bond of this order can be written with `/` or `\` (a
    /// zero-order or dative bond cannot; it would print as single).
    fn can_carry_marker(order: BondOrder) -> bool {
        matches!(
            order,
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Aromatic
        )
    }

    /// The substituent an alkene end prefers as its marker carrier: the
    /// lower-ranked one among bonds that can carry a marker, except that a
    /// heavy atom is preferred over an
    /// explicit hydrogen atom. RDKit's `RemoveHs` moves a marker off a
    /// removed H onto the end's other single bond; when that bond also
    /// flanks a conjugated double bond, the moved marker can contradict
    /// that double bond's own markers, and RDKit then drops or inverts its
    /// E/Z (#734: 22 of 572 E/Z molecules written with explicit H were read
    /// back differently).
    ///
    /// A heavy substituent that is itself the end of a double bond whose
    /// other substituent is an H atom is preferred too: that end carries
    /// its marker on this shared bond, and both ends choosing it keeps the
    /// choice free of an H carrier.
    fn preferred_carrier(&self, pair: &[(AtomIdx, BondIdx); 2]) -> usize {
        let is_h = |a: AtomIdx| self.mol.atom(a).element == chematic_core::Element::H;
        let far_end_has_h = |far: AtomIdx| {
            self.mol
                .neighbors(far)
                .any(|(_, b)| self.mol.bond(b).order == BondOrder::Double)
                && self.mol.neighbors(far).any(|(nb, _)| is_h(nb))
        };
        let key = |i: usize| {
            let far = pair[i].0;
            (
                !Self::can_carry_marker(self.mol.bond(pair[i].1).order),
                is_h(far),
                !far_end_has_h(far),
                self.ranks[pair[i].0.0 as usize],
            )
        };
        if key(1) < key(0) { 1 } else { 0 }
    }

    /// Jointly resolve one coupling component (see
    /// [`Self::coupling_components`]) — a cluster of ambiguous alkene ends,
    /// possibly of size 1 (no coupling at all, the common case: this method
    /// is the *only* code path that assigns marker carriers, unified rather
    /// than kept as a separate fast path, since a singleton component can
    /// never conflict with anything and this method's enumeration collapses
    /// to exactly the old single-end algorithm's own choice at negligible
    /// extra cost — 2 evaluations instead of a direct pick).
    ///
    /// Enumerates every combination of "which of its 2 candidate bonds does
    /// end `i` use as its own marker carrier" (bounded by
    /// [`Self::MAX_JOINT_COMPONENT_SIZE`]), keeps only combinations where
    /// every physical bond in the component receives a single, mutually
    /// consistent literal value from every end that votes on it (see
    /// [`Self::evaluate_choice`]), and — among those — picks the one
    /// deviating from each end's own canonical-rank-preferred choice the
    /// *fewest* times (primary key: total deviation count, ascending), with
    /// a **molecule-intrinsic** tie-break among equal-count combinations:
    /// ends sorted by ascending `self.ranks` (never `AtomIdx`/`BondIdx`/
    /// traversal order), comparing each combination's per-end "did I
    /// deviate" bit vector in that order, lexicographically smallest wins
    /// (deviations are pushed onto the higher-ranked end first).
    ///
    /// That two-part key is a bijection of `choice` (distinct combinations
    /// always produce distinct keys), so it alone can never report a tie —
    /// a real one still exists when two *different* minimal-count
    /// combinations differ only by which of two `self.ranks`-equal (hence
    /// provably automorphic) ends carries a given deviation: no signal this
    /// molecule offers distinguishes them, so [`Self::coupling_components`]'s
    /// deterministic-but-arbitrary traversal order must not be allowed to
    /// pick a "winner" by accident. Detected explicitly below rather than
    /// left to fall out of the key comparison. If no consistent combination
    /// exists at all, or a genuine tie is detected, the whole component is
    /// left exactly as the input spelled it: always safe, never guessed.
    fn resolve_component_jointly(&mut self, component: &[AtomIdx]) {
        let mut ordered: Vec<AtomIdx> = component.to_vec();
        ordered.sort_by_key(|&a| self.ranks[a.0 as usize]);
        let k = ordered.len();

        if k > Self::MAX_JOINT_COMPONENT_SIZE {
            #[cfg(test)]
            self.ez_shared_bond_abstains.extend_from_slice(&ordered);
            return; // too large to jointly resolve -- leave untouched, never partially applied
        }

        let subs: Vec<[(AtomIdx, BondIdx); 2]> = ordered
            .iter()
            .map(|&end| {
                let s = Self::substituents(self.mol, end);
                debug_assert_eq!(s.len(), 2, "component member must be a 2-substituent end");
                [s[0], s[1]]
            })
            .collect();

        if k == 1 {
            if self.resolve_singleton_component(ordered[0], subs[0]) {
                return;
            }
            #[cfg(test)]
            self.ez_shared_bond_abstains.push(ordered[0]);
            return;
        }
        let Some(plan) = self.best_joint_plan(&ordered, &subs) else {
            #[cfg(test)]
            self.ez_shared_bond_abstains.extend_from_slice(&ordered);
            return;
        };
        for (bidx, order) in plan {
            self.ez_marker.insert(bidx, order);
        }
    }

    fn resolve_singleton_component(
        &mut self,
        end: AtomIdx,
        pair: [(AtomIdx, BondIdx); 2],
    ) -> bool {
        let preferred = self.preferred_carrier(&pair);
        let Some(reference_up) = self.geometry_reference_up(end, &pair, preferred) else {
            return false;
        };
        for chosen in [preferred, 1 - preferred] {
            let demoted = pair[1 - chosen].1;
            if self.is_load_bearing_elsewhere(demoted, end) {
                continue;
            }
            let chosen_up = chosen == preferred && reference_up || chosen != preferred && !reference_up;
            let chosen_bond = pair[chosen].1;
            let order = Self::direction_for_up(self.mol.bond(chosen_bond).atom1, end, chosen_up);
            self.ez_marker.insert(chosen_bond, order);
            self.ez_marker
                .insert(demoted, Self::plain_order(self.mol.bond(demoted).order));
            return true;
        }
        false
    }

    fn best_joint_plan(
        &self,
        ordered: &[AtomIdx],
        substituents: &[[(AtomIdx, BondIdx); 2]],
    ) -> Option<HashMap<BondIdx, BondOrder>> {
        let preferred: Vec<usize> = substituents
            .iter()
            .map(|pair| self.preferred_carrier(pair))
            .collect();
        let reference_up: Vec<Option<bool>> = ordered
            .iter()
            .zip(substituents)
            .zip(&preferred)
            .map(|((&end, pair), &choice)| self.geometry_reference_up(end, pair, choice))
            .collect();
        let mut valid = Vec::new();
        for mask in 0u32..(1u32 << ordered.len()) {
            let choices: Vec<usize> = (0..ordered.len())
                .map(|index| ((mask >> index) & 1) as usize)
                .collect();
            if let Some(plan) = self.evaluate_choice(
                ordered,
                substituents,
                &preferred,
                &reference_up,
                &choices,
            ) {
                let deviations = choices
                    .iter()
                    .zip(&preferred)
                    .map(|(choice, preferred)| u8::from(choice != preferred))
                    .collect::<Vec<_>>();
                valid.push((deviations, plan));
            }
        }
        let minimum = valid.iter().map(|(bits, _)| deviation_count(bits)).min()?;
        let mut winners: Vec<_> = valid
            .iter()
            .filter(|(bits, _)| deviation_count(bits) == minimum)
            .collect();
        winners.sort_by(|left, right| left.0.cmp(&right.0));
        let (winner_bits, winner_plan) = winners[0];
        (!self.has_automorphic_joint_tie(ordered, winner_bits, &winners))
            .then(|| winner_plan.clone())
    }

    fn has_automorphic_joint_tie(
        &self,
        ordered: &[AtomIdx],
        winner: &[u8],
        candidates: &[&(Vec<u8>, HashMap<BondIdx, BondOrder>)],
    ) -> bool {
        (0..ordered.len()).any(|left| {
            (left + 1..ordered.len()).any(|right| {
                self.ranks[ordered[left].0 as usize] == self.ranks[ordered[right].0 as usize]
                    && winner[left] != winner[right]
                    && candidates.iter().any(|(bits, _)| {
                        bits[left] == winner[right]
                            && bits[right] == winner[left]
                            && (0..ordered.len()).all(|index| {
                                index == left || index == right || bits[index] == winner[index]
                            })
                    })
            })
        })
    }

    /// One end's own geometry fact: whether its rank-preferred ("reference")
    /// candidate substituent (`subs[pref_idx]`) reads as "up", read via a
    /// FIXED reference — the reference bond's own raw direction if it has
    /// one, else its sibling's raw direction inverted (the trigonal-carbon
    /// sibling-complement identity: the two substituents of a stereogenic
    /// alkene end always sit on opposite sides of the double-bond axis).
    /// Returns `None` only when *neither* candidate carries any raw
    /// direction at all (no geometry info at this end).
    ///
    /// This is deliberately never "whichever candidate happens to carry a
    /// raw mark" — always the SAME (rank-based) candidate — so its result
    /// depends only on molecule topology and encoded geometry, never on
    /// which of the two candidates the CURRENT input text happened to mark.
    /// Mirrors the test-only `up_of_reference` oracle's own logic
    /// (kept as two copies, not shared, because that one intentionally
    /// stays test-only scaffolding — see its doc comment).
    fn reference_up(
        &self,
        alkene_end: AtomIdx,
        subs: &[(AtomIdx, BondIdx); 2],
        pref_idx: usize,
    ) -> Option<bool> {
        let reference = subs[pref_idx];
        let sibling = subs[1 - pref_idx];
        if let Some(dir) = self.raw_input_direction(reference.1) {
            Some(Self::direction_is_up(
                dir,
                self.raw_direction_anchor(reference.1),
                alkene_end,
            ))
        } else {
            let dir = self.raw_input_direction(sibling.1)?;
            Some(!Self::direction_is_up(
                dir,
                self.raw_direction_anchor(sibling.1),
                alkene_end,
            ))
        }
    }

    /// Extract the E/Z facts encoded by the input before choosing a canonical
    /// marker carrier or a DFS/ring representation.  This is deliberately a
    /// read-only phase: it gives the eventual component-plan solver chemical
    /// relations to preserve, rather than parse-time marker locations to
    /// imitate.
    fn extract_ez_geometry_facts(&self) -> Vec<EzGeometryFact> {
        let mut facts = Vec::new();
        for (double_bond, bond) in self.mol.bonds() {
            if bond.order != BondOrder::Double {
                continue;
            }
            let Some(left) = self.extract_ez_geometry_end(bond.atom1) else {
                continue;
            };
            let Some(right) = self.extract_ez_geometry_end(bond.atom2) else {
                continue;
            };
            let mut double_key = (
                self.ranks[bond.atom1.0 as usize],
                self.ranks[bond.atom2.0 as usize],
            );
            if double_key.1 < double_key.0 {
                double_key = (double_key.1, double_key.0);
            }
            facts.push(EzGeometryFact {
                double_bond,
                double_key,
                ends: [left, right],
                same_side: left.up == right.up,
            });
        }
        facts.sort_by_key(|fact| fact.double_key);
        facts
    }

    /// Every physical bond whose token may need an explicit choice for an
    /// aromatic direction-stash E/Z system.  Coupled alkene substituents
    /// provide the structural candidates; raw carriers are included as well
    /// because a stash can legally be emitted from a bond outside that local
    /// component.  The latter is the essential #503 boundary: excluding it
    /// leaves no common semantic output for the held-out residuals.
    fn complete_ez_plan_slots(&self) -> Vec<BondIdx> {
        let ends = Self::compute_stereo_alkene_ends(self.mol);
        let mut slots: Vec<_> = Self::coupling_components(self.mol, &ends)
            .into_iter()
            .flatten()
            .flat_map(|end| Self::substituents(self.mol, end))
            .map(|(_, bond)| bond)
            .collect();
        slots.extend((0..self.mol.bond_count()).filter_map(|index| {
            let bond = BondIdx(index as u32);
            matches!(
                self.raw_input_direction(bond),
                Some(BondOrder::Up | BondOrder::Down)
            )
            .then_some(bond)
        }));
        slots.sort_unstable_by_key(|bond| bond.0);
        slots.dedup();
        slots
    }

    /// Enumerate the finite directional-token degrees of freedom left after
    /// canonical DFS/ring discovery, retain only reparsable candidates with
    /// exactly the original E/Z semantics, and select their lexicographic
    /// minimum. Enumeration order is intentionally irrelevant: no atom or
    /// bond index is a tie-breaker in the resulting public output.
    fn complete_aromatic_ez_plan(&self) -> Option<String> {
        const MAX_SLOTS: usize = 8;
        const MAX_CANDIDATES: usize = 65_536;

        let slots = self.complete_ez_plan_slots();
        if slots.is_empty() || slots.len() > MAX_SLOTS {
            return None;
        }
        let ring_slots: Vec<_> = slots
            .iter()
            .copied()
            .filter(|bond| self.ring_bonds[bond.0 as usize])
            .collect();
        let polarity_plans = 3usize.checked_pow(slots.len() as u32)?;
        let ring_plans = 1usize.checked_shl(ring_slots.len() as u32)?;
        if polarity_plans.checked_mul(ring_plans)? > MAX_CANDIDATES {
            return None;
        }

        let expected = ez_semantic_signature(self.mol);
        if expected.is_empty() {
            return None;
        }
        // Every candidate is written first; they are then checked in
        // ascending order, so the first that passes is the minimum the
        // exhaustive check would select, and the costly re-parse runs only
        // until it is found (row 4903 of the exposed 10k: 2,916 candidates).
        // A double bond keeps its E/Z on reparse only when each end writes a
        // directional token on one of its substituent bonds; an end whose
        // substituent bonds are all slots and all left plain cannot, so such
        // plans are not written at all.
        let slot_index: HashMap<BondIdx, usize> =
            slots.iter().enumerate().map(|(i, &b)| (b, i)).collect();
        let marker_sets: Vec<Vec<usize>> = self
            .extract_ez_geometry_facts()
            .iter()
            .flat_map(|fact| fact.ends.iter().map(|end| end.atom))
            .filter_map(|end| {
                Self::substituents(self.mol, end)
                    .iter()
                    .map(|&(_, bond)| slot_index.get(&bond).copied())
                    .collect::<Option<Vec<usize>>>()
            })
            .collect();
        let mut outputs: Vec<String> = Vec::new();
        for polarity_code in 0..polarity_plans {
            let digit = |slot: usize| (polarity_code / 3usize.pow(slot as u32)) % 3;
            if marker_sets
                .iter()
                .any(|set| set.iter().all(|&slot| digit(slot) == 0))
            {
                continue;
            }
            let mut code = polarity_code;
            let orders: HashMap<BondIdx, BondOrder> = slots
                .iter()
                .map(|&bond| {
                    let choice = code % 3;
                    code /= 3;
                    let order = match choice {
                        0 => Self::plain_order(self.mol.bond(bond).order),
                        1 => BondOrder::Up,
                        2 => BondOrder::Down,
                        _ => unreachable!(),
                    };
                    (bond, order)
                })
                .collect();
            for close_mask in 0..ring_plans {
                let mut candidate = self.clone();
                candidate.direct_ez_orders = Some(orders.clone());
                candidate.ring_marker_on_close = ring_slots
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &bond)| ((close_mask & (1 << index)) != 0).then_some(bond))
                    .collect();
                outputs.push(candidate.serialize_prepared());
            }
        }
        outputs.sort_unstable();
        outputs.dedup();
        outputs.into_iter().find(|output| {
            parse(output).is_ok_and(|reparsed| {
                ez_direction_assignments_consistent(&reparsed)
                    && ez_semantic_signature(&reparsed) == expected
            })
        })
    }

    /// Read one rank-fixed reference substituent's side.  A monosubstituted
    /// end has no carrier freedom; a disubstituted end uses the same
    /// rank-fixed/sibling-complement rule as [`Self::reference_up`].
    fn extract_ez_geometry_end(&self, alkene_end: AtomIdx) -> Option<EzGeometryEnd> {
        let subs = Self::substituents(self.mol, alkene_end);
        let (reference, up) = match subs.as_slice() {
            [reference] => {
                let direction = self.raw_input_direction(reference.1)?;
                (
                    *reference,
                    Self::direction_is_up(
                        direction,
                        self.raw_direction_anchor(reference.1),
                        alkene_end,
                    ),
                )
            }
            [first, second] => {
                let pair = [*first, *second];
                let preferred = usize::from(
                    self.ranks[pair[1].0.0 as usize] < self.ranks[pair[0].0.0 as usize],
                );
                (
                    pair[preferred],
                    self.reference_up(alkene_end, &pair, preferred)?,
                )
            }
            _ => return None,
        };
        Some(EzGeometryEnd {
            atom: alkene_end,
            reference_bond: reference.1,
            reference_rank: self.ranks[reference.0.0 as usize],
            up,
        })
    }

    /// Return a component end's already-extracted, rank-fixed geometry.  The
    /// carrier solver may choose a different output bond, but it must not
    /// derive its chemical side from the raw token on that candidate.
    fn geometry_reference_up(
        &self,
        alkene_end: AtomIdx,
        subs: &[(AtomIdx, BondIdx); 2],
        pref_idx: usize,
    ) -> Option<bool> {
        let reference_bond = subs[pref_idx].1;
        self.ez_geometry_facts
            .iter()
            .find_map(|fact| {
                debug_assert!(self.mol.bond(fact.double_bond).order == BondOrder::Double);
                debug_assert!(fact.double_key.0 <= fact.double_key.1);
                debug_assert_eq!(fact.same_side, fact.ends[0].up == fact.ends[1].up);
                fact.ends
                    .into_iter()
                    .find(|end| end.atom == alkene_end && end.reference_bond == reference_bond)
                    .map(|end| {
                        debug_assert_eq!(
                            end.reference_rank,
                            self.ranks[subs[pref_idx].0.0 as usize]
                        );
                        end.up
                    })
            })
            .or_else(|| self.reference_up(alkene_end, subs, pref_idx))
    }

    /// One candidate global assignment for a coupled component:
    /// `choice[i]` selects which of `subs[i]`'s two candidate bonds end
    /// `ordered[i]` uses as its own marker carrier. Returns the combined
    /// bond-value plan if every physical bond this component touches
    /// receives a single, consistent literal value from every end that
    /// votes on it (via [`Self::end_votes`]); `None` if any bond gets two
    /// disagreeing votes — the corruption [`Self::resolve_component_
    /// jointly`] exists to prevent.
    fn evaluate_choice(
        &self,
        ordered: &[AtomIdx],
        subs: &[[(AtomIdx, BondIdx); 2]],
        pref: &[usize],
        ref_up: &[Option<bool>],
        choice: &[usize],
    ) -> Option<HashMap<BondIdx, BondOrder>> {
        let mut votes: HashMap<BondIdx, BondOrder> = HashMap::new();
        for (i, &end) in ordered.iter().enumerate() {
            // A SMILES directional marker on a ring-closure bond is emitted
            // only at the DFS open side.  The close side is deliberately
            // suppressed because a ring digit cannot carry a second,
            // conflicting direction token.  Never select that close-side
            // occurrence as an E/Z carrier: the marker would be written at
            // the remote endpoint (which is not adjacent to this alkene),
            // and a re-parser would silently lose this stereochemistry.
            let chosen_bidx = subs[i][choice[i]].1;
            if !Self::can_carry_marker(self.mol.bond(chosen_bidx).order) {
                return None;
            }
            // The one exception: when the other candidate is an explicit H
            // atom, the marker goes on the ring bond and is written at its
            // opening (both readers take the direction from there), since
            // RDKit's `RemoveHs` misreads a marker on an H atom next to a
            // conjugated double bond (#734: an explicit-H macrolide).
            let sibling_is_h =
                self.mol.atom(subs[i][1 - choice[i]].0).element == chematic_core::Element::H;
            if self.ring_closure_is_close_side(chosen_bidx, end)
                && !self.ring_marker_is_permitted_on_close(chosen_bidx)
                && !sibling_is_h
            {
                return None;
            }
            // A carrier election is not geometry-neutral when the losing
            // (demoted) candidate is itself load-bearing for a *different*
            // stereo double bond -- demoting it would silently strip that
            // other double bond's own explicit geometry, and electing the
            // winner would (via `build_ez_groups`) hand that other double
            // bond's side a spurious value it never had (issue #390's actual
            // root cause). See `is_load_bearing_elsewhere`.
            let other_bidx = subs[i][1 - choice[i]].1;
            if self.is_load_bearing_elsewhere(other_bidx, end) {
                return None;
            }
            for (bidx, value) in self.end_votes(end, &subs[i], pref[i], ref_up[i], choice[i]) {
                match votes.get(&bidx) {
                    None => {
                        votes.insert(bidx, value);
                    }
                    Some(&existing) if existing == value => {}
                    Some(_) => return None,
                }
            }
        }
        Some(votes)
    }

}

fn deviation_count(bits: &[u8]) -> u32 {
    bits.iter().map(|&bit| u32::from(bit)).sum()
}
