impl<'a> CanonicalWriter<'a> {
    /// Whether `bidx` is written at `end` as the suppressed side of a ring
    /// closure. `find_ring_closures` populates both endpoint records before
    /// E/Z carrier resolution, so this is a topology-only fact and remains
    /// independent of the input atom ordering.
    fn ring_closure_is_close_side(&self, bidx: BondIdx, end: AtomIdx) -> bool {
        self.atom_ring_nums[end.0 as usize]
            .iter()
            .any(|&(_, is_open, _, ring_bidx)| ring_bidx == bidx && !is_open)
    }

    /// A complete E/Z plan may emit a directional ring marker at either legal
    /// ring-token occurrence. The set is empty on the ordinary path, which
    /// retains the opening-side default.
    fn ring_marker_is_permitted_on_close(&self, bidx: BondIdx) -> bool {
        self.ring_marker_on_close.contains(&bidx)
    }

    /// True if `bidx` (one of `owning_end`'s own two candidate substituent
    /// bonds) carries a genuine raw input mark *and* also flanks a
    /// *different* stereo double bond's end that is itself *not* ambiguous
    /// (exactly one substituent, so `bidx` is its ONLY possible source of
    /// geometric information) -- i.e. electing `owning_end`'s *other*
    /// candidate as carrier would demote `bidx` to a plain bond, silently
    /// stripping that other, unrelated double bond of the one explicit mark
    /// it has no alternative way to recover.
    ///
    /// Deliberately excludes the case where the far end is itself ambiguous
    /// (2 substituents, i.e. also in [`Self::compute_stereo_alkene_ends`]):
    /// such an end has its own resolution machinery (possibly jointly, via
    /// [`Self::coupling_components`]/[`Self::resolve_component_jointly`])
    /// and is not solely dependent on this one candidate's raw mark --
    /// forbidding its demotion here would block genuinely coupled,
    /// resolvable systems (`EZ_SHARED_CARRIER_FULLY_RESOLVED`) for no
    /// reason, since a real conflict there is already caught by
    /// `evaluate_choice`'s own vote-consistency check.
    ///
    /// This is deliberately narrower than "does `bidx` have a raw mark":
    /// demoting a raw-marked candidate that flanks nothing else (the common,
    /// intended case `resolve_ez_markers` exists for -- picking whichever of
    /// two candidates the canonical rank prefers, regardless of which one
    /// the input happened to mark) must stay allowed. Only a candidate that
    /// is *also* someone else's sole, non-ambiguous load-bearing mark is
    /// protected.
    fn is_load_bearing_elsewhere(&self, bidx: BondIdx, owning_end: AtomIdx) -> bool {
        if self.raw_input_direction(bidx).is_none() {
            return false;
        }
        let bond = self.mol.bond(bidx);
        let Some(other_end) = bond.other(owning_end) else {
            return false;
        };
        if Self::substituents(self.mol, other_end).len() != 1 {
            return false; // other_end is itself ambiguous -- has its own resolution path
        }
        // Only a double bond that is specified -- its far end carries a
        // marker too, and it is not held cis by a small ring -- depends on
        // this marker; a lone marker means nothing and is dropped on output.
        self.mol.neighbors(other_end).any(|(_, nb_bidx)| {
            nb_bidx != bidx
                && self.mol.bond(nb_bidx).order == BondOrder::Double
                && !crate::writer::in_ring_smaller_than_eight(self.mol, nb_bidx)
                && Self::end_has_substituent(self.mol, other_end)
                && self.mol.bond(nb_bidx).other(other_end).is_some_and(|far| {
                    Self::end_has_substituent(self.mol, far)
                        && self
                            .mol
                            .neighbors(far)
                            .any(|(_, sb)| sb != nb_bidx && self.raw_input_direction(sb).is_some())
                })
        })
    }

    /// One end's own votes on its two candidate bonds, given its
    /// precomputed [`Self::reference_up`] fact — generalizes the original
    /// single-end carrier/sibling branching to an arbitrary (not
    /// necessarily rank-preferred) `chosen` index, so a lone (uncoupled)
    /// end and a jointly-resolved coupled end share exactly the same
    /// per-end semantics. Always votes on *both* of the end's candidate
    /// bonds (the non-chosen one demoted to its plain order) — not
    /// conditionally, so the vote set depends only on `ref_up` (itself
    /// input-mark-placement-invariant), never on whether the non-chosen
    /// bond happened to carry a raw mark in this particular spelling.
    fn end_votes(
        &self,
        alkene_end: AtomIdx,
        subs: &[(AtomIdx, BondIdx); 2],
        pref_idx: usize,
        ref_up: Option<bool>,
        chosen_idx: usize,
    ) -> Vec<(BondIdx, BondOrder)> {
        let Some(ref_up) = ref_up else {
            return Vec::new(); // no direction info at this end at all
        };
        let chosen = subs[chosen_idx];
        let other = subs[1 - chosen_idx];
        let chosen_up = if chosen_idx == pref_idx {
            ref_up
        } else {
            !ref_up
        };
        // The input-side anchor is used only to recover `ref_up`. Once a
        // carrier has been selected, emit its direction in the stable
        // molecule-relative atom1→atom2 frame. Reusing a parser stash anchor
        // here would make the result depend on whether this physical bond
        // happened to be the original carrier in the input spelling.
        let picked = Self::direction_for_up(self.mol.bond(chosen.1).atom1, alkene_end, chosen_up);
        vec![
            (chosen.1, picked),
            (other.1, Self::plain_order(self.mol.bond(other.1).order)),
        ]
    }

    /// Union all directional single bonds flanking each stereo double bond
    /// into one group per connected E/Z system (order-independent — depends
    /// only on molecule topology, not on canonical ranks or write order).
    fn build_ez_groups(&mut self) {
        fn find(group: &mut HashMap<BondIdx, BondIdx>, x: BondIdx) -> BondIdx {
            let parent = *group.get(&x).unwrap_or(&x);
            if parent == x {
                x
            } else {
                let root = find(group, parent);
                group.insert(x, root);
                root
            }
        }

        fn union(group: &mut HashMap<BondIdx, BondIdx>, a: BondIdx, b: BondIdx) {
            let ra = find(group, a);
            let rb = find(group, b);
            if ra != rb {
                group.insert(ra, rb);
            }
        }

        for bidx in 0..self.mol.bond_count() {
            let bidx = BondIdx(bidx as u32);
            let bond = self.mol.bond(bidx);
            if bond.order != BondOrder::Double {
                continue;
            }
            let mut side_bonds = Vec::new();
            for endpoint in [bond.atom1, bond.atom2] {
                for (_, nb_bidx) in self.mol.neighbors(endpoint) {
                    if nb_bidx == bidx {
                        continue;
                    }
                    // `effective_order` (not the raw bond order/stash) so a
                    // marker `resolve_ez_markers` moved onto a different
                    // substituent is grouped by its NEW carrier bond, not
                    // its old (now-demoted) one.
                    if matches!(
                        self.effective_order(nb_bidx),
                        BondOrder::Up | BondOrder::Down
                    ) {
                        side_bonds.push(nb_bidx);
                    }
                }
            }
            let Some(&first) = side_bonds.first() else {
                continue;
            };
            // A double bond held cis by a small ring has no E/Z: its
            // flanking markers belong to the double bonds next to it, whose
            // spellings flip independently (two exocyclic imines on a
            // squaraine ring were written in either relative orientation).
            if side_bonds.len() > 1
                && *self
                    .small_ring_double
                    .borrow_mut()
                    .entry(bidx)
                    .or_insert_with(|| crate::writer::in_ring_smaller_than_eight(self.mol, bidx))
            {
                continue;
            }
            self.ez_group.entry(first).or_insert(first);
            for &b in &side_bonds[1..] {
                self.ez_group.entry(b).or_insert(b);
                union(&mut self.ez_group, first, b);
            }
        }
    }

    /// Normalize a directional bond so every member of its E/Z system is
    /// flipped consistently -- in mol-relative (`atom1`->`atom2`) terms -- to
    /// preserve geometry, with the group's shared sign pinned so the first
    /// occurrence in canonical write order prints as `Up` (`/`).
    ///
    /// This function does two genuinely different jobs, on two different
    /// frames, and conflating them is exactly how issue #390 happened:
    ///
    /// 1. **Propagation** (every call): flip [`Self::effective_order`] --
    ///    mol-relative, topology-fixed -- by the group's shared bit.
    ///    `build_ez_groups` unions bonds by molecule topology alone, so this
    ///    MUST run in that same topology-fixed frame. A bond's `atom1`/
    ///    `atom2` struct fields are fixed at parse/build time, but which
    ///    endpoint a *specific* DFS write visits first varies per occurrence
    ///    and per candidate canonical numbering; flipping an
    ///    already-write-oriented value mixes bonds visited "forward" with
    ///    bonds visited "backward" relative to their own storage, silently
    ///    inverting the very same/different relationship this function
    ///    exists to preserve whenever a group spans bonds with different
    ///    write directions (confirmed against RDKit: a coupled system with
    ///    one member's DFS direction reversed relative to another's produced
    ///    a genuinely wrong, not just differently spelled, configuration).
    /// 2. **Anchoring** (first call per group only): seed the group's shared
    ///    bit from whether *this specific occurrence*, re-oriented for
    ///    `from_atom` via [`Self::reorient_for_write`], would print as
    ///    `Down` -- because the point of the bit is to pin a *printed*
    ///    character, which only exists in write-perspective terms. Seeding
    ///    it from the mol-relative value instead pins an artifact of
    ///    whichever input text happened to be parsed (`atom1`/`atom2`
    ///    ordering is parse-order, not a canonical property): a different
    ///    spelling of the same molecule can swap `atom1`/`atom2`, seed a
    ///    different bit, and sign-flip the entire group -- geometry still
    ///    internally consistent (propagation was already fixed), but no
    ///    longer canonical or idempotent. `from_atom` is used **only** here,
    ///    never mixed into propagation.
    ///
    /// Every call site MUST separately apply [`Self::reorient_for_write`] to
    /// this function's *return value* for its own occurrence -- passing
    /// `from_atom` here does not do that; it only seeds the anchor.
    fn normalize_ez(&mut self, bidx: BondIdx, from_atom: AtomIdx) -> BondOrder {
        if self.direct_ez_orders.is_some() {
            let _ = from_atom;
            return self.effective_order(bidx);
        }
        let order = self.effective_order(bidx);
        if !matches!(order, BondOrder::Up | BondOrder::Down) {
            return order;
        }
        let root = {
            let mut x = *self.ez_group.get(&bidx).unwrap_or(&bidx);
            while let Some(&p) = self.ez_group.get(&x) {
                if p == x {
                    break;
                }
                x = p;
            }
            x
        };
        let flip = if let Some(&f) = self.ez_flip.get(&root) {
            f
        } else if let Some(forced) = &self.forced_ez_flip {
            let f = forced.get(&root).copied().unwrap_or(false);
            self.ez_flip.insert(root, f);
            f
        } else {
            let direction_anchor = self.raw_direction_anchor(bidx);
            let printed = Self::reorient_for_write(direction_anchor, from_atom, order);
            let seed = printed == BondOrder::Down;
            self.ez_flip.insert(root, seed);
            seed
        };
        if flip {
            match order {
                BondOrder::Up => BondOrder::Down,
                BondOrder::Down => BondOrder::Up,
                other => other,
            }
        } else {
            order
        }
    }

    /// Re-orient a directional order relative to its stored direction anchor
    /// for reading "from `from_atom` toward the bond's other endpoint" -- the
    /// per-occurrence step every [`Self::normalize_ez`] caller must apply to
    /// its return value. Mirrors `crate::writer::direction_from` exactly
    /// (kept as a second copy since this one only ever runs after
    /// `normalize_ez`, which has no equivalent in the plain writer).
    fn reorient_for_write(
        direction_anchor: AtomIdx,
        from_atom: AtomIdx,
        order: BondOrder,
    ) -> BondOrder {
        match order {
            BondOrder::Up => {
                if direction_anchor == from_atom {
                    BondOrder::Up
                } else {
                    BondOrder::Down
                }
            }
            BondOrder::Down => {
                if direction_anchor == from_atom {
                    BondOrder::Down
                } else {
                    BondOrder::Up
                }
            }
            other => other,
        }
    }

    /// Whether a directional bond can participate in an E/Z marker choice.
    ///
    /// Wedge/dash bonds are also represented as `Up`/`Down`, but they do not
    /// need the alkene marker solver.  Keep this check allocation-free so the
    /// common tetrahedral-stereo-only path avoids the SSSR and coupling graph
    /// work in [`Self::resolve_ez_markers`].
    fn has_directional_alkene_candidate(&self) -> bool {
        (0..self.mol.bond_count()).any(|i| {
            let bidx = BondIdx(i as u32);
            let bond = self.mol.bond(bidx);
            if bond.order == BondOrder::Double
                || !matches!(
                    crate::writer::raw_bond_direction(self.mol, bidx),
                    Some(BondOrder::Up | BondOrder::Down)
                )
            {
                return false;
            }
            [bond.atom1, bond.atom2].into_iter().any(|end| {
                self.mol
                    .neighbors(end)
                    .any(|(_, adjacent)| self.mol.bond(adjacent).order == BondOrder::Double)
            })
        })
    }

    pub(crate) fn write_all(self) -> String {
        self.write_all_impl(false).0
    }

    /// [`Self::write_all`], also returning the atoms in the order they appear
    /// in the string. The DFS visit order depends only on the ranks, the
    /// graph and the atom properties `canonical_cmp` reads — never on which
    /// E/Z carrier spelling is chosen — so every internal candidate spelling
    /// shares it, and it is read from one serialization of the prepared state.
    pub(crate) fn write_all_with_order(self) -> (String, Vec<AtomIdx>) {
        self.write_all_impl(true)
    }

    fn write_all_impl(mut self, want_order: bool) -> (String, Vec<AtomIdx>) {
        // Ring-closure side selection is independent of E/Z carrier choice,
        // but carrier choice must know which endpoint is the suppressed close
        // side. Discover the canonical DFS back-edges first so a directional
        // marker is never assigned to an occurrence that cannot be emitted.
        let starts = self.canonical_atom_list();
        self.find_ring_closures(&starts);

        // Phase -1/0 only matter when the molecule actually carries a
        // directional bond (literal Up/Down or an aromatic-direction stash).
        // Avoid scanning for stereo alkene ends and building union-find
        // groups for the common non-stereo case. `normalize_ez` remains safe
        // for ordinary bonds with an empty group map.
        let has_direct_plan = self.direct_ez_orders.is_some();

        if has_direct_plan {
            self.ez_geometry_facts = self.extract_ez_geometry_facts();
        } else if self.has_directional_alkene_candidate() {
            // Pick, for every stereo alkene end, which substituent carries the
            // marker, then group connected E/Z systems.
            self.resolve_ez_markers();
            self.build_ez_groups();
        }

        // Aromatic direction stashes can encode a valid marker on a carrier
        // outside the local coupled component. Search the complete bounded
        // slot universe first; a semantic reparse gate makes the selected
        // lexicographic minimum independent of input marker placement.
        let order_probe = want_order.then(|| self.clone());
        let order_of = |probe: Option<Self>| {
            probe.map_or_else(Vec::new, |p| p.serialize_prepared_with_order().1)
        };

        #[cfg(test)]
        let is_traversal_probe = self.traversal_bond_preference.is_some();
        #[cfg(not(test))]
        let is_traversal_probe = false;
        if !has_direct_plan
            && !is_traversal_probe
            && self.has_aromatic_direction_stash()
            && let Some(output) = self.complete_aromatic_ez_plan()
        {
            return (output, order_of(order_probe));
        }

        // Keep the established component-global polarity normalization as a
        // bounded fallback when the complete planner is inapplicable.
        if !has_direct_plan && self.has_aromatic_direction_stash() && !self.ez_group.is_empty() {
            let roots = self.ez_group.values().copied().collect::<HashSet<_>>();
            if roots.len() <= 8 {
                let roots = roots.into_iter().collect::<Vec<_>>();
                let mut best = None;
                for mask in 0..(1usize << roots.len()) {
                    let mut candidate = self.clone();
                    candidate.forced_ez_flip = Some(
                        roots
                            .iter()
                            .enumerate()
                            .map(|(index, &root)| (root, (mask & (1 << index)) != 0))
                            .collect(),
                    );
                    let output = candidate.serialize_prepared();
                    if best.as_ref().is_none_or(|current| output < *current) {
                        best = Some(output);
                    }
                }
                return (best.unwrap_or_default(), order_of(order_probe));
            }
        }

        let (out, order) = self.serialize_prepared_with_order();
        (out, if want_order { order } else { Vec::new() })
    }

    fn has_aromatic_direction_stash(&self) -> bool {
        self.mol.bonds().any(|(bidx, bond)| {
            bond.order == BondOrder::Aromatic && self.mol.bond_direction(bidx).is_some()
        })
    }

    fn serialize_prepared(self) -> String {
        self.serialize_prepared_with_order().0
    }

    fn serialize_prepared_with_order(mut self) -> (String, Vec<AtomIdx>) {
        // Phase 2: canonical DFS serialization.
        let mut first = true;
        for start in self.canonical_atom_list() {
            if self.written[start.0 as usize] {
                continue;
            }
            if !first {
                self.out.push('.');
            }
            first = false;
            self.write_chain(start, None, None);
        }

        (self.out, self.order)
    }

    /// Return all atoms sorted in canonical order: highest rank first, ties
    /// broken by chemical properties invariant across re-parses.
    fn canonical_atom_list(&self) -> Vec<AtomIdx> {
        let n = self.mol.atom_count();
        let sentinel = AtomIdx(u32::MAX);
        let mut atoms = vec![sentinel; n];
        if self.ranks.len() == n {
            let mut discrete = true;
            for (i, &rank) in self.ranks.iter().enumerate() {
                let rank = rank as usize;
                if rank >= n || atoms[n - 1 - rank] != sentinel {
                    discrete = false;
                    break;
                }
                atoms[n - 1 - rank] = AtomIdx(i as u32);
            }
            if discrete {
                // `n` unique ranks in 0..n necessarily fill every slot.
                return atoms;
            }
        }

        // Defensive fallback for tests and the legacy capped search, which
        // may intentionally hand the writer a rank vector containing ties.
        atoms.clear();
        atoms.extend((0..n).map(|i| AtomIdx(i as u32)));
        atoms.sort_by(|&a, &b| self.canonical_cmp(b, a)); // descending
        atoms
    }

    /// Canonical ordering comparator (ascending; negate for descending).
    /// Tie-breaking uses chemical properties only (not atom indices),
    /// so the order is invariant between runs on chemically identical molecules.
    fn canonical_cmp(&self, a: AtomIdx, b: AtomIdx) -> std::cmp::Ordering {
        let ra = self.ranks[a.0 as usize];
        let rb = self.ranks[b.0 as usize];
        if ra != rb {
            return ra.cmp(&rb);
        }

        let atom_a = self.mol.atom(a);
        let atom_b = self.mol.atom(b);

        // Break ties with: atomic_number → isotope → charge → aromatic → degree
        atom_a
            .element
            .atomic_number()
            .cmp(&atom_b.element.atomic_number())
            .then_with(|| {
                atom_a
                    .isotope
                    .unwrap_or(0)
                    .cmp(&atom_b.isotope.unwrap_or(0))
            })
            .then_with(|| atom_a.charge.cmp(&atom_b.charge))
            .then_with(|| (atom_a.aromatic as u8).cmp(&(atom_b.aromatic as u8)))
            .then_with(|| self.mol.degree(a).cmp(&self.mol.degree(b)))
    }

    /// Discover back-edges by running the same canonical DFS as the writer.
    /// Using identical traversal order ensures ring-closure numbers are stable.
    fn find_ring_closures(&mut self, starts: &[AtomIdx]) {
        let n = self.mol.atom_count();
        let mut visited = vec![false; n];
        let mut in_stack = vec![false; n];

        // Iterate in canonical order (the same allocation write_all reuses
        // for the serialization phase).
        for &start in starts {
            if !visited[start.0 as usize] {
                self.dfs_mark(start, None, &mut visited, &mut in_stack);
            }
        }
    }

    fn dfs_mark(
        &mut self,
        atom: AtomIdx,
        from_bond: Option<BondIdx>,
        visited: &mut Vec<bool>,
        in_stack: &mut Vec<bool>,
    ) {
        visited[atom.0 as usize] = true;
        in_stack[atom.0 as usize] = true;

        let mut neighbors: SmallVec<[(AtomIdx, BondIdx); 4]> = self.mol.neighbors(atom).collect();
        self.sort_neighbors_canonical(&mut neighbors);

        for (neighbor, bidx) in neighbors {
            if Some(bidx) == from_bond {
                continue;
            }
            if self.ring_bonds[bidx.0 as usize] {
                continue;
            }

            if !visited[neighbor.0 as usize] {
                if !self.canonical_tree_parent.is_empty() {
                    self.canonical_tree_parent[neighbor.0 as usize] = Some((atom, bidx));
                }
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                self.ring_bonds[bidx.0 as usize] = true;
                let rn = self.next_ring;
                self.next_ring += 1;
                // Discovery only decides ring numbering and which endpoint is
                // "open" (carries the marker, written when `neighbor` is
                // reached) vs "close" (suppressed, written when `atom` is
                // reached) -- NOT the actual `/`/`\` token. That is computed
                // later, in `write_chain`, from `normalize_ez`'s mol-relative
                // result re-oriented for whichever atom is actually being
                // written at consumption time (see `normalize_ez`'s doc
                // comment and `atom_ring_nums`'s field comment).
                self.atom_ring_nums[neighbor.0 as usize].push((rn, true, atom, bidx)); // open side; partner = close atom
                self.atom_ring_nums[atom.0 as usize].push((rn, false, neighbor, bidx)); // close side; partner = open atom
            }
        }

        in_stack[atom.0 as usize] = false;
    }

    /// `incoming_bond` is the already-oriented token for the edge leading to
    /// `atom` (`None` for the root or an implicit bond), not a `BondOrder`:
    /// a dative arrow's direction depends on which endpoint is written
    /// first, and `BondOrder::smiles_char` would truncate `"->"` to `'-'`
    /// anyway (issue #194). See `crate::writer::bond_token_from`.
    fn write_chain(
        &mut self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        incoming_bond: Option<&'static str>,
    ) {
        self.written[atom.0 as usize] = true;

        if let Some(token) = incoming_bond {
            self.out.push_str(token);
        }

        // Compute parity-corrected chirality before ring data is consumed.
        let corrected_chirality = self.corrected_chirality(atom, from_atom);
        self.order.push(atom);
        self.emit_atom(atom, corrected_chirality);

        self.write_ring_closures(atom);
        self.write_tree_children(atom, from_atom);
    }

    fn write_ring_closures(&mut self, atom: AtomIdx) {
        let rings = std::mem::take(&mut self.atom_ring_nums[atom.0 as usize]);
        for (number, is_open, partner, bond) in rings {
            let close_marker = self.ring_marker_is_permitted_on_close(bond);
            let emit_direction = is_open != close_marker;
            let order = if emit_direction {
                let normalized = self.normalize_ez(bond, atom);
                Self::reorient_for_write(self.raw_direction_anchor(bond), atom, normalized)
            } else {
                match self.effective_order(bond) {
                    BondOrder::Up | BondOrder::Down => {
                        Self::plain_order(self.mol.bond(bond).order)
                    }
                    other => other,
                }
            };
            let order = self.suppress_lone_marker(bond, order);
            if !self.bond_is_implicit(atom, partner, order) {
                self.out
                    .push_str(bond_token_from(self.mol, bond, order, atom));
            }
            self.write_ring_number(number);
        }
    }

    fn write_ring_number(&mut self, number: u32) {
        if number >= 100 {
            self.out.push_str("%(");
            self.out.push_str(&number.to_string());
            self.out.push(')');
        } else if number >= 10 {
            self.out.push('%');
            self.out.push(char::from_digit(number / 10, 10).unwrap());
            self.out.push(char::from_digit(number % 10, 10).unwrap());
        } else {
            self.out.push(char::from_digit(number, 10).unwrap());
        }
    }

    fn bond_is_implicit(&self, left: AtomIdx, right: AtomIdx, order: BondOrder) -> bool {
        let both_aromatic = self.mol.atom(left).aromatic && self.mol.atom(right).aromatic;
        match order {
            BondOrder::Single => !both_aromatic,
            BondOrder::Aromatic => both_aromatic,
            _ => false,
        }
    }

    fn write_tree_children(&mut self, atom: AtomIdx, from_atom: Option<AtomIdx>) {
        let mut children: SmallVec<[(AtomIdx, BondIdx); 4]> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds[bidx.0 as usize]
            })
            .collect();

        children.sort_by(|&(a, ..), &(b, ..)| self.canonical_cmp(a, b));
        let count = children.len();
        for (i, (child, bidx)) in children.into_iter().enumerate() {
            let normalized = self.normalize_ez(bidx, atom);
            let bond_order =
                Self::reorient_for_write(self.raw_direction_anchor(bidx), atom, normalized);
            let bond_order = self.suppress_lone_marker(bidx, bond_order);
            let written_bond = if self.bond_is_implicit(atom, child, bond_order) {
                None
            } else {
                Some(bond_token_from(self.mol, bidx, bond_order, atom))
            };
            if i + 1 < count {
                self.out.push('(');
                self.write_chain(child, Some(atom), written_bond);
                self.out.push(')');
            } else {
                self.write_chain(child, Some(atom), written_bond);
            }
        }
    }

    /// Sort a neighbor list in canonical order (for consistent DFS traversal).
    fn sort_neighbors_canonical(&self, neighbors: &mut [(AtomIdx, BondIdx)]) {
        neighbors.sort_by(|&(a, _a_bond), &(b, _b_bond)| {
            #[cfg(test)]
            if let Some(preferences) = &self.traversal_bond_preference {
                let a_preferred = preferences.get(&_a_bond).copied().unwrap_or(false);
                let b_preferred = preferences.get(&_b_bond).copied().unwrap_or(false);
                if a_preferred != b_preferred {
                    // A preferred edge is traversed first.  This comparator
                    // is used identically by ring discovery and emission.
                    return b_preferred.cmp(&a_preferred);
                }
            }
            self.canonical_cmp(b, a) // descending
        });
    }

    fn emit_atom(&mut self, idx: AtomIdx, chirality: Chirality) {
        let atom = self.mol.atom(idx);

        if atom.wildcard {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            self.out.push('*');
            if let Some(h) = atom.hydrogen_count.filter(|&h| h > 0) {
                self.out.push('H');
                if h > 1 {
                    self.out.push_str(&h.to_string());
                }
            }
            match atom.charge {
                0 => {}
                1 => self.out.push('+'),
                -1 => self.out.push('-'),
                c if c > 0 => self.out.push_str(&format!("+{c}")),
                c => self.out.push_str(&c.to_string()),
            }
            if let Some(map) = atom.atom_map {
                self.out.push(':');
                self.out.push_str(&map.to_string());
            }
            self.out.push(']');
            return;
        }

        // An explicit `hydrogen_count` only forces bracket notation here if it
        // carries genuine disambiguating information -- i.e. differs from
        // what organic-subset valence inference would produce anyway if the
        // atom were unbracketed. An atom whose stored H count merely repeats
        // that inferred value (e.g. `[Cl]` vs organic-subset `Cl`, both H=0)
        // is the same chemical species and must canonicalize to the same
        // (unbracketed) spelling regardless of which notation it was parsed
        // from -- unlike `writer.rs`'s plain (non-canonical) `emit_atom`,
        // which intentionally preserves original notation and is unaffected
        // by this distinction.
        let h_is_disambiguating = atom
            .hydrogen_count
            .is_some_and(|h| h != valence_inferred_hcount(self.mol, idx));

        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || h_is_disambiguating
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some()
            || chirality != Chirality::None;

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

            match chirality {
                Chirality::CounterClockwise => self.out.push('@'),
                Chirality::Clockwise => self.out.push_str("@@"),
                Chirality::None => {}
                Chirality::SquarePlanar(p) => self.out.push_str(square_planar_token(p)),
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

    /// Compute the parity-corrected chirality for `atom` when it is written
    /// with `from_atom` as the predecessor in the canonical DFS.
    ///
    /// Tetrahedral: returns the stored chirality unchanged when no stereo
    /// neighbor order is recorded (e.g. programmatically constructed
    /// molecules) -- with only 2 possible states, "unchanged, no better
    /// information" is a safe no-op. Square-planar: drops to
    /// [`Chirality::None`] instead whenever the neighbor order can't be
    /// verified (no recorded order, size mismatch, duplicate/foreign ids) --
    /// for a 3-state tag, passing the *original* tag through against a
    /// *reordered* neighbor list can silently describe a different,
    /// plausible-but-wrong stereoisomer, exactly the failure category this
    /// mechanism exists to eliminate. See
    /// [`chematic_core::remap_square_planar_tag`] (the generalized
    /// stereo-geometry module, `docs/rfcs/generalized_stereo_geometry_rfc.md`).
    fn corrected_chirality(&self, atom: AtomIdx, from_atom: Option<AtomIdx>) -> Chirality {
        let stored = self.mol.atom(atom).chirality;
        if stored == Chirality::None {
            return Chirality::None;
        }

        let Some(original) = self.mol.stereo_neighbor_order(atom) else {
            return if stored.is_tetrahedral() {
                stored
            } else {
                Chirality::None
            };
        };

        if stored.is_allene() {
            let Some(canonical) = self.canonical_allene_virtual_order(atom) else {
                return Chirality::None;
            };
            let Some(inverted) = crate::allene::allene_order_is_inverted(original, &canonical)
            else {
                return Chirality::None;
            };
            return if inverted {
                match stored {
                    Chirality::Allene1 => Chirality::Allene2,
                    Chirality::Allene2 => Chirality::Allene1,
                    _ => unreachable!(),
                }
            } else {
                stored
            };
        }

        let canonical = self.canonical_stereo_neighbor_order(atom, from_atom, original.len());
        if canonical.len() != original.len() {
            return if stored.is_tetrahedral() {
                stored // size mismatch → tetrahedral fallback (unchanged)
            } else {
                Chirality::None // square-planar: unverifiable → drop, don't guess
            };
        }

        self.remap_stored_chirality(stored, original, &canonical)
    }

    fn canonical_allene_virtual_order(&self, center: AtomIdx) -> Option<Vec<u32>> {
        let center_order = self.canonical_local_neighbor_order(center);
        let arms = crate::allene::allene_arms(self.mol, center, &center_order)?;
        let mut virtual_order = Vec::with_capacity(6);
        for arm in arms {
            let mut local = self.canonical_local_neighbor_order(arm.terminal);
            local.retain(|&atom| atom != arm.axis_neighbor.0);
            if local.len() != 2 {
                return None;
            }
            virtual_order.push(arm.terminal.0);
            virtual_order.extend(local);
        }
        Some(virtual_order)
    }

    /// Neighbor order produced by the already-discovered canonical DFS,
    /// independent of which atoms happen to have been emitted so far.
    fn canonical_local_neighbor_order(&self, atom: AtomIdx) -> Vec<u32> {
        let parent = self.canonical_tree_parent[atom.0 as usize].map(|(parent, _)| parent);
        let mut order = Vec::with_capacity(self.mol.degree(atom) + 1);
        if let Some(parent) = parent {
            order.push(parent.0);
        }
        if implicit_hcount(self.mol, atom) > 0 {
            order.push(STEREO_H_SENTINEL);
        }
        for &(_, _, partner, _) in &self.atom_ring_nums[atom.0 as usize] {
            order.push(partner.0);
        }
        let mut children: SmallVec<[AtomIdx; 4]> = self
            .mol
            .neighbors(atom)
            .filter(|(neighbor, bond)| {
                Some(*neighbor) != parent
                    && !self.ring_bonds[bond.0 as usize]
                    && self.canonical_tree_parent[neighbor.0 as usize]
                        .is_some_and(|(candidate, _)| candidate == atom)
            })
            .map(|(neighbor, _)| neighbor)
            .collect();
        children.sort_by(|&left, &right| self.canonical_cmp(left, right));
        order.extend(children.into_iter().map(|child| child.0));
        order
    }

    fn canonical_stereo_neighbor_order(
        &self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        capacity: usize,
    ) -> SmallVec<[u32; 4]> {
        let has_h = implicit_hcount(self.mol, atom) > 0;
        let mut canonical = SmallVec::with_capacity(capacity);
        if let Some(previous) = from_atom {
            canonical.push(previous.0);
        }
        if has_h {
            canonical.push(STEREO_H_SENTINEL);
        }
        for &(_, _, partner, _) in &self.atom_ring_nums[atom.0 as usize] {
            canonical.push(partner.0);
        }
        let mut children: SmallVec<[AtomIdx; 4]> = self
            .mol
            .neighbors(atom)
            .filter(|(neighbor, bond)| {
                Some(*neighbor) != from_atom
                    && !self.written[neighbor.0 as usize]
                    && !self.ring_bonds[bond.0 as usize]
            })
            .map(|(neighbor, _)| neighbor)
            .collect();
        children.sort_by(|&left, &right| self.canonical_cmp(left, right));
        canonical.extend(children.into_iter().map(|child| child.0));
        canonical
    }

    fn remap_stored_chirality(
        &self,
        stored: Chirality,
        original: &[u32],
        canonical: &[u32],
    ) -> Chirality {
        let original_arr = <[u32; 4]>::try_from(original).ok();
        let canonical_arr = <[u32; 4]>::try_from(canonical).ok();
        match stored {
            Chirality::Allene1 | Chirality::Allene2 => Chirality::None,
            Chirality::CounterClockwise | Chirality::Clockwise => {
                let is_odd = match (original_arr, canonical_arr) {
                    (Some(o), Some(c)) => remap_tetrahedral_parity(o, c).unwrap_or(false),
                    _ => permutation_is_odd(original, canonical),
                };
                if is_odd {
                    if stored == Chirality::CounterClockwise {
                        Chirality::Clockwise
                    } else {
                        Chirality::CounterClockwise
                    }
                } else {
                    stored
                }
            }
            Chirality::SquarePlanar(tag) => match (original_arr, canonical_arr) {
                (Some(o), Some(c)) => remap_square_planar_tag(tag, o, c)
                    .map(Chirality::SquarePlanar)
                    .unwrap_or(Chirality::None),
                _ => Chirality::None,
            },
            Chirality::None => Chirality::None,
            other => crate::writer::remap_nontetrahedral(other, original, canonical),
        }
    }
}
