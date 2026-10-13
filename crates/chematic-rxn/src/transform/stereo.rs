use super::*;

/// Copy each reactant bond's stashed direction ([`Molecule::bond_direction`]
/// with its anchor) onto the product bond between the same two atoms, when
/// that bond is still single or aromatic and one of its atoms keeps a double
/// bond to a third atom (#734: `[Cl:1]>>[Cl:1]` dropped the E/Z of a
/// squaraine imine).
pub(super) fn carry_direction_stash(
    product: &mut Molecule,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) {
    if !input_mols.iter().any(|m| m.has_bond_directions()) {
        return;
    }
    let mut source_of: FxHashMap<AtomIdx, (usize, AtomIdx)> =
        FxHashMap::with_capacity_and_hasher(src_to_new.len(), Default::default());
    for (&key, &new_idx) in src_to_new {
        source_of.insert(new_idx, key);
    }
    let stashed: Vec<(BondIdx, BondOrder, AtomIdx)> = product
        .bonds()
        .filter_map(|(pb, bond)| {
            if !matches!(bond.order, BondOrder::Single | BondOrder::Aromatic)
                || product.bond_direction(pb).is_some()
            {
                return None;
            }
            let &(ma, sa) = source_of.get(&bond.atom1)?;
            let &(mb, sb) = source_of.get(&bond.atom2)?;
            if ma != mb {
                return None;
            }
            let src = input_mols[ma];
            let (sbid, sbond) = src.bond_between(sa, sb)?;
            let dir = src.bond_direction(sbid)?;
            let anchor = src.bond_direction_anchor(sbid).unwrap_or(sbond.atom1);
            let new_anchor = *src_to_new.get(&(ma, anchor))?;
            let flanks_double = [bond.atom1, bond.atom2].into_iter().any(|end| {
                product
                    .neighbors(end)
                    .any(|(_, b)| b != pb && product.bond(b).order == BondOrder::Double)
            });
            flanks_double.then_some((pb, dir, new_anchor))
        })
        .collect();
    for (pb, dir, anchor) in stashed {
        product.set_bond_direction(pb, dir);
        product.set_bond_direction_anchor(pb, anchor);
    }
}

/// Per atom, whether it lies on a ring (has an incident non-bridge bond).
pub(super) fn atoms_in_rings(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut order = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut in_ring = vec![false; n];
    let mut counter = 0usize;
    // (atom, bond used to reach it, index of its next neighbour)
    let mut stack: Vec<(usize, Option<BondIdx>, usize)> = Vec::with_capacity(n);
    for root in 0..n {
        if order[root] != usize::MAX {
            continue;
        }
        stack.clear();
        stack.push((root, None, 0));
        order[root] = counter;
        low[root] = counter;
        counter += 1;
        while let Some(top) = stack.last_mut() {
            let (v, via, k) = *top;
            if let Some((w, b)) = mol.neighbors(AtomIdx(v as u32)).nth(k) {
                top.2 += 1;
                if Some(b) == via {
                    continue;
                }
                let w = w.0 as usize;
                if order[w] == usize::MAX {
                    order[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push((w, Some(b), 0));
                } else {
                    low[v] = low[v].min(order[w]);
                    // A back edge closes a cycle through both ends.
                    in_ring[v] = true;
                    in_ring[w] = true;
                }
            } else {
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[v]);
                    if low[v] <= order[p] {
                        // The tree edge p–v is not a bridge.
                        in_ring[v] = true;
                        in_ring[p] = true;
                    }
                }
            }
        }
    }
    in_ring
}

// ---------------------------------------------------------------------------
// SMIRKS chirality post-check (parity-aware)
// ---------------------------------------------------------------------------

/// Returns the parity of the permutation P that maps `from_seq` to `to_seq`:
/// `Some(true)` = even (same chirality sense), `Some(false)` = odd (inverted).
/// Returns `None` when the sequences differ in length or contain elements that
/// cannot be aligned (e.g. an unmapped neighbour).
pub(super) fn permutation_parity(from_seq: &[u32], to_seq: &[u32]) -> Option<bool> {
    let n = from_seq.len();
    if n != to_seq.len() {
        return None;
    }
    // Build perm[j] = index i in from_seq where from_seq[i] == to_seq[j].
    let mut perm = Vec::with_capacity(n);
    for &t in to_seq {
        let pos = from_seq.iter().position(|&f| f == t)?;
        perm.push(pos);
    }
    // Count inversions to determine parity.
    let mut inv = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if perm[i] > perm[j] {
                inv += 1;
            }
        }
    }
    Some(inv.is_multiple_of(2)) // true = even = chirality flags must agree for same config
}

/// Parity-aware stereo check for one (template, reactant, mapping) triple.
///
/// For each chiral atom in `tmpl`, maps the template's recorded SMILES stereo
/// neighbor order through the VF2 `match_map` into the reactant atom-index space,
/// then computes the parity of the permutation relative to the reactant's recorded
/// SMILES stereo neighbor order.
///
/// - Even parity → chirality flags must agree for same absolute configuration.
/// - Odd parity  → chirality flags must differ for same absolute configuration.
///
/// Returns `true` if all chiral centres are consistent, `false` if any mismatch.
pub(super) fn smirks_chirality_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &FxHashMap<usize, AtomIdx>,
) -> bool {
    for i in 0..tmpl.atom_count() {
        let tmpl_atom = tmpl.atom(AtomIdx(i as u32));
        if tmpl_atom.chirality == Chirality::None {
            continue;
        }

        // Template atom's SMILES stereo neighbour order (template atom indices).
        let Some(tmpl_order) = tmpl.stereo_neighbor_order(AtomIdx(i as u32)) else {
            continue; // No recorded order — skip this centre.
        };

        // Corresponding matched reactant atom.
        let Some(&react_idx) = match_map.get(&i) else {
            continue; // Template atom not in mapping (shouldn't happen for complete match).
        };

        let react_atom = reactant.atom(react_idx);
        if react_atom.chirality == Chirality::None {
            return false; // Template requires stereo; reactant atom has none.
        }

        // Map each template stereo-neighbour index to the corresponding reactant atom index.
        let mut mapped: Vec<u32> = Vec::with_capacity(tmpl_order.len());
        let mut all_mapped = true;
        for &t in tmpl_order {
            if t == STEREO_H_SENTINEL {
                mapped.push(STEREO_H_SENTINEL);
            } else {
                match match_map.get(&(t as usize)) {
                    Some(ri) => mapped.push(ri.0),
                    None => {
                        all_mapped = false;
                        break;
                    }
                }
            }
        }
        if !all_mapped {
            continue; // Partial substructure match — cannot verify chirality.
        }

        // Reactant atom's SMILES stereo neighbour order (reactant atom indices).
        let Some(react_order) = reactant.stereo_neighbor_order(react_idx) else {
            // No recorded order in reactant — fall back to raw flag comparison.
            if react_atom.chirality != tmpl_atom.chirality {
                return false;
            }
            continue;
        };

        let Some(even_parity) = permutation_parity(&mapped, react_order) else {
            continue; // Alignment failed — skip.
        };

        // Check: even parity → flags must agree; odd parity → flags must differ.
        let same_flag = tmpl_atom.chirality == react_atom.chirality;
        if same_flag != even_parity {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// SMIRKS E/Z double-bond stereo post-check
// ---------------------------------------------------------------------------

/// Returns the "outward" stereo-bond direction from `atom` (one endpoint of a
/// double bond) toward its Up/Down-annotated substituent, or `None` if no
/// such bond exists.
///
/// E/Z stereo is encoded on the *substituent* bonds adjacent to a C=C, not on
/// the double bond itself.  Whether the stored bond goes *into* or *out of*
/// `atom` determines how to read its direction:
/// - Outgoing (`atom` is `bond.atom1`): direction is as stored.
/// - Incoming (`atom` is `bond.atom2`): direction is flipped (Up ↔ Down).
///
/// The `other` parameter is the opposite endpoint of the double bond; that
/// bond is skipped so we look only at substituents.
pub(super) fn ez_stereo_outward(
    mol: &Molecule,
    atom: AtomIdx,
    other: AtomIdx,
) -> Option<BondOrder> {
    for (nb, bidx) in mol.neighbors(atom) {
        if nb == other {
            continue; // skip the double bond itself
        }
        let bond = mol.bond(bidx);
        match bond.order {
            BondOrder::Up | BondOrder::Down => {
                let outward = if bond.atom1 == atom {
                    // bond goes FROM atom outward → direction as stored
                    bond.order
                } else {
                    // bond comes INTO atom → flip to get outward direction
                    match bond.order {
                        BondOrder::Up => BondOrder::Down,
                        _ => BondOrder::Up,
                    }
                };
                return Some(outward);
            }
            _ => {}
        }
    }
    None
}

/// E/Z stereo post-check for one (template, reactant, mapping) triple.
///
/// For each double bond in `tmpl` that has Up/Down substituent bonds on **both**
/// sides, verify that the corresponding double bond in `reactant` (found via
/// the VF2 `match_map`) encodes the same E/Z parity.
///
/// Parity is determined by comparing the "outward direction" from each end of
/// the double bond (see [`ez_stereo_outward`]).  Two outward directions that are
/// equal → same-side (Z/cis); directions that differ → opposite-side (E/trans).
///
/// If only one side of a template double bond has a stereo bond (or the
/// reactant doesn't annotate stereo at all), the constraint is skipped — this
/// matches the behaviour of SMIRKS templates extracted from rdchiral where a
/// single-sided annotation is common.
///
/// Returns `true` if all constrained double bonds are consistent.
pub(super) fn smirks_ez_stereo_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &FxHashMap<usize, AtomIdx>,
) -> bool {
    for (_, bond) in tmpl.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let ta = bond.atom1;
        let tb = bond.atom2;

        // Outward directions from each end of the template double bond.
        let sa = ez_stereo_outward(tmpl, ta, tb);
        let sb = ez_stereo_outward(tmpl, tb, ta);

        // Both sides must be specified to establish an E/Z constraint.
        let (sa, sb) = match (sa, sb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue, // no constraint on this double bond
        };

        // Map template atoms to reactant atoms.
        let Some(&ra) = match_map.get(&(ta.0 as usize)) else {
            continue;
        };
        let Some(&rb) = match_map.get(&(tb.0 as usize)) else {
            continue;
        };

        // Outward directions from the corresponding reactant double bond ends.
        let ma = ez_stereo_outward(reactant, ra, rb);
        let mb = ez_stereo_outward(reactant, rb, ra);

        // If the reactant doesn't encode stereo on either end, skip (don't reject).
        let (ma, mb) = match (ma, mb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        // Compare E/Z parity: same outward direction = Z, different = E.
        // If template and reactant disagree, reject this mapping.
        if (sa == sb) != (ma == mb) {
            return false;
        }
    }
    true
}

/// Convert a SMIRKS reactant-template `Molecule` to a `QueryMolecule` for VF2.
///
/// Constraints included:
/// - `AtomicNum` and `Aromatic` (always)
/// - `Charge` when non-zero
/// - `HCount` when a bracket atom specifies H > 0 (e.g. `[NH2:1]`)
///   Zero-H bracket atoms (`[N:1]`) are treated as "any H count" because
///   the parser returns 0 for both unspecified and explicit-zero H.
///
/// Chirality (`@`/`@@`) is NOT encoded into the query here; it is checked after
/// VF2 matching via `smirks_chirality_ok`, which uses a permutation-parity
/// comparison so that the same absolute configuration is recognised regardless
/// of how the reactant molecule was written as SMILES.
pub(super) fn mol_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let mut q = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(
                atom.element.atomic_number(),
            ))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(atom.aromatic))),
        );

        if atom.charge != 0 {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge))),
            );
        }

        if let Some(mass) = atom.isotope {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Isotope(mass))),
            );
        }

        if let Some(h) = atom.hydrogen_count
            && h > 0
        {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::HCount(h))),
            );
        }

        qmol.add_atom_with_map(q, atom.atom_map);
    }

    for (_bidx, bond) in mol.bonds() {
        let bq = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => {
                BondQuery::Primitive(BondPrimitive::Single)
            }
            BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
            BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
            BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
            BondOrder::QuerySingleOrDouble => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
            ),
            BondOrder::QuerySingleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::QueryDoubleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::Quadruple => BondQuery::Primitive(BondPrimitive::Quadruple),
            BondOrder::Zero | BondOrder::QueryAny => BondQuery::Primitive(BondPrimitive::Any),
        };
        qmol.add_bond(bond.atom1.0 as usize, bond.atom2.0 as usize, bq);
    }

    qmol
}

/// Clear Up/Down stereo markers from bonds that have no adjacent double bond.
///
/// After a SMIRKS reaction, C=C → C=O conversions can leave stale Up/Down
/// markers (E/Z direction indicators) on single bonds that are no longer
/// adjacent to any double bond.  Such orphaned markers produce invalid SMILES
/// (`/C=O` is nonsensical) and must be demoted to plain Single bonds (RDKit #9339).
pub(super) fn clear_orphaned_stereo_bonds(mol: Molecule) -> Molecule {
    let orphaned: Vec<BondIdx> = mol
        .bonds()
        .filter_map(|(bidx, bond)| {
            if bond.order != BondOrder::Up && bond.order != BondOrder::Down {
                return None;
            }
            // Up/Down is valid only when at least one endpoint has an adjacent
            // double bond (the one that the Up/Down bond specifies direction for).
            let has_double = [bond.atom1, bond.atom2].iter().any(|&a| {
                mol.neighbors(a)
                    .any(|(_, nb_bidx)| mol.bond(nb_bidx).order == BondOrder::Double)
            });
            if has_double { None } else { Some(bidx) }
        })
        .collect();

    if orphaned.is_empty() {
        return mol;
    }

    let mut builder = chematic_core::MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        builder.add_atom(atom.clone());
    }
    for (bidx, bond) in mol.bonds() {
        let order = if orphaned.contains(&bidx) {
            BondOrder::Single
        } else {
            bond.order
        };
        let _ = builder.add_bond(bond.atom1, bond.atom2, order);
    }
    // copy_stereo_from copies stereo_neighbor_order but NOT stereo_groups.
    // Preserve both by applying each separately.
    builder.copy_stereo_from(&mol);
    let mut result = builder.build();
    // Restore enhanced stereo groups (ABS/OR/AND) that copy_stereo_from omits.
    result.copy_stereo_groups_from(&mol);
    result
}

// ---------------------------------------------------------------------------
// Product-side chirality correction (parity-aware)
// ---------------------------------------------------------------------------
//
// `build_product`'s Step 1 leaves every mapped core atom's chirality as
// whatever `src_atom.clone()` produced -- the reactant's own flag, if any.
// Two things can make that flag wrong or meaningless by the time the
// product molecule's bonds are fully built (Steps 3-4):
//
//   - An EXPLICIT product-template flag (`[C@:1]`/`[C@@:1]`) describes an
//     absolute configuration relative to the TEMPLATE's own neighbour
//     write-order, not the reactant's. Naively copying the symbol without
//     accounting for a reordered template ignores exactly the kind of
//     reorder that inverts configuration (mirroring the parity math
//     `smirks_chirality_ok`/`permutation_parity` already does for
//     REACTANT-side match validation, just never applied to the product).
//   - An INHERITED flag (no explicit template chirality) can survive on an
//     atom whose real bonding topology changed within the reaction's own
//     mapped/core region (e.g. a ring bond broken by the reaction) --
//     `build_product`'s old invalidation heuristic only ever compared
//     *unmapped* substituent element sets, never mapped-neighbour identity.
//
// Both are handled here, in a single post-build pass (run after all of
// `build_product`'s bonds exist, since validating either case requires the
// atom's REAL final neighbour set): an explicit template flag is
// transcribed together with its own neighbour order (validated against the
// atom's real adjacency, not re-derived); an inherited flag is kept only if
// both an order exists to trust it by and the core neighbourhood provably
// didn't change, and cleared otherwise. Two mechanisms, not one unified
// permutation-parity re-expression, because an unmapped/template-literal
// substituent has no reactant identity to map an order back to at all.

/// Map a template-space stereo order into product-index space via
/// `template_idx_to_new`. [`STEREO_H_SENTINEL`] passes through unchanged.
/// `None` if any real (non-sentinel) template atom didn't make it into the
/// product (shouldn't happen for a well-formed template, but fails closed).
pub(super) fn map_stereo_order(
    order: &[u32],
    template_idx_to_new: &[Option<AtomIdx>],
) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                template_idx_to_new
                    .get(t as usize)
                    .copied()
                    .flatten()
                    .map(|a| a.0)
            }
        })
        .collect()
}

/// True when `order` (already in product-index space) matches `new_idx`'s
/// REAL final neighbour set in `product`: at most one H-sentinel, no
/// duplicate real entries, right length, and the real entries are exactly
/// `new_idx`'s actual adjacency. Deliberately a set/degree check, not a
/// permutation-parity one -- `order` is stored verbatim once validated;
/// `corrected_chirality` (the SMILES writer) and `chematic-cip` already do
/// the write-order-independent parity math whenever they consume it.
pub(super) fn order_matches_final_topology(
    product: &Molecule,
    new_idx: AtomIdx,
    order: &[u32],
) -> bool {
    let sentinels = order.iter().filter(|&&x| x == STEREO_H_SENTINEL).count();
    if sentinels > 1 {
        return false;
    }
    let real: Vec<u32> = order
        .iter()
        .copied()
        .filter(|&x| x != STEREO_H_SENTINEL)
        .collect();
    let real_set: FxHashSet<u32> = real.iter().copied().collect();
    if real.len() != real_set.len() {
        return false; // Duplicate entry -- malformed order.
    }
    if order.len() != product.degree(new_idx) + sentinels {
        return false;
    }
    let actual: FxHashSet<u32> = product.neighbors(new_idx).map(|(nb, _)| nb.0).collect();
    real_set == actual
}

/// Bug-B validity gate for an atom that inherited its chirality flag from
/// the reactant (no explicit product-template chirality): the original
/// unmapped-substituent element-multiset comparison, PLUS a symmetric
/// mapped-neighbour atom-map-number-set comparison -- closing the gap where
/// a core/mapped neighbour's topology changes (e.g. a ring bond broken by
/// the reaction) invisibly to the element-multiset-only check.
pub(super) fn bug_b_topology_unchanged(
    product_template: &Molecule,
    ti: AtomIdx,
    reactant: &Molecule,
    src_idx: AtomIdx,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    mol_idx: usize,
    atom_to_map: &FxHashMap<(usize, AtomIdx), u16>,
) -> bool {
    let mut prod_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut prod_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in product_template.neighbors(ti) {
        match product_template.atom(nb).atom_map {
            None => {
                *prod_elems
                    .entry(product_template.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
            Some(am) => {
                prod_mapped.insert(am);
            }
        }
    }

    let mut rxn_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut rxn_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in reactant.neighbors(src_idx) {
        if !all_template_atoms.contains(&(mol_idx, nb)) {
            continue;
        }
        match atom_to_map.get(&(mol_idx, nb)) {
            Some(&am) => {
                rxn_mapped.insert(am);
            }
            None => {
                *rxn_elems
                    .entry(reactant.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
        }
    }

    // NOTE: `prod_elems.is_empty()` is deliberately NOT special-cased here.
    // An empty `prod_elems` legitimately means "product-template atom has no
    // unmapped neighbours" -- which is also true when the reactant's
    // template-matched unmapped substituents (F/Cl/Br etc., matched by the
    // reactant SMARTS itself, not carried-through remote substituents --
    // those never enter `rxn_elems` at all, since `rxn_elems` only counts
    // neighbours inside `all_template_atoms`) were silently DELETED by the
    // product template. `prod_elems == rxn_elems` alone correctly requires
    // both to be empty together, or both to hold the same multiset --
    // exactly what "unchanged" means; special-casing empty-on-one-side
    // would let a genuine deletion through undetected.
    prod_elems == rxn_elems && prod_mapped == rxn_mapped
}

/// Map a REACTANT-space stereo order (atom indices into
/// `input_mols[mol_idx]`) into product-index space via `src_to_new`.
/// [`STEREO_H_SENTINEL`] passes through unchanged. `None` if any real
/// (non-sentinel) reactant neighbour never made it into the product (e.g. it
/// was part of a leaving group the reaction consumed) -- the correspondence
/// is not derivable, so the caller must fail closed rather than guess.
pub(super) fn remap_reactant_stereo_order(
    order: &[u32],
    mol_idx: usize,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                src_to_new.get(&(mol_idx, AtomIdx(t))).map(|a| a.0)
            }
        })
        .collect()
}

/// The reactant stereo order with its implicit-H slot given to the one
/// neighbour the product atom gained, when that is the only change: every
/// other recorded neighbour maps into the product and stays bonded, and the
/// product atom has no H left. `None` otherwise.
pub(super) fn gained_neighbour_order(
    order: &[u32],
    mol_idx: usize,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    product: &Molecule,
    new_idx: AtomIdx,
) -> Option<Vec<u32>> {
    if order.iter().filter(|&&t| t == STEREO_H_SENTINEL).count() != 1 {
        return None;
    }
    let remapped = remap_reactant_stereo_order(order, mol_idx, src_to_new)?;
    let kept: FxHashSet<u32> = remapped
        .iter()
        .copied()
        .filter(|&t| t != STEREO_H_SENTINEL)
        .collect();
    let neighbours: Vec<AtomIdx> = product.neighbors(new_idx).map(|(nb, _)| nb).collect();
    if neighbours.len() != kept.len() + 1
        || !kept.iter().all(|&t| neighbours.contains(&AtomIdx(t)))
        || chematic_core::implicit_hcount(product, new_idx) != 0
    {
        return None;
    }
    let gained = neighbours.into_iter().find(|nb| !kept.contains(&nb.0))?;
    if product.atom(gained).element == chematic_core::Element::H {
        return None;
    }
    Some(
        remapped
            .into_iter()
            .map(|t| if t == STEREO_H_SENTINEL { gained.0 } else { t })
            .collect(),
    )
}

/// Post-build chirality correction. Must run after `build_product`'s Steps
/// 3-4 (all bonds added), since validation needs each atom's real final
/// degree/neighbour set. See the module-level doc above for why this is two
/// mechanisms (transcribe-and-validate for an explicit template flag,
/// gate-and-preserve-or-clear for an inherited one), not one.
pub(super) fn correct_product_stereo(
    mut product: Molecule,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Molecule {
    let atom_to_map: FxHashMap<(usize, AtomIdx), u16> =
        global_map.iter().map(|(&am, &k)| (k, am)).collect();

    for i in 0..product_template.atom_count() {
        let ti = AtomIdx(i as u32);
        let Some(new_idx) = template_idx_to_new[i] else {
            continue;
        };
        let tmpl_atom = product_template.atom(ti);

        if tmpl_atom.chirality != Chirality::None {
            // Explicit product-template chirality: the template is the sole
            // source of truth here regardless of whether this atom is a
            // matched core atom, an unmatched map number, or fully new --
            // no global_map lookup needed on this branch.
            let resolved = product_template
                .stereo_neighbor_order(ti)
                .and_then(|order| map_stereo_order(order, template_idx_to_new))
                .filter(|order| order_matches_final_topology(&product, new_idx, order));
            match resolved {
                Some(order) => {
                    product.set_chirality(new_idx, tmpl_atom.chirality);
                    product.set_stereo_neighbor_order(new_idx, order);
                }
                None => product.set_chirality(new_idx, Chirality::None),
            }
        } else if let Some(am) = tmpl_atom.atom_map
            && let Some(&(mol_idx, src_idx)) = global_map.get(&am)
        {
            // Matched core atom, no explicit template chirality: keep the
            // flag `src_atom.clone()` gave it in Step 1 only if it's
            // actually usable and provably still valid. Three-way rule,
            // all conditions required:
            //   1. bug_b_topology_unchanged: the TEMPLATE-level neighbour
            //      composition (unmapped element multiset + mapped atom-map
            //      set) is unchanged reactant-template -> product-template.
            //   2. remap_reactant_stereo_order succeeds: every atom the
            //      reactant's own recorded order references is uniquely
            //      identifiable in the product (via src_to_new) -- fails
            //      closed (None) for a leaving-group neighbour the reaction
            //      consumed, where no correspondence is derivable at all.
            //   3. order_matches_final_topology: the remapped order's real
            //      entries are EXACTLY this atom's real final adjacency in
            //      the built product molecule -- catches a neighbour that
            //      survived into the product (so (2) succeeds) but whose
            //      specific BOND to this atom did not (e.g. reattached
            //      elsewhere by the reaction), which a template-level
            //      heuristic alone cannot see.
            // Never keep a raw flag without a validated order alongside it.
            let src_atom = input_mols[mol_idx].atom(src_idx);
            if src_atom.chirality != Chirality::None {
                let topology_unchanged = bug_b_topology_unchanged(
                    product_template,
                    ti,
                    input_mols[mol_idx],
                    src_idx,
                    all_template_atoms,
                    mol_idx,
                    &atom_to_map,
                );
                let remapped_order = topology_unchanged
                    .then(|| input_mols[mol_idx].stereo_neighbor_order(src_idx))
                    .flatten()
                    .and_then(|order| remap_reactant_stereo_order(order, mol_idx, src_to_new))
                    .filter(|order| order_matches_final_topology(&product, new_idx, order));
                // Case B: the centre lost its one implicit H and gained one
                // new neighbour, every other neighbour kept. The new bond
                // takes the H's place, as RDKit's reaction runner copies the
                // tag (C-H functionalization keeps the configuration).
                let remapped_order = remapped_order.or_else(|| {
                    input_mols[mol_idx]
                        .stereo_neighbor_order(src_idx)
                        .and_then(|order| {
                            gained_neighbour_order(order, mol_idx, src_to_new, &product, new_idx)
                        })
                });
                match remapped_order {
                    Some(order) => {
                        product.set_chirality(new_idx, src_atom.chirality);
                        product.set_stereo_neighbor_order(new_idx, order);
                    }
                    None => product.set_chirality(new_idx, Chirality::None),
                }
            }
        }
        // else: no template chirality, and not a matched core atom --
        // new_atom.chirality is already Chirality::None from Step 1's clone
        // of a template-literal atom.
    }

    product
}
