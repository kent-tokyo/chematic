/// Find the dual-flank N-H ambiguity of aromatic lactam systems such as
/// cytosine, guanine, and hypoxanthine.  The carbonyl carbon is bonded to two
/// aromatic nitrogens; moving H between those nitrogens changes neither the
/// graph nor any bond order, but is required to make the two keto spellings
/// canonicalize identically.
fn find_dual_flank_matches(mol: &Molecule, rank: &[u32]) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut out = Vec::new();
    for (bridge, atom) in mol.atoms() {
        if !atom.aromatic || atom.element.atomic_number() != 6 || atom.charge != 0 {
            continue;
        }
        let has_carbonyl_oxygen = mol.neighbors(bridge).any(|(neighbor, bidx)| {
            mol.atom(neighbor).element.atomic_number() == 8
                && mol.bond(bidx).order == BondOrder::Double
        });
        if !has_carbonyl_oxygen {
            continue;
        }
        for ring in chematic_perception::find_sssr(mol).rings().iter() {
            if !ring.contains(&bridge) {
                continue;
            }
            let flank_n: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|&candidate| {
                    let a = mol.atom(candidate);
                    a.aromatic
                        && a.element.atomic_number() == 7
                        && a.charge == 0
                        && ring_distance(ring, bridge, candidate).is_some_and(|d| d % 2 == 1)
                })
                .collect();
            for &a in &flank_n {
                for &b in &flank_n {
                    if a == b {
                        continue;
                    }
                    let (donor, acceptor) = if mol.atom(a).hydrogen_count == Some(1)
                        && mol.atom(b).hydrogen_count.unwrap_or(0) == 0
                    {
                        (a, b)
                    } else {
                        continue;
                    };
                    out.push((donor, bridge, acceptor));
                }
            }
        }
    }
    out.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    out.dedup();
    out
}

/// Move H between the two aromatic N atoms flanking a carbonyl carbon.
/// Unlike the lactam/lactim shift, no bond order changes are needed.
fn transfer_hydrogen_dual_flank(
    mol: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    if mol.atom(donor).hydrogen_count != Some(1)
        || mol.atom(acceptor).hydrogen_count.unwrap_or(0) != 0
    {
        return None;
    }
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            atom.hydrogen_count = Some(0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(1);
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let b = mol.bond(BondIdx(i as u32));
        builder.add_bond(b.atom1, b.atom2, b.order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();
    chematic_core::kekulize(&next).ok().map(|_| next)
}

fn apply_dual_flank_tracked(
    mol: &Molecule,
    config: &TautomerConfig,
    _rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    // H placement is the ambiguity being normalized here.  Ranking the
    // input molecule directly lets that H change the root of a symmetric
    // pyrimidinone ring, so keto and enol spellings can choose opposite
    // flanks.  Use a rooted structural key made from the same graph with
    // explicit H annotations removed; atom indices remain unchanged.
    let rooted_rank = tautomer_skeleton_rank(mol);
    find_dual_flank_matches(mol, &rooted_rank)
        .into_iter()
        .find_map(|(donor, bridge, acceptor)| {
            let next = transfer_hydrogen_dual_flank(mol, donor, acceptor, &config.blocked_atoms)?;
            (rooted_structural_key(mol, donor) < rooted_structural_key(mol, acceptor))
                .then_some((next, donor, bridge, acceptor))
        })
}

/// Everything this transform must hold, checked on the *built* molecule, not
/// assumed from the field copy above: every atom except `donor`/`acceptor`
/// is byte-for-byte identical (element, charge, isotope, aromaticity,
/// chirality, H count); `donor`/`acceptor` keep everything except H count,
/// which must move by exactly 1 in the expected direction; every bond except
/// `bond_idx` keeps its endpoints and order; `bond_idx` keeps its endpoints
/// and its order must have actually flipped `Single` -> `Double`.
fn exocyclic_lactam_shift_preserves_invariants(
    before: &Molecule,
    after: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    bond_idx: BondIdx,
) -> bool {
    if before.atom_count() != after.atom_count() || before.bond_count() != after.bond_count() {
        return false;
    }
    for i in 0..before.atom_count() {
        let idx = AtomIdx(i as u32);
        let a = before.atom(idx);
        let b = after.atom(idx);
        if a.element != b.element
            || a.charge != b.charge
            || a.isotope != b.isotope
            || a.aromatic != b.aromatic
            || a.chirality != b.chirality
        {
            return false;
        }
        let h_before = implicit_hcount(before, idx);
        let h_after = implicit_hcount(after, idx);
        if idx == donor {
            if h_after != h_before.wrapping_sub(1) {
                return false;
            }
        } else if idx == acceptor {
            if h_after != h_before + 1 {
                return false;
            }
        } else if h_after != h_before {
            return false;
        }
    }
    for i in 0..before.bond_count() {
        let bidx = BondIdx(i as u32);
        let bb = before.bond(bidx);
        let ba = after.bond(bidx);
        if bb.atom1 != ba.atom1 || bb.atom2 != ba.atom2 {
            return false;
        }
        if bidx == bond_idx {
            if bb.order != BondOrder::Single || ba.order != BondOrder::Double {
                return false;
            }
        } else if bb.order != ba.order {
            return false;
        }
    }
    true
}

/// Apply one exocyclic lactam/lactim shift, chosen deterministically among
/// every eligible (donor, bridge, acceptor) triple by minimal canonical
/// SMILES of the result -- never `.first()`/index order, which would make
/// the choice track input atom order rather than structure (uracil has two
/// independent sites; cytosine's/guanine's carbonyl-bearing ring atom is
/// flanked by two ring nitrogens, only one a valid acceptor at a time, but
/// both may qualify simultaneously from a fully-shifted starting form). See
/// RFC section 4.4a's order-invariance requirement.
fn apply_exocyclic_lactam_shift_tracked(
    mol: &Molecule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    let mut candidates: Vec<(Molecule, AtomIdx, AtomIdx, AtomIdx)> =
        find_exocyclic_lactam_shift_matches(mol, rank)
            .into_iter()
            .filter_map(|(donor, bridge, acceptor)| {
                transfer_hydrogen_exocyclic_lactam(
                    mol,
                    donor,
                    bridge,
                    acceptor,
                    &config.blocked_atoms,
                    &config.blocked_bonds,
                )
                .map(|next| (next, donor, bridge, acceptor))
            })
            .collect();
    if candidates.is_empty() {
        return None;
    }
    candidates.sort_by(|x, y| {
        carbonyl_rooted_key(&y.0)
            .cmp(&carbonyl_rooted_key(&x.0))
            .then_with(|| {
                chematic_smiles::canonical_smiles(&x.0)
                    .cmp(&chematic_smiles::canonical_smiles(&y.0))
            })
    });
    // Candidates sharing the minimal canonical SMILES are meant to be the
    // same molecule (safely interchangeable -- pick the first, deterministic
    // given the sort above). Don't trust canonical_smiles alone for that:
    // this project has hit real canonical-SMILES bugs before. Cross-check
    // with mol_fingerprint, an independently-computed (sorted per-atom
    // element/charge/bond-order-sum) hash already used elsewhere in this
    // file for exactly this "is this really the same molecule" question. A
    // disagreement means the tie is not actually resolved -- fail closed
    // (report no move this iteration) rather than guess.
    let min_smiles = chematic_smiles::canonical_smiles(&candidates[0].0);
    let tied: Vec<&(Molecule, AtomIdx, AtomIdx, AtomIdx)> = candidates
        .iter()
        .take_while(|c| chematic_smiles::canonical_smiles(&c.0) == min_smiles)
        .collect();
    if tied.len() > 1 {
        let first_fp = mol_fingerprint(&tied[0].0);
        if !tied.iter().all(|c| mol_fingerprint(&c.0) == first_fp) {
            return None;
        }
    }
    candidates.into_iter().next()
}

use crate::hash::{FNV1A_OFFSET, FNV1A_PRIME};

/// Tautomer form score for canonical selection: prefers O-H > N-H > S-H and aromatic rings.
fn tautomer_score(mol: &Molecule) -> i32 {
    score_breakdown(mol).iter().map(|c| c.value).sum()
}

/// Same scoring rule as [`tautomer_score`], decomposed into named
/// contributions for [`TautomerAuditRecord::score_breakdown`]. This is the
/// single source of truth for the scoring rule; `tautomer_score` sums it.
fn score_breakdown(mol: &Molecule) -> Vec<ScoreContribution> {
    let mut contributions = Vec::new();
    let mut has_aromatic = false;

    for (_, atom) in mol.atoms() {
        if atom.aromatic {
            has_aromatic = true;
        }

        // Heteroatom hydrogen: O-H > N-H > S-H (explicit or implicit)
        let h_count = atom.hydrogen_count.unwrap_or(0) as i32;
        if h_count > 0 {
            let weight = match atom.element.atomic_number() {
                8 => Some(100), // O-H: highest priority
                7 => Some(50),  // N-H: medium priority
                16 => Some(25), // S-H: low priority
                _ => None,
            };
            if let Some(weight) = weight {
                contributions.push(ScoreContribution {
                    term: TautomerScoreTerm::HeteroatomHydrogen {
                        element: atom.element.atomic_number(),
                    },
                    value: h_count * weight,
                });
            }
        }
    }

    if has_aromatic {
        contributions.push(ScoreContribution {
            term: TautomerScoreTerm::AromaticRing,
            value: 1000,
        });
    }

    contributions
}

/// Order-independent structural hash for convergence detection.
fn mol_fingerprint(mol: &Molecule) -> u64 {
    let mut atoms: Vec<(u8, i8, u32, u8)> = (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let a = mol.atom(idx);
            let bos: u32 = mol
                .neighbors(idx)
                .map(|(_, bidx)| mol.bond(bidx).order.order_int() as u32)
                .sum();
            (
                a.element.atomic_number(),
                a.charge,
                bos,
                a.hydrogen_count.unwrap_or(0),
            )
        })
        .collect();
    atoms.sort();
    let mut hash = FNV1A_OFFSET;
    for (an, ch, bos, h) in atoms {
        hash ^= an as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= (ch as u8 as u64).wrapping_add(128);
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= bos as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= h as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
    }

    // BondOrder::Up/Down and the bond-direction side channel encode 2D E/Z
    // geometry.  The compact atom fingerprint above intentionally ignores
    // those directions, but using it alone makes E and Z forms collide in
    // tautomer search and causes one form to be discarded as "seen".  Add
    // the canonical serialization only for molecules that actually carry a
    // directional bond; this preserves the cheap/order-independent fast path
    // for the overwhelmingly common non-E/Z case.
    let has_direction = (0..mol.bond_count()).any(|i| {
        let bidx = BondIdx(i as u32);
        matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
            || mol.bond_direction(bidx).is_some()
    });
    // Explicit H placement is also position-sensitive: two equivalent
    // aromatic nitrogens can have the same sorted atom multiset while the H
    // sits on different atoms (the cytosine/guanine dual-flank case).
    let has_explicit_h =
        (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).hydrogen_count.is_some());
    if has_direction || has_explicit_h {
        hash ^= 0xff;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        for &byte in chematic_smiles::canonical_smiles(mol).as_bytes() {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(FNV1A_PRIME);
        }
    }
    hash
}

/// Whether `mol` carries an OpenSMILES directional marker that contributes to
/// E/Z stereochemistry.  A tautomer transform can move a double bond away
/// from the marked single bond; until the direction is explicitly remapped to
/// the new stereochemical carrier, transforming such a molecule would silently
/// erase its E/Z identity.
fn has_directional_bond(mol: &Molecule) -> bool {
    (0..mol.bond_count()).any(|i| {
        let bidx = BondIdx(i as u32);
        matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
            || mol.bond_direction(bidx).is_some()
    })
}

/// Find all (donor, bridge, acceptor) triples matching the rule in `mol`.
/// For path_len=3: donor-bridge-acceptor (standard 1,3-shift)
/// For path_len=5: donor-b1-b2-b3-acceptor (1,5-shift) where b2 is stored in "bridge"
fn find_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    rank: &[u32],
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut matches = if rule.path_len == 5 {
        find_15_shift_matches(mol, rule)
    } else {
        find_13_shift_matches(mol, rule)
    };

    // Sort by canonical rank, not raw AtomIdx (parse/insertion order) --
    // whichever match a caller treats as "first" must depend on the
    // molecule's structure alone, not on how its atoms happened to be
    // labeled. See the "Atom-order canonicalization" section above.
    matches.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    matches
}

fn is_rule_donor(mol: &Molecule, atom: AtomIdx, rule: &TautomerRule) -> bool {
    mol.atom(atom).element.atomic_number() == rule.donor_elem && implicit_hcount(mol, atom) > 0
}

fn bridge_element_matches(mol: &Molecule, atom: AtomIdx, expected: Option<u8>) -> bool {
    expected.is_none_or(|element| mol.atom(atom).element.atomic_number() == element)
}

fn find_13_shift_matches(
    mol: &Molecule,
    rule: &TautomerRule,
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut matches = Vec::new();
    for i in 0..mol.atom_count() {
        let donor = AtomIdx(i as u32);
        if !is_rule_donor(mol, donor, rule) {
            continue;
        }
        for (bridge, donor_bond) in mol.neighbors(donor) {
            if !rule.donor_bridge_order.matches(mol.bond(donor_bond).order)
                || !bridge_element_matches(mol, bridge, rule.bridge_elem)
            {
                continue;
            }
            for (acceptor, acceptor_bond) in mol.neighbors(bridge) {
                if acceptor == donor
                    || !rule
                        .bridge_acceptor_order
                        .matches(mol.bond(acceptor_bond).order)
                    || mol.atom(acceptor).element.atomic_number() != rule.acceptor_elem
                {
                    continue;
                }
                matches.push((donor, bridge, acceptor));
            }
        }
    }
    matches
}

fn find_15_shift_matches(
    mol: &Molecule,
    rule: &TautomerRule,
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut matches = Vec::new();
    for i in 0..mol.atom_count() {
        let donor = AtomIdx(i as u32);
        if !is_rule_donor(mol, donor, rule) {
            continue;
        }
        for (first, donor_bond) in mol.neighbors(donor) {
            if !rule.donor_bridge_order.matches(mol.bond(donor_bond).order)
                || mol.atom(first).element.atomic_number() != 6
            {
                continue;
            }
            collect_15_shift_from_first(mol, rule, donor, first, &mut matches);
        }
    }
    matches
}

fn collect_15_shift_from_first(
    mol: &Molecule,
    rule: &TautomerRule,
    donor: AtomIdx,
    first: AtomIdx,
    matches: &mut Vec<(AtomIdx, AtomIdx, AtomIdx)>,
) {
    for (central, _) in mol.neighbors(first) {
        if central == donor || !bridge_element_matches(mol, central, rule.bridge_elem) {
            continue;
        }
        for (third, _) in mol.neighbors(central) {
            if third == first || mol.atom(third).element.atomic_number() != 6 {
                continue;
            }
            for (acceptor, acceptor_bond) in mol.neighbors(third) {
                if acceptor != central
                    && rule
                        .bridge_acceptor_order
                        .matches(mol.bond(acceptor_bond).order)
                    && mol.atom(acceptor).element.atomic_number() == rule.acceptor_elem
                {
                    matches.push((donor, central, acceptor));
                }
            }
        }
    }
}

/// Apply a single tautomer transformation: donor-bridge bond Single → Double
/// and bridge-acceptor bond Double → Single. Returns `None` if the transform
/// would be invalid (e.g. donor has no explicit H to give up on a bracket atom).
///
/// For organic-subset atoms, implicit H counts adjust automatically through the
/// valence model; for bracket atoms with an explicit `hydrogen_count`, we
/// decrement the donor and increment the acceptor manually.
///
/// Donor/bridge/acceptor participate in a 1,3- or 1,5-shift, never a tetrahedral
/// stereocenter directly -- but the passthrough rebuild below still needs to copy
/// `stereo_neighbor_order`/`stereo_groups`/`bond_directions` verbatim for every
/// *other* atom, same as `transfer_hydrogen_aromatic` above, or any pre-existing
/// stereocenter elsewhere in the molecule silently loses the reference order its
/// `@`/`@@` is defined relative to on every tautomer step.
fn transfer_hydrogen(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor)
        || blocked_atoms.contains(&bridge)
        || blocked_atoms.contains(&acceptor)
    {
        return None;
    }
    let (db_bidx, _) = mol.bond_between(donor, bridge)?;
    let (ba_bidx, _) = mol.bond_between(bridge, acceptor)?;
    if blocked_bonds.contains(&db_bidx) || blocked_bonds.contains(&ba_bidx) {
        return None;
    }

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if let Some(h) = atom.hydrogen_count {
            if idx == donor {
                if h == 0 {
                    return None;
                }
                atom.hydrogen_count = Some(h - 1);
            } else if idx == acceptor {
                atom.hydrogen_count = Some(h.saturating_add(1));
            }
        }
        builder.add_atom(atom);
    }

    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = match bidx {
            x if x == db_bidx => BondOrder::Double,
            x if x == ba_bidx => BondOrder::Single,
            _ => b.order,
        };
        builder.add_bond(b.atom1, b.atom2, order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    Some(builder.build())
}

/// Like the first-match helper, but also returns the (donor, bridge, acceptor)
/// triple that was moved -- needed by [`tautomer_parent`] to record which
/// atoms/bonds each applied transform touched.
fn apply_first_match_tracked(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    find_matches(mol, rule, rank)
        .into_iter()
        .find_map(|(d, b, a)| {
            transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
                .map(|next| (next, d, b, a))
        })
}

/// Apply every matching transformation for `rule`; return all resulting molecules.
fn apply_all_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Vec<Molecule> {
    find_matches(mol, rule, rank)
        .into_iter()
        .filter_map(|(d, b, a)| {
            transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
        })
        .collect()
}

/// Apply all non-overlapping matches of one rule in a single bounded pass.
///
/// A rule match changes the donor/bridge/acceptor and the two bonds joining
/// them.  Matches that touch any of those atoms are therefore conservatively
/// treated as conflicting and deferred to a later pass.  This preserves the
/// old rule priority while allowing independent sites to consume one
/// `max_iter` unit together.  In particular, the result no longer depends on
/// which of many independent sites happened to receive the first raw atom
/// index when the transform budget is reached.
fn apply_nonconflicting_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> (Molecule, bool) {
    let (next, applied) = apply_nonconflicting_matches_tracked(mol, rule, config, rank, usize::MAX);
    (next, !applied.is_empty())
}

/// Tracked form of [`apply_nonconflicting_matches`]. `limit` is the maximum
/// number of transformations to apply in this pass, allowing the audit-aware
/// parent API to preserve its `max_transforms` bound while still making its
/// selection order canonical.
fn apply_nonconflicting_matches_tracked(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
    limit: usize,
) -> (Molecule, Vec<(AtomIdx, AtomIdx, AtomIdx)>) {
    let matches = find_matches(mol, rule, rank);
    if matches.is_empty() || limit == 0 {
        return (mol.clone(), Vec::new());
    }

    let mut current = mol.clone();
    let mut occupied: HashSet<AtomIdx> = HashSet::new();
    let mut applied = Vec::new();

    for (donor, bridge, acceptor) in matches {
        if applied.len() >= limit {
            break;
        }
        let Some((db, _)) = mol.bond_between(donor, bridge) else {
            continue;
        };
        let Some((ba, _)) = mol.bond_between(bridge, acceptor) else {
            continue;
        };
        let touched = [
            donor,
            bridge,
            acceptor,
            mol.bond(db).atom1,
            mol.bond(db).atom2,
            mol.bond(ba).atom1,
            mol.bond(ba).atom2,
        ];
        if touched.iter().any(|atom| occupied.contains(atom)) {
            continue;
        }
        let Some(next) = transfer_hydrogen(
            &current,
            donor,
            bridge,
            acceptor,
            &config.blocked_atoms,
            &config.blocked_bonds,
        ) else {
            continue;
        };
        current = next;
        occupied.extend(touched);
        applied.push((donor, bridge, acceptor));
    }

    (current, applied)
}
