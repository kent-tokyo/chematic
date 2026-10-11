/// Per-atom explicit hydrogen count vector (position-sensitive, for 1,2-shift dedup).
fn h_assignment(mol: &Molecule) -> Vec<Option<u32>> {
    (0..mol.atom_count())
        .map(|i| mol.atom(AtomIdx(i as u32)).hydrogen_count.map(|h| h as u32))
        .collect()
}

/// Enumerate all forms reachable from `start` by chained direct aromatic 1,2-shifts.
///
/// Returns every distinct form found via BFS (excluding `start` itself).
/// This covers the full ring orbit for systems like tetrazole where one H can
/// reach any of the N positions in a 5-membered ring through a chain of shifts.
fn enumerate_direct_aromatic_forms(
    start: &Molecule,
    blocked: &HashSet<AtomIdx>,
    max: usize,
    rank: &[u32],
) -> Vec<Molecule> {
    enumerate_direct_aromatic_forms_tracked(start, blocked, max, rank).0
}

/// Like [`enumerate_direct_aromatic_forms`], but also reports whether the
/// search was cut short by `max` while more candidates were still reachable
/// (`true`) versus exhausting the frontier naturally (`false`) -- needed by
/// [`tautomer_parent`] to distinguish `Completed` from `MaxTautomersReached`.
///
/// `rank` (see the "Atom-order canonicalization" section above) makes
/// `find_direct_aromatic_matches`'s per-step candidate order a function of
/// the molecule's structure, not raw `AtomIdx` -- relevant when `max` cuts
/// the BFS short: which subset of candidates got explored first must not
/// depend on the caller's atom labeling.
fn enumerate_direct_aromatic_forms_tracked(
    start: &Molecule,
    blocked: &HashSet<AtomIdx>,
    max: usize,
    rank: &[u32],
) -> (Vec<Molecule>, bool) {
    let mut result = Vec::new();
    let mut seen: HashSet<Vec<Option<u32>>> = HashSet::new();
    seen.insert(h_assignment(start));
    let mut frontier = vec![start.clone()];
    let mut truncated = false;

    'outer: while !frontier.is_empty() {
        if result.len() >= max {
            truncated = true;
            break;
        }
        let current = frontier.remove(0);
        for (d, a) in find_direct_aromatic_matches(&current, rank) {
            if result.len() >= max {
                truncated = true;
                break 'outer;
            }
            if let Some(next) = transfer_hydrogen_aromatic(&current, d, a, blocked) {
                let ha = h_assignment(&next);
                if !seen.contains(&ha) {
                    seen.insert(ha);
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
    }
    (result, truncated)
}

/// Find adjacent pairs of aromatic N atoms where the first (donor) carries an explicit H.
///
/// Returns (donor, acceptor) pairs for direct 1,2-shift (no bridge atom),
/// sorted by canonical rank (see the "Atom-order canonicalization" section
/// above) rather than raw `AtomIdx`/adjacency-list order.
fn find_direct_aromatic_matches(mol: &Molecule, rank: &[u32]) -> Vec<(AtomIdx, AtomIdx)> {
    let mut pairs = Vec::new();
    for (d, _) in mol.atoms() {
        let da = mol.atom(d);
        if !da.aromatic || da.element.atomic_number() != 7 {
            continue;
        }
        if da.hydrogen_count.is_none_or(|h| h == 0) {
            continue;
        }
        for (a, _) in mol.neighbors(d) {
            let aa = mol.atom(a);
            if aa.aromatic && (aa.element.atomic_number() == 7 || aa.element.atomic_number() == 8) {
                pairs.push((d, a));
            }
        }
    }
    pairs.sort_by_key(|&(d, a)| (rank[d.0 as usize], rank[a.0 as usize]));
    pairs
}

/// Transfer one H from an aromatic donor to an aromatic acceptor without changing bond orders.
///
/// Only handles atoms with an explicit `hydrogen_count` on the donor.
/// Returns `None` if donor or acceptor is in `blocked_atoms`.
///
/// Donor/acceptor are aromatic (sp2) atoms, never tetrahedral stereocenters, so this
/// function itself never touches `chirality`/`stereo_neighbor_order` -- but the
/// passthrough rebuild below still needs to copy those side channels (plus
/// `stereo_groups`/`bond_directions`) verbatim for every *other* atom, or any
/// pre-existing stereocenter elsewhere in the molecule silently loses the reference
/// order its `@`/`@@` is defined relative to on every tautomer step.
fn transfer_hydrogen_aromatic(
    mol: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    let donor_h = mol.atom(donor).hydrogen_count?;
    if donor_h == 0 {
        return None;
    }
    let acceptor_h = mol.atom(acceptor).hydrogen_count.unwrap_or(0);

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            // Use None (implicit valence) rather than Some(0) when the donor
            // reaches 0 H; this keeps the canonical SMILES clean (n, not [n]).
            atom.hydrogen_count = donor_h.checked_sub(1).filter(|&h| h > 0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(acceptor_h.saturating_add(1));
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        builder
            .add_bond(b.atom1, b.atom2, b.order)
            .expect("transfer_hydrogen_aromatic: bond from a valid molecule must be re-addable");
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();

    // Issue #415: `acceptor` being individually valence-legal (has an
    // available H slot in isolation) isn't sufficient in a fused/bridged
    // ring -- it can be ring-adjacent to *another* atom that already carries
    // an "extra" H, and two such atoms next to each other in one aromatic
    // ring can't both correctly contribute a lone pair to the same ring's pi
    // system. Confirmed via a real repro this exact function reaches
    // (`Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2`'s two ring nitrogens end up
    // adjacent-both-protonated, an over-valent state RDKit itself rejects).
    // Same general-purpose oracle used for the same reason in
    // `transfer_hydrogen_exocyclic_lactam` below -- see that function's own
    // comment for why `kekulize` (not `validate_valence`, which doesn't cap
    // aromatic atoms to their primary valence and misses this) is the right
    // check. Bounded cost: this function is only ever called on the (`max`,
    // default 16) candidates `enumerate_direct_aromatic_forms_tracked`'s BFS
    // frontier actually visits, not on every possible pair in the molecule.
    if chematic_core::kekulize(&next).is_err() {
        return None;
    }
    Some(next)
}

// ---------------------------------------------------------------------------
// Aromatic lactam/lactim shift (ROADMAP.md Phase 2, round 2C-2)
//
// Generalizes the ring-internal-only mechanism above to admit one more mobile-H
// position: an exocyclic oxygen singly bonded to an aromatic ring atom. See
// docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.4a for the full
// mechanism table and the empirical reason this must be a directional step
// (like the 41 `TautomerRule`s' own `prefer_forward` direction) rather than a
// new candidate generator feeding `enumerate_direct_aromatic_forms`'s
// score-ranked pool -- reusing that pool's `tautomer_score` would select the
// wrong (lactim) tautomer for 2-pyridone (measured: enol scores 1100 vs keto
// 1050).
// ---------------------------------------------------------------------------

/// Ring-path distance, in bond hops along the one ring `ring` (as returned by
/// [`chematic_perception::find_sssr`]), between two members of that ring.
/// `None` if either atom is not a member of `ring`.
fn ring_distance(ring: &[AtomIdx], a: AtomIdx, b: AtomIdx) -> Option<usize> {
    let ia = ring.iter().position(|&x| x == a)?;
    let ib = ring.iter().position(|&x| x == b)?;
    let d = ia.abs_diff(ib);
    Some(d.min(ring.len() - d))
}

/// Find every (donor, bridge, acceptor) triple eligible for the exocyclic
/// lactam/lactim shift. All three sides are deliberately narrow and
/// fail-closed -- restricted to exactly what the 5 confirmed-broken design
/// molecules evidence, not "aromatic ring + heteroatom" in general:
///
/// - `bridge`: a **neutral aromatic carbon**. Not any aromatic element --
///   an aromatic N/S/P/B bridge (e.g. an N-hydroxy heterocycle) is a
///   structurally different system this mechanism was never evidenced
///   against, and a charged bridge changes the valence arithmetic below.
/// - `donor`: a **neutral, exocyclic, non-aromatic oxygen** with exactly
///   one transferable H (checked via `implicit_hcount`, not
///   `.hydrogen_count` -- an unbracketed hydroxyl like `Oc1...` stores
///   `None` and relies on valence inference), bonded to `bridge` by a
///   `Single` bond that is its **only** bond (degree 1 -- an ether/bridging
///   oxygen with a second heavy-atom connection is not a lactam/lactim
///   hydroxyl and must never match).
/// - `acceptor`: a **neutral, pyridine-type aromatic nitrogen** -- no H,
///   not charged, and with **exactly 2** heavy-atom connections. Degree 2
///   is the valence-compatibility condition itself: an organic-subset N's
///   normal valence is 3, so a 2-connected aromatic N always has exactly
///   one free valence slot for the incoming H; a 3-connected aromatic N
///   (pyrrole-type, or a bridgehead/fusion position) has none, and adding
///   an H there would require a charge change this mechanism must never
///   make.
/// - sharing an SSSR ring with `bridge` at **odd** ring distance (the
///   condition for a real alternating single/double path between them to
///   exist at all -- RFC section 4.4a).
///
/// O-only, not O-or-N: every confirmed-broken molecule is this type (RFC
/// section 1.1); the analogous N-acceptor (amino/imino) case is a distinct,
/// unevidenced defect, deliberately out of scope (RFC section 4.4a).
///
/// Deduplicated and returned in a fixed (donor, bridge, acceptor) order,
/// never hash-iteration order: `bridge` can belong to more than one SSSR
/// ring in a fused system (hypoxanthine), which would otherwise report the
/// same triple once per ring it's found through.
fn find_exocyclic_lactam_shift_matches(
    mol: &Molecule,
    rank: &[u32],
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let rings = chematic_perception::find_sssr(mol);
    let mut out: HashSet<(AtomIdx, AtomIdx, AtomIdx)> = HashSet::new();
    for (bridge, bridge_atom) in mol.atoms() {
        if !bridge_atom.aromatic
            || bridge_atom.element.atomic_number() != 6
            || bridge_atom.charge != 0
        {
            continue;
        }
        for (donor, bond_idx) in mol.neighbors(bridge) {
            let donor_atom = mol.atom(donor);
            if donor_atom.aromatic
                || donor_atom.element.atomic_number() != 8
                || donor_atom.charge != 0
            {
                continue;
            }
            if implicit_hcount(mol, donor) != 1 {
                continue;
            }
            if mol.neighbors(donor).count() != 1 {
                continue;
            }
            if mol.bond(bond_idx).order != BondOrder::Single {
                continue;
            }
            for ring in rings.rings() {
                if !ring.contains(&bridge) {
                    continue;
                }
                for &acceptor in ring {
                    let acceptor_atom = mol.atom(acceptor);
                    if !acceptor_atom.aromatic
                        || acceptor_atom.element.atomic_number() != 7
                        || acceptor_atom.charge != 0
                    {
                        continue;
                    }
                    if implicit_hcount(mol, acceptor) != 0 {
                        continue;
                    }
                    if mol.neighbors(acceptor).count() != 2 {
                        continue;
                    }
                    if ring_distance(ring, bridge, acceptor).is_some_and(|d| d % 2 == 1) {
                        out.insert((donor, bridge, acceptor));
                    }
                }
            }
        }
    }
    let mut out: Vec<_> = out.into_iter().collect();
    out.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    out
}

/// Move the H from `donor` (exocyclic O) to `acceptor` (ring N), flipping
/// only the `donor`-`bridge` bond order (`Single` -> `Double`); every
/// ring-internal bond, including `bridge`-`acceptor` when they're directly
/// bonded, stays exactly as it was (`Aromatic`) -- confirmed structurally
/// true for every design molecule (RFC section 4.4a).
fn transfer_hydrogen_exocyclic_lactam(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    // implicit_hcount, not `.hydrogen_count`: an unbracketed hydroxyl's H is
    // stored as `None` (see find_exocyclic_lactam_shift_matches's own note).
    let donor_h = implicit_hcount(mol, donor);
    if donor_h == 0 {
        return None;
    }
    let (bond_idx, bond) = mol.bond_between(donor, bridge)?;
    if bond.order != BondOrder::Single || blocked_bonds.contains(&bond_idx) {
        return None;
    }
    let acceptor_h = mol.atom(acceptor).hydrogen_count.unwrap_or(0);

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            atom.hydrogen_count = donor_h.checked_sub(1).filter(|&h| h > 0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(acceptor_h.saturating_add(1));
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = if bidx == bond_idx {
            BondOrder::Double
        } else {
            b.order
        };
        builder.add_bond(b.atom1, b.atom2, order).expect(
            "transfer_hydrogen_exocyclic_lactam: bond from a valid molecule must be re-addable",
        );
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();

    if !exocyclic_lactam_shift_preserves_invariants(mol, &next, donor, acceptor, bond_idx) {
        // Defense in depth: the narrow preconditions above are meant to
        // guarantee this by construction, but `MoleculeBuilder::build()` can
        // re-perceive/normalize state beyond what the explicit field copy
        // above touched (the class of bug this project has hit before --
        // see Phase 1's fragment-extraction stereo-corruption fix). Fail
        // closed rather than return a molecule that silently changed
        // something this transform must never touch.
        return None;
    }

    // Second defense-in-depth layer (issue #415): the per-atom preconditions
    // above (`implicit_hcount(mol, acceptor) == 0`, degree 2) are necessary
    // but not sufficient in a fused/bridged ring system -- `find_sssr`'s own
    // ring choice for such systems isn't always unique (a known, separate
    // residual: the *same* physical molecule can get a different SSSR
    // decomposition depending on atom-labeling order, which can change which
    // ring this function's `ring_distance` check sees), and an acceptor that
    // looks individually valid (0 H, degree 2) can still be ring-adjacent to
    // *another* aromatic atom that already carries an "extra" H. Two such
    // atoms next to each other in one aromatic ring cannot both correctly
    // contribute a lone pair to the same ring's pi system -- confirmed via a
    // real repro (`Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2`) that this function
    // previously accepted and which produced an over-valent `[nH2]` nitrogen
    // RDKit itself rejects outright. `kekulize` is a general, already-used-
    // elsewhere oracle for "is this aromatic system actually valid" that
    // catches this (and would catch other, not-yet-evidenced pathological
    // shapes too) without needing to hand-enumerate every such shape here.
    if chematic_core::kekulize(&next).is_err() {
        return None;
    }
    Some(next)
}

/// Reverse the aromatic lactam/lactim edge: move H from a ring N back to a
/// carbonyl oxygen and reduce only the exocyclic C=O bond. This is used by
/// exhaustive tautomer enumeration so keto inputs can reach the same
/// candidate component as hydroxyl inputs; canonical selection remains
/// directional and continues to prefer the lactam side.
fn transfer_hydrogen_exocyclic_lactim(
    mol: &Molecule,
    bridge: AtomIdx,
    oxygen: AtomIdx,
    donor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&oxygen)
        || blocked_atoms.contains(&donor)
        || mol.atom(donor).hydrogen_count != Some(1)
    {
        return None;
    }
    let (bond_idx, bond) = mol.bond_between(oxygen, bridge)?;
    if bond.order != BondOrder::Double || blocked_bonds.contains(&bond_idx) {
        return None;
    }
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == oxygen {
            atom.hydrogen_count = Some(1);
        } else if idx == donor {
            atom.hydrogen_count = None;
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = if bidx == bond_idx {
            BondOrder::Single
        } else {
            b.order
        };
        builder.add_bond(b.atom1, b.atom2, order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();
    chematic_core::kekulize(&next).ok().map(|_| next)
}

fn reverse_exocyclic_lactim_candidates(mol: &Molecule, config: &TautomerConfig) -> Vec<Molecule> {
    let rings = chematic_perception::find_sssr(mol);
    let mut candidates = Vec::new();
    for (bridge, atom) in mol.atoms() {
        if !atom.aromatic || atom.element.atomic_number() != 6 || atom.charge != 0 {
            continue;
        }
        for (oxygen, bond_idx) in mol.neighbors(bridge) {
            if mol.atom(oxygen).element.atomic_number() != 8
                || mol.atom(oxygen).aromatic
                || mol.bond(bond_idx).order != BondOrder::Double
            {
                continue;
            }
            for ring in rings.rings() {
                if !ring.contains(&bridge) {
                    continue;
                }
                for &donor in ring {
                    let donor_atom = mol.atom(donor);
                    if !donor_atom.aromatic
                        || donor_atom.element.atomic_number() != 7
                        || donor_atom.charge != 0
                        || donor_atom.hydrogen_count != Some(1)
                        || mol.neighbors(donor).count() != 2
                        || ring_distance(ring, bridge, donor).is_none_or(|d| d % 2 != 1)
                    {
                        continue;
                    }
                    if let Some(next) = transfer_hydrogen_exocyclic_lactim(
                        mol,
                        bridge,
                        oxygen,
                        donor,
                        &config.blocked_atoms,
                        &config.blocked_bonds,
                    ) {
                        candidates.push(next);
                    }
                }
            }
        }
    }
    candidates.sort_by_key(|mol| (mol_fingerprint(mol), chematic_smiles::canonical_smiles(mol)));
    candidates.dedup_by(|a, b| mol_fingerprint(a) == mol_fingerprint(b));
    candidates
}

fn exocyclic_lactam_candidates(
    mol: &Molecule,
    rank: &[u32],
    config: &TautomerConfig,
) -> Vec<Molecule> {
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
        })
        .collect()
}

/// Canonical key rooted at the aromatic carbonyl. Explicit H annotations are
/// removed so equivalent lactam spellings are compared by the carbonyl's
/// ring environment rather than by the input proton placement.
fn carbonyl_rooted_key(mol: &Molecule) -> String {
    let root = mol.atoms().find_map(|(idx, atom)| {
        (atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(idx).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && mol.bond(bond).order == BondOrder::Double
            }))
        .then_some(idx)
    });
    let Some(root) = root else {
        return chematic_smiles::canonical_smiles(mol);
    };
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.hydrogen_count = None;
        if idx == root {
            atom.isotope = Some(65535);
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        builder.add_bond(bond.atom1, bond.atom2, bond.order).ok();
    }
    chematic_smiles::canonical_smiles(&builder.build())
}

fn has_aromatic_exocyclic_carbonyl(mol: &Molecule) -> bool {
    mol.atoms().any(|(bridge, atom)| {
        atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(bridge).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && mol.bond(bond).order == BondOrder::Double
            })
    })
}

fn aromatic_exocyclic_carbonyl_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(bridge, atom)| {
            atom.aromatic
                && atom.element.atomic_number() == 6
                && mol.neighbors(*bridge).any(|(neighbor, bond)| {
                    mol.atom(neighbor).element.atomic_number() == 8
                        && mol.bond(bond).order == BondOrder::Double
                })
        })
        .count()
}

/// An aromatic carbon with an exocyclic O that can take part in a
/// lactam/lactim shift: `c=O` or `c-[OH]`. An ether O (`c-OC`, methoxy) cannot;
/// counting it sent every methoxyarene through the global O/C selection,
/// whose O-H preference then turned unrelated amides into iminols.
fn has_aromatic_exocyclic_oxygen(mol: &Molecule) -> bool {
    mol.atoms().any(|(bridge, atom)| {
        atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(bridge).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && match mol.bond(bond).order {
                        BondOrder::Double => true,
                        BondOrder::Single => implicit_hcount(mol, neighbor) > 0,
                        _ => false,
                    }
            })
    })
}

/// Non-aromatic iminol groups, `C(-[OH])=N`.
fn iminol_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|&(c, atom)| {
            !atom.aromatic
                && atom.element.atomic_number() == 6
                && mol.neighbors(c).any(|(o, b)| {
                    mol.atom(o).element.atomic_number() == 8
                        && mol.bond(b).order == BondOrder::Single
                        && implicit_hcount(mol, o) > 0
                })
                && mol.neighbors(c).any(|(n, b)| {
                    mol.atom(n).element.atomic_number() == 7
                        && mol.bond(b).order == BondOrder::Double
                })
        })
        .count()
}
