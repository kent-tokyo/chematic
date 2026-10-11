// ---------------------------------------------------------------------------
// TautomerConfig
// ---------------------------------------------------------------------------

/// Configuration for tautomer enumeration and canonicalization.
///
/// # Rule indices
/// Rules are numbered 0-based in the order they appear in the built-in set.
/// Use [`TautomerConfig::rule_count`] to know the total and
/// [`TautomerConfig::rule_names`] to see what each index represents.
///
/// An empty `enabled_rules` (the default) activates **all** rules.
///
/// # Zone blocking
/// `blocked_atoms` and `blocked_bonds` prevent H-transfer through specific atoms/bonds.
/// Any tautomer move whose donor, bridge, or acceptor is in `blocked_atoms`,
/// or whose altered bond is in `blocked_bonds`, is suppressed.
/// Note: for 1,5-shift rules only the donor, bridge-central, and acceptor atoms are
/// checked; intermediate path atoms cannot be blocked without refactoring `find_matches`.
#[derive(Debug, Clone)]
pub struct TautomerConfig {
    /// Maximum iterations in [`canonical_tautomer_with_config`] (default 16).
    ///
    /// Each outer iteration applies all non-overlapping matches of the first
    /// active forward rule. Overlapping matches are deferred and reconsidered
    /// in the next iteration. This keeps the bound meaningful for chained or
    /// conflicting transformations while making independent sites consume one
    /// iteration together, independent of input atom insertion order.
    pub max_iter: usize,
    /// Maximum tautomers returned by [`enumerate_tautomers_with_config`] (default 32).
    pub max_tautomers: usize,
    /// 0-based indices of rules to activate.  Empty = all rules active.
    pub enabled_rules: Vec<usize>,
    /// Atoms that must not participate in any H-transfer (donor, bridge, or acceptor).
    /// Empty = no atom blocking (default).
    pub blocked_atoms: HashSet<AtomIdx>,
    /// Bonds that must not be altered by any H-transfer.
    /// Empty = no bond blocking (default).
    pub blocked_bonds: HashSet<BondIdx>,
}

impl Default for TautomerConfig {
    fn default() -> Self {
        Self {
            max_iter: 16,
            max_tautomers: 32,
            enabled_rules: Vec::new(),
            blocked_atoms: HashSet::new(),
            blocked_bonds: HashSet::new(),
        }
    }
}

impl TautomerConfig {
    /// Number of built-in tautomer rules available.
    pub fn rule_count() -> usize {
        RULES.len()
    }

    /// Names of all built-in rules, in index order.
    pub fn rule_names() -> Vec<&'static str> {
        RULES.iter().map(|r| r.name).collect()
    }

    /// Convenience: config with only the keto-enol rule (index 0) enabled.
    pub fn keto_enol_only() -> Self {
        Self {
            enabled_rules: vec![0],
            ..Self::default()
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helper: iterate only the active rules
// ---------------------------------------------------------------------------

/// Collect the rules that are active under `config`.
fn active_rules(config: &TautomerConfig) -> Vec<&'static TautomerRule> {
    RULES
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if config.enabled_rules.is_empty() || config.enabled_rules.contains(&i) {
                Some(r)
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Atom-order canonicalization (issue #415 residual)
//
// `find_sssr`'s ring choice (consumed directly by
// `find_exocyclic_lactam_shift_matches`) and `find_matches`'/
// `find_direct_aromatic_matches`'s "first match in raw AtomIdx order" policy
// are both dependent on the *caller's* parse/atom-insertion order, not just
// the molecule's structure -- so which tautomer this module reaches for the
// same molecule could depend on how its atoms happened to be labeled.
//
// Fixed by computing a canonical *rank* per atom once
// (`chematic_smiles::canonical_atom_order`, inverted into an `AtomIdx ->
// rank` lookup) and using it to break every such tie -- sort candidates by
// rank instead of by raw `AtomIdx`, and prefer the lowest-rank ring when a
// ring choice is needed. `AtomIdx` identity never changes across a
// tautomer-search step (every transform function rebuilds atom-for-atom via
// `for i in 0..mol.atom_count()`, same order in and out), so one rank
// vector computed from the *entry* molecule stays valid for every
// intermediate form the search visits.
//
// An earlier version of this fix instead rebuilt the molecule with atoms
// physically reordered into `canonical_atom_order`'s own order before
// searching. That was correct (verified: full order-invariance across 6
// permutations of issue #415's residual) but is NOT what's implemented here
// -- it uncovered a separate, pre-existing `chematic-smiles` defect where
// `canonical_smiles`/`canonical_atom_order`'s own individualize-refine
// search degenerates into a multi-minute hang specifically when its input's
// atom order already coincides with `canonical_atom_order`'s own output
// order, for at least one real highly-symmetric molecule
// (`chembl_accuracy_corpus_4999.smi` line 2741, a 94-atom molecule with 3
// near-equivalent Boc-benzylamine arms -- reproduces with zero
// `chematic-chem` involvement: `parse` -> `canonical_atom_order` -> rebuild
// in that order -> `canonical_smiles` hangs). Rank-based tie-breaking never
// constructs that pathological input shape (the molecule's own storage/atom
// order is never touched), so it sidesteps the landmine entirely rather
// than working around it. That `chematic-smiles` defect is filed separately
// (issue TBD) as discovered-but-deferred, out of scope for this fix.
// ---------------------------------------------------------------------------

/// Canonical rank per atom (`result[idx.0] = canonical position`), used to
/// break every atom-order-sensitive tie in this module's search
/// deterministically, as a function of `mol`'s structure alone.
fn atom_rank(mol: &Molecule) -> Vec<u32> {
    let order = chematic_smiles::canonical_atom_order(mol);
    let mut rank = vec![0u32; order.len()];
    for (pos, &orig) in order.iter().enumerate() {
        rank[orig] = pos as u32;
    }
    rank
}

/// Canonical rank for a tautomer skeleton, independent of explicit H
/// placement.  This is the rooted structural key used by transformations
/// whose only difference is which equivalent aromatic N carries H.
fn tautomer_skeleton_rank(mol: &Molecule) -> Vec<u32> {
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.hydrogen_count = None;
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        builder
            .add_bond(bond.atom1, bond.atom2, bond.order)
            .expect("tautomer skeleton preserves valid bonds");
    }
    atom_rank(&builder.build())
}

/// Return a canonical key for the H-free skeleton rooted at `root`.
/// Isotopically marking the root makes the comparison independent of input
/// atom order even when the unrooted skeleton has automorphisms.
fn rooted_structural_key(mol: &Molecule, root: AtomIdx) -> String {
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
        builder
            .add_bond(bond.atom1, bond.atom2, bond.order)
            .expect("rooted skeleton preserves valid bonds");
    }
    chematic_smiles::canonical_smiles(&builder.build())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Return the canonical (preferred) tautomer of `mol`.
///
/// Applies forward-preferred rules iteratively until no new form is found
/// or the iteration limit is reached; each iteration also tries the
/// directional exocyclic lactam/lactim shift (round 2C-2, RFC section
/// 4.4a) as a fallback when no rule matched. After that loop converges,
/// direct aromatic 1,2-shift tautomers (ring-internal only, e.g.
/// imidazole/pyrazole/tetrazole N-position choice) are compared and ranked
/// by `tautomer_score` (O-H > N-H > S-H, aromatic-ring bonus), with
/// canonical SMILES as the final tiebreaker.
///
/// **A stereocenter's CIP (R/S) label can legitimately change across this
/// shift, even at an atom this function never touches** (issue #402):
/// `stereo_neighbor_order`/`chirality` are always carried forward verbatim
/// for every atom not directly part of the transform (the real spatial
/// configuration never moves), but a CIP label is *computed* from
/// substituent priorities across the whole molecular graph, and a keto-enol
/// (or similar) shift a few bonds away can change those priorities. Compare
/// CIP labels only on the *same* tautomeric form on both sides -- comparing
/// a pre-shift label to a post-shift label for "did stereo survive
/// standardization" is not a valid check, and a difference there does not by
/// itself indicate a bug. Confirmed via independent RDKit CIP oracle on two
/// molecules where this shift changes a label: chematic's CIP on its own
/// shifted structure matches RDKit's CIP on that identical structure exactly
/// in both cases.
///
/// Uses [`TautomerConfig::default`] (all rules, max_iter=16).
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    canonical_tautomer_with_config(mol, &TautomerConfig::default())
}

/// Like [`canonical_tautomer`] but with explicit configuration.
///
/// Every internal candidate/ring tie is broken via [`atom_rank`] (see the
/// "Atom-order canonicalization" section above), computed once from `mol`
/// and reused for every intermediate form the search visits (atom identity
/// never changes across a tautomer-search step), so the result is
/// independent of the caller's parse/atom-insertion order.
pub fn canonical_tautomer_with_config(mol: &Molecule, config: &TautomerConfig) -> Molecule {
    // Preserve E/Z-bearing inputs verbatim until a transform can carry their
    // directional markers onto the new double bond.  Returning the original
    // form is preferable to silently collapsing distinct stereoisomers.
    if has_directional_bond(mol) {
        return mol.clone();
    }
    let rank = atom_rank(mol);

    // canonical_tautomer_search's own final direct-aromatic tie-break is a
    // one-shot step, not iterated to convergence -- confirmed real (not an
    // atom-order artifact) via issue #415's own residual, which needs
    // exactly one extra application to reach a stable form, order-
    // invariantly, once rank-based tie-breaking is used above. Re-run until
    // stable so this function is actually idempotent on its own output, not
    // just closer to it. The `canonical_smiles`-based equality check is
    // skipped entirely when a pass changes nothing (the common case for
    // molecules with no tautomerizable groups at all).
    let (mut result, changed) = canonical_tautomer_search(mol, config, &rank);
    if changed {
        let mut seen_smiles: HashSet<String> = HashSet::new();
        seen_smiles.insert(chematic_smiles::canonical_smiles(&result));
        for _ in 0..3 {
            let (next, next_changed) = canonical_tautomer_search(&result, config, &rank);
            if !next_changed {
                break;
            }
            let is_new = seen_smiles.insert(chematic_smiles::canonical_smiles(&next));
            result = next;
            if !is_new {
                break;
            }
        }
    }
    // The directed search above is intentionally cheap, but its greedy path
    // can make keto spellings of an asymmetric dual-flank ring land in
    // different components. For aromatic O/C tautomer systems, enumerate the
    // bounded connected component and select one representative globally.
    // Prefer the lactam (aromatic C=O) subset whenever it exists; this keeps
    // the established 2-pyridone preference while making keto inputs and
    // hydroxyl inputs share the same canonical result.
    if config.max_tautomers > 0 {
        let forms = enumerate_tautomers_with_config(mol, config);
        if forms.iter().any(has_aromatic_exocyclic_oxygen)
            && forms
                .iter()
                .all(|candidate| aromatic_exocyclic_carbonyl_count(candidate) <= 1)
        {
            let prefer_lactam = forms.iter().any(has_aromatic_exocyclic_carbonyl);
            // The selection is about the aromatic O/C system: a form that
            // also turns an amide elsewhere into an iminol is not a
            // candidate (the rule-based search prefers the amide).
            let max_iminols = iminol_count(&result);
            let mut candidates: Vec<Molecule> = forms
                .into_iter()
                .filter(|candidate| !prefer_lactam || has_aromatic_exocyclic_carbonyl(candidate))
                .filter(|candidate| iminol_count(candidate) <= max_iminols)
                .filter(|candidate| changes_confined_to_aromatic_system(&result, candidate))
                .collect();
            candidates.sort_by(|a, b| {
                tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                    carbonyl_rooted_key(a)
                        .cmp(&carbonyl_rooted_key(b))
                        .then_with(|| {
                            chematic_smiles::canonical_smiles(a)
                                .cmp(&chematic_smiles::canonical_smiles(b))
                        })
                })
            });
            if let Some(best) = candidates.into_iter().next() {
                result = best;
            }
        }
    }
    clear_moved_stereo(mol, &mut result);
    result
}

/// Whether `candidate` differs from `base` only inside aromatic O/C
/// systems: every atom whose H count or bonds changed is aromatic in one of
/// the two forms, or is a heteroatom directly on such an atom. Keeps the
/// global aromatic O/C selection from also enolizing carbonyls or turning
/// amides into enamine/enol forms elsewhere (a tyrosine phenol made the whole
/// peptide backbone a candidate, and the O-H-favouring score picked enols).
fn changes_confined_to_aromatic_system(base: &Molecule, candidate: &Molecule) -> bool {
    if base.atom_count() != candidate.atom_count() {
        return false;
    }
    let bond_orders = |m: &Molecule, idx: AtomIdx| {
        let mut v: Vec<(u32, BondOrder)> = m
            .neighbors(idx)
            .map(|(nb, b)| (nb.0, m.bond(b).order))
            .collect();
        v.sort_by_key(|&(nb, _)| nb);
        v
    };
    let aromatic_any = |idx: AtomIdx| base.atom(idx).aromatic || candidate.atom(idx).aromatic;
    (0..base.atom_count()).all(|i| {
        let idx = AtomIdx(i as u32);
        let changed = chematic_core::implicit_hcount(base, idx)
            != chematic_core::implicit_hcount(candidate, idx)
            || base.atom(idx).hydrogen_count != candidate.atom(idx).hydrogen_count
            || bond_orders(base, idx) != bond_orders(candidate, idx);
        if !changed || aromatic_any(idx) {
            return true;
        }
        let hetero = matches!(base.atom(idx).element.atomic_number(), 7 | 8 | 16);
        hetero && base.neighbors(idx).any(|(nb, _)| aromatic_any(nb))
    })
}

/// Drop tetrahedral tags the tautomer shift invalidated: an atom that gained
/// a multiple bond, changed aromaticity or changed its H count is no longer the
/// stereocentre it was (`N[C@H](C)C=O` → enol `NC(C)=CO`). RDKit's tautomer
/// canonicalizer removes the same stereo. Atoms the shift did not touch keep
/// theirs.
fn clear_moved_stereo(input: &Molecule, result: &mut Molecule) {
    let n = input.atom_count().min(result.atom_count());
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        if result.atom(idx).chirality == Chirality::None {
            continue;
        }
        let multiple = |m: &Molecule| {
            m.neighbors(idx).any(|(_, b)| {
                matches!(
                    m.bond(b).order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                )
            })
        };
        // An allene end keeps its double bond and its tag.
        let moved = (multiple(result) && !multiple(input))
            || result.atom(idx).aromatic != input.atom(idx).aromatic
            || chematic_core::implicit_hcount(result, idx)
                != chematic_core::implicit_hcount(input, idx);
        if moved {
            result.set_chirality(idx, Chirality::None);
        }
    }
}

/// The actual rule-based/exocyclic-lactam/direct-aromatic search. Returns
/// the result plus whether anything actually changed from `mol` -- lets
/// [`canonical_tautomer_with_config`]'s convergence loop skip re-verifying
/// (and the `canonical_smiles` calls that involves) when a pass was a no-op.
fn canonical_tautomer_search(
    mol: &Molecule,
    config: &TautomerConfig,
    rank: &[u32],
) -> (Molecule, bool) {
    let mut current = mol.clone();
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));
    let mut changed = false;
    let mut dual_flank_changed = false;
    let mut carbonyl_rooted = false;

    for _ in 0..config.max_iter {
        let mut iter_changed = false;
        for rule in active_rules(config)
            .into_iter()
            .filter(|r| r.prefer_forward)
        {
            let (next, applied) = apply_nonconflicting_matches(&current, rule, config, rank);
            if applied {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    break;
                }
            }
        }
        if !iter_changed {
            // The 41 rules above never match across an aromatic ring bond
            // (BondOrderMatch::Double never matches BondOrder::Aromatic).
            // This is the round 2C-2 generalization: an exocyclic O donor
            // across an aromatic ring atom -- see the mechanism block above
            // `transfer_hydrogen_exocyclic_lactam`.
            if let Some((next, ..)) = apply_exocyclic_lactam_shift_tracked(&current, config, rank) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    carbonyl_rooted = true;
                }
            }
            if !iter_changed
                && let Some((next, ..)) = apply_dual_flank_tracked(&current, config, rank)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    dual_flank_changed = true;
                }
            }
        }
        if !iter_changed {
            break;
        }
    }

    // The dual-flank step already selected a deterministic rooted structural
    // key.  Do not feed that result back into the broader direct-aromatic
    // score pool, which can otherwise move the proton to a third ring N and
    // undo the cytosine/guanine normalization.
    if dual_flank_changed || carbonyl_rooted {
        return (current, changed);
    }

    // Among direct aromatic 1,2-shift tautomers, pick by tautomer score
    // (O-H > N-H > S-H, aromatic rings), with canonical SMILES as tiebreaker.
    //
    // BFS via enumerate_direct_aromatic_forms ensures all ring positions are
    // explored (e.g. tetrazole has 4 N atoms; a single step would miss N1↔N3).
    // The canonical SMILES tiebreaker makes selection independent of input
    // SMILES write order (the previous h_assignment tiebreaker was not).
    let mut candidates: Vec<Molecule> = vec![current.clone()];
    candidates.extend(enumerate_direct_aromatic_forms(
        &current,
        &config.blocked_atoms,
        16,
        rank,
    ));
    if candidates.len() > 1 {
        candidates.sort_by(|a, b| {
            tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                chematic_smiles::canonical_smiles(a).cmp(&chematic_smiles::canonical_smiles(b))
            })
        });
        current = candidates.into_iter().next().unwrap();
        changed = true;
    }
    (current, changed)
}

/// Enumerate all reachable tautomers of `mol`, capped at 32.
///
/// Returns a `Vec<Molecule>` where the first element is the original molecule.
/// Includes both 1,3-shift (rule-based) and direct aromatic 1,2-shift tautomers.
///
/// Uses [`TautomerConfig::default`] (all rules, max_tautomers=32).
pub fn enumerate_tautomers(mol: &Molecule) -> Vec<Molecule> {
    enumerate_tautomers_with_config(mol, &TautomerConfig::default())
}

/// Like [`enumerate_tautomers`] but with explicit configuration.
pub fn enumerate_tautomers_with_config(mol: &Molecule, config: &TautomerConfig) -> Vec<Molecule> {
    let rank = atom_rank(mol);
    let mut result = vec![mol.clone()];
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(mol));
    // Separate seen-set for 1,2-shift (mol_fingerprint can't distinguish positional H isomers).
    let mut h_seen: HashSet<Vec<Option<u32>>> = HashSet::new();
    h_seen.insert(h_assignment(mol));
    let mut frontier = vec![mol.clone()];

    while !frontier.is_empty() && result.len() < config.max_tautomers {
        let current = frontier.remove(0);
        for rule in active_rules(config).into_iter() {
            for next in apply_all_matches(&current, rule, config, &rank) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                    if result.len() >= config.max_tautomers {
                        break;
                    }
                }
            }
            if result.len() >= config.max_tautomers {
                break;
            }
        }
        // Keep enumeration in lockstep with canonical search for aromatic
        // lactam/lactim systems. These shifts are outside the generic rule
        // table because aromatic ring bonds do not match an explicit Double
        // pattern, but they are still reachable tautomer edges.
        for next in [
            apply_exocyclic_lactam_shift_tracked(&current, config, &rank).map(|(next, ..)| next),
            apply_dual_flank_tracked(&current, config, &rank).map(|(next, ..)| next),
        ] {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) = next {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
        for next in exocyclic_lactam_candidates(&current, &rank, config) {
            if result.len() >= config.max_tautomers {
                break;
            }
            let fp = mol_fingerprint(&next);
            if !seen.contains(&fp) {
                seen.insert(fp);
                h_seen.insert(h_assignment(&next));
                frontier.push(next.clone());
                result.push(next);
            }
        }
        // Enumeration must retain both orientations of the dual-flank edge;
        // the canonical search helper intentionally returns only the rooted
        // preferred direction.
        let skeleton_rank = tautomer_skeleton_rank(&current);
        for (donor, _, acceptor) in find_dual_flank_matches(&current, &skeleton_rank) {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) =
                transfer_hydrogen_dual_flank(&current, donor, acceptor, &config.blocked_atoms)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
        for next in reverse_exocyclic_lactim_candidates(&current, config) {
            if result.len() >= config.max_tautomers {
                break;
            }
            let fp = mol_fingerprint(&next);
            if !seen.contains(&fp) {
                seen.insert(fp);
                h_seen.insert(h_assignment(&next));
                frontier.push(next.clone());
                result.push(next);
            }
        }
        // Direct aromatic 1,2-shift (e.g. pyrazole N1H ↔ N2H).
        for (d, a) in find_direct_aromatic_matches(&current, &rank) {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) = transfer_hydrogen_aromatic(&current, d, a, &config.blocked_atoms) {
                let ha = h_assignment(&next);
                if !h_seen.contains(&ha) {
                    h_seen.insert(ha);
                    seen.insert(mol_fingerprint(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tautomer Parent (ROADMAP.md Phase 2 round 2B) -- see
// docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.
// ---------------------------------------------------------------------------
