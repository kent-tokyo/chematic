//! Tetrahedral CIP assignment from constitutional and auxiliary priorities.
//!
//! After constitutional ranking, Rules 4a–4c compare auxiliary descriptors
//! in the same rooted digraph. Rule 5 can resolve further ties. The rank of
//! an auxiliary center's back-to-root ligand is based on Rule 1a.
//!
//! # Positions come from `stereo_neighbor_order`, ranks come from the new comparator
//!
//! [`crate::digraph::CipDigraph`]'s root children are built by iterating
//! `Molecule::neighbors()` -- raw adjacency order, which reflects bond-*creation*
//! time, not SMILES textual encounter order. That is exactly the wrong order for
//! interpreting a stereocenter's `@`/`@@` marker (its meaning is defined relative to
//! encounter order) -- precisely the bug `d0e726b` fixed in the older engine
//! (`crates/chematic-chem/src/cip.rs`) by switching to
//! `Molecule::stereo_neighbor_order`. This module sources the four substituent
//! *positions* from `stereo_neighbor_order` (mapping `STEREO_H_SENTINEL` to the
//! digraph's `ImplicitHydrogen` child), and only the *priority ranking* of those
//! positions from [`crate::compare::rank_children`]. The parity of the four physical
//! ranks determines whether this order agrees with the parsed chirality tag.
//!
//! # Physical ligands vs. duplicate nodes
//!
//! The digraph root's children are not always exactly the stereocenter's 4 physical
//! neighbors: a multiple bond *at* the stereocenter itself (e.g. a P=N phosphazene
//! center) adds one or more [`CipNodeKind::MultipleBondDuplicate`] siblings alongside the
//! real neighbor -- 5 root children for one double bond, not 4. `stereo_neighbor_order`
//! only ever names the 4 real physical neighbors (never a duplicate), so
//! [`position_node_ids`] already only ever resolves to real `Atom`/`ImplicitHydrogen`
//! nodes. What must NOT happen is treating a duplicate as if it were competing for one of
//! those 4 slots: [`assign_one`] ranks the *entire* root-children set (duplicates
//! included, since a duplicate's presence is real information for ranking a real
//! neighbor's own priority), but only ever treats a tie as unresolvable when it's between
//! two of the 4 *physical* positions -- a duplicate tying with anything doesn't block
//! assignment. A duplicate can occupy a rank slot between two physical positions
//! without changing their relative order or parity.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, BondOrder, Chirality, CipCode, Element, Molecule, STEREO_H_SENTINEL};

use crate::CipError;
use crate::budget::CipBudget;
use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::mancude::{PreparedKekule, prepare_kekule_view};
use crate::node::{CipNodeKind, NodeId};

/// Why a candidate stereocenter got no assignment -- distinct from "assigned but
/// mismatched," so a caller (the corpus report, in particular) can tell "we don't know"
/// from "we got it wrong."
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkipReason {
    /// Not exactly 4 resolvable substituent positions (not a plain tetrahedral center,
    /// or `stereo_neighbor_order` unavailable).
    NotFourSubstituents,
    /// Two or more substituents are mutually tied under Rules 1a/1b/2 alone -- a
    /// genuine tie CIP can't resolve without Rule 3+ (out of scope this milestone), not
    /// a guess.
    Tied,
    /// The underlying digraph or comparator exceeded its budget for this atom.
    BudgetExceeded,
    /// No chemically stable label is available for this stereocenter. A
    /// wildcard/dummy ligand has no specified elemental priority; a more
    /// distant wildcard can also change the label when instantiated.
    OracleUnstable,
    /// The centre carries a tetrahedral stereo tag but has only three explicit
    /// ligands (e.g. a bridgehead amine `[N@@]`, a sulfoxide `[S@](=O)`): its
    /// fourth ligand would be a lone pair, which this engine does not model as a
    /// phantom ligand. Reported instead of silently skipped (issue #634).
    LonePairCenter,
}

/// Result of the experimental tetrahedral-only assignment pass.
#[derive(Debug, Clone, Default)]
pub struct AccurateCipAssignment {
    pub assignments: Vec<(AtomIdx, CipCode)>,
    pub skipped: Vec<(AtomIdx, SkipReason)>,
}

/// Assign R/S to every tetrahedral stereocenter in `mol` that Rules 1a/1b/2 alone can
/// resolve. See module docs for scope and the positions-vs-ranks distinction.
///
/// Computes `mol`'s Kekulé-form clone and [`MancudeContext`] **once**, before the
/// per-atom loop when a tetrahedral annotation exists, and shares both across
/// every stereocenter's digraph -- never
/// recomputed per atom or per subtree expansion (a whole-molecule quantity). If Kekulé
/// form can't be computed at all for `mol` (rare -- e.g. a non-bipartite aromatic system
/// `chematic_core::kekulization::kekulize` can't resolve), falls back to the plain,
/// pre-Milestone-3B-1 digraph path (`CipDigraph::new` on the original aromatic-notation
/// `mol`, no `MancudeContext`) for that molecule rather than failing the whole
/// assignment -- exactly today's behavior for such a molecule, since it never had a
/// MANCUDE-fractional path to lose.
pub fn assign_cip_accurate_experimental(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    let result = assign_cip_without_wildcard_probe(mol, budget)?;
    Ok(guard_wildcard_uncertainty(
        mol,
        budget,
        result,
        assign_cip_without_wildcard_probe,
    ))
}

fn guard_wildcard_uncertainty(
    mol: &Molecule,
    budget: CipBudget,
    mut result: AccurateCipAssignment,
    probe_assignment: impl Fn(&Molecule, CipBudget) -> Result<AccurateCipAssignment, CipCompareError>,
) -> AccurateCipAssignment {
    if result.assignments.is_empty() {
        return result;
    }

    let wildcards: Vec<_> = mol
        .atoms()
        .filter_map(|(idx, atom)| atom.wildcard.then_some(idx))
        .collect();
    if wildcards.is_empty() {
        return result;
    }

    // CIP priorities cannot cross disconnected components. This also keeps
    // a separate wildcard reagent from affecting a fully specified center.
    let mut connected_to_wildcard = vec![false; mol.atom_count()];
    let mut stack = wildcards.clone();
    while let Some(idx) = stack.pop() {
        if connected_to_wildcard[idx.0 as usize] {
            continue;
        }
        connected_to_wildcard[idx.0 as usize] = true;
        stack.extend(mol.neighbors(idx).map(|(neighbor, _)| neighbor));
    }
    if !result
        .assignments
        .iter()
        .any(|(center, _)| connected_to_wildcard[center.0 as usize])
    {
        return result;
    }

    // A wildcard can affect the order of branches several bonds away. Detect
    // a concrete counterexample to an apparently stable label by replacing
    // each wildcard, one at a time, with carbon and fluorine. Carbon alone
    // misses cases where two carbon branches remain tied but a more distant
    // high-priority atom reverses their order. Compare only labels obtained
    // at the same center; an unresolved probe is not a reversal witness.
    // These probes do not establish stability for every possible element.
    const MAX_WILDCARDS: usize = 16;
    let mut unstable = HashSet::new();
    if wildcards.len() > MAX_WILDCARDS {
        // Bound the work for query-like inputs with many placeholders. Their
        // connected centers have unspecified priorities, so abstain rather
        // than publish an untested label.
        unstable.extend(
            result
                .assignments
                .iter()
                .filter(|(center, _)| connected_to_wildcard[center.0 as usize])
                .map(|(center, _)| *center),
        );
    } else {
        for wildcard in wildcards {
            for element in [Element::C, Element::F] {
                let concrete = mol.with_atom_element(wildcard, element);
                let Ok(probe) = probe_assignment(&concrete, budget) else {
                    continue;
                };
                let probe_labels: HashMap<_, _> = probe.assignments.into_iter().collect();
                for &(center, label) in &result.assignments {
                    if probe_labels
                        .get(&center)
                        .is_some_and(|&concrete_label| concrete_label != label)
                    {
                        unstable.insert(center);
                    }
                }
            }
        }
    }
    if !unstable.is_empty() {
        result
            .assignments
            .retain(|(center, _)| !unstable.contains(center));
        result.skipped.extend(
            unstable
                .into_iter()
                .map(|center| (center, SkipReason::OracleUnstable)),
        );
        result.skipped.sort_by_key(|(center, _)| center.0);
    }
    result
}

fn assign_cip_without_wildcard_probe(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    // Every assignment/skipped entry begins with a tetrahedral annotation.
    // With no such atom all three passes are empty, independent of budget.
    if !mol.atoms().any(|(_, atom)| atom.chirality.is_tetrahedral()) {
        return Ok(AccurateCipAssignment::default());
    }
    let kekule = prepare_kekule_view(mol).ok();
    let pass1 = assign_all(mol, budget, kekule.as_ref())?;
    // Refine constitutional ties by Rules 3 and 4 before trying Rule 5.
    // The Rule 5 pass starts from the same earlier rankings on unresolved centres.
    let pass2 = crate::resolver::apply_rule4b_pass(mol, budget, kekule.as_ref(), pass1);
    Ok(apply_rule5_pass(mol, budget, kekule.as_ref(), pass2))
}

/// Identical to [`assign_cip_accurate_experimental`], but never attaches a
/// [`MancudeContext`] -- reproduces exactly the pre-Milestone-3B-1b digraph
/// construction (plain `CipDigraph::new`, aromatic bonds contribute no
/// `MultipleBondDuplicate` nodes). Exists as a stable reference point for regression
/// tooling and tests that need to classify a stereocenter's wrong-vs-tied outcome
/// independent of whatever the live, MANCUDE-aware engine currently does -- see
/// `tests/common/mod.rs::is_bucket_misclassified`'s module docs for why that
/// independence matters (gating a *structural corpus scope* on the live engine's current
/// correctness makes the scope shrink every time the engine improves).
pub fn assign_cip_accurate_experimental_without_mancude(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    assign_all(mol, budget, None)
}

/// Whether the accurate label of tetrahedral centre `atom` can change
/// between neutral Kekulé spellings of one molecule: a phosphorus on a ring
/// (of at most eight atoms) carrying a double bond, such as a
/// cyclophosphazene P. RDKit's CIPLabeler flips such a label the same way
/// (220/220 labels agree over both spellings and 10 atom orders each), so
/// the label is reported, but it is not a property of the molecule.
pub fn label_depends_on_kekule_spelling(mol: &Molecule, atom: AtomIdx) -> bool {
    mol.atom(atom).element == Element::P && phosphorus_on_unsaturated_ring(mol, atom)
}

fn phosphorus_on_unsaturated_ring(mol: &Molecule, p: AtomIdx) -> bool {
    use chematic_core::BondOrder;
    // Breadth-first from each neighbour back to `p` without the first bond,
    // remembering whether the path used a double or aromatic bond.
    for (start, first) in mol.neighbors(p) {
        let first_unsat = matches!(
            mol.bond(first).order,
            BondOrder::Double | BondOrder::Aromatic
        );
        let mut seen = vec![false; mol.atom_count()];
        seen[start.0 as usize] = true;
        let mut frontier = vec![(start, first_unsat)];
        for _ in 0..7 {
            let mut next = Vec::new();
            for &(a, unsat) in &frontier {
                for (nb, b) in mol.neighbors(a) {
                    if b == first {
                        continue;
                    }
                    let u = unsat
                        || matches!(mol.bond(b).order, BondOrder::Double | BondOrder::Aromatic);
                    if nb == p {
                        if u {
                            return true;
                        }
                        continue;
                    }
                    if !seen[nb.0 as usize] {
                        seen[nb.0 as usize] = true;
                        next.push((nb, u));
                    }
                }
            }
            frontier = next;
        }
    }
    false
}

fn assign_all(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
) -> Result<AccurateCipAssignment, CipCompareError> {
    let mut result = AccurateCipAssignment::default();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        if !atom.chirality.is_tetrahedral() {
            continue;
        }

        // A P centre on a ring that carries a double bond (cyclophosphazene
        // P=N rings) is labelled as RDKit's CIPLabeler labels it. Both engines
        // flip that label between the ring's neutral Kekulé spellings, which
        // InChI reads as one molecule, so identity checks treat it as
        // spelling-dependent ([`label_depends_on_kekule_spelling`]).
        let Some(raw_stereo_order) = mol.stereo_neighbor_order(idx) else {
            result.skipped.push((idx, SkipReason::NotFourSubstituents));
            continue;
        };
        let stereo_order = rdkit_cip_position_order(mol, idx, raw_stereo_order);
        if stereo_order
            .iter()
            .any(|&position| position != STEREO_H_SENTINEL && mol.atom(AtomIdx(position)).wildcard)
        {
            result.skipped.push((idx, SkipReason::OracleUnstable));
            continue;
        }
        if stereo_order.len() == 3 {
            match assign_lone_pair_centre(mol, idx, atom.chirality, &stereo_order, budget, kekule) {
                Ok(Some(code)) => result.assignments.push((idx, code)),
                Ok(None) => result.skipped.push((idx, SkipReason::LonePairCenter)),
                Err(reason) => result.skipped.push((idx, reason)),
            }
            continue;
        }
        if stereo_order.len() != 4 {
            result.skipped.push((idx, SkipReason::NotFourSubstituents));
            continue;
        }

        match assign_one(mol, idx, atom.chirality, &stereo_order, budget, kekule) {
            Ok(Some(code)) => result.assignments.push((idx, code)),
            Ok(None) => result.skipped.push((idx, SkipReason::NotFourSubstituents)),
            Err(SkipReason::Tied) => result.skipped.push((idx, SkipReason::Tied)),
            Err(SkipReason::BudgetExceeded) => {
                result.skipped.push((idx, SkipReason::BudgetExceeded))
            }
            Err(other) => result.skipped.push((idx, other)),
        }
    }

    Ok(result)
}

/// Rule 5 refinement -- see module docs for scope. Only ever touches atoms Pass 1
/// (`assign_all`, above) and the Rule 4b pass (`crate::resolver::apply_rule4b_pass`) left
/// [`SkipReason::Tied`]; every other atom (resolved or skipped for any other reason) is
/// carried through unchanged. Not called by
/// [`assign_cip_accurate_experimental_without_mancude`] -- that function's whole purpose
/// is to stay a frozen, Rule-5-independent reference point (see its own doc comment), so
/// Rule 5 is deliberately only wired into the main, live entry point.
fn apply_rule5_pass(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
    pass1: AccurateCipAssignment,
) -> AccurateCipAssignment {
    let mut assignments = pass1.assignments;
    let mut skipped = Vec::with_capacity(pass1.skipped.len());

    for (idx, reason) in pass1.skipped {
        if reason != SkipReason::Tied {
            skipped.push((idx, reason));
            continue;
        }

        let atom = mol.atom(idx);
        let Some(raw_stereo_order) = mol.stereo_neighbor_order(idx) else {
            skipped.push((idx, reason));
            continue;
        };
        let stereo_order = rdkit_cip_position_order(mol, idx, raw_stereo_order);
        if stereo_order.len() != 4 && stereo_order.len() != 3 {
            skipped.push((idx, reason));
            continue;
        }

        match assign_one_with_rule5(mol, idx, atom.chirality, &stereo_order, budget, kekule) {
            Ok(Some(code)) => assignments.push((idx, code)),
            Err(SkipReason::BudgetExceeded) => skipped.push((idx, SkipReason::BudgetExceeded)),
            _ => skipped.push((idx, reason)),
        }
    }

    AccurateCipAssignment {
        assignments,
        skipped,
    }
}

/// Physical substituent positions used by RDKit's CIP labeler after
/// organometallic cleanup.  Cleanup can turn one bond at a nominally
/// five-coordinate tagged carbon into a dative bond.  RDKit excludes that
/// coordinate bond from the tetrahedral positions while retaining the
/// original stereo tag and its other four positions.
fn rdkit_cip_position_order<'a>(mol: &Molecule, center: AtomIdx, raw: &'a [u32]) -> Cow<'a, [u32]> {
    if raw.len() <= 4 {
        return Cow::Borrowed(raw);
    }
    let mut filtered: Vec<u32> = raw
        .iter()
        .copied()
        .filter(|&position| {
            position == STEREO_H_SENTINEL
                || mol
                    .bond_between(center, AtomIdx(position))
                    .is_none_or(|(_, bond)| bond.order != BondOrder::Dative)
        })
        .collect();
    if filtered.len() == 4 {
        // The SMILES tag was normalized while the original single bond was
        // still one of five encountered entries. RDKit's cleanup removes the
        // new dative ligand and applies the corresponding odd permutation.
        // Keep that parity correction explicit in the position order.
        filtered.swap(0, 1);
        Cow::Owned(filtered)
    } else {
        Cow::Borrowed(raw)
    }
}

/// Resolve a remaining tie with the auxiliary rules through Rule 5.
fn assign_one_with_rule5(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
) -> Result<Option<CipCode>, SkipReason> {
    // Unsupported unsaturated phosphazene auxiliaries can depend on Kekulé spelling.
    let lone_pair_centre = stereo_order.len() == 3;
    if mol.atom(idx).element == Element::P && !is_saturated_ring_lone_pair_centre(mol, idx) {
        return Err(SkipReason::Tied);
    }

    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).map_err(map_digraph_err)?
        }
        None => CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?,
    };
    let root = graph.root();
    let root_children = graph.expand_child_ids(root).map_err(map_digraph_err)?;

    let Some(position_nodes) = position_node_ids(&graph, stereo_order, &root_children) else {
        return Ok(None);
    };

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;

    let groups = crate::rule3::refine(&mut graph, groups, &mut ctx).map_err(map_compare_err)?;
    let mut cache = HashMap::new();
    let (resolved_groups, pseudo_count) =
        crate::rule4b::refine_stereo_groups(&mut graph, mol, groups, &position_nodes, &mut cache)
            .map_err(map_compare_err)?;
    if !lone_pair_centre
        && crate::rule6::has_reference_partition(&resolved_groups)
        && crate::rule6::eligible(mol, idx)
    {
        return crate::rule6::assign(
            &mut graph,
            &resolved_groups,
            &position_nodes,
            chirality,
            &mut ctx,
            &mut cache,
            pseudo_count % 2 == 1,
        )
        .map_err(map_compare_err)?
        .map(Some)
        .ok_or(SkipReason::Tied);
    }
    if resolved_groups.iter().any(|group| {
        group
            .iter()
            .filter(|node| position_nodes.contains(node))
            .count()
            > 1
    }) {
        return Err(SkipReason::Tied);
    }
    let pseudoasymmetric = pseudo_count % 2 == 1;

    let is_r = if lone_pair_centre {
        lone_pair_is_r(&resolved_groups, &position_nodes, chirality)?
    } else {
        resolve_is_r_from_groups(&resolved_groups, &position_nodes, chirality)
    };
    let Some(is_r) = is_r else {
        return Ok(None);
    };
    Ok(Some(match (is_r, pseudoasymmetric) {
        (true, true) => CipCode::LowerR,
        (false, true) => CipCode::LowerS,
        (true, false) => CipCode::R,
        (false, false) => CipCode::S,
    }))
}

fn assign_one(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
) -> Result<Option<CipCode>, SkipReason> {
    // `apply_kekule` preserves `AtomIdx` values exactly (verified Milestone 3B-0), so
    // `idx`/`stereo_order` (sourced from the original `mol`, above) name the same
    // physical atoms in `kekule_mol` -- no remapping needed either way.
    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).map_err(map_digraph_err)?
        }
        None => CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?,
    };
    let root = graph.root();
    let root_children = graph.expand_child_ids(root).map_err(map_digraph_err)?;

    let Some(position_nodes) = position_node_ids(&graph, stereo_order, &root_children) else {
        return Ok(None);
    };

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;

    // A tie only blocks resolution when two of the 4 *physical* positions land in the
    // same group -- a duplicate node tying with anything (another duplicate, or even a
    // physical position) doesn't compete for a stereo_neighbor_order slot. See module
    // docs ("Physical ligands vs. duplicate nodes").
    for group in &groups {
        if group.iter().filter(|n| position_nodes.contains(n)).count() > 1 {
            return Err(SkipReason::Tied);
        }
    }

    let Some(is_r) = resolve_is_r_from_groups(&groups, &position_nodes, chirality) else {
        return Ok(None);
    };
    Ok(Some(if is_r { CipCode::R } else { CipCode::S }))
}

/// A three-ligand centre uses the recorded ligand order and a lowest-priority
/// phantom ligand. Saturated ring centres and ring centres with a double bond
/// are supported. An aromatic ring centre or a centre with an H ligand
/// returns `None` (reported as `LonePairCenter`).
fn assign_lone_pair_centre(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
) -> Result<Option<CipCode>, SkipReason> {
    if stereo_order.len() != 3
        || stereo_order.contains(&STEREO_H_SENTINEL)
        || chematic_core::implicit_hcount(mol, idx) > 0
        || stereo_order
            .iter()
            .any(|&a| mol.atom(AtomIdx(a)).element == Element::H)
    {
        return Ok(None);
    }
    let saturated_ring_centre = is_saturated_ring_lone_pair_centre(mol, idx);
    if atom_in_ring(mol, idx)
        && !saturated_ring_centre
        && !mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == chematic_core::BondOrder::Double)
    {
        // An aromatic ring centre.
        return Ok(None);
    }
    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).map_err(map_digraph_err)?
        }
        None => CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?,
    };
    let root = graph.root();
    let root_children = graph.expand_child_ids(root).map_err(map_digraph_err)?;
    let Some(position_nodes) = position_node_ids(&graph, stereo_order, &root_children) else {
        return Ok(None);
    };
    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;
    for group in &groups {
        if group.iter().filter(|n| position_nodes.contains(n)).count() > 1 {
            return Err(SkipReason::Tied);
        }
    }
    Ok(lone_pair_is_r(&groups, &position_nodes, chirality)?
        .map(|r| if r { CipCode::R } else { CipCode::S }))
}

/// Whether a three-ligand centre with a lone pair is R, given its ligands'
/// priority groups (highest first, each physical ligand in its own group).
/// Reader-specific parity adjustments are already encoded in the graph's
/// chirality and ligand order.
pub(crate) fn lone_pair_is_r(
    groups: &[Vec<NodeId>],
    position_nodes: &[NodeId],
    chirality: Chirality,
) -> Result<Option<bool>, SkipReason> {
    // Highest group first: rank the three ligands 4, 3, 2; the lone pair is 1.
    let group_of = |node: NodeId| groups.iter().position(|g| g.contains(&node));
    let mut by_priority: Vec<(usize, usize)> = position_nodes
        .iter()
        .enumerate()
        .map(|(slot, &node)| group_of(node).map(|g| (g, slot)))
        .collect::<Option<_>>()
        .ok_or(SkipReason::NotFourSubstituents)?;
    by_priority.sort_unstable();
    let mut ranks3 = [0u8; 3];
    for (k, &(_, slot)) in by_priority.iter().enumerate() {
        ranks3[slot] = 4 - k as u8;
    }
    // Matched against RDKit 2026.03.6 `rdCIPLabeler` on acyclic sulfoxides,
    // sulfinamides and phosphines over random atom orders: the phantom takes
    // the second slot whether or not the centre has a preceding atom.
    let ranks = [ranks3[0], 1, ranks3[1], ranks3[2]];
    Ok(is_r_from_ranks(&ranks, chirality))
}

/// Whether a tagged ring centre has three single-bonded ligands and a
/// lone pair. RDKit's SMILES reader applies a special parity convention
/// to some spellings of these centres.
pub fn is_saturated_ring_lone_pair_centre(mol: &Molecule, atom: AtomIdx) -> bool {
    use chematic_core::BondOrder;
    mol.atom(atom).chirality.is_tetrahedral()
        && mol.degree(atom) == 3
        && chematic_core::implicit_hcount(mol, atom) == 0
        && mol.neighbors(atom).all(|(nb, b)| {
            mol.atom(nb).element != Element::H
                && matches!(
                    mol.bond(b).order,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down
                )
        })
        && atom_in_ring(mol, atom)
}

/// Whether `idx` lies on a ring: some neighbour reaches another without
/// passing through `idx`.
fn atom_in_ring(mol: &Molecule, idx: AtomIdx) -> bool {
    let nbrs: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();
    for (k, &start) in nbrs.iter().enumerate() {
        let mut seen = vec![false; mol.atom_count()];
        seen[idx.0 as usize] = true;
        seen[start.0 as usize] = true;
        let mut stack = vec![start];
        while let Some(a) = stack.pop() {
            for (nb, _) in mol.neighbors(a) {
                if nbrs[k + 1..].contains(&nb) {
                    return true;
                }
                if !seen[nb.0 as usize] {
                    seen[nb.0 as usize] = true;
                    stack.push(nb);
                }
            }
        }
    }
    false
}

/// Given a fully-ranked group partition (highest-priority group first, matching
/// [`rank_children`]'s convention) and the 4 physical positions in `stereo_order`
/// order, compute whether the center is *rectus* (R). Shared by [`assign_one`] (Rules
/// 1a/1b/2 only), [`assign_one_with_rule5`] (which passes a `groups` partition with
/// one tied group already split by a Rule 5 tiebreak), and (Milestone 4B-2)
/// `crate::resolver::assign_one_with_rule4b`/`resolve_chirality`. The same parity
/// calculation applies regardless of the rules used to order the ligands.
pub(crate) fn resolve_is_r_from_groups(
    groups: &[Vec<NodeId>],
    position_nodes: &[NodeId],
    chirality: Chirality,
) -> Option<bool> {
    let positions: &[NodeId; 4] = position_nodes.try_into().ok()?;
    let n = u8::try_from(groups.len()).ok()?;
    let mut ranks = [0; 4];
    for (node, rank) in positions.iter().zip(&mut ranks) {
        let group = groups.iter().rposition(|group| group.contains(node))?;
        *rank = n - group as u8;
    }
    is_r_from_ranks(&ranks, chirality)
}

/// R/S from four position ranks in stereo-neighbor order, high number first in priority.
fn is_r_from_ranks(ranks: &[u8; 4], chirality: Chirality) -> Option<bool> {
    let mut odd = false;
    for (i, rank) in ranks.iter().enumerate() {
        for next in &ranks[i + 1..] {
            if rank == next {
                return None;
            }
            odd ^= rank > next;
        }
    }
    Some((chirality == Chirality::Clockwise) ^ odd)
}

pub(crate) fn map_digraph_err(e: CipError) -> SkipReason {
    let CipError::BudgetExceeded { .. } = e;
    SkipReason::BudgetExceeded
}

pub(crate) fn map_compare_err(e: CipCompareError) -> SkipReason {
    match e {
        CipCompareError::BudgetExceeded { .. } | CipCompareError::Digraph(_) => {
            SkipReason::BudgetExceeded
        }
        CipCompareError::InvalidDigraph(_) => SkipReason::NotFourSubstituents,
    }
}

/// Map each `stereo_neighbor_order` position to the digraph node representing it.
///
/// For a tetrahedral stereocenter's own root children (the only case before Milestone
/// 4B-2), substituents are always single-bonded, so `candidates` only ever contains
/// `Atom`/`ImplicitHydrogen` kinds -- never a duplicate. Milestone 4B-2's Rule 4b
/// resolver additionally calls this for an *embedded* stereocenter's own 4 positions
/// (its forward children plus the back-to-root parent node,
/// `crate::resolver::resolve_chirality`), where a physical position can be a
/// `RingDuplicate` (the stereo_neighbor_order-named atom is reached via a ring
/// closure at that embedded position, not as a fresh real node) -- hence the
/// additional match arm below. Provably inert for every pre-existing caller
/// ([`assign_one`]/[`assign_one_with_rule5`]), which only ever pass a root's own
/// direct children.
pub(crate) fn position_node_ids(
    graph: &CipDigraph,
    stereo_order: &[u32],
    candidates: &[NodeId],
) -> Option<Vec<NodeId>> {
    let mut result = Vec::with_capacity(stereo_order.len());
    for &pos_val in stereo_order {
        let node_id = if pos_val == STEREO_H_SENTINEL {
            candidates
                .iter()
                .copied()
                .find(|&id| matches!(graph.node(id).kind, CipNodeKind::ImplicitHydrogen))?
        } else {
            let atom_idx = AtomIdx(pos_val);
            candidates
                .iter()
                .copied()
                .find(|&id| match graph.node(id).kind {
                    CipNodeKind::Atom { atom_idx: a } => a == atom_idx,
                    CipNodeKind::RingDuplicate { closure_atom, .. } => closure_atom == atom_idx,
                    _ => false,
                })?
        };
        result.push(node_id);
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::resolve_is_r_from_groups;
    use chematic_core::Chirality;

    use crate::node::NodeId;

    #[test]
    fn resolve_is_r_declines_malformed_position_state() {
        let four_positions = [NodeId(1), NodeId(2), NodeId(3), NodeId(4)];
        let complete_groups = four_positions
            .iter()
            .map(|&node| vec![node])
            .collect::<Vec<_>>();

        assert!(
            resolve_is_r_from_groups(&complete_groups, &four_positions, Chirality::Clockwise,)
                .is_some()
        );
        assert_eq!(
            resolve_is_r_from_groups(&complete_groups, &four_positions[..3], Chirality::Clockwise,),
            None
        );
        assert_eq!(
            resolve_is_r_from_groups(
                &complete_groups,
                &[NodeId(1), NodeId(2), NodeId(3), NodeId(99)],
                Chirality::Clockwise,
            ),
            None
        );
        assert_eq!(
            resolve_is_r_from_groups(&[], &four_positions, Chirality::Clockwise),
            None
        );
    }

    #[test]
    fn nonphysical_nodes_do_not_change_tetrahedral_parity() {
        let order = [NodeId(1), NodeId(2), NodeId(3), NodeId(4)];
        let groups = [
            vec![order[0], NodeId(15)],
            vec![order[1]],
            vec![NodeId(90)],
            vec![order[2]],
            vec![order[3]],
        ];
        assert_eq!(
            resolve_is_r_from_groups(&groups, &order, Chirality::Clockwise),
            Some(true)
        );
        assert_eq!(
            resolve_is_r_from_groups(
                &groups,
                &[order[1], order[0], order[2], order[3]],
                Chirality::Clockwise,
            ),
            Some(false)
        );
    }
}

#[cfg(test)]
#[path = "assignment_preparation_tests.rs"]
mod preparation_tests;
