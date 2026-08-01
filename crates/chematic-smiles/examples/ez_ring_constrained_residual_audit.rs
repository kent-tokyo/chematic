//! Wave 2C (issue #149) audit-only follow-up to Wave 2B (PR #229): tests the
//! hypothesis documented right above `EZ_SHARED_CARRIER_RING_CONSTRAINED_
//! RESIDUALS` in `crates/chematic-smiles/src/canonical.rs` -- that every one
//! of the 8 still-residual coupled components involves an alkene end whose
//! own C=X double bond is *endocyclic* in a 5- or 6-membered ring (fixed
//! cis/trans by ring topology, not a free choice), which
//! `compute_stereo_alkene_ends` has no gate for.
//!
//! Read-only: reuses `ez_shared_carrier_component_audit.rs`'s
//! `stereo_alkene_end_nodes`/`substituents`/coupling-graph logic verbatim as
//! its starting point, using only public `chematic_core`/`chematic_smiles`/
//! `chematic_perception` APIs -- no production code touched.
//!
//! For every stereo-alkene end in the corpus (coupled AND singleton -- the
//! blast-radius measurement downstream needs the full population, not just
//! the 31 coupled atoms), emits one JSON object per line to stdout with:
//!   - the input SMILES and the end atom's index/element
//!   - its C=X double-bond partner's index/element
//!   - coupling-component membership (singleton vs coupled, other members)
//!   - ring membership of the end atom itself (any SSSR ring, and sizes)
//!   - **endocyclic vs exocyclic**: does the C=X double bond itself lie
//!     within a single SSSR ring (both atoms share ring membership) -- the
//!     specific predicate the hypothesis is about, not just "atom is in a
//!     ring somewhere" (an end atom can sit in a small ring while its own
//!     double bond is still exocyclic and genuinely free -- see atom 16 in
//!     the worked fixture-1 example in the companion doc)
//!   - whether canonicalization actually placed a directional (`Up`/`Down`)
//!     marker on this end's own candidate bond(s) in the winning canonical
//!     output
//!
//! ## How `marker_placed` is determined without touching private API
//!
//! `canonical.rs`'s `ez_marker` map is `pub(crate)`, unreachable from an
//! example. Instead: `chematic_smiles::canonical_atom_order` is public and
//! is a pure, relabeling-invariant function of molecule structure (same
//! `winning_individualized_ranks` `canonical_smiles` itself uses) -- for the
//! SAME physical atom, `canonical_atom_order` assigns the SAME canonical
//! *position* whether called on the original molecule or on a reparse of
//! its own canonical-SMILES output, because `initial_invariant` collapses
//! `Single`/`Up`/`Down`/`Dative` to the identical bond-order-class `1`
//! (canonical.rs line ~451) -- rank computation does not see *which*
//! candidate bond carries a mark, only that a bond exists. So: compute
//! `canonical_atom_order` on both the original molecule and a reparse of its
//! own canonical output, use each as an atom-index -> canonical-position
//! map, and compose position(original) -> position(reparsed) to recover
//! atom correspondence -- all through public API, no molecule mutation
//! needed (an earlier design tried tagging atoms via `atom_map` for
//! correspondence and was rejected: `canonical_partition.rs` folds
//! `atom_map` into the initial invariant, so tagging would have perturbed
//! the very ranks/marker-choice this audit is trying to observe).
//!
//! This composition is NOT assumed safe -- it is verified per molecule by
//! replaying every one of the original molecule's bonds through the mapping
//! and checking the reparsed molecule has the same bond (same atom pair,
//! same bond-order *class*, treating `Single`/`Up`/`Down` as one class)
//! at the mapped positions. `correspondence_ok` is recorded on every
//! emitted row; when `false` (rare -- only possible if a genuine rank tie
//! among automorphic atoms picks a different winning branch between the
//! original and the reparse), `marker_placed`/`marker_placed_count` are
//! `null` rather than a guessed value.
//!
//! Run against the 18 pinned fixtures (no argument) or a full corpus (one
//! SMILES per line -- `scripts/descriptor_census_corpus.smi`, committed to
//! this repo, no external download needed):
//!
//! ```text
//! cargo run -p chematic-smiles --release --example ez_ring_constrained_residual_audit
//! cargo run -p chematic-smiles --release --example ez_ring_constrained_residual_audit -- scripts/descriptor_census_corpus.smi
//! ```

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

/// The 10 now-fully-resolved fixtures (`EZ_SHARED_CARRIER_FULLY_RESOLVED` in
/// `canonical.rs`) -- included for contrast, should NOT trip the
/// ring-endocyclic hypothesis.
const EZ_SHARED_CARRIER_FULLY_RESOLVED: &[&str] = &[
    r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
];

/// The 8 still-residual fixtures (`EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS`).
const EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS: &[&str] = &[
    r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
    r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
];

// ---------------------------------------------------------------------------
// Reused verbatim (topological reimplementation, public-API only) from
// `ez_shared_carrier_component_audit.rs` -- mirrors `canonical.rs`'s private
// `compute_stereo_alkene_ends`/`substituents` exactly.
// ---------------------------------------------------------------------------

fn end_has_substituent(mol: &Molecule, end: AtomIdx) -> bool {
    mol.neighbors(end)
        .any(|(_, b)| mol.bond(b).order != BondOrder::Double)
}

fn substituents(mol: &Molecule, end: AtomIdx) -> Vec<(AtomIdx, BondIdx)> {
    mol.neighbors(end)
        .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
        .collect()
}

fn stereo_alkene_end_nodes(mol: &Molecule) -> HashSet<AtomIdx> {
    let mut nodes = HashSet::new();
    for (_, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        if !end_has_substituent(mol, bond.atom1) || !end_has_substituent(mol, bond.atom2) {
            continue;
        }
        for end in [bond.atom1, bond.atom2] {
            if substituents(mol, end).len() == 2 {
                nodes.insert(end);
            }
        }
    }
    nodes
}

/// The double-bond partner of a stereo-alkene end atom (the other end of the
/// SAME double bond -- `end`'s only `Double`-order neighbor, guaranteed
/// unique since `end` is itself a double-bond terminus).
fn double_bond_partner(mol: &Molecule, end: AtomIdx) -> AtomIdx {
    mol.neighbors(end)
        .find(|&(_, b)| mol.bond(b).order == BondOrder::Double)
        .map(|(nb, _)| nb)
        .expect("stereo-alkene end atom must have exactly one Double-order neighbor")
}

/// Full connected-component partition of the shared-candidate-bond coupling
/// graph (nodes = `ends`, edges = a direct non-double bond between two
/// ends), INCLUDING singleton (size-1) components -- unlike the size>=2-only
/// filter in `ez_shared_carrier_component_audit.rs`'s `coupling_components`,
/// this audit needs every end's membership, singleton or not. Returns, per
/// end atom, the sorted list of every OTHER atom in its component (empty for
/// a singleton).
fn coupling_component_other_members(
    mol: &Molecule,
    ends: &HashSet<AtomIdx>,
) -> HashMap<AtomIdx, Vec<AtomIdx>> {
    let mut adjacency: HashMap<AtomIdx, Vec<AtomIdx>> = HashMap::new();
    for &n in ends {
        for (sub_atom, _) in substituents(mol, n) {
            if ends.contains(&sub_atom) {
                adjacency.entry(n).or_default().push(sub_atom);
            }
        }
    }

    let mut starts: Vec<AtomIdx> = ends.iter().copied().collect();
    starts.sort_by_key(|a| a.0);
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut result: HashMap<AtomIdx, Vec<AtomIdx>> = HashMap::new();
    for start in starts {
        if visited.contains(&start) {
            continue;
        }
        let mut queue = vec![start];
        let mut component = Vec::new();
        visited.insert(start);
        while let Some(cur) = queue.pop() {
            component.push(cur);
            let mut nbs: Vec<AtomIdx> = adjacency.get(&cur).cloned().unwrap_or_default();
            nbs.sort_by_key(|a| a.0);
            for nb in nbs {
                if visited.insert(nb) {
                    queue.push(nb);
                }
            }
        }
        component.sort_by_key(|a| a.0);
        for &member in &component {
            let others: Vec<AtomIdx> = component.iter().copied().filter(|&m| m != member).collect();
            result.insert(member, others);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Ring predicates (SSSR-based, `chematic_perception::find_sssr`).
// ---------------------------------------------------------------------------

fn ring_sizes_containing(rings: &[Vec<AtomIdx>], atom: AtomIdx) -> Vec<usize> {
    let mut sizes: Vec<usize> = rings
        .iter()
        .filter(|r| r.contains(&atom))
        .map(|r| r.len())
        .collect();
    sizes.sort_unstable();
    sizes
}

/// Sizes of every SSSR ring containing BOTH `a` and `b` -- the key predicate
/// distinguishing an *endocyclic* double bond (both termini share a ring)
/// from an atom that merely happens to sit in a ring while its own double
/// bond points outward (exocyclic).
fn endocyclic_ring_sizes(rings: &[Vec<AtomIdx>], a: AtomIdx, b: AtomIdx) -> Vec<usize> {
    let mut sizes: Vec<usize> = rings
        .iter()
        .filter(|r| r.contains(&a) && r.contains(&b))
        .map(|r| r.len())
        .collect();
    sizes.sort_unstable();
    sizes
}

// ---------------------------------------------------------------------------
// Atom correspondence across a canonicalize -> reparse round trip, via
// composition of `canonical_atom_order` (public, relabeling-invariant) --
// see the module doc comment for why this is safe and how it is verified.
// ---------------------------------------------------------------------------

/// Invert `canonical_atom_order`'s `Vec<usize>` (position -> atom index) into
/// atom index -> position.
fn position_of(order: &[usize]) -> Vec<usize> {
    let mut pos = vec![0usize; order.len()];
    for (position, &atom_idx) in order.iter().enumerate() {
        pos[atom_idx] = position;
    }
    pos
}

/// `Single`/`Up`/`Down` are the same bond CLASS for correspondence purposes
/// (the only axis this audit expects to legitimately differ between `mol`
/// and its own canonical-output reparse); every other `BondOrder` variant
/// must match exactly.
fn bond_order_class_eq(a: BondOrder, b: BondOrder) -> bool {
    let plain = |o: BondOrder| matches!(o, BondOrder::Single | BondOrder::Up | BondOrder::Down);
    if plain(a) && plain(b) {
        return true;
    }
    a == b
}

/// Verifies the `pos1 -> order2` composition actually reproduces `mol`'s
/// full bond structure inside `mol2` (see module doc comment). Returns
/// `false` (rather than a guessed mapping) if any bond fails to correspond --
/// e.g. a genuine automorphism-orbit rank tie that broke differently between
/// `mol` and the reparse of its own canonical output.
fn verify_correspondence(
    mol: &Molecule,
    pos1: &[usize],
    mol2: &Molecule,
    order2: &[usize],
) -> bool {
    if mol.atom_count() != mol2.atom_count() {
        return false;
    }
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1.0 as usize;
        let a2 = bond.atom2.0 as usize;
        if a1 >= pos1.len() || a2 >= pos1.len() {
            return false;
        }
        let p1 = pos1[a1];
        let p2 = pos1[a2];
        if p1 >= order2.len() || p2 >= order2.len() {
            return false;
        }
        let mapped_a = AtomIdx(order2[p1] as u32);
        let mapped_b = AtomIdx(order2[p2] as u32);
        match mol2.bond_between(mapped_a, mapped_b) {
            Some((_, mapped_bond)) if bond_order_class_eq(bond.order, mapped_bond.order) => {}
            _ => return false,
        }
    }
    true
}

/// Maps original atom `idx` to its counterpart in `mol2`, via
/// `pos1[idx] -> order2[pos1[idx]]`. Caller must have already confirmed
/// `verify_correspondence` for the pair.
fn map_atom(idx: AtomIdx, pos1: &[usize], order2: &[usize]) -> AtomIdx {
    AtomIdx(order2[pos1[idx.0 as usize]] as u32)
}

// ---------------------------------------------------------------------------
// JSON emission (hand-rolled -- no serde_json dependency, matching
// `kekulize_corpus_scan.rs`'s existing convention in this crate).
// ---------------------------------------------------------------------------

fn esc(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn usize_vec_json(v: &[usize]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn atomidx_vec_json(v: &[AtomIdx]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|a| a.0.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn opt_bool_json(v: Option<bool>) -> &'static str {
    match v {
        Some(true) => "true",
        Some(false) => "false",
        None => "null",
    }
}

fn opt_usize_json(v: Option<usize>) -> String {
    match v {
        Some(n) => n.to_string(),
        None => "null".to_string(),
    }
}

// ---------------------------------------------------------------------------
// Main per-molecule analysis.
// ---------------------------------------------------------------------------

fn analyze_molecule(smiles: &str) {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("parse failure, skipped: {smiles:?}: {e}");
            return;
        }
    };

    let ends = stereo_alkene_end_nodes(&mol);
    if ends.is_empty() {
        return;
    }

    let component_others = coupling_component_other_members(&mol, &ends);
    let rings = find_sssr(&mol).rings().to_vec();

    let canon = chematic_smiles::canonical_smiles(&mol);
    let order1 = chematic_smiles::canonical_atom_order(&mol);
    let pos1 = position_of(&order1);

    let (correspondence_ok, mol2_and_order2) = match chematic_smiles::parse(&canon) {
        Ok(mol2) => {
            let order2 = chematic_smiles::canonical_atom_order(&mol2);
            let ok = verify_correspondence(&mol, &pos1, &mol2, &order2);
            (ok, Some((mol2, order2)))
        }
        Err(_) => (false, None),
    };

    let mut sorted_ends: Vec<AtomIdx> = ends.iter().copied().collect();
    sorted_ends.sort_by_key(|a| a.0);

    for end in sorted_ends {
        let element = mol.atom(end).element.symbol();
        let partner = double_bond_partner(&mol, end);
        let partner_element = mol.atom(partner).element.symbol();

        let others = component_others.get(&end).cloned().unwrap_or_default();
        let coupled = !others.is_empty();

        let end_ring_sizes = ring_sizes_containing(&rings, end);
        let endo_sizes = endocyclic_ring_sizes(&rings, end, partner);
        let endocyclic = !endo_sizes.is_empty();

        let subs = substituents(&mol, end);

        let (marker_placed, marker_placed_count) = if correspondence_ok {
            if let Some((mol2, order2)) = &mol2_and_order2 {
                let mapped_end = map_atom(end, &pos1, order2);
                let mut marked = 0usize;
                for &(sub_atom, _) in &subs {
                    let mapped_sub = map_atom(sub_atom, &pos1, order2);
                    if let Some((_, b)) = mol2.bond_between(mapped_end, mapped_sub)
                        && matches!(b.order, BondOrder::Up | BondOrder::Down)
                    {
                        marked += 1;
                    }
                }
                (Some(marked > 0), Some(marked))
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };

        println!(
            "{{\"smiles\":\"{smiles}\",\"canonical_smiles\":\"{canon}\",\
             \"end_atom_idx\":{end_idx},\"end_element\":\"{element}\",\
             \"partner_atom_idx\":{partner_idx},\"partner_element\":\"{partner_element}\",\
             \"is_stereo_alkene_end\":true,\"candidate_substituent_count\":{sub_count},\
             \"coupled\":{coupled},\"component_size\":{component_size},\
             \"component_other_members\":{others_json},\
             \"end_atom_in_ring\":{end_in_ring},\"end_atom_ring_sizes\":{end_ring_json},\
             \"double_bond_endocyclic\":{endocyclic},\
             \"double_bond_endocyclic_ring_sizes\":{endo_ring_json},\
             \"correspondence_ok\":{correspondence_ok},\
             \"marker_placed\":{marker_placed_json},\
             \"marker_placed_count\":{marker_count_json}}}",
            smiles = esc(smiles),
            canon = esc(&canon),
            end_idx = end.0,
            element = element,
            partner_idx = partner.0,
            partner_element = partner_element,
            sub_count = subs.len(),
            coupled = coupled,
            component_size = others.len() + 1,
            others_json = atomidx_vec_json(&others),
            end_in_ring = !end_ring_sizes.is_empty(),
            end_ring_json = usize_vec_json(&end_ring_sizes),
            endocyclic = endocyclic,
            endo_ring_json = usize_vec_json(&endo_sizes),
            correspondence_ok = correspondence_ok,
            marker_placed_json = opt_bool_json(marker_placed),
            marker_count_json = opt_usize_json(marker_placed_count),
        );
    }
}

fn main() {
    let corpus_path = env::args().nth(1);
    let smiles_list: Vec<String> = match corpus_path {
        Some(path) => {
            let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
            text.lines()
                .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
                .filter(|l| !l.is_empty())
                .collect()
        }
        None => EZ_SHARED_CARRIER_FULLY_RESOLVED
            .iter()
            .chain(EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS.iter())
            .map(|s| s.to_string())
            .collect(),
    };

    for smiles in &smiles_list {
        analyze_molecule(smiles);
    }
}
