//! Wave 3 (issue #149) diagnosis-only audit of the shared-carrier coupling
//! population NOT explained by the ring-endocyclic mechanism PR #351 fixed
//! (Wave 2D). PR #351's own body states it closes only ~10% (3 of 31) of the
//! corpus's general coupling-component population -- the other ~90% (28 of
//! 31) is "a separate, still-unidentified mechanism." This tool measures
//! that remaining population directly, rather than assuming its size or
//! shape from the ~90% figure (which was a topological-presence count, not a
//! confirmed-permutation-invariance-failure count).
//!
//! Read-only, public-API-only (matching the established convention of both
//! prior audit examples in this file's family): every private
//! `canonical.rs` helper this tool needs (`compute_stereo_alkene_ends`
//! *including* PR #351's ring gate, `coupling_components`, `reference_up`,
//! `direction_is_up`/`direction_for_up`, `raw_bond_direction`,
//! `alternate_ez_markings`) is reimplemented verbatim from its production
//! counterpart using only `chematic_core`/`chematic_perception`/
//! `chematic_smiles` public API -- confirmed reimplementable via
//! `Molecule::bond_direction` (public) and `Molecule::with_bond_order`
//! (public). No production code is touched or generalized.
//!
//! ## Why this file, not either existing example
//!
//! `ez_shared_carrier_component_audit.rs` (Wave 2B) and
//! `ez_ring_constrained_residual_audit.rs` (Wave 2C) are each frozen
//! snapshots scoped to an already-answered hypothesis -- neither includes
//! PR #351's ring gate, so running either today still reports the *pre-fix*
//! topology (31 components on `scripts/descriptor_census_corpus.smi`,
//! reconfirmed unchanged by re-running the Wave 2B example directly).
//! This file measures the *current* population instead.
//!
//! ## Subcommands
//!
//! ```text
//! # Current (ring-gate-aware) topology + per-end/per-component facts.
//! # Default (no corpus arg): the 18 pinned EZ_SHARED_CARRIER_FULLY_RESOLVED
//! # fixtures + the 2 never-corrupts SMILES.
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- scan
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- scan scripts/descriptor_census_corpus.smi
//!
//! # Axis 2 (mark relocation, no RDKit): one SMILES per line in <file>.
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- axis2 <file>
//!
//! # Axis 1 (RDKit relabeling): <file> has "original_smiles\trelabeled_smiles"
//! # lines, relabeled_smiles supplied by the Python driver via RDKit.
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- axis1 <file>
//! ```

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

/// The 18 pinned `EZ_SHARED_CARRIER_FULLY_RESOLVED` fixtures (canonical.rs)
/// -- all already fully resolved post-PR-#351, included here only as a
/// negative-control population (Wave 3 expects 0 residual on these).
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
    r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
    r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
];

/// The permanent regression fixture (`ez_carrier_shared_bond_between_two_
/// stereo_systems_never_corrupts`, canonical.rs) -- two spellings of a real
/// corpus molecule (two exocyclic imines on a 4-membered ring sharing a
/// ring-closure bond) that do NOT converge to one canonical string, though
/// geometry is provably preserved either way. Both alkenes are exocyclic --
/// NOT the ring-endocyclic shape PR #351 fixed. This is the calibration
/// anchor: any new probe must reproduce "these two spellings give different
/// canonical output" before it is trusted on unlabeled corpus rows.
const EZ_NEVER_CORRUPTS_RESIDUAL: &[&str] = &[
    r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC",
    r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\1c(/c(c1O)O)=N/CCCCC",
];

// ---------------------------------------------------------------------------
// Topology, reimplemented verbatim (public-API only) from `canonical.rs`'s
// private `CanonicalWriter` methods -- see module doc comment for why this
// is safe (every helper below is a direct, line-for-line port confirmed
// against the production source, not a reinterpretation).
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

/// Mirrors `canonical.rs`'s `double_bond_endocyclic_in_small_ring` exactly
/// (ring size < 8, both endpoints of the double bond in the same SSSR ring).
fn double_bond_endocyclic_in_small_ring(rings: &[Vec<AtomIdx>], a: AtomIdx, b: AtomIdx) -> bool {
    rings
        .iter()
        .any(|r| r.len() < 8 && r.contains(&a) && r.contains(&b))
}

/// Mirrors `canonical.rs`'s CURRENT (post-PR-#351) `compute_stereo_alkene_
/// ends` exactly, including the ring-endocyclic exclusion -- unlike both
/// existing example files, which predate that gate.
fn stereo_alkene_end_nodes_current(mol: &Molecule) -> HashSet<AtomIdx> {
    let candidates: Vec<_> = mol
        .bonds()
        .filter(|(_, bond)| bond.order == BondOrder::Double)
        .filter(|(_, bond)| {
            end_has_substituent(mol, bond.atom1) && end_has_substituent(mol, bond.atom2)
        })
        .collect();
    if candidates.is_empty() {
        return HashSet::new();
    }

    let rings = find_sssr(mol);
    let rings = rings.rings();
    let mut ends = HashSet::new();
    for (_, bond) in candidates {
        if double_bond_endocyclic_in_small_ring(rings, bond.atom1, bond.atom2) {
            continue;
        }
        for end in [bond.atom1, bond.atom2] {
            if substituents(mol, end).len() == 2 {
                ends.insert(end);
            }
        }
    }
    ends
}

/// One coupling-component: nodes = ambiguous ends, edges = a direct
/// non-double bond between two such ends (the shared candidate bond).
/// Unlike `ez_shared_carrier_component_audit.rs`'s `ComponentReport`, this
/// keeps the actual edge list -- field 3 of the requested JSON output (the
/// "carrier conflict graph"), not just an aggregate size/cycle count.
struct Component {
    members: Vec<AtomIdx>,
    edges: Vec<(AtomIdx, AtomIdx, BondIdx)>,
}

fn coupling_components_with_edges(mol: &Molecule, ends: &HashSet<AtomIdx>) -> Vec<Component> {
    let mut adjacency: HashMap<AtomIdx, Vec<(AtomIdx, BondIdx)>> = HashMap::new();
    for &n in ends {
        for (sub_atom, bidx) in substituents(mol, n) {
            if ends.contains(&sub_atom) {
                adjacency.entry(n).or_default().push((sub_atom, bidx));
            }
        }
    }

    let mut starts: Vec<AtomIdx> = ends.iter().copied().collect();
    starts.sort_by_key(|a| a.0);
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut components = Vec::new();
    for start in starts {
        if visited.contains(&start) {
            continue;
        }
        let mut queue = vec![start];
        let mut members = Vec::new();
        visited.insert(start);
        while let Some(cur) = queue.pop() {
            members.push(cur);
            let mut nbs = adjacency.get(&cur).cloned().unwrap_or_default();
            nbs.sort_by_key(|&(a, _)| a.0);
            for (nb, _) in &nbs {
                if visited.insert(*nb) {
                    queue.push(*nb);
                }
            }
        }
        members.sort_by_key(|a| a.0);
        let member_set: HashSet<AtomIdx> = members.iter().copied().collect();
        let mut edges: Vec<(AtomIdx, AtomIdx, BondIdx)> = Vec::new();
        let mut seen_bonds: HashSet<BondIdx> = HashSet::new();
        for &m in &members {
            if let Some(nbs) = adjacency.get(&m) {
                for &(nb, bidx) in nbs {
                    if member_set.contains(&nb) && seen_bonds.insert(bidx) {
                        let (a, b) = if m.0 <= nb.0 { (m, nb) } else { (nb, m) };
                        edges.push((a, b, bidx));
                    }
                }
            }
        }
        components.push(Component { members, edges });
    }
    components
}

/// Mirrors `canonical.rs::CanonicalWriter::direction_is_up` exactly.
fn direction_is_up(dir: BondOrder, bond_atom1: AtomIdx, alkene_end: AtomIdx) -> bool {
    match dir {
        BondOrder::Up => bond_atom1 == alkene_end,
        BondOrder::Down => bond_atom1 != alkene_end,
        _ => false,
    }
}

/// Mirrors `canonical.rs::CanonicalWriter::direction_for_up` exactly.
fn direction_for_up(bond_atom1: AtomIdx, alkene_end: AtomIdx, want_up: bool) -> BondOrder {
    if (bond_atom1 == alkene_end) == want_up {
        BondOrder::Up
    } else {
        BondOrder::Down
    }
}

/// Mirrors `writer::raw_bond_direction` exactly, via the public
/// `Molecule::bond_direction` accessor (the aromatic-bond-direction-stash
/// reader) instead of the `pub(crate)` original.
fn raw_input_direction(mol: &Molecule, bidx: BondIdx) -> Option<BondOrder> {
    let order = mol.bond(bidx).order;
    if matches!(order, BondOrder::Up | BondOrder::Down) {
        return Some(order);
    }
    mol.bond_direction(bidx)
}

/// Mirrors `canonical.rs::CanonicalWriter::reference_up` exactly (the
/// rank-based, input-mark-placement-invariant geometry fact one end's own
/// rank-preferred candidate encodes).
fn reference_up(
    mol: &Molecule,
    alkene_end: AtomIdx,
    subs: &[(AtomIdx, BondIdx); 2],
    pref_idx: usize,
) -> Option<bool> {
    let reference = subs[pref_idx];
    let sibling = subs[1 - pref_idx];
    if let Some(dir) = raw_input_direction(mol, reference.1) {
        Some(direction_is_up(
            dir,
            mol.bond(reference.1).atom1,
            alkene_end,
        ))
    } else {
        let dir = raw_input_direction(mol, sibling.1)?;
        Some(!direction_is_up(dir, mol.bond(sibling.1).atom1, alkene_end))
    }
}

/// Which of `subs[0]`/`subs[1]` is rank-preferred (lower rank value =
/// higher canonical priority, matching `resolve_component_jointly`'s own
/// `pref` computation: index 1 wins when strictly lower-ranked).
fn pref_index(ranks: &[u64], subs: &[(AtomIdx, BondIdx); 2]) -> usize {
    if ranks[subs[1].0.0 as usize] < ranks[subs[0].0.0 as usize] {
        1
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Atom correspondence across a canonicalize -> reparse round trip (or,
// generally, across any two same-molecule spellings), via composition of
// `canonical_atom_order` (public, relabeling-invariant) -- verbatim from
// `ez_ring_constrained_residual_audit.rs`, see its module doc comment for
// why this is safe and how it is verified.
// ---------------------------------------------------------------------------

fn position_of(order: &[usize]) -> Vec<usize> {
    let mut pos = vec![0usize; order.len()];
    for (position, &atom_idx) in order.iter().enumerate() {
        pos[atom_idx] = position;
    }
    pos
}

fn bond_order_class_eq(a: BondOrder, b: BondOrder) -> bool {
    let plain = |o: BondOrder| matches!(o, BondOrder::Single | BondOrder::Up | BondOrder::Down);
    if plain(a) && plain(b) {
        return true;
    }
    a == b
}

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

fn map_atom(idx: AtomIdx, pos1: &[usize], order2: &[usize]) -> AtomIdx {
    AtomIdx(order2[pos1[idx.0 as usize]] as u32)
}

// ---------------------------------------------------------------------------
// JSON emission (hand-rolled, no serde_json -- matches this crate's
// existing example convention).
// ---------------------------------------------------------------------------

fn esc(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
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

fn atomidx_vec_json(v: &[AtomIdx]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|a| a.0.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn component_id(members: &[AtomIdx]) -> String {
    members
        .iter()
        .map(|a| a.0.to_string())
        .collect::<Vec<_>>()
        .join("-")
}

// ---------------------------------------------------------------------------
// Correspondence-trick marker-placement recovery (verbatim technique from
// `ez_ring_constrained_residual_audit.rs`): canonicalize `mol`, reparse its
// own output, verify correspondence, and check whether the reparsed
// molecule marks `end`'s candidate bonds `Up`/`Down`.
// ---------------------------------------------------------------------------

struct MarkerReadout {
    canonical_smiles: String,
    correspondence_ok: bool,
    /// `end atom idx -> (marker_placed, marker_placed_count)`, only for ends
    /// passed in; `None` for both fields when correspondence failed.
    per_end: HashMap<AtomIdx, (Option<bool>, Option<usize>)>,
}

fn read_marker_placement(mol: &Molecule, ends: &[AtomIdx]) -> MarkerReadout {
    let canon = chematic_smiles::canonical_smiles(mol);
    let order1 = chematic_smiles::canonical_atom_order(mol);
    let pos1 = position_of(&order1);

    let (correspondence_ok, mol2_and_order2) = match chematic_smiles::parse(&canon) {
        Ok(mol2) => {
            let order2 = chematic_smiles::canonical_atom_order(&mol2);
            let ok = verify_correspondence(mol, &pos1, &mol2, &order2);
            (ok, Some((mol2, order2)))
        }
        Err(_) => (false, None),
    };

    let mut per_end = HashMap::new();
    for &end in ends {
        let value = if correspondence_ok {
            if let Some((mol2, order2)) = &mol2_and_order2 {
                let mapped_end = map_atom(end, &pos1, order2);
                let subs = substituents(mol, end);
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
        per_end.insert(end, value);
    }

    MarkerReadout {
        canonical_smiles: canon,
        correspondence_ok,
        per_end,
    }
}

// ---------------------------------------------------------------------------
// `scan`: current (ring-gate-aware) topology + per-end/per-component facts.
// ---------------------------------------------------------------------------

fn scan_molecule(source_tag: &str, smiles: &str) {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("parse failure, skipped ({source_tag}): {smiles:?}: {e}");
            return;
        }
    };

    let ends = stereo_alkene_end_nodes_current(&mol);
    if ends.is_empty() {
        return;
    }

    let ranks = chematic_smiles::morgan_ranks(&mol);
    let rings = find_sssr(&mol).rings().to_vec();
    let components = coupling_components_with_edges(&mol, &ends);

    let mut ends_sorted: Vec<AtomIdx> = ends.iter().copied().collect();
    ends_sorted.sort_by_key(|a| a.0);
    let readout = read_marker_placement(&mol, &ends_sorted);

    // One "component" row per component (including size-1 singletons, so
    // the branch-point/shape claim is auditable against the full
    // population, not just the coupled subset).
    for comp in &components {
        let shape = if comp.members.len() < 2 {
            "singleton"
        } else if comp.edges.len() == comp.members.len() {
            "cycle"
        } else {
            "path"
        };
        let mut tie_break_order = comp.members.clone();
        tie_break_order.sort_by_key(|&a| ranks[a.0 as usize]);
        let edges_json = format!(
            "[{}]",
            comp.edges
                .iter()
                .map(|(a, b, bidx)| format!(
                    "{{\"a\":{},\"b\":{},\"bond_idx\":{}}}",
                    a.0, b.0, bidx.0
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
        println!(
            "{{\"kind\":\"component\",\"source\":\"{source}\",\"smiles\":\"{smiles}\",\
             \"canonical_smiles\":\"{canon}\",\"component_id\":\"{cid}\",\
             \"members\":{members},\"size\":{size},\"shape\":\"{shape}\",\
             \"edges\":{edges},\"tie_break_order\":{tbo},\
             \"propagation_order_applicable\":false}}",
            source = esc(source_tag),
            smiles = esc(smiles),
            canon = esc(&readout.canonical_smiles),
            cid = component_id(&comp.members),
            members = atomidx_vec_json(&comp.members),
            size = comp.members.len(),
            shape = shape,
            edges = edges_json,
            tbo = atomidx_vec_json(&tie_break_order),
        );
    }

    let mut member_of: HashMap<AtomIdx, (String, usize)> = HashMap::new();
    for comp in &components {
        for &m in &comp.members {
            member_of.insert(m, (component_id(&comp.members), comp.members.len()));
        }
    }

    for &end in &ends_sorted {
        let element = mol.atom(end).element.symbol();
        let subs = substituents(&mol, end);
        debug_assert_eq!(subs.len(), 2);
        let subs2: [(AtomIdx, BondIdx); 2] = [subs[0], subs[1]];
        let pref = pref_index(&ranks, &subs2);
        let ref_up = reference_up(&mol, end, &subs2, pref);

        let candidate_bonds = format!(
            "[{}]",
            subs2
                .iter()
                .map(|&(other, bidx)| format!(
                    "{{\"other_atom_idx\":{},\"bond_idx\":{},\"current_bond_order\":\"{:?}\"}}",
                    other.0,
                    bidx.0,
                    mol.bond(bidx).order
                ))
                .collect::<Vec<_>>()
                .join(",")
        );

        let end_ring_sizes: Vec<usize> = {
            let mut v: Vec<usize> = rings
                .iter()
                .filter(|r| r.contains(&end))
                .map(|r| r.len())
                .collect();
            v.sort_unstable();
            v
        };
        let end_ring_json = format!(
            "[{}]",
            end_ring_sizes
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );

        let (component_id_str, component_size) = member_of
            .get(&end)
            .cloned()
            .unwrap_or_else(|| (component_id(&[end]), 1));
        let coupled = component_size >= 2;

        let (marker_placed, marker_placed_count) =
            readout.per_end.get(&end).copied().unwrap_or((None, None));

        println!(
            "{{\"kind\":\"end\",\"source\":\"{source}\",\"smiles\":\"{smiles}\",\
             \"canonical_smiles\":\"{canon}\",\"end_atom_idx\":{end_idx},\
             \"end_element\":\"{element}\",\"end_rank\":{end_rank},\
             \"candidate_bonds\":{candidate_bonds},\"pref_idx\":{pref},\
             \"reference_up\":{ref_up_json},\
             \"end_atom_in_ring\":{end_in_ring},\"end_atom_ring_sizes\":{end_ring_json},\
             \"component_id\":\"{cid}\",\"component_size\":{component_size},\
             \"coupled\":{coupled},\"correspondence_ok\":{correspondence_ok},\
             \"marker_placed\":{marker_placed_json},\
             \"marker_placed_count\":{marker_placed_count_json}}}",
            source = esc(source_tag),
            smiles = esc(smiles),
            canon = esc(&readout.canonical_smiles),
            end_idx = end.0,
            element = element,
            end_rank = ranks[end.0 as usize],
            candidate_bonds = candidate_bonds,
            pref = pref,
            ref_up_json = opt_bool_json(ref_up),
            end_in_ring = !end_ring_sizes.is_empty(),
            end_ring_json = end_ring_json,
            cid = component_id_str,
            component_size = component_size,
            coupled = coupled,
            correspondence_ok = readout.correspondence_ok,
            marker_placed_json = opt_bool_json(marker_placed),
            marker_placed_count_json = opt_usize_json(marker_placed_count),
        );
    }
}

fn run_scan(corpus_path: Option<String>) {
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
            .chain(EZ_NEVER_CORRUPTS_RESIDUAL.iter())
            .map(|s| s.to_string())
            .collect(),
    };
    for smiles in &smiles_list {
        scan_molecule("scan", smiles);
    }
}

// ---------------------------------------------------------------------------
// `axis2`: mark-relocation alternates, no RDKit. Reimplements
// `canonical.rs`'s private `alternate_ez_markings` test helper verbatim via
// public `Molecule::with_bond_order`.
// ---------------------------------------------------------------------------

/// Mirrors `canonical.rs`'s test-only `geometry_fingerprint` exactly: each
/// stereogenic double bond's E/Z parity, read via the rank-based
/// `reference_up`-equivalent oracle (marker-*placement*-invariant), keyed
/// in an intrinsic (rank-based) bond order so it is comparable across two
/// spellings of the same molecule.
fn geometry_fingerprint(mol: &Molecule, ranks: &[u64]) -> Vec<Option<bool>> {
    let mut doubles: Vec<(u64, BondIdx)> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .filter(|(_, b)| end_has_substituent(mol, b.atom1) && end_has_substituent(mol, b.atom2))
        .map(|(bidx, b)| {
            let key = ranks[b.atom1.0 as usize].min(ranks[b.atom2.0 as usize]);
            (key, bidx)
        })
        .collect();
    doubles.sort_by_key(|&(k, _)| k);

    doubles
        .into_iter()
        .map(|(_, bidx)| {
            let bond = mol.bond(bidx);
            let ua = up_of_reference(mol, ranks, bond.atom1)?;
            let ub = up_of_reference(mol, ranks, bond.atom2)?;
            Some(ua != ub)
        })
        .collect()
}

/// Mirrors `canonical.rs`'s test-only `up_of_reference` oracle exactly.
fn up_of_reference(mol: &Molecule, ranks: &[u64], end: AtomIdx) -> Option<bool> {
    let subs = substituents(mol, end);
    match subs.len() {
        1 => {
            let (_, bidx) = subs[0];
            let dir = raw_input_direction(mol, bidx)?;
            Some(direction_is_up(dir, mol.bond(bidx).atom1, end))
        }
        2 => {
            let reference = *subs
                .iter()
                .min_by_key(|&&(a, _)| ranks[a.0 as usize])
                .expect("subs has 2 elements");
            let sibling = if reference.0 == subs[0].0 {
                subs[1]
            } else {
                subs[0]
            };
            if let Some(dir) = raw_input_direction(mol, reference.1) {
                Some(direction_is_up(dir, mol.bond(reference.1).atom1, end))
            } else {
                let dir = raw_input_direction(mol, sibling.1)?;
                Some(!direction_is_up(dir, mol.bond(sibling.1).atom1, end))
            }
        }
        _ => None,
    }
}

struct Alternate {
    moved_end: AtomIdx,
    moved_from_bond: BondIdx,
    moved_to_bond: BondIdx,
    mol: Molecule,
}

/// Mirrors `canonical.rs`'s test-only `alternate_ez_markings` exactly:
/// every geometry-preserving relocation of an explicit `/`/`\` marker from
/// one of an end's two candidate bonds onto the other.
fn alternate_ez_markings(mol: &Molecule) -> Vec<Alternate> {
    let ranks = chematic_smiles::morgan_ranks(mol);
    let baseline_geo = geometry_fingerprint(mol, &ranks);
    let ends = stereo_alkene_end_nodes_current(mol);
    let mut alternates = Vec::new();
    let mut sorted_ends: Vec<AtomIdx> = ends.iter().copied().collect();
    sorted_ends.sort_by_key(|a| a.0);

    for end in sorted_ends {
        let subs = substituents(mol, end);
        if subs.len() != 2 {
            continue;
        }
        for i in 0..2 {
            let chosen = subs[i];
            let other = subs[1 - i];
            let chosen_order = mol.bond(chosen.1).order;
            let other_order = mol.bond(other.1).order;
            let chosen_dir = match chosen_order {
                BondOrder::Up | BondOrder::Down => chosen_order,
                _ => continue,
            };
            if other_order != BondOrder::Single {
                continue;
            }
            let up = direction_is_up(chosen_dir, mol.bond(chosen.1).atom1, end);
            let new_other_dir = direction_for_up(mol.bond(other.1).atom1, end, !up);

            let alt = mol
                .with_bond_order(chosen.1, BondOrder::Single)
                .with_bond_order(other.1, new_other_dir);
            if geometry_fingerprint(&alt, &ranks) != baseline_geo {
                continue;
            }
            alternates.push(Alternate {
                moved_end: end,
                moved_from_bond: chosen.1,
                moved_to_bond: other.1,
                mol: alt,
            });
        }
    }
    alternates
}

fn run_axis2(path: &str) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    for line in text.lines() {
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        let mol = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("parse failure, skipped: {smiles:?}: {e}");
                continue;
            }
        };
        let ends = stereo_alkene_end_nodes_current(&mol);
        if ends.is_empty() {
            continue;
        }
        let mut ends_sorted: Vec<AtomIdx> = ends.iter().copied().collect();
        ends_sorted.sort_by_key(|a| a.0);

        let baseline = read_marker_placement(&mol, &ends_sorted);
        println!(
            "{{\"kind\":\"axis2_variant\",\"source_smiles\":\"{smiles}\",\
             \"variant\":\"baseline\",\"canonical_smiles\":\"{canon}\",\
             \"correspondence_ok\":{ok}}}",
            smiles = esc(smiles),
            canon = esc(&baseline.canonical_smiles),
            ok = baseline.correspondence_ok,
        );

        for (idx, alt) in alternate_ez_markings(&mol).into_iter().enumerate() {
            let readout = read_marker_placement(&alt.mol, &ends_sorted);
            println!(
                "{{\"kind\":\"axis2_variant\",\"source_smiles\":\"{smiles}\",\
                 \"variant\":\"alt-{idx}\",\"moved_end_atom_idx\":{end},\
                 \"moved_from_bond_idx\":{from},\"moved_to_bond_idx\":{to},\
                 \"canonical_smiles\":\"{canon}\",\"correspondence_ok\":{ok},\
                 \"differs_from_baseline\":{differs}}}",
                smiles = esc(smiles),
                idx = idx,
                end = alt.moved_end.0,
                from = alt.moved_from_bond.0,
                to = alt.moved_to_bond.0,
                canon = esc(&readout.canonical_smiles),
                ok = readout.correspondence_ok,
                differs = readout.canonical_smiles != baseline.canonical_smiles,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// `axis1`: RDKit-relabeled variants. Correspondence is two-hop: original ->
// relabeled (atom-identity across respellings, via the same
// canonical_atom_order composition trick) and relabeled -> reparse of
// relabeled's own canonical output (to read where the marker landed).
// ---------------------------------------------------------------------------

fn run_axis1(path: &str) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let (Some(original), Some(relabeled)) = (parts.next(), parts.next()) else {
            eprintln!("malformed line (expected original<TAB>relabeled), skipped: {line:?}");
            continue;
        };

        let mol1 = match chematic_smiles::parse(original) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("original parse failure, skipped: {original:?}: {e}");
                continue;
            }
        };
        let ends1 = stereo_alkene_end_nodes_current(&mol1);
        if ends1.is_empty() {
            continue;
        }
        let mut ends1_sorted: Vec<AtomIdx> = ends1.iter().copied().collect();
        ends1_sorted.sort_by_key(|a| a.0);

        let order1 = chematic_smiles::canonical_atom_order(&mol1);
        let pos1 = position_of(&order1);

        let mol2 = match chematic_smiles::parse(relabeled) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("relabeled parse failure, skipped: {relabeled:?}: {e}");
                continue;
            }
        };
        let order2 = chematic_smiles::canonical_atom_order(&mol2);
        let cross_ok = verify_correspondence(&mol1, &pos1, &mol2, &order2);

        let readout2 = read_marker_placement(&mol2, &[]);
        println!(
            "{{\"kind\":\"axis1_variant\",\"original_smiles\":\"{orig}\",\
             \"relabeled_smiles\":\"{relab}\",\"cross_correspondence_ok\":{cross_ok},\
             \"canonical_smiles\":\"{canon}\",\"correspondence_ok\":{ok}}}",
            orig = esc(original),
            relab = esc(relabeled),
            cross_ok = cross_ok,
            canon = esc(&readout2.canonical_smiles),
            ok = readout2.correspondence_ok,
        );

        if !cross_ok {
            continue;
        }
        for &end in &ends1_sorted {
            let mapped_end = map_atom(end, &pos1, &order2);
            let subs2 = substituents(&mol2, mapped_end);
            if subs2.len() != 2 {
                // ring gate or topology differs at this mapped atom in mol2 --
                // report, don't guess.
                println!(
                    "{{\"kind\":\"axis1_end\",\"original_smiles\":\"{orig}\",\
                     \"relabeled_smiles\":\"{relab}\",\"original_end_atom_idx\":{end},\
                     \"mapped_end_atom_idx\":{mapped},\"mapped_is_stereo_alkene_end\":false}}",
                    orig = esc(original),
                    relab = esc(relabeled),
                    end = end.0,
                    mapped = mapped_end.0,
                );
                continue;
            }
            let readout_for_end = read_marker_placement(&mol2, &[mapped_end]);
            let (marker_placed, marker_placed_count) = readout_for_end
                .per_end
                .get(&mapped_end)
                .copied()
                .unwrap_or((None, None));
            println!(
                "{{\"kind\":\"axis1_end\",\"original_smiles\":\"{orig}\",\
                 \"relabeled_smiles\":\"{relab}\",\"original_end_atom_idx\":{end},\
                 \"mapped_end_atom_idx\":{mapped},\"mapped_is_stereo_alkene_end\":true,\
                 \"marker_placed\":{marker_placed_json},\
                 \"marker_placed_count\":{marker_placed_count_json}}}",
                orig = esc(original),
                relab = esc(relabeled),
                end = end.0,
                mapped = mapped_end.0,
                marker_placed_json = opt_bool_json(marker_placed),
                marker_placed_count_json = opt_usize_json(marker_placed_count),
            );
        }
    }
}

fn main() {
    let mut args = env::args().skip(1);
    let subcommand = args.next().unwrap_or_else(|| "scan".to_string());
    let file_arg = args.next();

    match subcommand.as_str() {
        "scan" => run_scan(file_arg),
        "axis2" => run_axis2(&file_arg.expect("axis2 requires a <file> argument")),
        "axis1" => run_axis1(&file_arg.expect("axis1 requires a <file> argument")),
        other => {
            eprintln!("unknown subcommand {other:?} -- expected scan|axis2|axis1");
            std::process::exit(1);
        }
    }
}
