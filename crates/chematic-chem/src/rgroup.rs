//! R-group decomposition: split molecules into a common scaffold core and variable
//! R-group substituents.
//!
//! Given a SMARTS scaffold pattern and a list of molecules, this module identifies
//! for each molecule which atoms belong to the scaffold core and extracts the
//! substituents at each attachment point as separate R-group SMILES strings.
//!
//! # Usage
//!
//! ```
//! # use chematic_smiles::parse;
//! # use chematic_chem::rgroup::{rgroup_decompose, RGroupError};
//! let mols = vec![
//!     parse("CCc1ccccc1").unwrap(),
//!     parse("CCCc1ccccc1").unwrap(),
//!     parse("Nc1ccccc1").unwrap(),
//! ];
//! let refs: Vec<&_> = mols.iter().collect();
//! let results = rgroup_decompose("c1ccccc1", &refs).unwrap();
//! for r in results.iter().flatten() {
//!     println!("R1: {:?}", r.r_groups.get(&1));
//! }
//! ```

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, MatchConfig, QueryMolecule, find_matches_with_config, parse_smarts,
};
use chematic_smiles::{canonical_smiles, canonical_smiles_with_atom_order, parse};

/// Error type for R-group decomposition.
#[derive(Debug, Clone, PartialEq)]
pub enum RGroupError {
    /// The scaffold SMARTS string could not be parsed.
    InvalidSmarts(String),
    /// An internally generated canonical SMILES could not be parsed again.
    Canonicalization(String),
    /// A labelled decomposition requires positive, unique R-group labels.
    InvalidLabel(u16),
    /// The same R-group label appears more than once in the core query.
    DuplicateLabel(u16),
    /// A mapped wildcard attachment must be terminal in the core query.
    NonTerminalAttachmentLabel(u16),
    /// More than one substituent would be assigned to the same output column.
    MultipleAttachments { mol_idx: usize, label: u16 },
    /// The bounded matcher reached its cap before a deterministic match could be selected.
    MatchLimitExceeded { mol_idx: usize, limit: usize },
}

impl std::fmt::Display for RGroupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSmarts(s) => write!(f, "invalid SMARTS: {s}"),
            Self::Canonicalization(s) => {
                write!(
                    f,
                    "failed to normalize molecule before R-group matching: {s}"
                )
            }
            Self::InvalidLabel(label) => {
                write!(f, "R-group labels must be positive; found {label}")
            }
            Self::DuplicateLabel(label) => {
                write!(f, "R-group label R{label} appears more than once")
            }
            Self::NonTerminalAttachmentLabel(label) => write!(
                f,
                "mapped wildcard R{label} must be terminal in the core query"
            ),
            Self::MultipleAttachments { mol_idx, label } => write!(
                f,
                "molecule {mol_idx} has more than one substituent for R{label}"
            ),
            Self::MatchLimitExceeded { mol_idx, limit } => write!(
                f,
                "molecule {mol_idx} reached the bounded R-group match limit ({limit})"
            ),
        }
    }
}

impl std::error::Error for RGroupError {}

/// Result of decomposing one molecule against a scaffold.
#[derive(Debug, Clone)]
pub struct RGroupResult {
    /// Index of the molecule in the input slice.
    pub mol_idx: usize,
    /// Canonical SMILES of the scaffold core.
    /// Attachment positions (where R-groups were cut) carry `[*]` wildcard atoms.
    pub core_smiles: String,
    /// R-groups keyed by 1-based attachment-point index.
    ///
    /// The attachment points are ordered by the scaffold-atom index (ascending)
    /// at which the cut was made.  Each value is a SMILES string where `[*]`
    /// represents the bond back to the scaffold core.
    pub r_groups: HashMap<u8, String>,
}

/// Label-stable result compatible with RDKit-style R-group tables.
///
/// Labels come from atom-map numbers on the core SMARTS. A terminal mapped
/// wildcard (for example ``[*:2]``) labels its attachment, while a mapped core
/// atom labels substituents attached to that atom. Unlabelled attachment sites
/// receive deterministic labels based on query-atom order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LabeledRGroupResult {
    /// Index of the molecule in the input slice.
    pub mol_idx: usize,
    /// Canonical core SMILES with mapped wildcard attachment points.
    pub core_smiles: String,
    /// R-groups keyed by the stable R-group label.
    pub r_groups: BTreeMap<u16, String>,
}

const LABELED_MATCH_LIMIT: usize = 1024;

struct LabeledCandidate {
    missing_explicit: usize,
    signature: Vec<(u16, String, usize)>,
    result: LabeledRGroupResult,
}

/// Decompose `mols` against `scaffold_smarts`.
///
/// For each molecule the algorithm:
/// 1. Finds the first SMARTS match of the scaffold.
/// 2. Treats the matched atoms as the **core**.
/// 3. Identifies every bond crossing from a core atom to a non-core atom —
///    each such bond defines one R-group attachment point.
/// 4. Extracts the non-core subgraph reachable from that bond's endpoint (BFS,
///    blocked at core atoms) and writes it as a SMILES fragment with `[*]` at
///    the attachment.
///
/// Attachment points are numbered 1, 2, … in ascending order of the core atom's
/// molecule-atom index.
///
/// Returns `Vec<Option<RGroupResult>>` parallel to `mols`:
/// - `Some(result)` — the scaffold matched and R-groups were extracted.
/// - `None` — the scaffold did not match this molecule.
///
/// # Errors
///
/// Returns `RGroupError::InvalidSmarts` if `scaffold_smarts` fails to parse.
pub fn rgroup_decompose(
    scaffold_smarts: &str,
    mols: &[&Molecule],
) -> Result<Vec<Option<RGroupResult>>, RGroupError> {
    let query =
        parse_smarts(scaffold_smarts).map_err(|e| RGroupError::InvalidSmarts(e.to_string()))?;

    let cfg = MatchConfig {
        max_matches: Some(1),
        ..Default::default()
    };

    let mut results = Vec::with_capacity(mols.len());

    for (mol_idx, &mol) in mols.iter().enumerate() {
        let matches = find_matches_with_config(&query, mol, &cfg);

        if matches.is_empty() {
            results.push(None);
            continue;
        }

        // Core atoms: molecule atom indices that were matched by the scaffold.
        // Note: find_matches returns HashMap<query_atom_idx, mol_atom_idx>.
        let core: HashSet<AtomIdx> = matches[0].values().copied().collect();

        // Find attachment points: pairs (core_atom, rgroup_root) for bonds
        // crossing the core boundary.  Collect them sorted by core atom index
        // for stable R-group numbering.
        let mut attachment_pairs: Vec<(AtomIdx, AtomIdx)> = Vec::new();
        {
            let mut seen_bonds: HashSet<(u32, u32)> = HashSet::new();
            let mut core_sorted: Vec<AtomIdx> = core.iter().copied().collect();
            core_sorted.sort_by_key(|a| a.0);

            for &ca in &core_sorted {
                for (nb, _) in mol.neighbors(ca) {
                    if !core.contains(&nb) {
                        let key = (ca.0.min(nb.0), ca.0.max(nb.0));
                        if seen_bonds.insert(key) {
                            attachment_pairs.push((ca, nb));
                        }
                    }
                }
            }
        }

        // For each attachment point, BFS the non-core subgraph and build a fragment.
        let mut r_groups: HashMap<u8, String> = HashMap::new();
        for (rg_idx, &(_core_atom, rg_root)) in attachment_pairs.iter().enumerate() {
            let rg_num = (rg_idx + 1) as u8;
            let fragment = extract_rgroup(mol, rg_root, &core, None);
            r_groups.insert(rg_num, canonical_smiles(&fragment));
        }

        // Build core SMILES: copy core atoms, replacing boundary bonds with [*].
        let core_mol = build_core_mol(mol, &core, &attachment_pairs, None);
        let core_smiles = canonical_smiles(&core_mol);

        results.push(Some(RGroupResult {
            mol_idx,
            core_smiles,
            r_groups,
        }));
    }

    Ok(results)
}

/// Decompose molecules into a label-stable R-group table.
///
/// This is the bounded counterpart of RDKit's labelled
/// ``RGroupDecomposition`` surface. The core query may label attachment sites
/// in either of these forms:
///
/// - terminal mapped wildcard: ``c1cc([*:1])ccc1[*:2]``;
/// - mapped core atom: ``[c:1]1ccccc1``.
///
/// Terminal wildcard atoms are attachment placeholders and are not retained as
/// part of the core. Every output core and R-group carries the corresponding
/// atom-map label on its wildcard. If the query contains symmetric matches, a
/// canonical target-atom order is used before comparing output signatures, so
/// reordering the input molecule does not change the selected row.
///
/// Unlabelled attachment sites are assigned stable labels from query-atom
/// order, skipping explicit labels. A molecule that does not match the core is
/// returned as ``None``. Ambiguous multiple substituents for one label and
/// queries that hit the bounded match limit are refused with typed errors.
pub fn rgroup_decompose_labeled(
    scaffold_smarts: &str,
    mols: &[&Molecule],
) -> Result<Vec<Option<LabeledRGroupResult>>, RGroupError> {
    let query =
        parse_smarts(scaffold_smarts).map_err(|e| RGroupError::InvalidSmarts(e.to_string()))?;
    let labels = query_labels(&query)?;
    let auto_labels = auto_labels(&query, &labels);
    let cfg = MatchConfig {
        max_matches: Some(LABELED_MATCH_LIMIT),
        // Labelled symmetric sites must retain their alternative embeddings;
        // deduplicating by target-atom set would make R1/R2 depend on input
        // atom order.
        uniquify: false,
        ..Default::default()
    };

    let mut rows = Vec::with_capacity(mols.len());
    let mut column_counts: BTreeMap<u16, HashMap<String, usize>> = BTreeMap::new();
    for (mol_idx, &input_mol) in mols.iter().enumerate() {
        // Match a canonicalized copy. The API returns strings rather than
        // source atom indices, so this removes atom-order dependence without
        // losing any result identity.
        let normalized_smiles = canonical_smiles(input_mol);
        let normalized = parse(&normalized_smiles)
            .map_err(|err| RGroupError::Canonicalization(err.to_string()))?;
        let matches = find_matches_with_config(&query, &normalized, &cfg);
        if matches.is_empty() {
            rows.push(None);
            continue;
        }
        if matches.len() == LABELED_MATCH_LIMIT {
            return Err(RGroupError::MatchLimitExceeded {
                mol_idx,
                limit: LABELED_MATCH_LIMIT,
            });
        }

        let canonical_rank = canonical_atom_ranks(&normalized);
        let mut candidates = Vec::with_capacity(matches.len());
        let mut first_error = None;
        for mapping in &matches {
            match labeled_candidate(
                mol_idx,
                &normalized,
                &query,
                mapping,
                &labels,
                &auto_labels,
                &canonical_rank,
            ) {
                Ok(candidate) => candidates.push(candidate),
                Err(err) => {
                    first_error.get_or_insert(err);
                }
            };
        }

        candidates.sort_by(|a, b| compare_labeled_candidates(a, b, &column_counts));
        if let Some(candidate) = candidates.into_iter().next() {
            for (&label, smiles) in &candidate.result.r_groups {
                *column_counts
                    .entry(label)
                    .or_default()
                    .entry(smiles.clone())
                    .or_default() += 1;
            }
            rows.push(Some(candidate.result));
        } else if let Some(err) = first_error {
            return Err(err);
        } else {
            rows.push(None);
        }
    }
    Ok(rows)
}

fn compare_labeled_candidates(
    a: &LabeledCandidate,
    b: &LabeledCandidate,
    column_counts: &BTreeMap<u16, HashMap<String, usize>>,
) -> std::cmp::Ordering {
    let alignment_score = |candidate: &LabeledCandidate| {
        candidate
            .result
            .r_groups
            .iter()
            .map(|(label, smiles)| {
                column_counts
                    .get(label)
                    .and_then(|counts| counts.get(smiles))
                    .copied()
                    .unwrap_or(0)
            })
            .sum::<usize>()
    };
    a.missing_explicit
        .cmp(&b.missing_explicit)
        .then_with(|| alignment_score(b).cmp(&alignment_score(a)))
        // RDKit's row scorer uses the first row to break a symmetric tie and
        // then favours agreement with established columns. Descending output
        // signatures reproduce that bounded behaviour without a GA.
        .then_with(|| b.signature.cmp(&a.signature))
}

fn query_labels(query: &QueryMolecule) -> Result<HashMap<usize, u16>, RGroupError> {
    let mut labels = HashMap::new();
    let mut seen = HashSet::new();
    for (query_idx, atom) in query.atoms.iter().enumerate() {
        let Some(label) = atom.atom_map else {
            continue;
        };
        if label == 0 {
            return Err(RGroupError::InvalidLabel(label));
        }
        if !seen.insert(label) {
            return Err(RGroupError::DuplicateLabel(label));
        }
        if is_wildcard(&atom.query) && query.adj[query_idx].len() != 1 {
            return Err(RGroupError::NonTerminalAttachmentLabel(label));
        }
        labels.insert(query_idx, label);
    }
    Ok(labels)
}

fn auto_labels(query: &QueryMolecule, explicit: &HashMap<usize, u16>) -> HashMap<usize, u16> {
    let mut used: HashSet<u16> = explicit.values().copied().collect();
    let mut labels = HashMap::new();
    let mut next = 1u16;
    for query_idx in 0..query.atoms.len() {
        if explicit.contains_key(&query_idx) || is_terminal_mapped_wildcard(query, query_idx) {
            continue;
        }
        while used.contains(&next) {
            next = next.saturating_add(1);
        }
        labels.insert(query_idx, next);
        used.insert(next);
        next = next.saturating_add(1);
    }
    labels
}

fn is_wildcard(query: &AtomQuery) -> bool {
    matches!(query, AtomQuery::Primitive(AtomPrimitive::Wildcard))
}

fn is_terminal_mapped_wildcard(query: &QueryMolecule, idx: usize) -> bool {
    query.atoms[idx].atom_map.is_some()
        && is_wildcard(&query.atoms[idx].query)
        && query.adj[idx].len() == 1
}

fn canonical_atom_ranks(mol: &Molecule) -> Vec<usize> {
    let (_, order) = canonical_smiles_with_atom_order(mol);
    let mut ranks = vec![usize::MAX; mol.atom_count()];
    for (rank, atom) in order.into_iter().enumerate() {
        ranks[atom.0 as usize] = rank;
    }
    ranks
}

#[allow(clippy::too_many_arguments)]
fn labeled_candidate(
    mol_idx: usize,
    mol: &Molecule,
    query: &QueryMolecule,
    mapping: &rustc_hash::FxHashMap<usize, AtomIdx>,
    explicit_labels: &HashMap<usize, u16>,
    auto_labels: &HashMap<usize, u16>,
    canonical_rank: &[usize],
) -> Result<LabeledCandidate, RGroupError> {
    let placeholder_queries: HashSet<usize> = (0..query.atoms.len())
        .filter(|&idx| is_terminal_mapped_wildcard(query, idx))
        .collect();
    let core: HashSet<AtomIdx> = mapping
        .iter()
        .filter_map(|(query_idx, atom_idx)| {
            (!placeholder_queries.contains(query_idx)).then_some(*atom_idx)
        })
        .collect();
    let target_to_query: HashMap<AtomIdx, usize> = mapping
        .iter()
        .filter_map(|(query_idx, atom_idx)| {
            (!placeholder_queries.contains(query_idx)).then_some((*atom_idx, *query_idx))
        })
        .collect();

    let mut attachments: BTreeMap<u16, (AtomIdx, AtomIdx)> = BTreeMap::new();
    for &placeholder_idx in &placeholder_queries {
        let label = explicit_labels[&placeholder_idx];
        let core_query_idx = query.adj[placeholder_idx][0].1;
        let pair = (mapping[&core_query_idx], mapping[&placeholder_idx]);
        insert_attachment(&mut attachments, mol_idx, label, pair)?;
    }

    let placeholder_pairs: HashSet<(AtomIdx, AtomIdx)> = attachments.values().copied().collect();
    for &core_atom in &core {
        let Some(&query_idx) = target_to_query.get(&core_atom) else {
            continue;
        };
        let label = explicit_labels
            .get(&query_idx)
            .or_else(|| auto_labels.get(&query_idx))
            .copied()
            .expect("each core query atom has an explicit or automatic label");
        for (neighbor, _) in mol.neighbors(core_atom) {
            if core.contains(&neighbor) || placeholder_pairs.contains(&(core_atom, neighbor)) {
                continue;
            }
            insert_attachment(&mut attachments, mol_idx, label, (core_atom, neighbor))?;
        }
    }

    let mut r_groups = BTreeMap::new();
    let mut signature = Vec::with_capacity(attachments.len());
    for (&label, &(core_atom, root)) in &attachments {
        let fragment = extract_rgroup_labeled(mol, root, &core, label);
        let smiles = canonical_smiles(&fragment);
        signature.push((label, smiles.clone(), canonical_rank[core_atom.0 as usize]));
        r_groups.insert(label, smiles);
    }
    signature.sort();
    let core_mol = build_core_mol_labeled(mol, &core, &attachments);
    let core_smiles = canonical_smiles(&core_mol);
    signature.push((0, core_smiles.clone(), 0));
    let explicit_values: HashSet<u16> = explicit_labels.values().copied().collect();
    let missing_explicit = explicit_values
        .iter()
        .filter(|label| !attachments.contains_key(label))
        .count();

    Ok(LabeledCandidate {
        missing_explicit,
        signature,
        result: LabeledRGroupResult {
            mol_idx,
            core_smiles,
            r_groups,
        },
    })
}

fn insert_attachment(
    attachments: &mut BTreeMap<u16, (AtomIdx, AtomIdx)>,
    mol_idx: usize,
    label: u16,
    pair: (AtomIdx, AtomIdx),
) -> Result<(), RGroupError> {
    if attachments.insert(label, pair).is_some() {
        return Err(RGroupError::MultipleAttachments { mol_idx, label });
    }
    Ok(())
}

fn extract_rgroup_labeled(
    mol: &Molecule,
    root: AtomIdx,
    core: &HashSet<AtomIdx>,
    label: u16,
) -> Molecule {
    extract_rgroup(mol, root, core, Some(label))
}

fn build_core_mol_labeled(
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    attachments: &BTreeMap<u16, (AtomIdx, AtomIdx)>,
) -> Molecule {
    let pairs: Vec<(AtomIdx, AtomIdx)> = attachments.values().copied().collect();
    let labels: Vec<u16> = attachments.keys().copied().collect();
    build_core_mol(mol, core, &pairs, Some(&labels))
}

/// Extract the R-group subgraph starting at `root`, blocked by the `core` atom set.
///
/// Returns a `Molecule` containing all atoms reachable from `root` through bonds
/// that do not enter the core, plus a `[*]` wildcard stub at the attachment point.
fn extract_rgroup(
    mol: &Molecule,
    root: AtomIdx,
    core: &HashSet<AtomIdx>,
    wildcard_label: Option<u16>,
) -> Molecule {
    // BFS to collect the R-group atom set.
    let mut rg_atoms: HashSet<AtomIdx> = HashSet::new();
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();
    rg_atoms.insert(root);
    queue.push_back(root);

    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if !core.contains(&nb) && !rg_atoms.contains(&nb) {
                rg_atoms.insert(nb);
                queue.push_back(nb);
            }
        }
    }

    // Build subgraph molecule, remapping indices.
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Add root first so its index is predictable (index 0).
    let new_root = builder.add_atom(mol.atom(root).clone());
    remap.insert(root, new_root);

    let mut sorted: Vec<AtomIdx> = rg_atoms.iter().copied().filter(|&a| a != root).collect();
    sorted.sort_by_key(|a| a.0);
    for &old_idx in &sorted {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }

    // Add bonds within the R-group.
    for bidx in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(bidx as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }

    // Add a [*] wildcard attached to root representing the attachment to the core.
    let mut wildcard = Atom::wildcard();
    wildcard.atom_map = wildcard_label;
    let wc = builder.add_atom(wildcard);
    let _ = builder.add_bond(new_root, wc, BondOrder::Single);

    builder.build()
}

/// Build the core molecule with `[*]` stubs at each attachment point.
fn build_core_mol(
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    attachment_pairs: &[(AtomIdx, AtomIdx)],
    labels: Option<&[u16]>,
) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Add all core atoms (stable order).
    let mut core_sorted: Vec<AtomIdx> = core.iter().copied().collect();
    core_sorted.sort_by_key(|a| a.0);
    for &old_idx in &core_sorted {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }

    // Add bonds between core atoms.
    for bidx in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(bidx as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    // Add [*] wildcards at each attachment point.
    for (attachment_idx, &(core_atom, _rg_root)) in attachment_pairs.iter().enumerate() {
        if let Some(&new_ca) = remap.get(&core_atom) {
            let mut wildcard = Atom::wildcard();
            wildcard.atom_map = labels.map(|values| values[attachment_idx]);
            let wc = builder.add_atom(wildcard);
            let _ = builder.add_bond(new_ca, wc, BondOrder::Single);
        }
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    fn decompose<'a>(
        scaffold: &str,
        smiles: impl IntoIterator<Item = &'a str>,
    ) -> Vec<Option<RGroupResult>> {
        let mols: Vec<Molecule> = smiles.into_iter().map(mol).collect();
        let refs: Vec<&Molecule> = mols.iter().collect();
        rgroup_decompose(scaffold, &refs).expect("decompose failed")
    }

    #[test]
    fn test_rgroup_monosubstituted_benzene() {
        // Three alkylbenzenes — scaffold = benzene ring (6 C, no attachment markers).
        let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1"]);
        assert_eq!(results.len(), 3);
        for (i, r) in results.iter().enumerate() {
            let r = r.as_ref().unwrap_or_else(|| panic!("mol {i} should match"));
            assert_eq!(r.mol_idx, i);
            // Must have exactly one R-group (one substituted position).
            assert_eq!(r.r_groups.len(), 1, "mol {i}: expected 1 R-group");
            let r1 = r.r_groups.get(&1).expect("R1 missing");
            // R1 must contain a wildcard marker and the appropriate alkyl chain.
            assert!(
                r1.contains('*'),
                "R1 must contain [*] attachment, got {r1:?}"
            );
        }
    }

    #[test]
    fn test_rgroup_mol_idx_preserved() {
        let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCCC", "CCCc1ccccc1"]);
        // CCCC does not match benzene scaffold → None
        assert!(results[0].is_some());
        assert!(results[1].is_none());
        assert!(results[2].is_some());
        assert_eq!(results[0].as_ref().unwrap().mol_idx, 0);
        assert_eq!(results[2].as_ref().unwrap().mol_idx, 2);
    }

    #[test]
    fn test_rgroup_no_substituents() {
        // Benzene itself: scaffold matches, zero R-groups.
        let results = decompose("c1ccccc1", ["c1ccccc1"]);
        let r = results[0]
            .as_ref()
            .expect("benzene should match its own scaffold");
        assert_eq!(r.r_groups.len(), 0, "benzene has no R-groups vs itself");
    }

    #[test]
    fn test_rgroup_two_substituents() {
        // para-xylene: two methyl groups → two R-groups.
        let results = decompose("c1ccccc1", ["Cc1ccc(C)cc1"]);
        let r = results[0].as_ref().expect("should match");
        assert_eq!(r.r_groups.len(), 2, "para-xylene has two R-groups");
    }

    #[test]
    fn test_rgroup_invalid_smarts() {
        let mols = [mol("c1ccccc1")];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let err = rgroup_decompose("((invalid", &refs);
        assert!(matches!(err, Err(RGroupError::InvalidSmarts(_))));
    }

    #[test]
    fn test_rgroup_different_heteroatom_substituents() {
        // Aniline and phenol vs benzene scaffold.
        let results = decompose("c1ccccc1", ["Nc1ccccc1", "Oc1ccccc1"]);
        let aniline = results[0].as_ref().expect("aniline matches");
        let phenol = results[1].as_ref().expect("phenol matches");
        assert_eq!(aniline.r_groups.len(), 1);
        assert_eq!(phenol.r_groups.len(), 1);
        // Both R1 values must start with N or O (amine / hydroxyl).
        let r1_a = aniline.r_groups.get(&1).unwrap();
        let r1_p = phenol.r_groups.get(&1).unwrap();
        assert!(
            r1_a.starts_with('N') || r1_a.contains('N'),
            "aniline R1 should contain N, got {r1_a}"
        );
        assert!(
            r1_p.starts_with('O') || r1_p.contains('O'),
            "phenol R1 should contain O, got {r1_p}"
        );
    }

    #[test]
    fn labelled_dummy_atoms_align_symmetric_columns_like_rdkit() {
        let mols = [mol("Cc1ccc(CC)cc1"), mol("CCc1ccc(N)cc1")];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let rows = rgroup_decompose_labeled("c1cc([*:1])ccc1[*:2]", &refs)
            .expect("labelled decomposition");
        let first = rows[0].as_ref().expect("first molecule matches");
        let second = rows[1].as_ref().expect("second molecule matches");

        assert_eq!(first.core_smiles, second.core_smiles);
        assert_eq!(first.r_groups.keys().copied().collect::<Vec<_>>(), [1, 2]);
        assert!(first.r_groups[&1].contains("CC"));
        assert!(first.r_groups[&2].contains('C'));
        assert!(second.r_groups[&1].contains("CC"));
        assert!(second.r_groups[&2].contains('N'));
        assert!(first.core_smiles.contains("*:1"));
        assert!(first.core_smiles.contains("*:2"));
        assert!(first.r_groups[&1].contains("*:1"));
        assert!(first.r_groups[&2].contains("*:2"));
    }

    #[test]
    fn mapped_core_atom_labels_its_substituent() {
        let mols = [mol("Cc1ccccc1")];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let rows =
            rgroup_decompose_labeled("[c:7]1ccccc1", &refs).expect("mapped-core decomposition");
        let row = rows[0].as_ref().expect("toluene matches");
        assert_eq!(row.r_groups.keys().copied().collect::<Vec<_>>(), [7]);
        assert!(row.core_smiles.contains("*:7"));
        assert!(row.r_groups[&7].contains("*:7"));
    }

    #[test]
    fn duplicate_and_nonterminal_labels_are_typed_errors() {
        let mols = [mol("Cc1ccccc1")];
        let refs: Vec<&Molecule> = mols.iter().collect();

        let duplicate = rgroup_decompose_labeled("[c:1]1cccc[c:1]1", &refs);
        assert!(matches!(duplicate, Err(RGroupError::DuplicateLabel(1))));

        let nonterminal = rgroup_decompose_labeled("c1cc([*:2]C)ccc1", &refs);
        assert!(matches!(
            nonterminal,
            Err(RGroupError::NonTerminalAttachmentLabel(2))
        ));
    }
}
