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

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
    STEREO_VACANT_SENTINEL, StereoGroup, implicit_hcount,
};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, MatchConfig, QueryMolecule,
    find_matches_with_config, parse_smarts,
};
use chematic_smiles::{canonical_smiles, canonical_smiles_with_atom_order};

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

/// Parsed, validated core query for repeated labelled R-group decomposition.
///
/// Constructing this object once avoids reparsing the same SMARTS core for
/// every batch.  The decomposition semantics are identical to
/// [`rgroup_decompose_labeled`]; the free function is a convenience wrapper
/// around this prepared form.
pub struct PreparedRGroupCore {
    query: QueryMolecule,
    explicit_labels: HashMap<usize, u16>,
    auto_labels: HashMap<usize, u16>,
}

impl PreparedRGroupCore {
    /// Parse and validate a labelled R-group core.
    pub fn new(scaffold_smarts: &str) -> Result<Self, RGroupError> {
        let query =
            parse_smarts(scaffold_smarts).map_err(|e| RGroupError::InvalidSmarts(e.to_string()))?;
        let explicit_labels = query_labels(&query)?;
        let auto_labels = auto_labels(&query, &explicit_labels);
        Ok(Self {
            query,
            explicit_labels,
            auto_labels,
        })
    }

    /// Decompose `mols` using the prepared core.
    pub fn decompose_labeled(
        &self,
        mols: &[&Molecule],
    ) -> Result<Vec<Option<LabeledRGroupResult>>, RGroupError> {
        decompose_labeled_prepared(&self.query, &self.explicit_labels, &self.auto_labels, mols)
    }
}

struct LabeledCandidate {
    missing_explicit: usize,
    signature: Vec<(u16, String, usize)>,
    attachment_elements: Vec<(u16, u8)>,
    single_attachment_site: bool,
    bridged_attachment_sites: bool,
    rgroup_fingerprints: BTreeMap<u16, chematic_fp::BitVec2048>,
    core: HashSet<AtomIdx>,
    attachments: BTreeMap<u16, (AtomIdx, AtomIdx)>,
    hydrogen_sites: BTreeMap<u16, AtomIdx>,
    core_stub_orders: BTreeMap<u16, BondOrder>,
    result: LabeledRGroupResult,
}

#[derive(Clone)]
struct PlaceholderSites {
    attachments: BTreeMap<u16, (AtomIdx, AtomIdx)>,
    hydrogen_sites: BTreeMap<u16, AtomIdx>,
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

        // A connected substituent can bridge two core sites. Each column must
        // carry the complete fragment with all of its boundary stubs.
        let boundaries: Vec<_> = attachment_pairs
            .iter()
            .map(|&(core_atom, root)| (core_atom, root, None))
            .collect();
        let mut r_groups: HashMap<u8, String> = HashMap::new();
        for (rg_idx, &(_, rg_root)) in attachment_pairs.iter().enumerate() {
            let rg_num = (rg_idx + 1) as u8;
            let fragment = extract_rgroup(mol, rg_root, &core, &boundaries);
            r_groups.insert(rg_num, canonical_smiles(&fragment));
        }

        // Build core SMILES: copy core atoms, replacing boundary bonds with [*].
        let core_mol = build_core_mol(mol, &core, &attachment_pairs, None, None);
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
    PreparedRGroupCore::new(scaffold_smarts)?.decompose_labeled(mols)
}

fn decompose_labeled_prepared(
    query: &QueryMolecule,
    labels: &HashMap<usize, u16>,
    auto_labels: &HashMap<usize, u16>,
    mols: &[&Molecule],
) -> Result<Vec<Option<LabeledRGroupResult>>, RGroupError> {
    let optional_query = without_terminal_mapped_wildcards(query);
    let cfg = MatchConfig {
        max_matches: Some(LABELED_MATCH_LIMIT),
        // Labelled symmetric sites must retain their alternative embeddings;
        // deduplicating by target-atom set would make R1/R2 depend on input
        // atom order.
        uniquify: false,
        ..Default::default()
    };

    let mut rows: Vec<Option<LabeledCandidate>> = Vec::with_capacity(mols.len());
    let mut column_counts: BTreeMap<u16, HashMap<String, usize>> = BTreeMap::new();
    let mut column_fingerprints: BTreeMap<u16, Vec<chematic_fp::BitVec2048>> = BTreeMap::new();
    for (mol_idx, &input_mol) in mols.iter().enumerate() {
        // Enumerate every bounded embedding on the original graph, then use
        // canonical atom ranks and canonical fragment strings below to pick a
        // stable row.  The former implementation wrote a canonical SMILES and
        // parsed it back before matching.  That round trip did not add any
        // information to the result identity, but duplicated canonicalization
        // and graph construction for every input molecule.
        let mut matches = find_matches_with_config(query, input_mol, &cfg);
        if matches.is_empty()
            && let Some((reduced, original_indices)) = &optional_query
        {
            matches = find_matches_with_config(reduced, input_mol, &cfg)
                .into_iter()
                .map(|mapping| {
                    mapping
                        .into_iter()
                        .map(|(idx, atom)| (original_indices[idx], atom))
                        .collect()
                })
                .collect();
        }
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

        let canonical_rank = canonical_atom_ranks(input_mol);
        let mut candidates = Vec::with_capacity(matches.len());
        let mut first_error = None;
        for mapping in &matches {
            match labeled_candidates(
                mol_idx,
                input_mol,
                query,
                mapping,
                labels,
                auto_labels,
                &canonical_rank,
            ) {
                Ok(mut found) => candidates.append(&mut found),
                Err(err) => {
                    first_error.get_or_insert(err);
                }
            };
        }

        candidates
            .sort_by(|a, b| compare_labeled_candidates(a, b, &column_counts, &column_fingerprints));
        if let Some(candidate) = candidates.into_iter().next() {
            for (&label, smiles) in &candidate.result.r_groups {
                *column_counts
                    .entry(label)
                    .or_default()
                    .entry(smiles.clone())
                    .or_default() += 1;
            }
            for (label, fingerprint) in &candidate.rgroup_fingerprints {
                column_fingerprints
                    .entry(*label)
                    .or_default()
                    .push(fingerprint.clone());
            }
            rows.push(Some(candidate));
        } else if let Some(err) = first_error {
            return Err(err);
        } else {
            rows.push(None);
        }
    }
    Ok(finalize_labeled_rows(rows, mols))
}

fn finalize_labeled_rows(
    rows: Vec<Option<LabeledCandidate>>,
    mols: &[&Molecule],
) -> Vec<Option<LabeledRGroupResult>> {
    // RDKit emits an R-group hydrogen only when another row establishes that
    // column. A lone unsubstituted core therefore has no empty R column.
    let active_labels: HashSet<u16> = rows
        .iter()
        .flatten()
        .flat_map(|candidate| candidate.result.r_groups.keys().copied())
        .collect();
    rows.into_iter()
        .map(|row| {
            row.map(|mut candidate| {
                let mut attachments = candidate.attachments;
                let mut filled_hydrogen = false;
                for (label, site) in candidate.hydrogen_sites {
                    if active_labels.contains(&label) {
                        filled_hydrogen = true;
                        candidate
                            .result
                            .r_groups
                            .insert(label, canonical_smiles(&hydrogen_rgroup(label)));
                        // The second tuple component is unused by the core
                        // builder for a hydrogen attachment.
                        attachments.insert(label, (site, AtomIdx(u32::MAX)));
                    }
                }
                if filled_hydrogen {
                    candidate.result.core_smiles = canonical_smiles(&build_core_mol_labeled(
                        mols[candidate.result.mol_idx],
                        &candidate.core,
                        &attachments,
                        &candidate.core_stub_orders,
                    ));
                }
                candidate.result
            })
        })
        .collect()
}

// A terminal mapped wildcard describes a substituent site. Removing it from
// the fallback query lets an implicit hydrogen occupy that site when no full
// embedding exists. Keep original indices so the candidate builder sees one
// consistent query regardless of which matcher produced the embedding.
fn without_terminal_mapped_wildcards(query: &QueryMolecule) -> Option<(QueryMolecule, Vec<usize>)> {
    for idx in 0..query.atom_count() {
        if is_terminal_mapped_wildcard(query, idx) {
            let (bond_idx, _) = query.adj[idx][0];
            if !matches!(
                query.bonds[bond_idx].query,
                BondQuery::Any | BondQuery::Primitive(BondPrimitive::Single)
            ) {
                // A missing multiple-bond attachment cannot be hydrogen.
                return None;
            }
        }
    }
    let mut reduced = QueryMolecule::new();
    let mut original_indices = Vec::new();
    let mut new_indices = vec![None; query.atom_count()];
    for (idx, atom) in query.atoms.iter().enumerate() {
        if is_terminal_mapped_wildcard(query, idx) {
            continue;
        }
        new_indices[idx] = Some(reduced.add_atom_with_map(atom.query.clone(), atom.atom_map));
        original_indices.push(idx);
    }
    if original_indices.len() == query.atom_count() || original_indices.is_empty() {
        return None;
    }
    for bond in &query.bonds {
        if let (Some(a), Some(b)) = (new_indices[bond.atom1], new_indices[bond.atom2]) {
            reduced.add_bond(a, b, bond.query.clone());
        }
    }
    Some((reduced, original_indices))
}

fn compare_labeled_candidates(
    a: &LabeledCandidate,
    b: &LabeledCandidate,
    column_counts: &BTreeMap<u16, HashMap<String, usize>>,
    column_fingerprints: &BTreeMap<u16, Vec<chematic_fp::BitVec2048>>,
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
    let similarity_score = |candidate: &LabeledCandidate| {
        candidate
            .rgroup_fingerprints
            .iter()
            .map(|(label, fingerprint)| {
                column_fingerprints
                    .get(label)
                    .map(|column| {
                        column
                            .iter()
                            .map(|known| fingerprint.tanimoto(known))
                            .sum::<f64>()
                    })
                    .unwrap_or(0.0)
            })
            .sum::<f64>()
    };
    a.missing_explicit
        .cmp(&b.missing_explicit)
        .then_with(|| alignment_score(b).cmp(&alignment_score(a)))
        .then_with(|| similarity_score(b).total_cmp(&similarity_score(a)))
        // When a symmetric, partially substituted core has no established
        // column yet, RDKit puts the real group in the lowest available R
        // label and leaves the later label for hydrogen.
        .then_with(|| a.attachments.keys().cmp(b.attachments.keys()))
        // With no established column, RDKit places a heteroatom-rooted
        // substituent before a carbon-rooted one at equivalent core sites.
        .then_with(|| {
            if a.single_attachment_site && b.single_attachment_site {
                b.attachment_elements.cmp(&a.attachment_elements)
            } else {
                std::cmp::Ordering::Equal
            }
        })
        // For a ring segment spanning two core sites, RDKit assigns the
        // lower-numbered column to the carbon-rooted end before a heteroatom.
        .then_with(|| {
            if a.bridged_attachment_sites && b.bridged_attachment_sites {
                a.attachment_elements.cmp(&b.attachment_elements)
            } else {
                std::cmp::Ordering::Equal
            }
        })
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

fn placeholder_sites(
    mol_idx: usize,
    mol: &Molecule,
    query: &QueryMolecule,
    mapping: &rustc_hash::FxHashMap<usize, AtomIdx>,
    explicit_labels: &HashMap<usize, u16>,
    placeholder_queries: &HashSet<usize>,
    core: &HashSet<AtomIdx>,
) -> Result<Option<Vec<PlaceholderSites>>, RGroupError> {
    let mut variants = vec![PlaceholderSites {
        attachments: BTreeMap::new(),
        hydrogen_sites: BTreeMap::new(),
    }];
    let mut placeholders_by_anchor: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &placeholder_idx in placeholder_queries {
        placeholders_by_anchor
            .entry(query.adj[placeholder_idx][0].1)
            .or_default()
            .push(placeholder_idx);
    }
    for (anchor_idx, mut placeholders) in placeholders_by_anchor {
        placeholders.sort_by_key(|idx| explicit_labels[idx]);
        let site = mapping[&anchor_idx];
        if placeholders.iter().all(|idx| mapping.contains_key(idx)) {
            for variant in &mut variants {
                for &placeholder_idx in &placeholders {
                    let label = explicit_labels[&placeholder_idx];
                    insert_attachment(
                        &mut variant.attachments,
                        mol_idx,
                        label,
                        (site, mapping[&placeholder_idx]),
                    )?;
                }
            }
            continue;
        }
        if placeholders.iter().any(|idx| mapping.contains_key(idx)) {
            return Ok(None);
        }
        let external: Vec<_> = mol
            .neighbors(site)
            .filter_map(|(neighbor, _)| (!core.contains(&neighbor)).then_some(neighbor))
            .collect();
        // With two external groups, both assignments to the first two labels
        // must be scored against established columns. Keep the enumeration
        // bounded; larger sites still require a full-query embedding.
        if external.len() > 2
            || external.len() > placeholders.len()
            || usize::from(implicit_hcount(mol, site)) < placeholders.len() - external.len()
        {
            return Ok(None);
        }
        let mut orders = vec![external.clone()];
        if external.len() == 2 {
            orders.push(vec![external[1], external[0]]);
        }
        let mut next = Vec::new();
        for variant in variants {
            for order in &orders {
                let mut candidate = variant.clone();
                for (index, &placeholder_idx) in placeholders.iter().enumerate() {
                    let label = explicit_labels[&placeholder_idx];
                    if let Some(&root) = order.get(index) {
                        insert_attachment(
                            &mut candidate.attachments,
                            mol_idx,
                            label,
                            (site, root),
                        )?;
                    } else {
                        candidate.hydrogen_sites.insert(label, site);
                    }
                }
                next.push(candidate);
            }
        }
        if next.len() > 16 {
            return Ok(None);
        }
        variants = next;
    }
    Ok(Some(variants))
}

#[allow(clippy::too_many_arguments)]
fn labeled_candidates(
    mol_idx: usize,
    mol: &Molecule,
    query: &QueryMolecule,
    mapping: &rustc_hash::FxHashMap<usize, AtomIdx>,
    explicit_labels: &HashMap<usize, u16>,
    auto_labels: &HashMap<usize, u16>,
    canonical_rank: &[usize],
) -> Result<Vec<LabeledCandidate>, RGroupError> {
    let placeholder_queries: HashSet<usize> = (0..query.atoms.len())
        .filter(|&idx| is_terminal_mapped_wildcard(query, idx))
        .collect();
    let core_stub_orders: BTreeMap<u16, BondOrder> = placeholder_queries
        .iter()
        .filter_map(|&idx| {
            let (bond_idx, _) = query.adj[idx][0];
            let order = match query.bonds[bond_idx].query {
                BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                BondQuery::Primitive(BondPrimitive::Quadruple) => BondOrder::Quadruple,
                _ => return None,
            };
            Some((explicit_labels[&idx], order))
        })
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

    let Some(sites_variants) = placeholder_sites(
        mol_idx,
        mol,
        query,
        mapping,
        explicit_labels,
        &placeholder_queries,
        &core,
    )?
    else {
        return Ok(Vec::new());
    };
    sites_variants
        .into_iter()
        .map(|sites| {
            labeled_candidate_for_sites(
                mol_idx,
                mol,
                &core,
                &target_to_query,
                explicit_labels,
                auto_labels,
                canonical_rank,
                sites,
                &core_stub_orders,
            )
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn labeled_candidate_for_sites(
    mol_idx: usize,
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    target_to_query: &HashMap<AtomIdx, usize>,
    explicit_labels: &HashMap<usize, u16>,
    auto_labels: &HashMap<usize, u16>,
    canonical_rank: &[usize],
    sites: PlaceholderSites,
    core_stub_orders: &BTreeMap<u16, BondOrder>,
) -> Result<LabeledCandidate, RGroupError> {
    let PlaceholderSites {
        mut attachments,
        mut hydrogen_sites,
    } = sites;
    let placeholder_pairs: HashSet<(AtomIdx, AtomIdx)> = attachments.values().copied().collect();
    for &core_atom in core {
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
    let mut rgroup_fingerprints = BTreeMap::new();
    let mut signature = Vec::with_capacity(attachments.len());
    let mut attachment_elements = Vec::with_capacity(attachments.len());
    let mut attachment_sites = attachments.values().map(|&(site, _)| site);
    let first_site = attachment_sites.next();
    let single_attachment_site = attachment_sites.all(|site| Some(site) == first_site);
    let boundaries: Vec<_> = attachments
        .iter()
        .map(|(&label, &(core_atom, root))| (core_atom, root, Some(label)))
        .collect();
    for (&label, &(core_atom, root)) in &attachments {
        let fragment = extract_rgroup(mol, root, core, &boundaries);
        let smiles = canonical_smiles(&fragment);
        rgroup_fingerprints.insert(label, chematic_fp::ecfp4(&fragment));
        signature.push((label, smiles.clone(), canonical_rank[core_atom.0 as usize]));
        attachment_elements.push((label, mol.atom(root).element.atomic_number()));
        r_groups.insert(label, smiles);
    }
    let bridged_attachment_sites = r_groups.values().collect::<HashSet<_>>().len() < r_groups.len();
    signature.sort();
    let core_mol = build_core_mol_labeled(mol, core, &attachments, core_stub_orders);
    let core_smiles = canonical_smiles(&core_mol);
    signature.push((0, core_smiles.clone(), 0));
    let explicit_values: HashSet<u16> = explicit_labels.values().copied().collect();
    let missing_explicit = explicit_values
        .iter()
        .filter(|label| !attachments.contains_key(label))
        .count();
    // An automatically labelled site needs the same deferred hydrogen fill as
    // an explicitly labelled one once another row establishes its column.
    for (&site, &query_idx) in target_to_query {
        let label = explicit_labels
            .get(&query_idx)
            .or_else(|| auto_labels.get(&query_idx))
            .copied()
            .expect("each core query atom has an explicit or automatic label");
        if !attachments.contains_key(&label) && implicit_hcount(mol, site) > 0 {
            hydrogen_sites.insert(label, site);
        }
    }

    Ok(LabeledCandidate {
        missing_explicit,
        signature,
        attachment_elements,
        single_attachment_site,
        bridged_attachment_sites,
        rgroup_fingerprints,
        core: core.clone(),
        attachments,
        hydrogen_sites,
        core_stub_orders: core_stub_orders.clone(),
        result: LabeledRGroupResult {
            mol_idx,
            core_smiles,
            r_groups,
        },
    })
}

fn hydrogen_rgroup(label: u16) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let hydrogen = builder.add_atom(Atom::new(Element::H));
    let mut wildcard = Atom::wildcard();
    wildcard.atom_map = Some(label);
    let attachment = builder.add_atom(wildcard);
    let _ = builder.add_bond(hydrogen, attachment, BondOrder::Single);
    builder.build()
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

fn build_core_mol_labeled(
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    attachments: &BTreeMap<u16, (AtomIdx, AtomIdx)>,
    core_stub_orders: &BTreeMap<u16, BondOrder>,
) -> Molecule {
    let pairs: Vec<(AtomIdx, AtomIdx)> = attachments.values().copied().collect();
    let labels: Vec<u16> = attachments.keys().copied().collect();
    build_core_mol(mol, core, &pairs, Some(&labels), Some(core_stub_orders))
}

/// Carry stereochemical reference frames across a subgraph cut. Every removed
/// neighbor of a retained stereocenter must be replaced by its own attachment
/// wildcard in the original neighbor position, not just by a new bond.
fn copy_cut_stereo(
    mol: &Molecule,
    builder: &mut MoleculeBuilder,
    remap: &HashMap<AtomIdx, AtomIdx>,
    bond_remap: &HashMap<BondIdx, BondIdx>,
    cut_neighbors: &HashMap<(AtomIdx, AtomIdx), AtomIdx>,
) {
    for (&old_center, &new_center) in remap {
        if let Some(order) = mol.stereo_neighbor_order(old_center) {
            let mapped = order
                .iter()
                .map(|&neighbor| {
                    if matches!(neighbor, STEREO_H_SENTINEL | STEREO_VACANT_SENTINEL) {
                        Some(neighbor)
                    } else {
                        remap
                            .get(&AtomIdx(neighbor))
                            .or_else(|| cut_neighbors.get(&(old_center, AtomIdx(neighbor))))
                            .map(|idx| idx.0)
                    }
                })
                .collect::<Option<Vec<_>>>();
            if let Some(mapped) = mapped {
                builder.set_stereo_neighbor_order(new_center, mapped);
            }
        }
        if let Some(label) = mol.r_group_label(old_center) {
            builder.set_r_group(new_center, label);
        }
        if let Some(tag) = mol.atom_tag(old_center) {
            builder.set_tag(new_center, Some(tag.get()));
        }
    }

    for (&old_bond, &new_bond) in bond_remap {
        if let Some(direction) = mol.bond_direction(old_bond) {
            builder.set_bond_direction(new_bond, direction);
        }
        if let Some(anchor) = mol.bond_direction_anchor(old_bond)
            && let Some(&new_anchor) = remap.get(&anchor)
        {
            builder.set_bond_direction_anchor(new_bond, new_anchor);
        }
    }

    for (group_idx, group) in mol.stereo_groups().iter().enumerate() {
        let atoms = group
            .atom_indices
            .iter()
            .filter_map(|old| remap.get(old).copied())
            .collect::<Vec<_>>();
        let bonds = mol
            .stereo_group_bond_indices(group_idx)
            .unwrap_or_default()
            .iter()
            .filter_map(|old| bond_remap.get(old).copied())
            .collect::<Vec<_>>();
        if !atoms.is_empty() || !bonds.is_empty() {
            builder.add_stereo_group_with_bonds(StereoGroup::new(group.kind.clone(), atoms), bonds);
        }
    }
}

/// A directional boundary bond still describes an alkene when its retained
/// endpoint remains adjacent to the double bond after cutting.
fn cut_bond_order(
    mol: &Molecule,
    bond_idx: BondIdx,
    retained: AtomIdx,
    retained_atoms: &HashSet<AtomIdx>,
    preserve_multiple: bool,
) -> BondOrder {
    let bond = mol.bond(bond_idx);
    if matches!(
        bond.order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple | BondOrder::Aromatic
    ) {
        return if preserve_multiple {
            bond.order
        } else {
            BondOrder::Single
        };
    }
    let direction = mol.bond_direction(bond_idx).unwrap_or(bond.order);
    if matches!(direction, BondOrder::Up | BondOrder::Down)
        && mol.neighbors(retained).any(|(neighbor, bidx)| {
            retained_atoms.contains(&neighbor)
                && bidx != bond_idx
                && mol.bond(bidx).order == BondOrder::Double
        })
    {
        let anchor = mol.bond_direction_anchor(bond_idx).unwrap_or(bond.atom1);
        if anchor == retained {
            direction
        } else if direction == BondOrder::Up {
            BondOrder::Down
        } else {
            BondOrder::Up
        }
    } else {
        BondOrder::Single
    }
}

/// Extract the R-group subgraph starting at `root`, blocked by the `core` atom set.
///
/// Returns a `Molecule` containing all atoms reachable from `root` through bonds
/// that do not enter the core, plus one wildcard for each core boundary bond
/// incident on that component. A bridging substituent is therefore represented
/// by the same multi-attachment fragment in all of its R-group columns.
fn extract_rgroup(
    mol: &Molecule,
    root: AtomIdx,
    core: &HashSet<AtomIdx>,
    boundaries: &[(AtomIdx, AtomIdx, Option<u16>)],
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
    let mut bond_remap = HashMap::new();
    let mut cut_neighbors = HashMap::new();

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
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2))
            && let Ok(new_bond) = builder.add_bond(new_a, new_b, bond.order)
        {
            bond_remap.insert(BondIdx(bidx as u32), new_bond);
        }
    }

    for &(core_atom, boundary_root, label) in boundaries {
        let Some(&new_boundary_root) = remap.get(&boundary_root) else {
            continue;
        };
        let (boundary_idx, _) = mol
            .bond_between(core_atom, boundary_root)
            .expect("R-group boundary must be a molecular bond");
        let stub_order = cut_bond_order(mol, boundary_idx, boundary_root, &rg_atoms, true);
        let mut wildcard = Atom::wildcard();
        wildcard.atom_map = label;
        let wc = builder.add_atom(wildcard);
        let _ = builder.add_bond(new_boundary_root, wc, stub_order);
        cut_neighbors.insert((boundary_root, core_atom), wc);
    }

    copy_cut_stereo(mol, &mut builder, &remap, &bond_remap, &cut_neighbors);
    builder.build()
}

/// Build the core molecule with `[*]` stubs at each attachment point.
fn build_core_mol(
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    attachment_pairs: &[(AtomIdx, AtomIdx)],
    labels: Option<&[u16]>,
    core_stub_orders: Option<&BTreeMap<u16, BondOrder>>,
) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut bond_remap = HashMap::new();
    let mut cut_neighbors = HashMap::new();

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
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2))
            && let Ok(new_bond) = builder.add_bond(na, nb, bond.order)
        {
            bond_remap.insert(BondIdx(bidx as u32), new_bond);
        }
    }

    // Add [*] wildcards at each attachment point.
    for (attachment_idx, &(core_atom, rg_root)) in attachment_pairs.iter().enumerate() {
        if let Some(&new_ca) = remap.get(&core_atom) {
            let mut wildcard = Atom::wildcard();
            wildcard.atom_map = labels.map(|values| values[attachment_idx]);
            let wc = builder.add_atom(wildcard);
            let order = labels
                .and_then(|values| core_stub_orders?.get(&values[attachment_idx]))
                .copied()
                .unwrap_or_else(|| {
                    mol.bond_between(core_atom, rg_root)
                        .map(|(bond_idx, _)| cut_bond_order(mol, bond_idx, core_atom, core, false))
                        .unwrap_or(BondOrder::Single)
                });
            let _ = builder.add_bond(new_ca, wc, order);
            if mol.bond_between(core_atom, rg_root).is_some() {
                cut_neighbors.insert((core_atom, rg_root), wc);
            }
        }
    }

    copy_cut_stereo(mol, &mut builder, &remap, &bond_remap, &cut_neighbors);
    builder.build()
}

#[cfg(test)]
#[path = "rgroup_tests.rs"]
mod tests;
