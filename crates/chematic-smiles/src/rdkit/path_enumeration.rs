//! Ordered paths following RDKit Release_2026_09_1 Subgraphs.cpp.
//! Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors.
//! BSD-3-Clause; see THIRD_PARTY_NOTICES.md. Uses the existing sanitized model.
use super::{RdkitSmilesError, graph_queries, mol::Mol};
use chematic_core::Molecule;
mod search;

/// Options for RDKit's FindAllPathsOfLengthN, in sanitized model numbering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkitPathParams {
    /// True: length and output indices count bonds; false: they count atoms.
    pub use_bonds: bool,
    pub use_hs: bool,
    /// None starts at every atom. Out-of-range roots return no paths.
    pub rooted_at_atom: Option<usize>,
    /// Retains only geodesics; upstream considers graph Hs regardless of use_hs.
    pub only_shortest_paths: bool,
}
impl Default for RdkitPathParams {
    fn default() -> Self {
        Self {
            use_bonds: true,
            use_hs: false,
            rooted_at_atom: None,
            only_shortest_paths: false,
        }
    }
}

/// Independent work/storage caps. Exhaustion returns an error, never partial paths.
#[derive(Debug, Clone, Copy)]
pub struct RdkitPathLimits {
    /// Neighbor extensions examined, including rejected repeats/shortcuts.
    pub max_visits: usize,
    /// Prefix count in either breadth-first frontier, before deduplication.
    pub max_frontier: usize,
    /// Number of unique returned paths.
    pub max_results: usize,
    /// Total indices in either frontier, or in the returned paths, separately.
    pub max_path_indices: usize,
}
impl Default for RdkitPathLimits {
    fn default() -> Self {
        Self {
            max_visits: 1_000_000,
            max_frontier: 100_000,
            max_results: 100_000,
            max_path_indices: 1_000_000,
        }
    }
}

fn unsupported(detail: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(detail.into())
}
fn reserve<T>(v: &mut Vec<T>, count: usize) -> Result<(), RdkitSmilesError> {
    v.try_reserve(count)
        .map_err(|e| unsupported(format!("path allocation failed: {e}")))
}
type Neighbors = Vec<Vec<(usize, usize)>>;

/// FindAllPathsOfLengthN with ordered atom-index expansion and bond-set dedup.
/// Unlike branched subgraphs, prefixes are linear. The last step can close a
/// ring to any earlier atom except the immediately previous one. Both zero
/// lengths return empty; atom length one includes isolated/graph-H atoms.
pub fn rdkit_paths_of_length(
    mol: &Molecule,
    length: u32,
    params: &RdkitPathParams,
) -> Result<Vec<Vec<usize>>, RdkitSmilesError> {
    rdkit_paths_of_length_with_limits(mol, length, params, &RdkitPathLimits::default())
}

/// Same sanitized model and <=1024 input atoms /16384 bonds as graph queries.
/// Very large lengths return empty once impossible, without looping to length.
pub fn rdkit_paths_of_length_with_limits(
    mol: &Molecule,
    length: u32,
    params: &RdkitPathParams,
    limits: &RdkitPathLimits,
) -> Result<Vec<Vec<usize>>, RdkitSmilesError> {
    let model = graph_queries::prepare(mol)?;
    let n = model.atoms.len();
    if length == 0 || params.rooted_at_atom.is_some_and(|root| root >= n) {
        return Ok(Vec::new());
    }
    // Compare before adding one so even a 32-bit target cannot overflow usize.
    if u64::from(length) > n as u64 + u64::from(!params.use_bonds) {
        return Ok(Vec::new());
    }
    let target = length as usize + usize::from(params.use_bonds);
    // With shortest-path selection RDKit builds adjacency from the unfiltered
    // distance matrix, which includes graph Hs even when useHs is false.
    let neighbors = neighbors(&model, params.use_hs || params.only_shortest_paths)?;
    search::enumerate(&neighbors, target, params, limits)
}

fn neighbors(model: &Mol, use_hs: bool) -> Result<Neighbors, RdkitSmilesError> {
    let mut result = Vec::new();
    reserve(&mut result, model.atoms.len())?;
    result.resize_with(model.atoms.len(), Vec::new);
    for (index, bond) in model.bonds.iter().enumerate() {
        if use_hs || (model.atoms[bond.begin].anum != 1 && model.atoms[bond.end].anum != 1) {
            reserve(&mut result[bond.begin], 1)?;
            result[bond.begin].push((bond.end, index));
            reserve(&mut result[bond.end], 1)?;
            result[bond.end].push((bond.begin, index));
        }
    }
    for row in &mut result {
        row.sort_unstable();
        // The upstream adjacency matrix has one edge per atom pair and bond
        // lookup selects the first bond. Keep that same first-index relation.
        row.dedup_by_key(|(atom, _)| *atom);
    }
    Ok(result)
}
