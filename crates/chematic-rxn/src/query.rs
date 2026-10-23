//! Reaction SMARTS querying for chemical reaction matching.

use chematic_smarts::{MatchConfig, QueryMolecule, find_matches_with_config, parse_smarts};
use rustc_hash::FxHashMap;
use std::collections::HashSet;

use crate::reaction::Reaction;

/// A reaction query consisting of reactant and product SMARTS patterns.
#[derive(Clone, Debug)]
pub struct ReactionQuery {
    /// SMARTS queries for reactant pattern matching.
    pub reactant_patterns: Vec<QueryMolecule>,
    /// Disconnected SMARTS components grouped by `|` alternative.
    pub reactant_components: Vec<Vec<QueryMolecule>>,
    /// SMARTS queries for agent pattern matching.
    pub agent_patterns: Vec<QueryMolecule>,
    /// Disconnected agent components grouped by `|` alternative.
    pub agent_components: Vec<Vec<QueryMolecule>>,
    /// SMARTS queries for product pattern matching.
    pub product_patterns: Vec<QueryMolecule>,
    /// Disconnected SMARTS components grouped by `|` alternative.
    pub product_components: Vec<Vec<QueryMolecule>>,
}

/// A reaction SMARTS pattern with atom mapping and agent support.
#[derive(Clone, Debug)]
pub struct ReactionSmartsPattern {
    /// SMARTS queries for reactant pattern matching.
    pub reactant_patterns: Vec<QueryMolecule>,
    /// Disconnected reactant components grouped by `|` alternative.
    pub reactant_components: Vec<Vec<QueryMolecule>>,
    /// SMARTS queries for agent pattern matching (optional middle section).
    pub agent_patterns: Vec<QueryMolecule>,
    /// Disconnected agent components grouped by `|` alternative.
    pub agent_components: Vec<Vec<QueryMolecule>>,
    /// SMARTS queries for product pattern matching.
    pub product_patterns: Vec<QueryMolecule>,
    /// Disconnected product components grouped by `|` alternative.
    pub product_components: Vec<Vec<QueryMolecule>>,
    /// Metadata about atom map numbers (e.g., :1, :2, etc.)
    pub map_number_info: MapNumberInfo,
}

/// Information about atom map numbers in a reaction SMARTS pattern.
#[derive(Clone, Debug)]
pub struct MapNumberInfo {
    /// All unique map numbers found in the pattern.
    pub all_map_numbers: HashSet<u16>,
    /// Map numbers found in reactants.
    pub reactant_maps: HashSet<u16>,
    /// Map numbers found in agents.
    pub agent_maps: HashSet<u16>,
    /// Map numbers found in products.
    pub product_maps: HashSet<u16>,
}

/// Detailed information about a single reactant or product pattern match.
#[derive(Clone, Debug)]
pub struct MoleculeMatch {
    /// Index of the molecule in the reaction (reactants or products list).
    pub molecule_index: usize,
    /// Index of the pattern that matched.
    pub pattern_index: usize,
    /// Indices of the matched atoms in the molecule.
    pub atom_indices: Vec<usize>,
}

/// Detailed information about reactant pattern matches in a reaction.
#[derive(Clone, Debug)]
pub struct ReactantMatches {
    /// Matches for each reactant pattern.
    pub pattern_matches: Vec<Vec<MoleculeMatch>>,
}

impl ReactantMatches {
    /// Get all molecules that matched a specific reactant pattern.
    pub fn get_pattern_matches(&self, pattern_index: usize) -> Option<&[MoleculeMatch]> {
        self.pattern_matches
            .get(pattern_index)
            .map(|v| v.as_slice())
    }

    /// Check if a specific reactant pattern matched any molecule.
    pub fn pattern_matched(&self, pattern_index: usize) -> bool {
        self.pattern_matches
            .get(pattern_index)
            .is_some_and(|matches| !matches.is_empty())
    }
}

/// Detailed information about product pattern matches in a reaction.
#[derive(Clone, Debug)]
pub struct ProductMatches {
    /// Matches for each product pattern.
    pub pattern_matches: Vec<Vec<MoleculeMatch>>,
}

impl ProductMatches {
    /// Get all molecules that matched a specific product pattern.
    pub fn get_pattern_matches(&self, pattern_index: usize) -> Option<&[MoleculeMatch]> {
        self.pattern_matches
            .get(pattern_index)
            .map(|v| v.as_slice())
    }

    /// Check if a specific product pattern matched any molecule.
    pub fn pattern_matched(&self, pattern_index: usize) -> bool {
        self.pattern_matches
            .get(pattern_index)
            .is_some_and(|matches| !matches.is_empty())
    }
}

/// Detailed information about agent pattern matches in a reaction.
///
/// Agent patterns separated by `|` are alternatives. A non-empty list is
/// therefore considered satisfied when at least one pattern has a match.
#[derive(Clone, Debug)]
pub struct AgentMatches {
    /// Matches for each agent alternative pattern.
    pub pattern_matches: Vec<Vec<MoleculeMatch>>,
}

impl AgentMatches {
    /// Get all molecules that matched a specific agent pattern.
    pub fn get_pattern_matches(&self, pattern_index: usize) -> Option<&[MoleculeMatch]> {
        self.pattern_matches
            .get(pattern_index)
            .map(|v| v.as_slice())
    }

    /// Check if at least one agent alternative matched.
    pub fn any_pattern_matched(&self) -> bool {
        self.pattern_matches
            .iter()
            .any(|matches| !matches.is_empty())
    }
}

/// Detailed match information for a reaction against a SMARTS pattern.
#[derive(Clone, Debug)]
pub struct ReactionSmartsMatch {
    /// Reactant pattern matches (all patterns must have at least one match for overall match).
    pub reactant_matches: ReactantMatches,
    /// Agent alternative matches (`|` uses OR semantics).
    pub agent_matches: AgentMatches,
    /// Product pattern matches (all patterns must have at least one match for overall match).
    pub product_matches: ProductMatches,
    /// Whether the reactant component groups were satisfied, including their
    /// dot-separated injective and pipe-separated alternative semantics.
    /// This is intentionally stored separately from the flattened detailed
    /// match list, which cannot recover which alternatives belong to one
    /// group.
    pub reactants_satisfied: bool,
    /// Whether the product component groups were satisfied.
    pub products_satisfied: bool,
    /// Whether the agent component groups were satisfied.
    pub agents_satisfied: bool,
    /// Whether all patterns matched (true if reaction is valid against the query).
    pub is_complete_match: bool,
}

impl ReactionSmartsMatch {
    /// Check if all reactant patterns matched.
    pub fn all_reactants_matched(&self) -> bool {
        self.reactants_satisfied
    }

    /// Check if all product patterns matched.
    pub fn all_products_matched(&self) -> bool {
        self.products_satisfied
    }
}

impl MapNumberInfo {
    /// Check if all map numbers are consistent across reactants and products.
    pub fn validate(&self) -> Result<(), String> {
        // All map numbers in reactants must appear in products
        let mut missing_in_products: Vec<u16> = self
            .reactant_maps
            .iter()
            .filter(|&m| !self.product_maps.contains(m))
            .copied()
            .collect();
        // HashSet iteration is intentionally unordered. Sort diagnostics so
        // the same malformed query produces the same error in every process
        // and binding (important for reproducible validation reports).
        missing_in_products.sort_unstable();

        if !missing_in_products.is_empty() {
            return Err(format!(
                "map numbers in reactants missing from products: {:?}",
                missing_in_products
            ));
        }

        // All map numbers in products should be in reactants
        let mut undefined_in_reactants: Vec<u16> = self
            .product_maps
            .iter()
            .filter(|&m| !self.reactant_maps.contains(m))
            .copied()
            .collect();
        undefined_in_reactants.sort_unstable();

        if !undefined_in_reactants.is_empty() {
            return Err(format!(
                "map numbers in products not found in reactants: {:?}",
                undefined_in_reactants
            ));
        }

        Ok(())
    }
}

/// Error type for reaction query operations.
#[derive(Debug)]
pub enum ReactionQueryError {
    /// A configured batch-query limit was exceeded.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    /// Failed to parse a SMARTS pattern.
    SmartsParseError { smarts: String, source: String },
    /// Missing arrow delimiter in reaction SMARTS.
    MissingArrowDelimiter,
    /// Invalid agents section (middle part between > and >).
    InvalidAgentsSection,
    /// Map number inconsistency (e.g., :1 in reactants but not products).
    MapNumberMismatch { map_num: u16, message: String },
}

impl core::fmt::Display for ReactionQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "reaction query {resource} exceeds limit {limit} (got {actual})"
            ),
            Self::SmartsParseError { smarts, source } => {
                write!(f, "failed to parse SMARTS '{smarts}': {source}")
            }
            Self::MissingArrowDelimiter => {
                write!(f, "reaction SMARTS must contain '>' or '>>' delimiter")
            }
            Self::InvalidAgentsSection => {
                write!(f, "invalid agents section in reaction SMARTS")
            }
            Self::MapNumberMismatch { map_num, message } => {
                write!(f, "map number :{} error: {}", map_num, message)
            }
        }
    }
}

impl std::error::Error for ReactionQueryError {}

/// RDKit compatibility information and validation.
pub mod rdkit_compat {
    use super::*;

    /// Result of RDKit compatibility checking.
    #[derive(Clone, Debug)]
    pub struct RDKitCompatReport {
        /// Whether the reaction SMARTS is fully RDKit-compatible.
        pub is_compatible: bool,
        /// Warnings about potential differences in behavior.
        pub warnings: Vec<String>,
        /// RDKit format (legacy ">>" vs new ">").
        pub rdkit_format: RDKitFormat,
    }

    /// RDKit reaction SMARTS format variants.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub enum RDKitFormat {
        /// Legacy 2-part format: "reactants>>products"
        Legacy,
        /// New 3-part format: "reactants>agents>products"
        WithAgents,
    }

    /// Check RDKit compatibility of a reaction SMARTS pattern.
    ///
    /// Returns a `RDKitCompatReport` with compatibility status and any warnings.
    /// Both chematic and RDKit support the same reaction SMARTS format, but
    /// there may be subtle differences in pattern matching behavior.
    pub fn check_rdkit_compatibility(
        smarts: &str,
    ) -> Result<RDKitCompatReport, ReactionQueryError> {
        let mut warnings = Vec::new();

        // Detect format
        let seps = crate::reaction::reaction_separators(smarts);
        let rdkit_format = if seps.windows(2).any(|w| w[1] == w[0] + 1) {
            RDKitFormat::Legacy
        } else if !seps.is_empty() {
            RDKitFormat::WithAgents
        } else {
            return Err(ReactionQueryError::MissingArrowDelimiter);
        };

        // Try to parse to validate syntax
        let _pattern = super::parse_reaction_smarts(smarts)?;

        // Check for RDKit-specific concerns
        if smarts.contains("[#") {
            // Atomic number queries are fully supported
        }

        if smarts.contains("[*") {
            // Wildcard atoms
            warnings.push("wildcard atoms [*] may match behavior differently".to_string());
        }

        if smarts.contains("$") {
            // Recursive SMARTS are complex
            warnings.push("recursive SMARTS ($(...)) have complex matching behavior".to_string());
        }

        // Check for pipe-separated OR patterns
        let section_count = seps.len();
        let has_or_patterns = smarts.contains('|');
        if has_or_patterns && section_count == 0 {
            warnings
                .push("pipe-separated patterns (|) require at least one > delimiter".to_string());
        }

        Ok(RDKitCompatReport {
            is_compatible: warnings.is_empty(),
            warnings,
            rdkit_format,
        })
    }

    /// Convert chematic ReactionSmartsMatch result to RDKit-compatible format.
    ///
    /// Returns a summary suitable for comparison with RDKit results.
    pub fn match_to_rdkit_format(smarts_match: &ReactionSmartsMatch) -> RDKitMatchSummary {
        RDKitMatchSummary {
            matched: smarts_match.is_complete_match,
            reactant_count: smarts_match.reactant_matches.pattern_matches.len(),
            product_count: smarts_match.product_matches.pattern_matches.len(),
            reactant_molecules_matched: smarts_match
                .reactant_matches
                .pattern_matches
                .iter()
                .filter(|m| !m.is_empty())
                .count(),
            product_molecules_matched: smarts_match
                .product_matches
                .pattern_matches
                .iter()
                .filter(|m| !m.is_empty())
                .count(),
        }
    }

    /// Summary of match result in RDKit-comparable format.
    #[derive(Clone, Debug)]
    pub struct RDKitMatchSummary {
        /// Whether the reaction matched the pattern.
        pub matched: bool,
        /// Number of reactant patterns in the query.
        pub reactant_count: usize,
        /// Number of product patterns in the query.
        pub product_count: usize,
        /// How many reactant patterns matched.
        pub reactant_molecules_matched: usize,
        /// How many product patterns matched.
        pub product_molecules_matched: usize,
    }

    /// Validate that a reaction SMARTS follows RDKit conventions.
    pub fn validate_rdkit_conventions(smarts: &str) -> Result<(), String> {
        // Check for balanced arrow delimiters
        let arrow_count = crate::reaction::reaction_separators(smarts).len();
        if arrow_count == 0 {
            return Err("must contain at least one '>' or '>>' delimiter".to_string());
        }

        // Check for empty sections (three > in a row would be >>>)
        if smarts.contains(">>>") {
            return Err("invalid delimiter sequence >>>".to_string());
        }

        // Map numbers must be between 1 and 999 in RDKit
        // Extract all map numbers
        let map_nums: Vec<u16> = smarts
            .match_indices(':')
            .filter_map(|(i, _)| {
                let remainder = &smarts[i + 1..];
                let digits: String = remainder
                    .chars()
                    .take_while(|c| c.is_ascii_digit())
                    .collect();
                if !digits.is_empty() {
                    digits.parse().ok()
                } else {
                    None
                }
            })
            .collect();

        for map_num in map_nums {
            if map_num == 0 || map_num > 999 {
                return Err(format!(
                    "map number :{} outside RDKit range [1-999]",
                    map_num
                ));
            }
        }

        Ok(())
    }
}

/// Parse a reaction SMARTS pattern with agents section and atom mapping support.
///
/// Format: `"reactants>[agents]>products"` (new) or `"reactants>>products"` (legacy)
/// where each section is a pipe-separated (|) list of SMARTS alternatives;
/// each alternative may contain dot-separated (.) components that are matched
/// to distinct reaction molecules.
///
/// Validates atom map numbers for consistency between reactants and products.
///
/// # Example
///
/// ```ignore
/// // With agents section
/// let pattern = parse_reaction_smarts("[C:1][C:2]>[Pd]>[C:2][C:1]").unwrap();
/// // Legacy format
/// let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
/// ```
pub fn parse_reaction_smarts(
    smarts_str: &str,
) -> Result<ReactionSmartsPattern, ReactionQueryError> {
    // Sections split at the arrow separators (a dative `->` is a bond):
    // `reactants>>products`, `reactants>agents>products`, or a single `>`
    // read as `reactants>products`.
    let seps = crate::reaction::reaction_separators(smarts_str);
    let (reactant_strs, agent_strs, product_strs) = match seps[..] {
        [] => return Err(ReactionQueryError::MissingArrowDelimiter),
        [a] => (&smarts_str[..a], "", &smarts_str[a + 1..]),
        [a, b, ..] => (
            &smarts_str[..a],
            &smarts_str[a + 1..b],
            &smarts_str[b + 1..],
        ),
    };

    let reactant_components = parse_pattern_components(reactant_strs)?;
    let product_components = parse_pattern_components(product_strs)?;
    let agent_components = parse_pattern_components(agent_strs)?;
    let reactant_patterns = flatten_pattern_components(&reactant_components);
    let product_patterns = flatten_pattern_components(&product_components);
    let agent_patterns = flatten_pattern_components(&agent_components);

    let map_number_info = extract_map_numbers(smarts_str)?;

    Ok(ReactionSmartsPattern {
        reactant_patterns,
        reactant_components,
        agent_patterns,
        agent_components,
        product_patterns,
        product_components,
        map_number_info,
    })
}

/// Extracts and validates atom map numbers from a reaction SMARTS string.
fn extract_map_numbers(smarts_str: &str) -> Result<MapNumberInfo, ReactionQueryError> {
    let mut all_map_numbers = HashSet::new();
    let mut reactant_maps = HashSet::new();
    let mut agent_maps = HashSet::new();
    let mut product_maps = HashSet::new();

    // Find section boundaries (a dative `->` is a bond, not an arrow).
    let seps = crate::reaction::reaction_separators(smarts_str);
    let (reactant_end, agent_start, agent_end, product_start) = match seps[..] {
        [] => {
            return Ok(MapNumberInfo {
                all_map_numbers,
                reactant_maps,
                agent_maps,
                product_maps,
            });
        }
        [idx] => (idx, idx + 1, idx + 1, idx + 1),
        [first, second, ..] => (first, first + 1, second, second + 1),
    };

    // Extract maps from reactants
    let reactant_section = &smarts_str[..reactant_end];
    for map_num in extract_map_numbers_from_section(reactant_section) {
        reactant_maps.insert(map_num);
        all_map_numbers.insert(map_num);
    }

    // Extract maps from agents (if present)
    if agent_start < agent_end && agent_end <= smarts_str.len() {
        let agent_section = &smarts_str[agent_start..agent_end];
        if !agent_section.is_empty() {
            for map_num in extract_map_numbers_from_section(agent_section) {
                agent_maps.insert(map_num);
                all_map_numbers.insert(map_num);
            }
        }
    }

    // Extract maps from products
    if product_start < smarts_str.len() {
        let product_section = &smarts_str[product_start..];
        for map_num in extract_map_numbers_from_section(product_section) {
            product_maps.insert(map_num);
            all_map_numbers.insert(map_num);
        }
    }

    let info = MapNumberInfo {
        all_map_numbers,
        reactant_maps,
        agent_maps,
        product_maps,
    };

    // Validate consistency (map numbers must exist in both reactants and products)
    info.validate().map_err(|msg| {
        // Select the smallest offending map rather than an arbitrary element
        // of the HashSet. The message already identifies the whole mismatch;
        // this field must nevertheless be stable for machine consumers.
        let map_num = info
            .reactant_maps
            .symmetric_difference(&info.product_maps)
            .min()
            .copied()
            .unwrap_or(0);
        ReactionQueryError::MapNumberMismatch {
            map_num,
            message: msg,
        }
    })?;

    Ok(info)
}

/// Extract map numbers (format `:123`) from a SMARTS section string.
fn extract_map_numbers_from_section(smarts: &str) -> Vec<u16> {
    let mut map_numbers = Vec::new();
    let bytes = smarts.as_bytes();
    let mut bracket_depth: u32 = 0;
    let mut i = 0;

    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                bracket_depth += 1;
                i += 1;
            }
            b']' => {
                bracket_depth = bracket_depth.saturating_sub(1);
                i += 1;
            }
            b':' if bracket_depth > 0 => {
                // Atom map number: only valid inside a bracket atom `[…:N]`.
                let mut j = i + 1;
                let mut num_str = String::new();
                while j < bytes.len() && bytes[j].is_ascii_digit() {
                    num_str.push(bytes[j] as char);
                    j += 1;
                }
                if !num_str.is_empty()
                    && let Ok(num) = num_str.parse::<u16>()
                {
                    map_numbers.push(num);
                }
                i = j;
            }
            _ => {
                i += 1;
            }
        }
    }

    map_numbers
}

/// Parse pipe-separated alternatives and dot-separated components.
fn parse_pattern_components(side: &str) -> Result<Vec<Vec<QueryMolecule>>, ReactionQueryError> {
    if side.is_empty() {
        return Ok(Vec::new());
    }
    side.split('|')
        .filter(|p| !p.is_empty())
        .map(|alternative| {
            alternative
                .split('.')
                .filter(|component| !component.is_empty())
                .map(|component| {
                    parse_smarts(component).map_err(|e| ReactionQueryError::SmartsParseError {
                        smarts: component.to_string(),
                        source: e.to_string(),
                    })
                })
                .collect()
        })
        .collect()
}

fn flatten_pattern_components(components: &[Vec<QueryMolecule>]) -> Vec<QueryMolecule> {
    components
        .iter()
        .flat_map(|group| group.iter().cloned())
        .collect()
}

/// Parse a reaction SMARTS query with reactant and product patterns.
///
/// Format: `"reactant_smarts>>product_smarts"`
/// where each side is a pipe-separated (|) list of SMARTS alternatives and
/// each alternative may contain dot-separated (.) components.
///
/// Example: `"[C:1]([#6])[C:2]>>[C:1][C:2]"` matches reactions that break and reform C-C bonds.
pub fn parse_reaction_query(s: &str) -> Result<ReactionQuery, ReactionQueryError> {
    let pattern = parse_reaction_smarts(s)?;
    Ok(ReactionQuery {
        reactant_patterns: pattern.reactant_patterns,
        reactant_components: pattern.reactant_components,
        agent_patterns: pattern.agent_patterns,
        agent_components: pattern.agent_components,
        product_patterns: pattern.product_patterns,
        product_components: pattern.product_components,
    })
}

// Query reporting only retains one embedding per molecule. Asking the VF2
// engine for one match avoids allocating every symmetric embedding just to
// inspect `.first()` below.
fn first_match_config() -> MatchConfig {
    MatchConfig {
        max_matches: Some(1),
        ..MatchConfig::default()
    }
}

/// Return an embedding that satisfies both the SMARTS query and, when the
/// target molecule carries maps, the query's atom-map metadata.
///
/// Unmapped target reactions intentionally retain the historical behavior:
/// map labels in a query are metadata and cannot be checked against absent
/// target labels. Once target maps are present, accepting an embedding with a
/// different map would make reaction-side correspondence meaningless.
fn reaction_match_embedding(
    pattern: &QueryMolecule,
    mol: &chematic_core::Molecule,
) -> Option<FxHashMap<usize, chematic_core::AtomIdx>> {
    let query_has_maps = pattern.atoms.iter().any(|atom| atom.atom_map.is_some());
    let target_has_maps = mol.atoms().any(|(_, atom)| atom.atom_map.is_some());
    let mut config = first_match_config();
    if query_has_maps && target_has_maps {
        // The first VF2 embedding may use the right topology but the wrong
        // map labels. Inspect all bounded embeddings before failing closed.
        config.max_matches = Some(10_000);
    }
    find_matches_with_config(pattern, mol, &config)
        .into_iter()
        .find(|embedding| {
            if !query_has_maps || !target_has_maps {
                return true;
            }
            pattern
                .atoms
                .iter()
                .enumerate()
                .all(|(query_index, query_atom)| {
                    query_atom.atom_map.is_none_or(|map| {
                        embedding.get(&query_index).is_some_and(|target_index| {
                            mol.atom(*target_index).atom_map == Some(map)
                        })
                    })
                })
        })
}

/// Convert an embedding to the public atom-index order promised by
/// `MoleculeMatch`: query atom 0, query atom 1, and so on.  The embedding is
/// stored in a hash map because VF2 lookup is keyed by query index; exposing
/// `.values()` directly would make the result order an implementation detail.
fn ordered_atom_indices(embedding: &FxHashMap<usize, chematic_core::AtomIdx>) -> Vec<usize> {
    let mut entries: Vec<_> = embedding.iter().collect();
    entries.sort_unstable_by_key(|(query_index, _)| **query_index);
    entries
        .into_iter()
        .map(|(_, atom_index)| atom_index.0 as usize)
        .collect()
}

/// Match every disconnected component in one alternative to a distinct target
/// molecule. A single-component group is intentionally equivalent to the
/// historical per-pattern check.
fn match_distinct_components(
    components: &[QueryMolecule],
    molecules: &[chematic_core::Molecule],
    component_index: usize,
    used: &mut [bool],
) -> bool {
    if component_index == components.len() {
        return true;
    }
    molecules
        .iter()
        .enumerate()
        .any(|(molecule_index, molecule)| {
            !used[molecule_index]
                && reaction_match_embedding(&components[component_index], molecule).is_some()
                && {
                    used[molecule_index] = true;
                    let matched =
                        match_distinct_components(components, molecules, component_index + 1, used);
                    used[molecule_index] = false;
                    matched
                }
        })
}

fn match_component_groups(
    groups: &[Vec<QueryMolecule>],
    molecules: &[chematic_core::Molecule],
) -> bool {
    // Pipe-separated groups are alternatives for every reaction section.
    // Within one group, dot-separated components still require an injective
    // assignment to distinct target molecules. Keeping these levels separate
    // means `C|N` requires either C or N, while `C.N` requires both.
    groups.is_empty()
        || groups.iter().any(|components| {
            let mut used = vec![false; molecules.len()];
            match_distinct_components(components, molecules, 0, &mut used)
        })
}

/// Check if a reaction matches the given query pattern.
///
/// Returns `true` if:
/// - At least one reactant alternative matches its required components, AND
/// - At least one agent pattern matches when agent alternatives are present, AND
/// - At least one product alternative matches its required components
///
/// If the query has no patterns (empty reaction query), returns `true` (trivial match).
pub fn has_reaction_substructure_match(rxn: &Reaction, query: &ReactionQuery) -> bool {
    if !match_component_groups(&query.reactant_components, &rxn.reactants) {
        return false;
    }

    // Agent patterns are matched against the middle reaction section. An
    // empty agent section remains a valid legacy query and requires no match.
    if !query.agent_components.is_empty()
        && !query.agent_components.iter().any(|components| {
            let mut used = vec![false; rxn.agents.len()];
            match_distinct_components(components, &rxn.agents, 0, &mut used)
        })
    {
        return false;
    }

    match_component_groups(&query.product_components, &rxn.products)
}

/// Get detailed match information for a reaction against a query pattern.
///
/// Returns a `ReactionSmartsMatch` containing:
/// - All matched atoms for each reactant/product pattern
/// - Which molecules matched which patterns
/// - Overall match status
///
/// This is useful for understanding *how* a reaction matches, not just *whether* it matches.
pub fn get_reaction_smarts_matches(rxn: &Reaction, query: &ReactionQuery) -> ReactionSmartsMatch {
    // Collect all reactant pattern matches
    let mut reactant_pattern_matches = Vec::new();
    for (pattern_idx, pattern) in query.reactant_patterns.iter().enumerate() {
        let mut matches_for_pattern = Vec::new();
        for (mol_idx, mol) in rxn.reactants.iter().enumerate() {
            let atom_matches = reaction_match_embedding(pattern, mol);
            // Use first match if any exist (we only care that it matched)
            if let Some(first_match) = atom_matches.as_ref() {
                let atom_indices = ordered_atom_indices(first_match);
                matches_for_pattern.push(MoleculeMatch {
                    molecule_index: mol_idx,
                    pattern_index: pattern_idx,
                    atom_indices,
                });
            }
        }
        reactant_pattern_matches.push(matches_for_pattern);
    }

    // Collect all agent alternative matches. Keep every matching alternative
    // in the detailed result even though only one is required for completion.
    let mut agent_pattern_matches = Vec::new();
    for (pattern_idx, pattern) in query.agent_patterns.iter().enumerate() {
        let mut matches_for_pattern = Vec::new();
        for (mol_idx, mol) in rxn.agents.iter().enumerate() {
            if let Some(first_match) = reaction_match_embedding(pattern, mol).as_ref() {
                let atom_indices = ordered_atom_indices(first_match);
                matches_for_pattern.push(MoleculeMatch {
                    molecule_index: mol_idx,
                    pattern_index: pattern_idx,
                    atom_indices,
                });
            }
        }
        agent_pattern_matches.push(matches_for_pattern);
    }

    // Collect all product pattern matches
    let mut product_pattern_matches = Vec::new();
    for (pattern_idx, pattern) in query.product_patterns.iter().enumerate() {
        let mut matches_for_pattern = Vec::new();
        for (mol_idx, mol) in rxn.products.iter().enumerate() {
            let atom_matches = reaction_match_embedding(pattern, mol);
            // Use first match if any exist (we only care that it matched)
            if let Some(first_match) = atom_matches.as_ref() {
                let atom_indices = ordered_atom_indices(first_match);
                matches_for_pattern.push(MoleculeMatch {
                    molecule_index: mol_idx,
                    pattern_index: pattern_idx,
                    atom_indices,
                });
            }
        }
        product_pattern_matches.push(matches_for_pattern);
    }

    // Check if all patterns matched
    let all_reactants_matched = match_component_groups(&query.reactant_components, &rxn.reactants);
    let all_products_matched = match_component_groups(&query.product_components, &rxn.products);
    let all_agents_matched = query.agent_components.is_empty()
        || query.agent_components.iter().any(|components| {
            let mut used = vec![false; rxn.agents.len()];
            match_distinct_components(components, &rxn.agents, 0, &mut used)
        });
    let is_complete_match = all_reactants_matched && all_agents_matched && all_products_matched;

    ReactionSmartsMatch {
        reactant_matches: ReactantMatches {
            pattern_matches: reactant_pattern_matches,
        },
        agent_matches: AgentMatches {
            pattern_matches: agent_pattern_matches,
        },
        product_matches: ProductMatches {
            pattern_matches: product_pattern_matches,
        },
        reactants_satisfied: all_reactants_matched,
        products_satisfied: all_products_matched,
        agents_satisfied: all_agents_matched,
        is_complete_match,
    }
}

/// Batch query results with statistics and filtering.
#[derive(Clone, Debug)]
pub struct BatchQueryResults {
    /// Total reactions queried
    pub total_reactions: usize,
    /// Reactions that matched the query
    pub matching_reactions: usize,
    /// Match percentage (0-100)
    pub match_percentage: f64,
    /// Individual match results
    pub matches: Vec<(usize, bool)>, // (reaction_index, is_match)
}

/// Resource limits for batch reaction SMARTS queries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchQueryLimits {
    /// Maximum number of reactions processed in one batch.
    pub max_reactions: usize,
    /// Maximum number of patterns evaluated in one library batch.
    pub max_patterns: usize,
}

impl Default for BatchQueryLimits {
    fn default() -> Self {
        Self {
            max_reactions: 100_000,
            max_patterns: 10_000,
        }
    }
}

impl BatchQueryResults {
    /// Get indices of all matching reactions
    pub fn matching_indices(&self) -> Vec<usize> {
        self.matches
            .iter()
            .filter_map(|(idx, matched)| if *matched { Some(*idx) } else { None })
            .collect()
    }

    /// Get indices of all non-matching reactions
    pub fn non_matching_indices(&self) -> Vec<usize> {
        self.matches
            .iter()
            .filter_map(|(idx, matched)| if !*matched { Some(*idx) } else { None })
            .collect()
    }
}

/// Pattern library for fast reaction querying against multiple patterns
#[derive(Clone, Debug)]
pub struct ReactionPatternLibrary {
    /// Named reaction patterns
    pub patterns: std::collections::HashMap<String, ReactionSmartsPattern>,
}

impl ReactionPatternLibrary {
    /// Create a new empty pattern library
    pub fn new() -> Self {
        ReactionPatternLibrary {
            patterns: std::collections::HashMap::new(),
        }
    }

    /// Add a named pattern to the library
    pub fn add_pattern(&mut self, name: String, pattern: ReactionSmartsPattern) {
        self.patterns.insert(name, pattern);
    }

    /// Add a pattern from SMARTS string
    pub fn add_pattern_from_smarts(
        &mut self,
        name: String,
        smarts: &str,
    ) -> Result<(), ReactionQueryError> {
        let pattern = parse_reaction_smarts(smarts)?;
        self.add_pattern(name, pattern);
        Ok(())
    }

    /// Get number of patterns in library
    pub fn len(&self) -> usize {
        self.patterns.len()
    }

    /// Check if library is empty
    pub fn is_empty(&self) -> bool {
        self.patterns.is_empty()
    }

    /// Get pattern by name
    pub fn get(&self, name: &str) -> Option<&ReactionSmartsPattern> {
        self.patterns.get(name)
    }

    /// List all pattern names
    pub fn pattern_names(&self) -> Vec<String> {
        self.patterns.keys().cloned().collect()
    }
}

impl Default for ReactionPatternLibrary {
    fn default() -> Self {
        Self::new()
    }
}

/// Query a single reaction with a SMARTS pattern
pub fn query_reaction(
    rxn: &Reaction,
    smarts: &str,
) -> Result<ReactionSmartsMatch, ReactionQueryError> {
    let pattern = parse_reaction_smarts(smarts)?;
    // Convert ReactionSmartsPattern to ReactionQuery for matching
    let query = ReactionQuery {
        reactant_patterns: pattern.reactant_patterns,
        reactant_components: pattern.reactant_components,
        agent_patterns: pattern.agent_patterns,
        agent_components: pattern.agent_components,
        product_patterns: pattern.product_patterns,
        product_components: pattern.product_components,
    };
    Ok(get_reaction_smarts_matches(rxn, &query))
}

/// Batch query reactions against a single SMARTS pattern
pub fn batch_query_reactions(
    reactions: &[Reaction],
    smarts: &str,
) -> Result<BatchQueryResults, ReactionQueryError> {
    batch_query_reactions_with_limits(reactions, smarts, &BatchQueryLimits::default())
}

/// Batch query reactions with an explicit reaction-count limit.
pub fn batch_query_reactions_with_limits(
    reactions: &[Reaction],
    smarts: &str,
    limits: &BatchQueryLimits,
) -> Result<BatchQueryResults, ReactionQueryError> {
    if reactions.len() > limits.max_reactions {
        return Err(ReactionQueryError::ResourceLimit {
            resource: "batch reactions",
            actual: reactions.len(),
            limit: limits.max_reactions,
        });
    }
    let pattern = parse_reaction_smarts(smarts)?;
    let query = ReactionQuery {
        reactant_patterns: pattern.reactant_patterns,
        reactant_components: pattern.reactant_components,
        agent_patterns: pattern.agent_patterns,
        agent_components: pattern.agent_components,
        product_patterns: pattern.product_patterns,
        product_components: pattern.product_components,
    };

    let mut matches = Vec::new();
    let mut matching_count = 0;

    for (idx, rxn) in reactions.iter().enumerate() {
        let result = get_reaction_smarts_matches(rxn, &query);
        let is_match = result.is_complete_match;
        if is_match {
            matching_count += 1;
        }
        matches.push((idx, is_match));
    }

    let match_percentage = if reactions.is_empty() {
        0.0
    } else {
        (matching_count as f64 / reactions.len() as f64) * 100.0
    };

    Ok(BatchQueryResults {
        total_reactions: reactions.len(),
        matching_reactions: matching_count,
        match_percentage,
        matches,
    })
}

/// Batch query reactions against a pattern library
pub fn batch_query_with_library(
    reactions: &[Reaction],
    library: &ReactionPatternLibrary,
) -> std::collections::HashMap<String, BatchQueryResults> {
    let mut results = std::collections::HashMap::new();

    for (pattern_name, pattern) in &library.patterns {
        let query = ReactionQuery {
            reactant_patterns: pattern.reactant_patterns.clone(),
            reactant_components: pattern.reactant_components.clone(),
            agent_patterns: pattern.agent_patterns.clone(),
            agent_components: pattern.agent_components.clone(),
            product_patterns: pattern.product_patterns.clone(),
            product_components: pattern.product_components.clone(),
        };

        let mut matches = Vec::new();
        let mut matching_count = 0;

        for (idx, rxn) in reactions.iter().enumerate() {
            let result = get_reaction_smarts_matches(rxn, &query);
            let is_match = result.is_complete_match;
            if is_match {
                matching_count += 1;
            }
            matches.push((idx, is_match));
        }

        let match_percentage = if reactions.is_empty() {
            0.0
        } else {
            (matching_count as f64 / reactions.len() as f64) * 100.0
        };

        results.insert(
            pattern_name.clone(),
            BatchQueryResults {
                total_reactions: reactions.len(),
                matching_reactions: matching_count,
                match_percentage,
                matches,
            },
        );
    }

    results
}

/// Batch query a pattern library with explicit reaction and pattern limits.
pub fn batch_query_with_library_with_limits(
    reactions: &[Reaction],
    library: &ReactionPatternLibrary,
    limits: &BatchQueryLimits,
) -> Result<std::collections::HashMap<String, BatchQueryResults>, ReactionQueryError> {
    if reactions.len() > limits.max_reactions {
        return Err(ReactionQueryError::ResourceLimit {
            resource: "batch reactions",
            actual: reactions.len(),
            limit: limits.max_reactions,
        });
    }
    if library.patterns.len() > limits.max_patterns {
        return Err(ReactionQueryError::ResourceLimit {
            resource: "query patterns",
            actual: library.patterns.len(),
            limit: limits.max_patterns,
        });
    }

    let mut results = std::collections::HashMap::with_capacity(library.patterns.len());
    for (pattern_name, pattern) in &library.patterns {
        let query = ReactionQuery {
            reactant_patterns: pattern.reactant_patterns.clone(),
            reactant_components: pattern.reactant_components.clone(),
            agent_patterns: pattern.agent_patterns.clone(),
            agent_components: pattern.agent_components.clone(),
            product_patterns: pattern.product_patterns.clone(),
            product_components: pattern.product_components.clone(),
        };
        let mut matches = Vec::with_capacity(reactions.len());
        let mut matching_count = 0;
        for (idx, rxn) in reactions.iter().enumerate() {
            let is_match = get_reaction_smarts_matches(rxn, &query).is_complete_match;
            if is_match {
                matching_count += 1;
            }
            matches.push((idx, is_match));
        }
        let match_percentage = if reactions.is_empty() {
            0.0
        } else {
            (matching_count as f64 / reactions.len() as f64) * 100.0
        };
        results.insert(
            pattern_name.clone(),
            BatchQueryResults {
                total_reactions: reactions.len(),
                matching_reactions: matching_count,
                match_percentage,
                matches,
            },
        );
    }
    Ok(results)
}

#[cfg(test)]
mod tests;
