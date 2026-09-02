//! Reaction SMARTS querying for chemical reaction matching.

use chematic_smarts::{MatchConfig, QueryMolecule, find_matches_with_config, parse_smarts};
use std::collections::HashSet;

use crate::reaction::Reaction;

/// A reaction query consisting of reactant and product SMARTS patterns.
#[derive(Clone, Debug)]
pub struct ReactionQuery {
    /// SMARTS queries for reactant pattern matching.
    pub reactant_patterns: Vec<QueryMolecule>,
    /// SMARTS queries for product pattern matching.
    pub product_patterns: Vec<QueryMolecule>,
}

/// A reaction SMARTS pattern with atom mapping and agent support.
#[derive(Clone, Debug)]
pub struct ReactionSmartsPattern {
    /// SMARTS queries for reactant pattern matching.
    pub reactant_patterns: Vec<QueryMolecule>,
    /// SMARTS queries for agent pattern matching (optional middle section).
    pub agent_patterns: Vec<QueryMolecule>,
    /// SMARTS queries for product pattern matching.
    pub product_patterns: Vec<QueryMolecule>,
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

/// Detailed match information for a reaction against a SMARTS pattern.
#[derive(Clone, Debug)]
pub struct ReactionSmartsMatch {
    /// Reactant pattern matches (all patterns must have at least one match for overall match).
    pub reactant_matches: ReactantMatches,
    /// Product pattern matches (all patterns must have at least one match for overall match).
    pub product_matches: ProductMatches,
    /// Whether all patterns matched (true if reaction is valid against the query).
    pub is_complete_match: bool,
}

impl ReactionSmartsMatch {
    /// Check if all reactant patterns matched.
    pub fn all_reactants_matched(&self) -> bool {
        self.reactant_matches
            .pattern_matches
            .iter()
            .all(|m| !m.is_empty())
    }

    /// Check if all product patterns matched.
    pub fn all_products_matched(&self) -> bool {
        self.product_matches
            .pattern_matches
            .iter()
            .all(|m| !m.is_empty())
    }
}

impl MapNumberInfo {
    /// Check if all map numbers are consistent across reactants and products.
    pub fn validate(&self) -> Result<(), String> {
        // All map numbers in reactants must appear in products
        let missing_in_products: Vec<u16> = self
            .reactant_maps
            .iter()
            .filter(|&m| !self.product_maps.contains(m))
            .copied()
            .collect();

        if !missing_in_products.is_empty() {
            return Err(format!(
                "map numbers in reactants missing from products: {:?}",
                missing_in_products
            ));
        }

        // All map numbers in products should be in reactants
        let undefined_in_reactants: Vec<u16> = self
            .product_maps
            .iter()
            .filter(|&m| !self.reactant_maps.contains(m))
            .copied()
            .collect();

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
        let rdkit_format = if smarts.contains(">>") {
            RDKitFormat::Legacy
        } else if smarts.contains('>') {
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
        let section_count = smarts.matches('>').count();
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
        let arrow_count = smarts.matches('>').count();
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
/// where each section is a pipe-separated (|) list of SMARTS patterns.
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
    // Detect format by looking for ">>" or ">"
    let has_double_arrow = smarts_str.contains(">>");

    let (reactant_strs, agent_strs, product_strs) = if has_double_arrow {
        // Legacy format: reactants >> products
        let parts: Vec<&str> = smarts_str.splitn(2, ">>").collect();
        if parts.len() != 2 {
            return Err(ReactionQueryError::MissingArrowDelimiter);
        }
        (parts[0], "", parts[1])
    } else {
        // New format: reactants > agents > products
        let parts: Vec<&str> = smarts_str.splitn(3, '>').collect();
        if parts.len() < 2 {
            return Err(ReactionQueryError::MissingArrowDelimiter);
        }
        if parts.len() == 2 {
            // Only one '>' found, treat as reactants > products (no agents)
            (parts[0], "", parts[1])
        } else {
            // Two '>' found: reactants > agents > products
            (parts[0], parts[1], parts[2])
        }
    };

    let reactant_patterns = parse_patterns(reactant_strs)?;
    let product_patterns = parse_patterns(product_strs)?;
    let agent_patterns = if agent_strs.is_empty() {
        Vec::new()
    } else {
        parse_patterns(agent_strs)?
    };

    let map_number_info = extract_map_numbers(smarts_str)?;

    Ok(ReactionSmartsPattern {
        reactant_patterns,
        agent_patterns,
        product_patterns,
        map_number_info,
    })
}

/// Extracts and validates atom map numbers from a reaction SMARTS string.
fn extract_map_numbers(smarts_str: &str) -> Result<MapNumberInfo, ReactionQueryError> {
    let mut all_map_numbers = HashSet::new();
    let mut reactant_maps = HashSet::new();
    let mut agent_maps = HashSet::new();
    let mut product_maps = HashSet::new();

    // Find section boundaries
    let has_double_arrow = smarts_str.contains(">>");

    let (reactant_end, agent_start, agent_end, product_start) = if has_double_arrow {
        let idx = smarts_str.find(">>").unwrap();
        (idx, idx, idx, idx + 2)
    } else {
        let arrow_positions: Vec<_> = smarts_str.match_indices('>').collect();
        match arrow_positions.len() {
            0 => {
                return Ok(MapNumberInfo {
                    all_map_numbers,
                    reactant_maps,
                    agent_maps,
                    product_maps,
                });
            }
            1 => {
                let idx = arrow_positions[0].0;
                (idx, idx + 1, idx + 1, idx + 1)
            }
            _ => {
                let first = arrow_positions[0].0;
                let second = arrow_positions[1].0;
                (first, first + 1, second, second + 1)
            }
        }
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
        let map_num = info.all_map_numbers.iter().next().copied().unwrap_or(0);
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

/// Helper to parse pipe-separated SMARTS patterns.
fn parse_patterns(side: &str) -> Result<Vec<QueryMolecule>, ReactionQueryError> {
    if side.is_empty() {
        return Ok(Vec::new());
    }
    side.split('|')
        .filter(|p| !p.is_empty())
        .map(|p| {
            parse_smarts(p).map_err(|e| ReactionQueryError::SmartsParseError {
                smarts: p.to_string(),
                source: e.to_string(),
            })
        })
        .collect()
}

/// Parse a reaction SMARTS query with reactant and product patterns.
///
/// Format: `"reactant_smarts>>product_smarts"`
/// where each side is a pipe-separated (|) list of SMARTS patterns.
///
/// Example: `"[C:1]([#6])[C:2]>>[C:1][C:2]"` matches reactions that break and reform C-C bonds.
pub fn parse_reaction_query(s: &str) -> Result<ReactionQuery, ReactionQueryError> {
    let parts: Vec<&str> = s.splitn(2, ">>").collect();
    if parts.len() != 2 {
        return Err(ReactionQueryError::SmartsParseError {
            smarts: s.to_string(),
            source: "reaction query must contain '>>'".to_string(),
        });
    }

    Ok(ReactionQuery {
        reactant_patterns: parse_patterns(parts[0])?,
        product_patterns: parse_patterns(parts[1])?,
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

/// Check if a reaction matches the given query pattern.
///
/// Returns `true` if:
/// - All reactant patterns match at least one reactant molecule, AND
/// - All product patterns match at least one product molecule
///
/// If the query has no patterns (empty reaction query), returns `true` (trivial match).
pub fn has_reaction_substructure_match(rxn: &Reaction, query: &ReactionQuery) -> bool {
    // Check if all reactant patterns are satisfied
    for pattern in &query.reactant_patterns {
        let mut matched = false;
        for mol in &rxn.reactants {
            if !find_matches_with_config(pattern, mol, &first_match_config()).is_empty() {
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }

    // Check if all product patterns are satisfied
    for pattern in &query.product_patterns {
        let mut matched = false;
        for mol in &rxn.products {
            if !find_matches_with_config(pattern, mol, &first_match_config()).is_empty() {
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }

    true
}

/// Get detailed match information for a reaction against a query pattern.
///
/// Returns a `ReactionSmartsMatch` containing:
/// - All matched atoms for each pattern
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
            let atom_matches = find_matches_with_config(pattern, mol, &first_match_config());
            // Use first match if any exist (we only care that it matched)
            if let Some(first_match) = atom_matches.first() {
                let atom_indices: Vec<usize> = first_match
                    .values()
                    .map(|atom_idx| atom_idx.0 as usize)
                    .collect();
                matches_for_pattern.push(MoleculeMatch {
                    molecule_index: mol_idx,
                    pattern_index: pattern_idx,
                    atom_indices,
                });
            }
        }
        reactant_pattern_matches.push(matches_for_pattern);
    }

    // Collect all product pattern matches
    let mut product_pattern_matches = Vec::new();
    for (pattern_idx, pattern) in query.product_patterns.iter().enumerate() {
        let mut matches_for_pattern = Vec::new();
        for (mol_idx, mol) in rxn.products.iter().enumerate() {
            let atom_matches = find_matches_with_config(pattern, mol, &first_match_config());
            // Use first match if any exist (we only care that it matched)
            if let Some(first_match) = atom_matches.first() {
                let atom_indices: Vec<usize> = first_match
                    .values()
                    .map(|atom_idx| atom_idx.0 as usize)
                    .collect();
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
    let all_reactants_matched = reactant_pattern_matches.iter().all(|m| !m.is_empty());
    let all_products_matched = product_pattern_matches.iter().all(|m| !m.is_empty());
    let is_complete_match = all_reactants_matched && all_products_matched;

    ReactionSmartsMatch {
        reactant_matches: ReactantMatches {
            pattern_matches: reactant_pattern_matches,
        },
        product_matches: ProductMatches {
            pattern_matches: product_pattern_matches,
        },
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
        product_patterns: pattern.product_patterns,
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
        product_patterns: pattern.product_patterns,
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
            product_patterns: pattern.product_patterns.clone(),
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
            product_patterns: pattern.product_patterns.clone(),
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
mod tests {
    use super::*;
    use crate::parse_reaction;

    fn rxn(s: &str) -> Reaction {
        crate::reaction::parse_reaction(s).unwrap()
    }

    #[test]
    fn test_parse_reaction_query_basic() {
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        assert_eq!(query.reactant_patterns.len(), 1);
        assert_eq!(query.product_patterns.len(), 1);
    }

    #[test]
    fn test_parse_reaction_query_multiple_patterns() {
        let query = parse_reaction_query("[#6]|[#7]>>[#8]|[#9]").unwrap();
        assert_eq!(query.reactant_patterns.len(), 2);
        assert_eq!(query.product_patterns.len(), 2);
    }

    #[test]
    fn test_has_reaction_substructure_match_simple() {
        // Reaction: ethane to ethane (trivial)
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        assert!(has_reaction_substructure_match(&rxn, &query));
    }

    #[test]
    fn test_has_reaction_substructure_match_no_match() {
        // Reaction: ethane to ethane
        let rxn = rxn("CC>>CC");
        // Query: looking for nitrogen (not present)
        let query = parse_reaction_query("[#7]>>[#7]").unwrap();
        assert!(!has_reaction_substructure_match(&rxn, &query));
    }

    #[test]
    fn test_has_reaction_substructure_match_product_mismatch() {
        // Reaction: ethane to methane
        let rxn = rxn("CC>>C");
        // Query: looking for ethane in products (not present)
        let query = parse_reaction_query("[#6]>>CC").unwrap();
        assert!(!has_reaction_substructure_match(&rxn, &query));
    }

    #[test]
    fn test_has_reaction_substructure_match_empty_query() {
        // Empty query should match any reaction (trivial)
        let rxn = rxn("CC>>C");
        let query = parse_reaction_query(">>").unwrap();
        assert!(has_reaction_substructure_match(&rxn, &query));
    }

    // ===== Phase 1: Agents Section Support =====

    #[test]
    fn test_parse_reaction_smarts_legacy_format() {
        // Legacy 2-part format should still work
        let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
        assert_eq!(pattern.reactant_patterns.len(), 1);
        assert_eq!(pattern.agent_patterns.len(), 0);
        assert_eq!(pattern.product_patterns.len(), 1);
    }

    #[test]
    fn test_parse_reaction_smarts_with_agents() {
        // New 3-part format with agents
        let pattern = parse_reaction_smarts("[C:1][C:2]>[Pd]>[C:2][C:1]").unwrap();
        assert_eq!(pattern.reactant_patterns.len(), 1);
        assert_eq!(pattern.agent_patterns.len(), 1);
        assert_eq!(pattern.product_patterns.len(), 1);
    }

    #[test]
    fn test_parse_reaction_smarts_empty_agents() {
        // Format: reactants > (empty) > products
        let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
        assert_eq!(pattern.agent_patterns.len(), 0);
    }

    #[test]
    fn test_parse_reaction_smarts_multiple_agent_patterns() {
        // Multiple agent patterns with pipe-separated OR logic
        let pattern = parse_reaction_smarts("[C:1]>[Pd]|[Ni]>[C:1]").unwrap();
        assert_eq!(pattern.agent_patterns.len(), 2);
    }

    #[test]
    fn test_extract_map_numbers_basic() {
        let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
        assert!(pattern.map_number_info.reactant_maps.contains(&1));
        assert!(pattern.map_number_info.reactant_maps.contains(&2));
        assert!(pattern.map_number_info.product_maps.contains(&1));
        assert!(pattern.map_number_info.product_maps.contains(&2));
    }

    #[test]
    fn test_extract_map_numbers_with_agents() {
        let pattern = parse_reaction_smarts("[C:1]>[Pd:3]>[C:1]").unwrap();
        assert!(pattern.map_number_info.reactant_maps.contains(&1));
        assert!(pattern.map_number_info.agent_maps.contains(&3));
        assert!(pattern.map_number_info.product_maps.contains(&1));
    }

    #[test]
    fn test_map_number_validation_missing_in_products() {
        // Map number :1 in reactants but not in products should fail
        let result = parse_reaction_smarts("[C:1][C:2]>>[C:2]");
        assert!(result.is_err());
        if let Err(ReactionQueryError::MapNumberMismatch { message, .. }) = result {
            assert!(message.contains("missing from products"));
        } else {
            panic!("expected MapNumberMismatch error");
        }
    }

    #[test]
    fn test_map_number_validation_undefined_in_reactants() {
        // Map number :1 in products but not in reactants should fail
        let result = parse_reaction_smarts("[C:2]>>[C:1][C:2]");
        assert!(result.is_err());
        if let Err(ReactionQueryError::MapNumberMismatch { message, .. }) = result {
            assert!(message.contains("not found in reactants"));
        } else {
            panic!("expected MapNumberMismatch error");
        }
    }

    #[test]
    fn test_map_number_validation_consistent() {
        // Valid: all map numbers present in both reactants and products
        let pattern = parse_reaction_smarts("[C:1][C:2][C:3]>>[C:3][C:1][C:2]").unwrap();
        assert_eq!(pattern.map_number_info.reactant_maps.len(), 3);
        assert_eq!(pattern.map_number_info.product_maps.len(), 3);
    }

    #[test]
    fn test_map_numbers_no_maps() {
        // Pattern with no map numbers is valid
        let pattern = parse_reaction_smarts("[#6]>>[#6]").unwrap();
        assert_eq!(pattern.map_number_info.all_map_numbers.len(), 0);
    }

    #[test]
    fn test_map_numbers_large_numbers() {
        // Handle large map numbers (e.g., :100, :999)
        let pattern = parse_reaction_smarts("[C:100][N:999]>>[N:999][C:100]").unwrap();
        assert!(pattern.map_number_info.reactant_maps.contains(&100));
        assert!(pattern.map_number_info.reactant_maps.contains(&999));
    }

    #[test]
    fn test_parse_reaction_smarts_missing_delimiter() {
        // No '>' or '>>' delimiter should error
        let result = parse_reaction_smarts("[C:1][C:2][C:1][C:2]");
        assert!(result.is_err());
    }

    #[test]
    fn test_agents_section_in_map_numbers() {
        // Agents can have their own map numbers (independent of reactant/product validation)
        let pattern = parse_reaction_smarts("[C:1]>[Pd:99]>[C:1]").unwrap();
        assert!(pattern.map_number_info.agent_maps.contains(&99));
        // Agent map :99 is not required in reactants/products
        assert!(!pattern.map_number_info.reactant_maps.contains(&99));
    }

    // ===== Phase 2: Detailed Match Information =====

    #[test]
    fn test_get_reaction_smarts_matches_basic() {
        // Simple match: ethane to ethane with carbon pattern
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Should have one reactant pattern and one product pattern
        assert_eq!(matches.reactant_matches.pattern_matches.len(), 1);
        assert_eq!(matches.product_matches.pattern_matches.len(), 1);

        // Both patterns should match
        assert!(matches.all_reactants_matched());
        assert!(matches.all_products_matched());
        assert!(matches.is_complete_match);
    }

    #[test]
    fn test_get_reaction_smarts_matches_multiple_reactants() {
        // Multiple reactants: C + C >> CC
        let rxn = rxn("C.C>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Reactant pattern should match both molecules
        assert_eq!(matches.reactant_matches.pattern_matches[0].len(), 2);
        // Product pattern should match the combined molecule
        assert_eq!(matches.product_matches.pattern_matches[0].len(), 1);
    }

    #[test]
    fn test_get_reaction_smarts_matches_incomplete() {
        // Query for nitrogen (not present) should have incomplete match
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#7]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Reactant pattern should not match
        assert!(!matches.all_reactants_matched());
        assert!(!matches.is_complete_match);
    }

    #[test]
    fn test_get_reaction_smarts_matches_empty_patterns() {
        // Empty query should match anything
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query(">>").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Empty patterns should be considered complete matches
        assert!(matches.is_complete_match);
    }

    #[test]
    fn test_molecule_match_has_correct_indices() {
        // Verify that matched atoms have correct indices
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Reactant carbon pattern matches (first match = first atom)
        let reactant_matches = &matches.reactant_matches.pattern_matches[0];
        assert_eq!(reactant_matches.len(), 1);
        assert_eq!(reactant_matches[0].molecule_index, 0);
        assert!(!reactant_matches[0].atom_indices.is_empty()); // At least one carbon matched

        // Product carbon pattern matches
        let product_matches = &matches.product_matches.pattern_matches[0];
        assert_eq!(product_matches.len(), 1);
        assert_eq!(product_matches[0].molecule_index, 0);
        assert!(!product_matches[0].atom_indices.is_empty()); // At least one carbon matched
    }

    #[test]
    fn test_reactant_matches_get_pattern_matches() {
        // Test the get_pattern_matches helper method
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Should be able to retrieve matches for a specific pattern
        let reactant_pattern_0 = matches.reactant_matches.get_pattern_matches(0);
        assert!(reactant_pattern_0.is_some());
        assert_eq!(reactant_pattern_0.unwrap().len(), 1);

        // Out of bounds pattern should return None
        let reactant_pattern_1 = matches.reactant_matches.get_pattern_matches(1);
        assert!(reactant_pattern_1.is_none());
    }

    #[test]
    fn test_product_matches_get_pattern_matches() {
        // Test the get_pattern_matches helper method for products
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        let product_pattern_0 = matches.product_matches.get_pattern_matches(0);
        assert!(product_pattern_0.is_some());
        assert_eq!(product_pattern_0.unwrap().len(), 1);
    }

    #[test]
    fn test_pattern_matched_helper() {
        // Test the pattern_matched helper method
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#7]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Reactant pattern should match
        assert!(matches.reactant_matches.pattern_matched(0));

        // Product pattern should not match (looking for nitrogen)
        assert!(!matches.product_matches.pattern_matched(0));
    }

    #[test]
    fn test_multiple_patterns_or_logic() {
        // Test with multiple patterns (OR logic)
        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]|[#7]>>[#6]|[#8]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Should have two patterns each (carbon|nitrogen and carbon|oxygen)
        assert_eq!(matches.reactant_matches.pattern_matches.len(), 2);
        assert_eq!(matches.product_matches.pattern_matches.len(), 2);

        // Both reactant patterns should match (we have carbons)
        assert!(matches.reactant_matches.pattern_matched(0)); // [#6]
        assert!(!matches.reactant_matches.pattern_matched(1)); // [#7] - no nitrogen

        // Both product patterns should match
        assert!(matches.product_matches.pattern_matched(0)); // [#6]
        assert!(!matches.product_matches.pattern_matched(1)); // [#8] - no oxygen

        // Overall match should be false (nitrogen pattern in reactants didn't match)
        assert!(!matches.is_complete_match);
    }

    #[test]
    fn test_complex_reaction_matches() {
        // Test with more complex reaction
        let rxn = rxn("CC(C)>>CC=C");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        // Both should match
        assert!(matches.all_reactants_matched());
        assert!(matches.all_products_matched());
        assert!(matches.is_complete_match);

        // Reactant has one molecule matching the pattern
        assert_eq!(matches.reactant_matches.pattern_matches[0].len(), 1);
        assert!(
            !matches.reactant_matches.pattern_matches[0][0]
                .atom_indices
                .is_empty()
        );

        // Product has one molecule matching the pattern
        assert_eq!(matches.product_matches.pattern_matches[0].len(), 1);
        assert!(
            !matches.product_matches.pattern_matches[0][0]
                .atom_indices
                .is_empty()
        );
    }

    // ===== Phase 3: RDKit Compatibility Tooling =====

    #[test]
    fn test_rdkit_compat_check_legacy_format() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[C:1][C:2]>>[C:2][C:1]").unwrap();
        assert!(result.is_compatible);
        assert_eq!(result.rdkit_format, RDKitFormat::Legacy);
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_rdkit_compat_check_agents_format() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[C:1][C:2]>[Pd]>[C:2][C:1]").unwrap();
        assert!(result.is_compatible);
        assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_rdkit_compat_wildcard_warning() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[*:1]>>[*:1]").unwrap();
        assert!(!result.is_compatible);
        assert!(result.warnings.iter().any(|w| w.contains("wildcard")));
    }

    #[test]
    fn test_rdkit_compat_recursive_smarts_warning() {
        use crate::query::rdkit_compat::*;

        // Use simpler recursive SMARTS pattern
        let result = check_rdkit_compatibility("[#6;$([#6]~[#8])]>>[#6]").unwrap();
        assert!(!result.is_compatible);
        assert!(result.warnings.iter().any(|w| w.contains("recursive")));
    }

    #[test]
    fn test_rdkit_compat_missing_delimiter() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[C:1][C:2][C:1][C:2]");
        assert!(result.is_err());
    }

    #[test]
    fn test_rdkit_format_detection_legacy() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("C>>C").unwrap();
        assert_eq!(result.rdkit_format, RDKitFormat::Legacy);
    }

    #[test]
    fn test_rdkit_format_detection_agents() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("C>[Pd]>C").unwrap();
        assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
    }

    #[test]
    fn test_rdkit_match_summary_complete() {
        use crate::query::rdkit_compat::*;

        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#6]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        let summary = match_to_rdkit_format(&matches);
        assert!(summary.matched);
        assert_eq!(summary.reactant_count, 1);
        assert_eq!(summary.product_count, 1);
        assert_eq!(summary.reactant_molecules_matched, 1);
        assert_eq!(summary.product_molecules_matched, 1);
    }

    #[test]
    fn test_rdkit_match_summary_incomplete() {
        use crate::query::rdkit_compat::*;

        let rxn = rxn("CC>>CC");
        let query = parse_reaction_query("[#7]>>[#6]").unwrap();
        let matches = get_reaction_smarts_matches(&rxn, &query);

        let summary = match_to_rdkit_format(&matches);
        assert!(!summary.matched);
        assert_eq!(summary.reactant_molecules_matched, 0);
    }

    #[test]
    fn test_rdkit_validate_conventions_valid() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:1][C:2]>>[C:2][C:1]");
        assert!(result.is_ok());
    }

    #[test]
    fn test_rdkit_validate_conventions_no_delimiter() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:1][C:2][C:1][C:2]");
        assert!(result.is_err());
    }

    #[test]
    fn test_rdkit_validate_conventions_triple_arrow() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:1]>>>[C:1]");
        assert!(result.is_err());
    }

    #[test]
    fn test_rdkit_validate_conventions_map_number_zero() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:0]>>[C:0]");
        assert!(result.is_err());
        if let Err(msg) = result {
            assert!(msg.contains("outside RDKit range"));
        }
    }

    #[test]
    fn test_rdkit_validate_conventions_map_number_out_of_range() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:1000]>>[C:1000]");
        assert!(result.is_err());
        if let Err(msg) = result {
            assert!(msg.contains("outside RDKit range"));
        }
    }

    #[test]
    fn test_rdkit_validate_conventions_map_number_valid() {
        use crate::query::rdkit_compat::*;

        let result = validate_rdkit_conventions("[C:1][N:999]>>[N:999][C:1]");
        assert!(result.is_ok());
    }

    #[test]
    fn test_rdkit_compat_agents_pipes() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[C:1]>[Pd]|[Ni]>[C:1]").unwrap();
        assert!(result.is_compatible);
        assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
    }

    #[test]
    fn test_rdkit_compat_multiple_patterns() {
        use crate::query::rdkit_compat::*;

        let result = check_rdkit_compatibility("[#6]|[#7]>>[#8]|[#9]").unwrap();
        assert!(result.is_compatible);
    }

    // B7 Enhancement Tests: Batch querying and pattern libraries

    #[test]
    fn test_batch_query_reactions_empty() {
        let reactions: Vec<Reaction> = vec![];
        let smarts = "[C:1]>>[C:1]";

        let result = batch_query_reactions(&reactions, smarts);
        assert!(result.is_ok());

        let batch = result.unwrap();
        assert_eq!(batch.total_reactions, 0);
        assert_eq!(batch.matching_reactions, 0);
        assert_eq!(batch.match_percentage, 0.0);
    }

    #[test]
    fn test_batch_query_limits_reject_oversized_batch() {
        let reactions = vec![parse_reaction("C>>C").unwrap()];
        let err = batch_query_reactions_with_limits(
            &reactions,
            "[C:1]>>[C:1]",
            &BatchQueryLimits {
                max_reactions: 0,
                ..Default::default()
            },
        );
        assert!(matches!(
            err,
            Err(ReactionQueryError::ResourceLimit {
                resource: "batch reactions",
                ..
            })
        ));
    }

    #[test]
    fn test_batch_query_reactions_single_match() {
        let rxn = parse_reaction("C>>C").unwrap();
        let reactions = vec![rxn];
        let smarts = "[C:1]>>[C:1]";

        let result = batch_query_reactions(&reactions, smarts);
        assert!(result.is_ok());

        let batch = result.unwrap();
        assert_eq!(batch.total_reactions, 1);
        assert!(batch.matching_reactions > 0);
    }

    #[test]
    fn test_batch_query_reactions_multiple() {
        let rxn1 = parse_reaction("C>>C").unwrap();
        let rxn2 = parse_reaction("CC>>C").unwrap();
        let rxn3 = parse_reaction("CCC>>CC").unwrap();

        let reactions = vec![rxn1, rxn2, rxn3];
        let smarts = "[C:1]>>[C:1]";

        let result = batch_query_reactions(&reactions, smarts);
        assert!(result.is_ok());

        let batch = result.unwrap();
        assert_eq!(batch.total_reactions, 3);
        assert_eq!(batch.matches.len(), 3);
    }

    #[test]
    fn test_batch_query_results_indices() {
        let rxn1 = parse_reaction("C>>C").unwrap();
        let rxn2 = parse_reaction("CC>>C").unwrap();
        let rxn3 = parse_reaction("CCC>>CC").unwrap();

        let reactions = vec![rxn1, rxn2, rxn3];
        let smarts = "[C:1]>>[C:1]";

        let batch = batch_query_reactions(&reactions, smarts).unwrap();

        // Both matching and non-matching indices should return valid results
        let matching = batch.matching_indices();
        let non_matching = batch.non_matching_indices();

        assert_eq!(matching.len() + non_matching.len(), batch.total_reactions);
    }

    #[test]
    fn test_reaction_pattern_library_new() {
        let library = ReactionPatternLibrary::new();
        assert!(library.is_empty());
        assert_eq!(library.len(), 0);
    }

    #[test]
    fn test_reaction_pattern_library_add_pattern() {
        let mut library = ReactionPatternLibrary::new();

        let pattern = parse_reaction_smarts("[C:1]>>[C:1]").unwrap();
        library.add_pattern("simple_carbon".to_string(), pattern);

        assert_eq!(library.len(), 1);
        assert!(!library.is_empty());
        assert!(library.get("simple_carbon").is_some());
    }

    #[test]
    fn test_reaction_pattern_library_add_from_smarts() {
        let mut library = ReactionPatternLibrary::new();

        let result = library.add_pattern_from_smarts(
            "acylation".to_string(),
            "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[N:3]",
        );
        assert!(result.is_ok());
        assert_eq!(library.len(), 1);
    }

    #[test]
    fn test_reaction_pattern_library_multiple_patterns() {
        let mut library = ReactionPatternLibrary::new();

        let _ = library.add_pattern_from_smarts("reaction1".to_string(), "[C:1]>>[C:1]");
        let _ = library.add_pattern_from_smarts("reaction2".to_string(), "[N:1]>>[N:1]");
        let _ = library.add_pattern_from_smarts("reaction3".to_string(), "[O:1]>>[O:1]");

        assert_eq!(library.len(), 3);

        let names = library.pattern_names();
        assert_eq!(names.len(), 3);
        assert!(names.contains(&"reaction1".to_string()));
    }

    #[test]
    fn test_batch_query_with_library() {
        let mut library = ReactionPatternLibrary::new();
        let _ = library.add_pattern_from_smarts("carbon_only".to_string(), "[C:1]>>[C:1]");

        let rxn1 = parse_reaction("C>>C").unwrap();
        let rxn2 = parse_reaction("CC>>C").unwrap();
        let reactions = vec![rxn1, rxn2];

        let results = batch_query_with_library(&reactions, &library);

        assert!(!results.is_empty());
        assert!(results.contains_key("carbon_only"));

        let batch = &results["carbon_only"];
        assert_eq!(batch.total_reactions, 2);
    }

    #[test]
    fn test_query_reaction_simple() {
        let rxn = parse_reaction("C>>C").unwrap();
        let smarts = "[C:1]>>[C:1]";

        let result = query_reaction(&rxn, smarts);
        assert!(result.is_ok());

        let match_result = result.unwrap();
        // Check that it produced a valid result structure
        assert!(!match_result.reactant_matches.pattern_matches.is_empty());
    }

    #[test]
    fn test_batch_query_match_percentage() {
        let rxn1 = parse_reaction("C>>C").unwrap();
        let rxn2 = parse_reaction("CC>>C").unwrap();
        let rxn3 = parse_reaction("CCC>>CC").unwrap();

        let reactions = vec![rxn1, rxn2, rxn3];
        let smarts = "[C:1]>>[C:1]";

        let batch = batch_query_reactions(&reactions, smarts).unwrap();

        // Match percentage should be between 0 and 100
        assert!(batch.match_percentage >= 0.0);
        assert!(batch.match_percentage <= 100.0);

        // Verify calculation
        let expected_percentage =
            (batch.matching_reactions as f64 / batch.total_reactions as f64) * 100.0;
        assert!((batch.match_percentage - expected_percentage).abs() < 0.01);
    }

    #[test]
    fn test_reaction_pattern_library_default() {
        let library = ReactionPatternLibrary::default();
        assert!(library.is_empty());
    }
}
