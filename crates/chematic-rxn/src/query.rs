//! Reaction SMARTS querying for chemical reaction matching.

use chematic_smarts::{find_matches, parse_smarts, QueryMolecule};

use crate::reaction::Reaction;

/// A reaction query consisting of reactant and product SMARTS patterns.
#[derive(Clone, Debug)]
pub struct ReactionQuery {
    /// SMARTS queries for reactant pattern matching.
    pub reactant_patterns: Vec<QueryMolecule>,
    /// SMARTS queries for product pattern matching.
    pub product_patterns: Vec<QueryMolecule>,
}

/// Error type for reaction query operations.
#[derive(Debug)]
pub enum ReactionQueryError {
    /// Failed to parse a SMARTS pattern.
    SmartsParseError { smarts: String, source: String },
}

impl core::fmt::Display for ReactionQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SmartsParseError { smarts, source } => {
                write!(f, "failed to parse SMARTS '{smarts}': {source}")
            }
        }
    }
}

impl std::error::Error for ReactionQueryError {}

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

    let parse_patterns = |side: &str| -> Result<Vec<QueryMolecule>, ReactionQueryError> {
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
    };

    Ok(ReactionQuery {
        reactant_patterns: parse_patterns(parts[0])?,
        product_patterns: parse_patterns(parts[1])?,
    })
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
            if !find_matches(pattern, mol).is_empty() {
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
            if !find_matches(pattern, mol).is_empty() {
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
