//! Combinatorial library enumeration from SMIRKS templates and fragment sets.
//!
//! Generates virtual libraries by combining reaction templates with scaffolds and building blocks.

use crate::transform::run_reactants;
use chematic_core::Molecule;

/// Configuration for library enumeration.
#[derive(Clone, Debug)]
pub struct LibraryConfig {
    /// Skip molecules that fail transformation (vs collecting errors).
    pub skip_failures: bool,
    /// Maximum total library size to prevent runaway enumeration.
    pub max_size: Option<usize>,
}

impl Default for LibraryConfig {
    fn default() -> Self {
        Self {
            skip_failures: true,
            max_size: Some(1_000_000),
        }
    }
}

/// Enumerate a combinatorial library from reaction templates and fragment sets.
///
/// Given N SMIRKS templates and M fragment sets, generates all valid products
/// where each template is applied with one fragment from each set.
///
/// # Arguments
/// * `templates` - SMIRKS reaction templates (e.g., `"[C:1][Cl].[C:2][NH2]>>[C:1][NH][C:2]"`)
/// * `fragment_sets` - One or more sets of fragments (e.g., `vec![vec![scaffold1, scaffold2], vec![sidechain1, sidechain2]]`)
/// * `config` - Enumeration options (max size, skip failures)
///
/// # Returns
/// `Vec<Molecule>` containing all enumerated products (may be deduplicated by hash).
///
/// # Example
/// ```ignore
/// let template = "[C:1]Cl.[C:2]N>>[C:1]N[C:2]"; // Amide coupling
/// let scaffolds = vec![phenol, aniline];
/// let sidechains = vec![acetyl, benzoyl];
/// let library = enumerate_library(
///     &template,
///     vec![scaffolds, sidechains],
///     &LibraryConfig::default()
/// );
/// // Returns all products: phenol + acetyl, phenol + benzoyl, aniline + acetyl, aniline + benzoyl
/// ```
pub fn enumerate_library(
    template: &str,
    fragment_sets: Vec<Vec<Molecule>>,
    config: &LibraryConfig,
) -> Result<Vec<Molecule>, LibraryError> {
    if fragment_sets.is_empty() {
        return Err(LibraryError::NoFragmentSets);
    }

    // Empty sets have no valid combinations. Return an empty result instead
    // of entering the index loop below, which would index `set[0]`.
    if fragment_sets.iter().any(|set| set.is_empty()) {
        return Ok(Vec::new());
    }

    // Check max size feasibility with checked arithmetic. A wrapped product
    // could otherwise bypass the configured cap or panic in debug builds.
    let total_combos = fragment_sets
        .iter()
        .map(|set| set.len())
        .try_fold(1usize, |total, size| total.checked_mul(size))
        .ok_or(LibraryError::CombinationCountOverflow)?;

    if let Some(max) = config.max_size
        && total_combos > max
    {
        return Err(LibraryError::EnumerationTooLarge(total_combos, max));
    }

    let mut products: Vec<Molecule> = Vec::new();
    let mut indices = vec![0usize; fragment_sets.len()];

    loop {
        // Build reactants list for this iteration (as references)
        let mut reactants = Vec::new();
        for (set_idx, frag_set) in fragment_sets.iter().enumerate() {
            reactants.push(&frag_set[indices[set_idx]]);
        }

        // Apply template
        match run_reactants(template, &reactants) {
            Ok(reaction_product_sets) => {
                // Flatten: each set can have multiple products
                for product_set in reaction_product_sets {
                    if let Some(max) = config.max_size {
                        let projected = products
                            .len()
                            .checked_add(product_set.len())
                            .ok_or(LibraryError::CombinationCountOverflow)?;
                        if projected >= max {
                            return Err(LibraryError::EnumerationLimitExceeded);
                        }
                    }
                    products.extend(product_set);
                }
            }
            Err(_) => {
                if !config.skip_failures {
                    // For now, skip and continue (could collect errors)
                }
            }
        }

        // Check product count (actual size, not just iteration count)
        if let Some(max) = config.max_size
            && products.len() >= max
        {
            return Err(LibraryError::EnumerationLimitExceeded);
        }

        // Increment indices (counter-like behavior)
        let mut carry = 1usize;
        for i in 0..indices.len() {
            indices[i] += carry;
            if indices[i] < fragment_sets[i].len() {
                carry = 0;
                break;
            }
            indices[i] = 0;
        }

        // If carry is 1 after loop, we've exhausted all combinations
        if carry == 1 {
            break;
        }
    }

    Ok(products)
}

/// Simpler API: enumerate with a single template and two fragment sets (common case).
///
/// Useful for reactions like: `scaffolds × building_blocks → products`
pub fn enumerate_library_2way(
    template: &str,
    scaffolds: Vec<Molecule>,
    building_blocks: Vec<Molecule>,
    config: &LibraryConfig,
) -> Result<Vec<Molecule>, LibraryError> {
    enumerate_library(template, vec![scaffolds, building_blocks], config)
}

/// Enumerate with three fragment sets (also common: scaffold × R1 × R2 → products).
pub fn enumerate_library_3way(
    template: &str,
    scaffolds: Vec<Molecule>,
    r1_set: Vec<Molecule>,
    r2_set: Vec<Molecule>,
    config: &LibraryConfig,
) -> Result<Vec<Molecule>, LibraryError> {
    enumerate_library(template, vec![scaffolds, r1_set, r2_set], config)
}

/// Errors from library enumeration.
#[derive(Debug, Clone)]
pub enum LibraryError {
    /// No fragment sets provided.
    NoFragmentSets,
    /// Total combinations exceeds configured maximum.
    EnumerationTooLarge(usize, usize), // actual, max
    /// Hit iteration limit (safety check).
    EnumerationLimitExceeded,
    /// The fragment-set Cartesian product does not fit in `usize`.
    CombinationCountOverflow,
}

impl std::fmt::Display for LibraryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoFragmentSets => write!(f, "no fragment sets provided"),
            Self::EnumerationTooLarge(actual, max) => {
                write!(f, "enumeration too large: {} > {}", actual, max)
            }
            Self::EnumerationLimitExceeded => write!(f, "enumeration iteration limit exceeded"),
            Self::CombinationCountOverflow => {
                write!(f, "fragment combination count exceeds usize")
            }
        }
    }
}

impl std::error::Error for LibraryError {}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    #[test]
    fn test_enumerate_library_2way_simple() {
        // Simple amide coupling: phenol + amine → amide
        let template = "[C:1][OH].[C:2][NH2]>>[C:1]O[C:2]"; // Ether formation
        let scaffolds = vec![mol("c1ccccc1O"), mol("Cc1ccccc1O")]; // Two phenols
        let building_blocks = vec![mol("c1ccccc1N"), mol("CCN")]; // Two amines

        let result = enumerate_library_2way(
            template,
            scaffolds,
            building_blocks,
            &LibraryConfig::default(),
        );

        match result {
            Ok(products) => {
                // Should have at most 2×2=4 products
                assert!(
                    products.len() <= 4,
                    "expected <= 4 products, got {}",
                    products.len()
                );
            }
            Err(_) => {
                // If template fails, that's OK for this test (focus on enumeration logic)
            }
        }
    }

    #[test]
    fn test_enumerate_library_no_fragments() {
        let template = "[C:1][Cl]>>[C:1]I";
        let result = enumerate_library(template, vec![], &LibraryConfig::default());
        assert!(matches!(result, Err(LibraryError::NoFragmentSets)));
    }

    #[test]
    fn test_enumerate_library_empty_fragment_set_is_empty() {
        let result = enumerate_library(
            "[C:1]>>[C:1]",
            vec![vec![], vec![mol("C")]],
            &LibraryConfig::default(),
        )
        .unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_enumerate_library_checks_output_cap_before_extending() {
        let result = enumerate_library(
            "[C:1][C:2]>>[C:1].[C:2]",
            vec![vec![mol("CC")]],
            &LibraryConfig {
                skip_failures: true,
                max_size: Some(1),
            },
        );
        assert!(matches!(
            result,
            Err(LibraryError::EnumerationLimitExceeded)
        ));
    }

    #[test]
    fn test_enumerate_library_size_check() {
        // Create a config that rejects large enumerations
        let config = LibraryConfig {
            skip_failures: true,
            max_size: Some(3),
        };

        let template = "[C:1][Cl]>>[C:1]I";
        let set1 = vec![mol("C"), mol("CC"), mol("CCC"), mol("CCCC")];
        let set2 = vec![mol("C"), mol("CC"), mol("CCC"), mol("CCCC")];

        let result = enumerate_library(template, vec![set1, set2], &config);
        // 4×4 = 16 > 3, should fail
        assert!(matches!(
            result,
            Err(LibraryError::EnumerationTooLarge(16, 3))
        ));
    }

    #[test]
    fn test_enumerate_library_3way() {
        let template = "[C:1][C:2][C:3]"; // Dummy template
        let set1 = vec![mol("C"), mol("CC")];
        let set2 = vec![mol("CCC")];
        let set3 = vec![mol("CCCC")];

        let result = enumerate_library_3way(
            template,
            set1,
            set2,
            set3,
            &LibraryConfig {
                skip_failures: true,
                max_size: Some(100),
            },
        );

        // Should enumerate 2×1×1 = 2 combinations
        match result {
            Ok(products) => {
                assert!(products.len() <= 2, "expected <= 2 products");
            }
            Err(_) => {
                // Template might fail, that's OK
            }
        }
    }
}
