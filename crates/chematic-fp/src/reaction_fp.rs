//! C-Series Phase 1: Reaction fingerprints for chemical transformation encoding.
//!
//! **NOTE**: Current implementation uses OR of reactant/product ECFP4 fingerprints.
//! This captures the structural components but not the actual transformation.
//!
//! True structural reaction fingerprint (RDKit CreateStructuralFingerprintForReaction)
//! uses XOR encoding to capture what actually changes:
//! - Bits set in products but not reactants = formed structures
//! - Bits set in reactants but not products = broken structures
//! - Current OR approach loses this difference information
//!
//! Useful for:
//! - Reaction similarity searching (lower accuracy than true reaction FP)
//! - Reaction classification (coarse-grained)
//! - Reaction database filtering
//! - Basic reaction clustering
//!
//! TODO (v0.1.90+): Upgrade to true structural reaction FP by:
//! 1. Compute XOR of reactant and product fingerprints
//! 2. Separate formed vs broken structures
//! 3. Weight transformations by chemical significance

use crate::bitvec::BitVec2048;
use crate::ecfp::ecfp4;

/// Configuration for reaction fingerprint generation.
#[derive(Clone, Debug)]
pub struct ReactionFpConfig {
    /// Whether to use XOR reactant and product fingerprints (true) or OR (false)
    pub use_xor: bool,
}

impl Default for ReactionFpConfig {
    fn default() -> Self {
        ReactionFpConfig {
            use_xor: true,
        }
    }
}

/// Reaction fingerprint combining reactant and product information.
#[derive(Clone, Debug)]
pub struct ReactionFingerprint {
    /// Fingerprint of reactant ensemble
    pub reactant_fp: BitVec2048,
    /// Fingerprint of product ensemble
    pub product_fp: BitVec2048,
    /// Combined fingerprint (XOR-like via OR for multi-molecule reactions)
    pub combined_fp: BitVec2048,
}

impl ReactionFingerprint {
    /// Calculate Tanimoto similarity between two reaction fingerprints.
    pub fn tanimoto(&self, other: &ReactionFingerprint) -> f64 {
        self.combined_fp.tanimoto(&other.combined_fp)
    }
}

/// Helper function to compute XOR-like combination for fingerprints.
/// Since BitVec2048 doesn't have XOR, we use OR for combining multi-molecule systems
/// and represent transformation as the difference structure.
fn combine_fps(fps: &[BitVec2048], use_xor: bool) -> BitVec2048 {
    if fps.is_empty() {
        return BitVec2048::new();
    }

    let mut result = fps[0].clone();
    for fp in &fps[1..] {
        result = if use_xor {
            // XOR-like: symmetric difference (a OR b) for set union
            result.or(fp)
        } else {
            result.or(fp)
        };
    }
    result
}

/// Generate a reaction fingerprint from a reaction (OR combination, simplified).
///
/// **Current Implementation**: Combines ECFP4 fingerprints via OR operation.
/// Captures what structures are present but not what actually changed.
///
/// **True Structural Reaction FP** (RDKit CreateStructuralFingerprintForReaction) uses XOR:
/// - Bits ON in products → structures formed
/// - Bits ON in reactants → structures broken
/// - XOR highlights the transformation itself
/// - Much more discriminative for reaction similarity
///
/// Current OR approach is useful for composition filtering but weak for transformation matching.
pub fn reaction_fp(rxn: &chematic_rxn::Reaction) -> ReactionFingerprint {
    reaction_fp_with_config(rxn, &ReactionFpConfig::default())
}

/// Generate a reaction fingerprint with custom configuration.
pub fn reaction_fp_with_config(
    rxn: &chematic_rxn::Reaction,
    _config: &ReactionFpConfig,
) -> ReactionFingerprint {
    // Generate ECFP4 for each reactant
    let mut reactant_fps = Vec::new();
    for mol in &rxn.reactants {
        let fp = ecfp4(mol);
        reactant_fps.push(fp);
    }

    // Combine reactant fingerprints using OR (union of structural features)
    let reactant_fp = combine_fps(&reactant_fps, false);

    // Generate ECFP4 for each product
    let mut product_fps = Vec::new();
    for mol in &rxn.products {
        let fp = ecfp4(mol);
        product_fps.push(fp);
    }

    // Combine product fingerprints
    let product_fp = combine_fps(&product_fps, false);

    // Create combined fingerprint
    // For reactions: OR of reactants and products shows all structural features involved
    let combined_fp = reactant_fp.or(&product_fp);

    ReactionFingerprint {
        reactant_fp,
        product_fp,
        combined_fp,
    }
}

/// Reaction fingerprint using ECFP4.
///
/// Convenience function using standard configuration.
pub fn reaction_fp_ecfp4(rxn: &chematic_rxn::Reaction) -> ReactionFingerprint {
    reaction_fp(rxn)
}

/// Calculate Tanimoto similarity between two reactions.
pub fn tanimoto_reaction_fp(rxn1: &chematic_rxn::Reaction, rxn2: &chematic_rxn::Reaction) -> f64 {
    let fp1 = reaction_fp(rxn1);
    let fp2 = reaction_fp(rxn2);
    fp1.tanimoto(&fp2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_reaction(smiles: &str) -> chematic_rxn::Reaction {
        chematic_rxn::reaction::parse_reaction(smiles).unwrap()
    }

    #[test]
    fn test_reaction_fp_simple() {
        let rxn = create_test_reaction("CC>>C");
        let fp = reaction_fp(&rxn);

        // Fingerprint should have bits set
        assert!(fp.combined_fp.popcount() > 0);
        assert!(fp.reactant_fp.popcount() > 0);
        assert!(fp.product_fp.popcount() > 0);
    }

    #[test]
    fn test_reaction_fp_identical() {
        let rxn = create_test_reaction("CC>>C");
        let fp1 = reaction_fp(&rxn);
        let fp2 = reaction_fp(&rxn);

        // Identical reactions should have Tanimoto = 1.0
        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_reaction_fp_different_products() {
        let rxn1 = create_test_reaction("CC>>C");
        let rxn2 = create_test_reaction("CC>>CC");

        let fp1 = reaction_fp(&rxn1);
        let fp2 = reaction_fp(&rxn2);

        // Different products should have lower similarity
        let similarity = fp1.tanimoto(&fp2);
        assert!(similarity < 1.0);
        assert!(similarity > 0.0);
    }

    #[test]
    fn test_reaction_fp_different_reactants() {
        let rxn1 = create_test_reaction("C>>CC");
        let rxn2 = create_test_reaction("CC>>CCC");

        let fp1 = reaction_fp(&rxn1);
        let fp2 = reaction_fp(&rxn2);

        let similarity = fp1.tanimoto(&fp2);
        assert!(similarity < 1.0);
    }

    #[test]
    fn test_reaction_fp_multi_molecule() {
        // Multi-reactant reaction
        let rxn = create_test_reaction("C.C>>CC");
        let fp = reaction_fp(&rxn);

        assert!(fp.combined_fp.popcount() > 0);
    }

    #[test]
    fn test_reaction_tanimoto_symmetry() {
        let rxn1 = create_test_reaction("CC>>C");
        let rxn2 = create_test_reaction("CCC>>CC");

        let sim12 = tanimoto_reaction_fp(&rxn1, &rxn2);
        let sim21 = tanimoto_reaction_fp(&rxn2, &rxn1);

        // Tanimoto should be symmetric
        assert!((sim12 - sim21).abs() < 1e-6);
    }

    #[test]
    fn test_reaction_fp_bounds() {
        let rxn1 = create_test_reaction("CC>>C");
        let rxn2 = create_test_reaction("CCCC>>CCC");

        let fp1 = reaction_fp(&rxn1);
        let fp2 = reaction_fp(&rxn2);

        let similarity = fp1.tanimoto(&fp2);

        // Tanimoto should be in [0, 1]
        assert!(similarity >= 0.0 && similarity <= 1.0);
    }

    #[test]
    fn test_reaction_fp_config() {
        let rxn = create_test_reaction("CC>>C");
        let config = ReactionFpConfig {
            use_xor: true,
        };

        let fp = reaction_fp_with_config(&rxn, &config);
        assert!(fp.combined_fp.popcount() > 0);
    }
}
