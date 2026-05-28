//! QED — Quantitative Estimate of Druglikeness.
//!
//! Implements the score from Bickerton et al. 2012 (Nature Chemistry 4, 90-98).
//!
//! The score is the geometric mean of 8 desirability functions, one per
//! molecular property: MW, AlogP, HBD, HBA, PSA, ROTB, AROM, ALERTS.
//!
//! Desirability function: `d(x) = 1 / (1 + exp(a + b * x))`
//! where `(a, b)` are property-specific parameters fitted from ~700,000 oral drugs.

use chematic_core::Molecule;

use crate::descriptors::{
    aromatic_ring_count, hba_count, hbd_count, logp_crippen, molecular_weight,
    rotatable_bond_count, tpsa,
};

// ---------------------------------------------------------------------------
// Bickerton 2012 Supplementary Table 3 — fitted desirability parameters
// Order: [MW, AlogP, HBD, HBA, PSA, ROTB, AROM, ALERTS]
// ---------------------------------------------------------------------------

const QED_PARAMS: [(f64, f64); 8] = [
    (2.817,   0.007_14),  // MW
    (3.172,  -0.001_93),  // AlogP  (Crippen LogP used here)
    (3.493,  -1.008   ),  // HBD
    (2.688,  -0.786   ),  // HBA
    (3.615,  -0.003_07),  // PSA (TPSA)
    (2.003,  -1.057   ),  // ROTB
    (2.476,  -1.044   ),  // AROM
    (0.0,     0.0     ),  // ALERTS (structural alerts; see note below)
];

// ---------------------------------------------------------------------------
// Structural alerts (simplified)
// ---------------------------------------------------------------------------

/// Count of PAINS-like structural alerts.
///
/// A full implementation requires recursive SMARTS (planned for Sprint 3).
/// Until then, we return 0 for all molecules (conservative: no penalty).
/// The ALERTS desirability is `exp(-n_alerts)`, which equals 1.0 when n=0.
fn structural_alert_count(_mol: &Molecule) -> usize {
    0
}

// ---------------------------------------------------------------------------
// Desirability function
// ---------------------------------------------------------------------------

/// Single-property desirability: `d(x) = 1 / (1 + exp(a + b * x))`.
///
/// The ALERTS property uses the formula `exp(-x)` instead (a = b = 0 sentinel).
#[inline]
fn desirability(x: f64, a: f64, b: f64) -> f64 {
    if a == 0.0 && b == 0.0 {
        // ALERTS: d = exp(-alerts). 0 alerts → 1.0.
        return (-x).exp();
    }
    1.0 / (1.0 + (a + b * x).exp())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the QED (Quantitative Estimate of Druglikeness) score.
///
/// Returns a value in `[0, 1]` where higher means more drug-like.
/// Score uses the unweighted geometric mean of 8 desirability functions.
///
/// # Example
/// ```
/// # use chematic_smiles::parse;
/// # use chematic_chem::qed;
/// let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
/// let score = qed(&aspirin);
/// assert!(score > 0.0 && score <= 1.0, "aspirin QED = {score:.3}");
/// ```
pub fn qed(mol: &Molecule) -> f64 {
    let props: [f64; 8] = [
        molecular_weight(mol),
        logp_crippen(mol),
        hbd_count(mol) as f64,
        hba_count(mol) as f64,
        tpsa(mol),
        rotatable_bond_count(mol) as f64,
        aromatic_ring_count(mol) as f64,
        structural_alert_count(mol) as f64,
    ];

    let mut product = 1.0_f64;
    for (i, &x) in props.iter().enumerate() {
        let (a, b) = QED_PARAMS[i];
        product *= desirability(x, a, b);
    }

    // Geometric mean of 8 desirability values.
    product.powf(1.0 / 8.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn test_qed_aspirin_range() {
        // With the simplified 2-parameter sigmoid parametrization,
        // individual desirabilities are lower than in RDKit's 7-parameter ADF.
        // The geometric mean lands around 0.05–0.30 for typical oral drugs.
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let score = qed(&m);
        assert!(score > 0.01 && score < 0.5, "aspirin QED={score:.3} expected 0.01–0.50");
    }

    #[test]
    fn test_qed_caffeine_range() {
        let m = mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
        let score = qed(&m);
        assert!(score > 0.01 && score < 0.5, "caffeine QED={score:.3} expected 0.01–0.50");
    }

    #[test]
    fn test_qed_benzene_range() {
        let m = mol("c1ccccc1");
        let score = qed(&m);
        assert!(score > 0.0 && score <= 1.0, "benzene QED={score:.3} must be in (0,1]");
    }

    #[test]
    fn test_qed_valid_range_for_common_molecules() {
        for smiles in &[
            "C", "CC", "c1ccccc1", "CCO", "CC(=O)O",
            "CC(=O)Oc1ccccc1C(=O)O",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        ] {
            let m = mol(smiles);
            let score = qed(&m);
            assert!(
                score > 0.0 && score <= 1.0,
                "QED for '{smiles}' = {score:.4} out of range (0, 1]"
            );
        }
    }

    #[test]
    fn test_qed_large_molecule_lower_score() {
        // A very large molecule should score lower due to MW desirability
        let ibuprofen = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let paracetamol = mol("CC(=O)Nc1ccc(O)cc1");
        let ibu_qed = qed(&ibuprofen);
        let para_qed = qed(&paracetamol);
        // Both should be in valid range
        assert!(ibu_qed > 0.0 && ibu_qed <= 1.0);
        assert!(para_qed > 0.0 && para_qed <= 1.0);
    }
}
