//! ESOL aqueous solubility prediction (Delaney 2004).
//!
//! Implements the Delaney linear regression model:
//!
//! ```text
//! logS = 0.16 − 0.63·cLogP − 0.0062·MW + 0.066·RotB − 0.74·AP
//! ```
//!
//! where:
//! - `cLogP` = Crippen LogP
//! - `MW`    = average molecular weight (Da)
//! - `RotB`  = rotatable bond count
//! - `AP`    = aromatic proportion (aromatic heavy atoms / total heavy atoms)
//!
//! Reference: Delaney, J.S. (2004). *ESOL: Estimating aqueous solubility
//! directly from molecular structure*. J. Chem. Inf. Comput. Sci. 44, 1000–1005.

use chematic_core::{Element, Molecule};

use crate::descriptors::logp_crippen;
use crate::descriptors::{heavy_atom_count, molecular_weight, rotatable_bond_count};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Predict aqueous solubility using the Delaney (ESOL) model.
///
/// Returns `log(S)` in log mol/L.  Typical values range from −10 (insoluble)
/// to 0 (fully miscible).
///
/// For molecules with no heavy atoms, returns `0.0`.
pub fn esol_solubility(mol: &Molecule) -> f64 {
    let hac = heavy_atom_count(mol);
    if hac == 0 {
        return 0.0;
    }

    let clogp = logp_crippen(mol);
    let mw = molecular_weight(mol);
    let rot_b = rotatable_bond_count(mol) as f64;
    let ap = aromatic_proportion(mol);

    0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rot_b - 0.74 * ap
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Aromatic proportion: fraction of heavy atoms that are aromatic.
fn aromatic_proportion(mol: &Molecule) -> f64 {
    let hac = heavy_atom_count(mol);
    if hac == 0 {
        return 0.0;
    }
    let aromatic = mol
        .atoms()
        .filter(|(_, atom)| atom.aromatic && atom.element != Element::H)
        .count();
    aromatic as f64 / hac as f64
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn esol(smiles: &str) -> f64 {
        esol_solubility(&parse(smiles).unwrap())
    }

    #[test]
    fn test_water_highly_soluble() {
        // Water should have high predicted solubility (logS close to 0 or positive).
        let s = esol("O");
        assert!(s > -2.0, "water logS should be > -2, got {s:.3}");
    }

    #[test]
    fn test_benzene_moderate_solubility() {
        // Benzene: logS ≈ -1.9 (Delaney model).
        let s = esol("c1ccccc1");
        assert!(
            s > -4.0 && s < 0.0,
            "benzene logS should be in -4..0, got {s:.3}"
        );
    }

    #[test]
    fn test_lipophilic_molecule_low_solubility() {
        // Long alkyl chain → very low solubility.
        let octane = esol("CCCCCCCC"); // octane
        let ethane = esol("CC"); // ethane
        assert!(octane < ethane, "octane should be less soluble than ethane");
    }

    #[test]
    fn test_aromatic_reduces_solubility() {
        // Aromatic compounds tend to be less soluble than aliphatic.
        let benzene = esol("c1ccccc1");
        let cyclohexane = esol("C1CCCCC1");
        assert!(
            benzene < cyclohexane,
            "benzene logS {benzene:.3} should be < cyclohexane {cyclohexane:.3}"
        );
    }

    #[test]
    fn test_returns_finite() {
        assert!(esol("CC(=O)Oc1ccccc1C(=O)O").is_finite()); // aspirin
    }
}
