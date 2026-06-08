//! Green chemistry metrics for reaction efficiency assessment.
//!
//! Implements the core atom-economy and mass-efficiency metrics defined by
//! Trost (1991) and Sheldon (1992):
//!
//! - [`atom_economy`] — percentage of reactant atoms incorporated into the product.
//! - [`e_factor`] — kg waste per kg product (lower is greener).
//! - [`pmi_rxn`] — Process Mass Intensity: total mass in / product mass out.
//!
//! All functions that need molecular weights compute them internally using
//! average atomic masses so that no external crate dependency is required.

use chematic_core::{Molecule, implicit_hcount};

use crate::reaction::Reaction;

// ---------------------------------------------------------------------------
// Internal MW calculation
// ---------------------------------------------------------------------------

/// Average atomic mass (Da) for the most common elements.
fn avg_mass(an: u8) -> f64 {
    match an {
        1 => 1.008,    // H
        2 => 4.003,    // He
        3 => 6.941,    // Li
        4 => 9.012,    // Be
        5 => 10.811,   // B
        6 => 12.011,   // C
        7 => 14.007,   // N
        8 => 15.999,   // O
        9 => 18.998,   // F
        11 => 22.990,  // Na
        12 => 24.305,  // Mg
        13 => 26.982,  // Al
        14 => 28.086,  // Si
        15 => 30.974,  // P
        16 => 32.065,  // S
        17 => 35.453,  // Cl
        19 => 39.098,  // K
        20 => 40.078,  // Ca
        26 => 55.845,  // Fe
        29 => 63.546,  // Cu
        30 => 65.38,   // Zn
        35 => 79.904,  // Br
        53 => 126.904, // I
        n => n as f64,
    }
}

fn molecular_weight(mol: &Molecule) -> f64 {
    let mut mw = 0.0f64;
    for (idx, atom) in mol.atoms() {
        mw += avg_mass(atom.element.atomic_number());
        mw += implicit_hcount(mol, idx) as f64 * 1.008;
    }
    mw
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the atom economy of a reaction (%).
///
/// Atom economy = Σ(MW of desired products) / Σ(MW of all reactants) × 100.
///
/// Only `rxn.products` count as "desired product"; `rxn.agents` (solvents,
/// catalysts) are excluded from both numerator and denominator, consistent
/// with Trost's original definition.
///
/// Returns `0.0` if there are no reactants.
pub fn atom_economy(rxn: &Reaction) -> f64 {
    let reactant_mw: f64 = rxn.reactants.iter().map(molecular_weight).sum();
    let product_mw: f64 = rxn.products.iter().map(molecular_weight).sum();
    if reactant_mw == 0.0 {
        return 0.0;
    }
    (product_mw / reactant_mw) * 100.0
}

/// Compute the E-factor (Environmental factor) of a reaction.
///
/// E-factor = waste_mass / product_mass.
///
/// A lower E-factor means less waste per unit product (greener process).
/// Typical values: bulk chemicals ~1–5, fine chemicals ~5–50,
/// pharmaceuticals ~25–100+.
///
/// `waste_mass` and `product_mass` are in consistent units (e.g. both kg).
///
/// Returns `f64::INFINITY` if `product_mass` is zero.
pub fn e_factor(waste_mass: f64, product_mass: f64) -> f64 {
    if product_mass == 0.0 {
        return f64::INFINITY;
    }
    waste_mass / product_mass
}

/// Compute the Process Mass Intensity (PMI) of a reaction.
///
/// PMI = Σ(all input masses) / product_mass.
///
/// `all_masses` should include all inputs: reactants, solvents, reagents,
/// water, catalysts.  The theoretical minimum PMI is 1.0 (100 % atom economy
/// with no solvent waste).
///
/// Returns `f64::INFINITY` if `product_mass` is zero.
pub fn pmi_rxn(all_masses: &[f64], product_mass: f64) -> f64 {
    if product_mass == 0.0 {
        return f64::INFINITY;
    }
    all_masses.iter().sum::<f64>() / product_mass
}

/// Compute the reaction mass efficiency (RME, %).
///
/// RME = product_mass / Σ(reactant_masses) × 100.
///
/// A higher RME means a greater fraction of reactant mass ends up in the
/// product.  Returns `0.0` if `reactant_masses` is empty or sums to zero.
pub fn reaction_mass_efficiency(reactant_masses: &[f64], product_mass: f64) -> f64 {
    let total: f64 = reactant_masses.iter().sum();
    if total == 0.0 {
        return 0.0;
    }
    (product_mass / total) * 100.0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reaction::parse_reaction;

    #[test]
    fn test_atom_economy_100_percent() {
        // A → A: all atoms end up in the product.
        let rxn = parse_reaction("CC>>CC").unwrap();
        let ae = atom_economy(&rxn);
        assert!(
            (ae - 100.0).abs() < 0.1,
            "AA→AA should be 100% atom economy, got {ae:.2}"
        );
    }

    #[test]
    fn test_atom_economy_less_than_100() {
        // Acylation with HCl byproduct: only the acyl product is "desired".
        // C.O=O>>O=C=O represents CH4 + O2 → CO2 (incomplete combustion, ignoring H2O).
        // Reactant MW ≈ CH4 (16) + O2 (32) = 48; Product MW = CO2 ≈ 44 → ~91.7%
        let rxn = parse_reaction("C.O=O>>O=C=O").unwrap();
        let ae = atom_economy(&rxn);
        assert!(ae > 0.0, "atom economy should be > 0%");
        assert!(
            ae < 100.0,
            "partial product should give < 100% atom economy, got {ae:.2}"
        );
    }

    #[test]
    fn test_atom_economy_no_reactants() {
        let rxn = parse_reaction(">>CC").unwrap();
        assert_eq!(atom_economy(&rxn), 0.0);
    }

    #[test]
    fn test_e_factor_basic() {
        // 10 kg waste, 2 kg product → E-factor = 5
        assert!((e_factor(10.0, 2.0) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_e_factor_zero_product() {
        assert!(e_factor(1.0, 0.0).is_infinite());
    }

    #[test]
    fn test_pmi_minimum() {
        // PMI = 1 means all input became product (theoretical minimum)
        assert!((pmi_rxn(&[5.0], 5.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_pmi_zero_product() {
        assert!(pmi_rxn(&[1.0, 2.0], 0.0).is_infinite());
    }

    #[test]
    fn test_rme_basic() {
        // 10 kg reactants → 8 kg product: 80% RME
        assert!((reaction_mass_efficiency(&[10.0], 8.0) - 80.0).abs() < 1e-9);
    }

    #[test]
    fn test_rme_empty() {
        assert_eq!(reaction_mass_efficiency(&[], 5.0), 0.0);
    }
}
