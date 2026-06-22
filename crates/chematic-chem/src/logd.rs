//! LogD — pH-dependent apparent lipophilicity (distribution coefficient).
//!
//! LogD accounts for ionisation: at physiological pH many drugs carry a
//! net charge that reduces their membrane permeability compared to the
//! neutral form captured by LogP.
//!
//! # Model
//! This implementation uses a simplified Henderson-Hasselbalch correction
//! based on the molecule's formal charge and estimated pKa.  Because full
//! pKa prediction is a research-grade problem, the pKa is approximated from
//! Crippen atom contributions and formal charge counts:
//!
//! - **Acids** (net negative charge at high pH): estimated pKa ≈ 4.0
//! - **Bases** (net positive charge at low pH): estimated pKa ≈ 9.5
//! - **Zwitterions** or **neutral**: no correction applied, LogD ≈ LogP
//!
//! For a single ionisable group the Henderson-Hasselbalch formula gives:
//! - Acid: `LogD = LogP + log(1 / (1 + 10^(pH − pKa)))`
//! - Base: `LogD = LogP + log(1 / (1 + 10^(pKa − pH)))`
//!
//! **Accuracy**: ±1 log unit vs experiment; suitable for rank-ordering and
//! qualitative rule-of-5 assessments.  For high-accuracy LogD a dedicated
//! pKa model is required.

use crate::descriptors::logp_crippen;
use crate::pka::{pka_acid, pka_base};
use chematic_core::{Molecule, implicit_hcount};

// ---------------------------------------------------------------------------
// Estimated pKa constants
// ---------------------------------------------------------------------------

/// Estimated pKa for typical carboxylic acid / phenol (acidic) groups.
const PKA_ACID: f64 = 4.0;
/// Estimated pKa for typical amine (basic) groups.
const PKA_BASE: f64 = 9.5;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Classify the molecule's ionisation character from its formal charge sum
/// and the presence of heteroatoms that commonly ionise.
///
/// Returns `(is_acid, is_base)`.  Both can be false for neutral molecules.
fn ionisation_class(mol: &Molecule) -> (bool, bool) {
    // Count atoms that are typical acid/base sites.
    let mut acid_sites = 0u32; // COOH, SH, phenol OH
    let mut base_sites = 0u32; // primary/secondary/tertiary amines

    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let h = implicit_hcount(mol, idx) as u32;
        match an {
            // Oxygen: O-H with a double-bonded O neighbor → COOH / phenol
            8 if h > 0 => {
                let has_carbonyl = mol.neighbors(idx).any(|(nb, _)| {
                    mol.atom(nb).element.atomic_number() == 6
                        && mol.neighbors(nb).any(|(nb2, bid2)| {
                            nb2 != idx
                                && mol.atom(nb2).element.atomic_number() == 8
                                && mol.bond(bid2).order == chematic_core::BondOrder::Double
                        })
                });
                if has_carbonyl {
                    acid_sites += 1;
                }
            }
            // Sulfur: S-H
            16 if h > 0 => {
                acid_sites += 1;
            }
            // Nitrogen: N-H or NR3 with sp3 character → base
            7 if atom.charge >= 0 => {
                let heavy_deg = mol.neighbors(idx).count() as u32;
                // sp3 N with at most 3 heavy-atom neighbors can be protonated.
                if heavy_deg <= 3 && !atom.aromatic {
                    base_sites += 1;
                }
            }
            _ => {}
        }
    }

    (acid_sites > 0, base_sites > 0)
}

/// Henderson-Hasselbalch correction to LogP.
fn hh_correction(logp: f64, ph: f64, pka: f64, is_acid: bool) -> f64 {
    if is_acid {
        // LogD_acid = LogP − log(1 + 10^(pH − pKa))
        logp - (1.0 + 10f64.powf(ph - pka)).log10()
    } else {
        // LogD_base = LogP − log(1 + 10^(pKa − pH))
        logp - (1.0 + 10f64.powf(pka - ph)).log10()
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Predict LogD at the given pH using a simplified Henderson-Hasselbalch model.
///
/// `ph` should be in the range 0–14; typical physiological pH is 7.4.
///
/// Returns LogD (dimensionless, same scale as LogP).
pub fn logd_simple(mol: &Molecule, ph: f64) -> f64 {
    let logp = logp_crippen(mol);
    logd_from_logp(logp, mol, ph)
}

/// Like [`logd_simple`] but accepts a pre-computed LogP value.
///
/// Use this when LogP has already been computed to avoid recomputing the
/// 117-pattern Crippen SMARTS matching a second time:
///
/// ```ignore
/// let logp = logp_crippen(&mol);
/// let logd = logd_from_logp(logp, &mol, 7.4);
/// ```
pub fn logd_from_logp(logp: f64, mol: &Molecule, ph: f64) -> f64 {
    let (is_acid_legacy, is_base_legacy) = ionisation_class(mol);

    // Use pKa module values when available; fall back to legacy constants.
    let pka_a = pka_acid(mol).unwrap_or(PKA_ACID);
    let pka_b = pka_base(mol).unwrap_or(PKA_BASE);
    let is_acid = is_acid_legacy || pka_acid(mol).is_some();
    let is_base = is_base_legacy || pka_base(mol).is_some();

    match (is_acid, is_base) {
        (true, false) => hh_correction(logp, ph, pka_a, true),
        (false, true) => hh_correction(logp, ph, pka_b, false),
        // Zwitterion: apply both corrections (approximate).
        (true, true) => {
            let acid_corr = hh_correction(logp, ph, pka_a, true);
            let base_corr = hh_correction(logp, ph, pka_b, false);
            (acid_corr + base_corr) / 2.0
        }
        // Neutral: LogD ≈ LogP
        (false, false) => logp,
    }
}

/// Compute LogD at a series of pH values.
///
/// Returns a `Vec<(ph, logd)>` with `steps` evenly-spaced points from
/// `ph_start` to `ph_end` (inclusive).
pub fn logd_profile(mol: &Molecule, ph_start: f64, ph_end: f64, steps: usize) -> Vec<(f64, f64)> {
    if steps == 0 {
        return Vec::new();
    }
    if steps == 1 {
        let ph = (ph_start + ph_end) / 2.0;
        return vec![(ph, logd_simple(mol, ph))];
    }
    let step_size = (ph_end - ph_start) / (steps - 1) as f64;
    (0..steps)
        .map(|i| {
            let ph = ph_start + i as f64 * step_size;
            (ph, logd_simple(mol, ph))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_neutral_logd_equals_logp() {
        // Benzene has no ionisable groups → LogD = LogP at any pH.
        let m = mol("c1ccccc1");
        let logp = logp_crippen(&m);
        let logd = logd_simple(&m, 7.4);
        assert!(
            (logd - logp).abs() < 0.01,
            "neutral: logD={logd:.3} logP={logp:.3}"
        );
    }

    #[test]
    fn test_acid_logd_lower_at_high_ph() {
        // Acetic acid: at pH >> pKa (4.0), logD << logP.
        let m = mol("CC(=O)O");
        let logp = logp_crippen(&m);
        let logd_acid = logd_simple(&m, 9.0); // well above pKa
        assert!(
            logd_acid < logp,
            "acid at high pH: logD={logd_acid:.3} should be < logP={logp:.3}"
        );
    }

    #[test]
    fn test_logd_profile_length() {
        let m = mol("CC");
        let profile = logd_profile(&m, 1.0, 14.0, 5);
        assert_eq!(profile.len(), 5);
        assert!((profile[0].0 - 1.0).abs() < 1e-9);
        assert!((profile[4].0 - 14.0).abs() < 1e-9);
    }

    #[test]
    fn test_logd_profile_zero_steps() {
        let m = mol("CC");
        assert!(logd_profile(&m, 0.0, 14.0, 0).is_empty());
    }
}
