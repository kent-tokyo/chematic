//! Multi-mode canonical SMILES — analogous to OCL's 5-mode IDcode canonicalization.

use chematic_core::Molecule;
use chematic_smiles::canonical_smiles;

use crate::standardize::{remove_isotopes, remove_stereo, uncharge};
use crate::tautomer::canonical_tautomer;

/// Canonicalization mode for [`canonical_smiles_mode`].
///
/// Inspired by the five modes in OCL's IDcode canonicalization, but outputs
/// canonical SMILES rather than a binary IDcode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalMode {
    /// Standard canonical SMILES — stereo, charges, and isotopes preserved.
    Normal,
    /// Canonical SMILES with all stereochemistry removed
    /// (`@`/`@@`/`/`/`\` and wedge bonds cleared).
    NoStereo,
    /// Backbone topology only: charges, isotopes, and stereo removed.
    /// Element symbols and bond orders are preserved.
    Backbone,
    /// Canonical tautomer, then canonical SMILES.
    /// Useful for tautomer-invariant matching.
    Tautomer,
    /// Canonical tautomer + stereo removal.
    NoStereoTautomer,
}

/// Return the canonical SMILES for `mol` under the given [`CanonicalMode`].
///
/// # Examples
///
/// ```ignore
/// use chematic_chem::canonical::{CanonicalMode, canonical_smiles_mode};
/// let mol = parse("[C@@H](N)(C)C(=O)O").unwrap(); // L-alanine
/// assert_eq!(canonical_smiles_mode(&mol, CanonicalMode::Normal),      "[C@@H](N)(C)C(=O)O");
/// assert_eq!(canonical_smiles_mode(&mol, CanonicalMode::NoStereo),    "CC(N)C(=O)O");
/// assert_eq!(canonical_smiles_mode(&mol, CanonicalMode::Backbone),    "CC(N)C(=O)O");
/// ```
pub fn canonical_smiles_mode(mol: &Molecule, mode: CanonicalMode) -> String {
    match mode {
        CanonicalMode::Normal => canonical_smiles(mol),
        CanonicalMode::NoStereo => canonical_smiles(&remove_stereo(mol)),
        CanonicalMode::Backbone => {
            let m = remove_stereo(mol);
            let m = remove_isotopes(&m);
            let m = uncharge(&m);
            canonical_smiles(&m)
        }
        CanonicalMode::Tautomer => canonical_smiles(&canonical_tautomer(mol)),
        CanonicalMode::NoStereoTautomer => {
            canonical_smiles(&remove_stereo(&canonical_tautomer(mol)))
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn normal_preserves_stereo() {
        let m = mol("[C@@H](N)(C)C(=O)O"); // L-alanine
        let smi = canonical_smiles_mode(&m, CanonicalMode::Normal);
        // Canonical SMILES should contain stereo marker
        assert!(
            smi.contains("@@") || smi.contains("@"),
            "Normal mode should preserve stereo: {smi}"
        );
    }

    #[test]
    fn nostereo_removes_stereo() {
        let m = mol("[C@@H](N)(C)C(=O)O");
        let with = canonical_smiles_mode(&m, CanonicalMode::Normal);
        let without = canonical_smiles_mode(&m, CanonicalMode::NoStereo);
        assert!(
            !without.contains('@') && !without.contains('/'),
            "NoStereo should strip stereo: {without}"
        );
        assert_ne!(
            with, without,
            "Normal and NoStereo should differ for chiral mol"
        );
    }

    #[test]
    fn nostereo_enantiomers_identical() {
        // L- and D-alanine differ only in stereo; NoStereo must give same string
        let l = mol("[C@@H](N)(C)C(=O)O");
        let d = mol("[C@H](N)(C)C(=O)O");
        let l_ns = canonical_smiles_mode(&l, CanonicalMode::NoStereo);
        let d_ns = canonical_smiles_mode(&d, CanonicalMode::NoStereo);
        assert_eq!(
            l_ns, d_ns,
            "enantiomers should be identical in NoStereo mode"
        );
    }

    #[test]
    fn backbone_strips_charge_and_isotope() {
        // [13C] isotope and N+ charge
        let m = mol("[13CH3][N+](C)(C)C"); // isotope-labelled and charged
        let backbone = canonical_smiles_mode(&m, CanonicalMode::Backbone);
        assert!(
            !backbone.contains("13"),
            "Backbone should strip isotope: {backbone}"
        );
        assert!(
            !backbone.contains('+'),
            "Backbone should strip charge: {backbone}"
        );
        assert!(
            !backbone.contains('@'),
            "Backbone should strip stereo: {backbone}"
        );
    }

    #[test]
    fn tautomer_mode_normalises() {
        // Keto form and enol form should give same Tautomer-mode SMILES
        let keto = mol("CC(=O)C"); // acetone
        let enol = mol("CC(O)=C"); // enol — unusual tautomer, but tests the mode
        let k_tau = canonical_smiles_mode(&keto, CanonicalMode::Tautomer);
        let e_tau = canonical_smiles_mode(&enol, CanonicalMode::Tautomer);
        // canonical_tautomer should prefer one form for both
        assert_eq!(
            k_tau, e_tau,
            "keto/enol should collapse to same tautomer: '{k_tau}' vs '{e_tau}'"
        );
    }

    #[test]
    fn nostereo_tautomer_combines_both() {
        // Chiral keto compound — NoStereoTautomer strips both
        let m = mol("[C@@H](N)(C)C(=O)O");
        let smi = canonical_smiles_mode(&m, CanonicalMode::NoStereoTautomer);
        assert!(
            !smi.contains('@'),
            "NoStereoTautomer should strip stereo: {smi}"
        );
        assert!(smi.is_ascii(), "result should be valid SMILES: {smi}");
    }

    #[test]
    fn all_modes_return_valid_smiles() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        for mode in [
            CanonicalMode::Normal,
            CanonicalMode::NoStereo,
            CanonicalMode::Backbone,
            CanonicalMode::Tautomer,
            CanonicalMode::NoStereoTautomer,
        ] {
            let smi = canonical_smiles_mode(&m, mode);
            assert!(!smi.is_empty(), "{mode:?} returned empty SMILES");
            // Should be re-parseable
            assert!(
                parse(&smi).is_ok(),
                "{mode:?} produced unparseable SMILES: {smi}"
            );
        }
    }
}
